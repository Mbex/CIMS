## FUNCTIONS

class Peaklist(object):

    def __init__(self):
        
        self.elements = {
            "A": 0, "D": 0, "E": 0, "G": 0, "H": 0, "J": 0, "L": 0, "M": 0,
            "N": 0, "O": 0, "Q": 0, "R": 0, "T": 0, "V": 0, "W": 0, "X": 0,
            "Z": 0, 'C': 12.0000, 'B': 11.009305, 'H': 1.007825, 'O': 15.994915,
            'N': 14.003074, 'I': 126.904477, 'Cl': 34.968853, 'Br': 78.918336,
            'S': 31.972072, 'F': 18.998403}
        
        self.kendrick_bases = ["CH2", "O", "OH", "H2", "CO", "CHO", "NO2", "NO3"]

           
    def remove_punctuation(self, string):
        
        """remove minus from string"""
        
        return string.replace("-","").replace("'","")
        
        
    def _Get_molecule_params(self, molecule):
        
        """Hit API to get molecule parameters"""
        
         import urllib, urllib2, json

         params = urllib.urlencode({'mf': molecule})
        
         response = urllib2.urlopen('http://chemcalc.org/chemcalc/mf',params)
         data = json.loads(response.read())
        
         return data

    
    def mass_calculator(self, moiety):

        """Calculates exact mass for a given formula (string)."""
    
        data = self._Get_molecule_params(moiety)
        
        return data['em']
    
    
    def count_elements(self, moiety):

        """Returns dict with keys of elements and values of number of elements."""

        import urllib, urllib2, json

        ans = collections.OrderedDict({})
    
        data = self._Get_molecule_params(moiety)

        contents = data['parts'][0]['ea']

        for block in contents:
            ans[block['element']] = block['number']
            
        return ans    

    
    def kendrick_bases_apart(self, amu1, amu2, kendrick_base):
    
        """How many kendrick_base amu1 and amu2 are different from each other."""

        ans = (np.round(amu1) - np.round(amu2)) / np.round(kendrick_base)
   
        if ans % 1 != 0:
            
            ans = "not integer bases apart"
        
        return ans 
    

    def fAmuApart(self, amu1, amu2, kendrick_base):

        """
        Returns True if difference between amu1 and
        amu2 when divided by moiety is equal to 0.
        """

        return abs(amu1 - amu2) % np.round(kendrick_base) == 0


    def fWithinError(self, KMD1, KMD1e, KMD2, KMD2e):

        """Returns True if KMD1 +/- KMD1e is within KMD2 +/- KMD2e."""

        KMD1plus = KMD1 + KMD1e
        KMD1minus = KMD1 - KMD1e

        KMD2plus = KMD2 + KMD2e
        KMD2minus = KMD2 - KMD2e

        #return (KMD2minus < KMD1minus < KMD2plus) or (KMD2minus < KMD1plus < KMD2plus)
        return KMD2minus < KMD1 < KMD2plus


    def load_peak_list(self, directory, fname, **kwargs):

        """Returns dataframe with peaklist information."""

        import pandas as pd

        df = pd.read_csv(directory+fname, header=0, **kwargs)
        df = df[['ion','sumFormula','x0']]
        #df.columns = ['name', 'formula', 'exactMass']
        df.sort_values('sumFormula').reset_index(drop=True)

        return df

    
    def rename_unknowns(self, df, col, pattern="unknown"):

        "used for renaming unknowns. reset numbers at end of unknown pattern"

        names = peakList[col]
        count = 1
        for i, name in enumerate(names):
            if pattern in name:
                newname = pattern + "%04d" % (count,)
                peakList.loc[i, col] = newname
                count = count + 1

        return df


    def calc_mass_defect(self, df, reagent_ion_exact_mass):

        """
        Makes integer and exact mass columns with
        and without reagent ion then
        calculates the mass defect from them.
        """

        import numpy as np

        regaentIonIntegerMass = np.round(reagent_ion_exact_mass)
        df['integerMass'] = np.round(df['x0'])
          
        # Duplicate - prep for 'removing reagent ion' columns
        df['exactMassNoRI'] = df['x0']
        df['integerMassNoRI'] = df['integerMass']
        # Remove reagent ion
        df.loc[df['integerMass'] > regaentIonIntegerMass, 'integerMassNoRI'] = df['integerMass'] - regaentIonIntegerMass
        df.loc[df['x0'] > regaentIonIntegerMass, 'exactMassNoRI'] = df['x0'] - regaentIonIntegerMass 
        # mass defect
        df['massDefect'] = df['exactMassNoRI'] - df['integerMassNoRI']

        return df  
    

    def calc_kendrick_mass_defect(self, df, kendrick_base, kendrick_base_mass):

        """
        Calculate the kendrick mass defect for the passed kendrick base.
        Requries that calc_mass_defect has been run first.
        """

        import numpy as np

        # Calculate Kendrick Exact Mass (KEM) and Kendrick Mass Defect (KMD) for speciesList
        kem = "KEM_"+kendrick_base
        kmd = "KMD_"+kendrick_base
        kmdError = "KMD_"+kendrick_base+"_error"

        normalise = (np.round(kendrick_base_mass) / kendrick_base_mass)
        speciesList[kem] = df['exactMassNoRI'] * normalise
        speciesList[kmd] = df['integerMassNoRI'] - df[kem]
        speciesList[kmdError] = df['exactMassNoRI'] * (20/1e6)

        speciesList.loc[:,"KMD_"+kendrick_base+"_matches"] = None

        return speciesList


    def match_peaks_on_kmd(self, df, kendrick_base, kendrick_base_mass):

        # amu apart
        amuApart1 = np.array(map(lambda y :
                                 map(lambda x:
                                     self.fAmuApart(x, y, np.round(kendrick_base_mass)), df['integerMass']
                                     ),
                                speciesList['integerMass']
                                )
                             )

        # within error
        withinError1 = np.array(map(lambda y, z:
                                    map(lambda x, w:
                                        self.fWithinError(x, w, y, z), df["KMD_"+kendrick_base], df["KMD_"+kendrick_base+"_error"]
                                        ),
                                    df["KMD_"+kendrick_base], df["KMD_"+kendrick_base+"_error"]
                                    )
                                )

        # 2d arrays for each entry we have a list of bools saying if the index of that list is a match or not
        isMatch = amuApart1 & withinError1    

        # put matches into one array
        matches = np.array(map(lambda y :
                           [i for i, x in enumerate(list(isMatch[y])) if x],
                           range(isMatch.shape[1])
                          ))

        # For every entry append the matching list into the matches column
        for i, entry in enumerate(matches): 

            name = df.ix[i,'ion']
            match = map(lambda x: df.ix[x,'ion'], entry)
            if name in match:
                match.remove(name)
            toAppend = str(match).replace('[',"").replace(']',"")
            df.loc[i,"KMD_"+kendrick_base+"_matches"] = toAppend  


        return df
    
    
    def calc_unknown_formula(self, unknown_mass, formula, kendrick_base):

        """
        Returns estimated formula for unknown_mass based on formula, and kendrick_base
        """
        
        # get mass and elements for provided formula
        formula_mass = self.mass_calculator(self.remove_punctuation(formula))
        formula_elements = self.count_elements(self.remove_punctuation(formula))
        
        # get mass and elements for provided kendrick base
        kendrick_base_mass = self.mass_calculator(self.remove_punctuation(kendrick_base))
        kendrick_base_elements = self.count_elements(kendrick_base)
               
        # integer multiple of kendrick bases
        kmd_units = self.kendrick_bases_apart(formula_mass, unknown_mass, kendrick_base_mass)
        
        # account for when kendrick base element is not in the provided formula
        if kendrick_base not in formula_elements:
            for element in self.count_elements(kendrick_base).keys():
                formula_elements[element] = 0
    
        # populate dictionary with keys of elements and values of element count
        unknown_formula_elements = collections.OrderedDict()
        for element in formula_elements:

            if element in kendrick_base_elements:
                unknown_formula_elements[element] = int(formula_elements[element] - (kmd_units * kendrick_base_elements[element]))
            else:
                unknown_formula_elements[element] = int(formula_elements[element])
        
        # order the dictionary
        unknown_formula_elements = collections.OrderedDict((k, v) for k, v in unknown_formula_elements.iteritems() if v)
        
        # if the kmd_units come back as a string then the known and unknown didnt match
        if isinstance(kmd_units, str):
            estimated = "not integer kmd_units away"
        elif sum(1 for number in unknown_formula_elements.values() if number < 0) > 0:
            estimated = "cant have negative subscript in formula"
        else:
            estimated = ''.join("%s%r" % (key,val) for (key,val) in unknown_formula_elements.iteritems())  
            
        return str(estimated)


    def get_kmd_suggestions(self, df, kendrick_base):

        same = ""
        count = 0
        exact_masses = df.loc[:,"x0"]
        match_column = df.loc[:,"KMD_"+kendrick_base+"_matches"]
        df[kendrick_base+"_suggested"] = None

        names = df.loc[:,"ion"]
        #for every entry
        for i, name in enumerate(names):

            matches = match_column[i]  
            #extract where the peak is unknown but matches a known
            if "unknown" in name and "-" in matches:
                # get mass of that unknown
                unknown_mass = float(exact_masses[i])
                # sort knowns into proper format (could be more than one)
                knowns = [entry for entry in matches.replace("'","").split(", ") if "unknown" not in entry]

                estimates = []
                # for every known
                for known in knowns:
                    # calculate what it thinks the formula of the unknown mass is and put in dataframe
                    estimated =  pl.calc_unknown_formula(unknown_mass, known, kendrick_base)
                    known_estimated = {"known": known, "estimated": estimated}
                    estimates.append(known_estimated)
                    # print count, known + " suggests <" + name + "> (" + str(unknown_mass) + ") is " + estimated 
                df.loc[i,kendrick_base+"_suggested"] = str(estimates)
                #print "\n"
                count += 1

        return df    
    
    
    
    

    
    
        
#initialise peak list object
pl = Peaklist()

