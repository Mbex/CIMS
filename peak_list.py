
class Peaklist(object):

    def __init__(self, reagent_ion, kendrick_bases):
        
        self.reagent_ion = reagent_ion
        self.kendrick_bases = kendrick_bases
        self.kendrick_base_mass = [self.mass_calculator(mass) for mass in self.kendrick_bases]
        self.ion = []
        self.sumFormula = []
        self.x0 = []     
        self._KEM = {}
        self._KMD = {}
        self._KmdError = {}
        self.KMD_matches = {}
    
        
        # Calculate mass defect
        self.calc_mass_defect(self.mass_calculator(self.reagent_ion))

        
    def calc_mass_defect(self, reagent_ion_exact_mass):

        """
        Makes integer and exact mass columns with
        and without reagent ion then
        calculates the mass defect from them.
        """

        print "calculating mass defects"

        regaentIonIntegerMass = np.round(reagent_ion_exact_mass).astype(int)   
        self._integerMass = [int(np.round(n)) for n in self.x0]
        self._integerMassNoRI = [x - regaentIonIntegerMass if x > regaentIonIntegerMass else x for x in self._integerMass]   
        self._exactMassNoRI = [x - regaentIonIntegerMass if x > regaentIonIntegerMass else x for x in self.x0]   
        self._massDefect = [x - y for x, y in zip(self._exactMassNoRI, self._integerMassNoRI)]
        
        
    def remove_punctuation(self, string):
        
        """Remove minus from string."""
        
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

        return pyteomics.mass.calculate_mass(formula=moiety)
 

    def count_elements(self, moiety):

        """Returns dict with keys of elements and values of number of elements."""
   
        #re.findall(r'([A-Z][a-z]*)(\d*)', '"H15C10O')
        
        ans = collections.OrderedDict({})
    
        data = self._Get_molecule_params(moiety)

        contents = data['parts'][0]['ea']

        for block in contents:
            ans[block['element']] = block['number']
            
        return ans 
    
    
    def is_hydrocarbon(self, moiety):
    
        """Returns true if moiety is a hydrocarbon."""
    
        moiety = moiety.replace("C","")
        moiety = moiety.replace("H","")
        try:
            int(moiety)
            ishydrocarbon = True
        except ValueError:
            ishydrocarbon = False

        return ishydrocarbon

    
    def kendrick_bases_apart(self, amu1, amu2, kendrick_base):
    
        """How many kendrick_base amu1 and amu2 are different from each other."""

        ans = (int(np.round(amu1)) - int(np.round(amu2))) / np.round(kendrick_base)
        
        if ans % 1 != 0.0: 
            ans = "not integer bases apart"
       
        return ans 
    

    def fAmuApart(self, amu1, amu2, kendrick_base_mass):

        """
        Returns True if difference between amu1 and
        amu2 when divided by moiety is equal to 0.
        """

        return abs(int(amu1) - int(amu2)) % int(kendrick_base_mass) == 0


    def fWithinError(self, KMD1, KMD1e, KMD2, KMD2e):

        """Returns True if KMD1 +/- KMD1e is within KMD2 +/- KMD2e."""

        KMD1plus = KMD1 + KMD1e
        KMD2plus = KMD2 + KMD2e
        KMD1minus = KMD1 - KMD1e
        KMD2minus = KMD2 - KMD2e

        return (KMD2minus < KMD1minus < KMD2plus) or (KMD2minus < KMD1plus < KMD2plus)


    def load_peak_list(self, directory, fname, **kwargs):

        """Returns dataframe with peaklist information."""

        #import pandas as pd

        df = pd.read_csv(directory+fname, header=0, **kwargs)
        self.ion = df["ion"].values.tolist()
        self.sumFormula = df['sumFormula'].values.tolist()
        self.x0 = df['x0'].values.tolist()
        
        self.main()
        
                                   
        return None

    
    def rename_unknowns(self, col, pattern="unknown"):

        """used for renaming unknowns. reset numbers at end of unknown pattern"""

        names = peakList[col]
        count = 1
        for i, name in enumerate(names):
            if pattern in name:
                newname = pattern + "%04d" % (count,)
                peakList.loc[i, col] = newname
                count = count + 1

        return None

   

    def calc_kendrick_mass_defect(self, kendrick_base, kendrick_base_mass):

        """
        Calculate the kendrick mass defect for the passed kendrick base.
        Requries that calc_mass_defect has been run first.
        """

        print "calculating %s kendrick mass defects" % kendrick_base
             
        error = (20/1e6)
        normalise = (np.round(kendrick_base_mass) / kendrick_base_mass)
        kem = [x * normalise for x in self._exactMassNoRI]
        kmd = [x - y for x, y in zip(self._integerMassNoRI, kem)]
        kmdError = [x * error for x in self._exactMassNoRI]
        
        self._KEM[kendrick_base] = kem
        self._KMD[kendrick_base] = kmd
        self._KmdError[kendrick_base] = kmdError

        return None


    def match_peaks_on_kmd(self, kendrick_base):
        
        """
        Generates new column of species that have the same 
        kendrick mass defect as that in the ion column.
        
        Must run self.calc_kendrick_mass_defect first to ensure
        the relevent columns are there:
            + ion
            + integerMass
            + KMD_<kendrick_base>
            + KMD_<kendrick_base>_error
        """

        print "matching %s kendrick mass defects" % kendrick_base

        peak_names = self.ion
        kbm_integer = int(np.round(self.mass_calculator(kendrick_base)))
        length = len(self.ion)
        amuApart1 = []
        amuApart1 = []
        withinError1 = []
        withinError1 = []
        self.KMD_matches[kendrick_base] = ["NaN" for x in peak_names]

        # amu apart
        for y in self._integerMass:
            for x in self._integerMass:
                amuApart1.append(self.fAmuApart(x, y, kbm_integer))
                
        # within error
        for y, z in zip(self._KMD[kendrick_base], self._KmdError[kendrick_base]):
            for x, w in zip(self._KMD[kendrick_base], self._KmdError[kendrick_base]):
                withinError1.append(self.fWithinError(x, w, y, z))

        # both criteria are met
        isMatch = [a and b for a, b in zip(amuApart1, withinError1)]   
   
        # loop step makes i the starting point of bool array for each row
        # then mask peaknames where isMatch is true, append into dataframe
        for i, name in enumerate(peak_names):
            
            mask = isMatch[i * length : (i * length) + length]        
            matches = [x for j, x in enumerate(peak_names) if mask[j]]
                        
            if len(matches) > 0:
                toAppend = str(list(matches)).replace('[',"").replace(']',"").replace("'","")
                self.KMD_matches[kendrick_base][i] = toAppend
            
        return None

            
    def calc_unknown_formula(self, known_formula, unknown_mass, kendrick_base):

        """
        Returns estimated formula for unknown_mass based on formula, and kendrick_base
        known_formula : string
        unknown_mass  : float
        kendrick_base : string
        """
        
        # get mass and elements for provided formula
        formula_mass = self.mass_calculator(self.remove_punctuation(known_formula))
        formula_elements = self.count_elements(self.remove_punctuation(known_formula))
     
        # get mass and elements for provided kendrick base
        kendrick_base_mass = self.mass_calculator(self.remove_punctuation(kendrick_base))
        kendrick_base_elements = self.count_elements(kendrick_base)
               
        # integer multiple of kendrick bases
        kmd_units = self.kendrick_bases_apart(formula_mass, unknown_mass, kendrick_base_mass)
         
        # list all elements in kendrick base and known formula
        all_elements = formula_elements.copy()
        all_elements.update(kendrick_base_elements)

        # if the kmd_units come back as a string then the known and unknown didnt match
        if isinstance(kmd_units, str):
            estimated = "not integer kmd_units away"
        else:

            # new dictionary that multiplies the kendrick_bases elements by how many repeating kmd bases there are
            kendrick_base_update = collections.Counter()
            for element in kendrick_base_elements:
                kendrick_base_update[element] = int(kendrick_base_elements[element] * -kmd_units)

            # update unknown_formula_elements to contain suggested formula
            unknown_formula_elements = collections.Counter()
            unknown_formula_elements.update(kendrick_base_update)
            unknown_formula_elements.update(formula_elements)

            # get rid of any 0 values
            for k,v in unknown_formula_elements.items():
                if v == 0:
                    del unknown_formula_elements[k]

            # if the suggested formula have negative subscript it cant be real        
            if sum(1 for number in unknown_formula_elements.values() if number < 0) > 0:
                estimated = "cant have negative subscript in formula"
            else:
                estimated = ''.join("%s%r" % (key,val) for (key,val) in unknown_formula_elements.iteritems())  

        return str(estimated)

        
    def most_likely_formula(self, suggested_formulas, unknown_mass):
        
        """
        Takes a list of estimated formulas and returns the most likely candidates.
        The results are collated and returned as a dictionary with unique formula keys
        and count as the values.
        
        This is the function that decides if what the solver has returned is rubbish or not.
        """
        
        # if chemspider has no result, get rid 
        suggested_formulas = [suggested_formula.replace('I1','').replace('I2','').replace('I3','') for suggested_formula in suggested_formulas]
        potential_formulas = [x for x in suggested_formulas if len(cs.simple_search_by_formula(x)) > 0]
        #if only a hydrocarbon, get rid
        potential_formulas = [x for x in potential_formulas if not self.is_hydrocarbon(x)]
        
        # count how many times a formula is guessed and make a dictionary, then apend to list
        instances = 1
        assoc_formula = []        
        for potential_formula in potential_formulas:
            assoc_formula.append({potential_formula: potential_formulas.count(potential_formula)})

        # Collect where formulas are repeated
        guesses = collections.defaultdict(int)
        for formula in assoc_formula:
            guesses[formula.keys()[0]] += formula.values()[0]
            
        return dict(guesses) #most_likely                     
               
        
    def collate_best_guesses(self, best_guesses):

        """
        Takes a list of dicts with str keys of suggested formulas and
        int values of the times they are chosen as a match.
        """

        unique_keys = (set().union(*(d.keys() for d in best_guesses)))
        answer = {el:0 for el in unique_keys}

        if len(unique_keys) < 1:
            pass
        else:

            for key in unique_keys:
                val = 0
                for item in best_guesses:
                    try:
                        val += item[key]
                    except KeyError:
                        pass
                    
                answer[key] = val
            
        return answer   
    
    
    def known_in_any_kmd_match(self, dict_of_array_vals):
    
        """
        Returns index mask array where the - (indicating a known species)character
        is found in the nth element of the array values for all keys in a dictionary.
        Specific useage for kmd_matches array dict.
        e.g. {"CH2":["NaN","C2H5O-","NaN","NaN"], "OH":["CH4O2-","NaN","NaN","NaN"]}

        """

        mask = []
        for i in xrange(len(dict_of_array_vals.values()[0])):
            for k in (dict_of_array_vals.keys()):
                item = dict_of_array_vals[k][i]
                if "-" in item:
                    mask.append(i)
                    break
        return mask

        
    def get_kmd_suggestions(self):

        """
        Works out what the matching peak could be based on the kendrick_base
        and returns a list of the unknown name, unknown mass and chemspider objects
        """
        
        print "Getting suggested structures. This may take some time"

             
        # no point iterating over every row. we are only interested where the
        # ion is unknown but at least one of the matches IS known.
        self._unknown_ions = [i for i, x in enumerate(self.ion) if "unknown" in x]
        self._known_matches = self.known_in_any_kmd_match(self.KMD_matches)
        self.common_rows = list(set(self._unknown_ions) & set(self._known_matches))     
        self._unknown_masses = [self.x0[x] for x in self.common_rows]
        self._unknown_names = [self.ion[x] for x in self.common_rows]
        
        results = []
     
        # for each unknown peak
        for i, row in enumerate(self.common_rows):
            
            unknown_mass = self._unknown_masses[i] 
            unknown_name = self._unknown_names[i] 
            kmd_matches = [self.KMD_matches[key][row].split(', ') for key in self.kendrick_bases]           
            kmd_matches = map(lambda x: [match for match in kmd_matches[x] if "-" in match], range(len(kmd_matches)))
            
            # for each known match in that KMD_XXX_MATCH column 
            suggested_formulas = []
            
            # iterate over the matches and their bases
            for pair in zip(self.kendrick_bases, kmd_matches):
                
                kendrick_base = pair[0]
                kmd_match = pair[1]
                #calculate the suggested formulas for each kmd match 
                suggested = []
                map(lambda x:suggested.append(self.calc_unknown_formula(str(x), float(unknown_mass), kendrick_base)), kmd_match)
                suggested_formulas.append(suggested)
                
            # get the suggested formulas based on logic applied in self.most_likely_formula
            best_guesses = map( lambda x: self.most_likely_formula(suggested_formulas[x], unknown_mass), range(len(suggested_formulas)))

            # Collapse best guesses list into dict that counts unique guesses and their frequency
            weighted_best_guess = self.collate_best_guesses(best_guesses)

            print "%d of %d" % (i, len(self.common_rows))
            
            # Get common names of weighted best_guesses
            compounds = []
            common_names = []           
            for formula in weighted_best_guess.keys():
                f = cs.simple_search_by_formula(formula)
                compounds.append(f)

            ans = (unknown_name, unknown_mass, weighted_best_guess, compounds)
            results.append(ans)
           
        return results
 

    def main(self):

	"""Called after peaklist is loaded."""
        
        # Calculate mass defect
        self.calc_mass_defect(self.mass_calculator(self.reagent_ion))

        # # Calculate kendrick mass defects
        [self.calc_kendrick_mass_defect(base, self.kendrick_base_mass[i]) for i, base in enumerate(self.kendrick_bases)]
  
        # Match kendrick mass defects
        [self.match_peaks_on_kmd(base) for base in self.kendrick_bases]
        
        # # Get suggested structures
        self.unknown_tables = pd.DataFrame().from_records(self.get_kmd_suggestions(),
                                                                   columns = ["unknown_name",
                                                                       "unknown_mass",
                                                                       "weighted_formulas",
                                                                       "chemspider_compounds"
                                                                    ])
        
 
    def condition_for_visible(self, compound):

        """
        Condition for whether the suggested molecule is visible by CIMS.
        + Can't be dimer.
        + Can't have + or - charge.
        + Cant have partial charge.
        + No wierd character in name that is odd encoding of a greek letter
            (indicates non typical oxidation state)
        """

        return ("." not in compound.smiles) and ("+" not in compound.smiles) and ("-" not in compound.smiles) and ("$" not in compound.common_name) and ("{" not in compound.common_name) and ("^" not in compound.common_name)

 
    def get_compound_names(self, list_of_csCompounds):

        """Takes a list of cs.Compounds and returns their
        names if they are visible by CIMS."""

        list_of_names = []    
        for compound in list_of_csCompounds:
            try:
                condition_met = self.condition_for_visible(compound)
                if condition_met:
                    list_of_names.append(compound.common_name)
                    if len(list_of_names) > 4:
                        break
            except KeyError:
                list_of_names.append("No Common Name")

        return list_of_names

