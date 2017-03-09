#!/usr/bin/python

class Peaklist(object):

    def __init__(self, reagent_ion, kendrick_bases):

        """
        self.reagent_ion: str, reagent ion e.g. I
        self.kendrick_bases: str list, base moeities
        self.kendrick_base_mass: float list, exact masses of kendrick bases
        self.ion: str list, names of peaks
        self.sumFormula: str list, formulas of peaks
        self.x0: float list, exact masses of formulas
        self.exact_mass_minus_reagent_ion: copy of x0 minus exact mass of reagent ion
        self.integer_mass: int list, rounded copy of x0
        self.integer_mass_minus_reagent_ion: copy of integer_mass minus exact mass of reagent ion
        self.mass_defect: float list, mass defect of peaks
        self._KEM: dict, keys=kendrick base, vals=float list, normalised mass to new kendrick base
        self._KMD: dict, keys=kendrick base, vals=float list, mass defect normalised to new kendrick base
        self._KMD_error: dict, keys=kendrick base, vals=float list, mass defect error to new kendrick base
        self.KMD_matches: dict, keys=kendrick base, vals=str list,
                          for ith self.ion, ith KMD_matches[kendrick_base] are peaks that meet matching c
                          criteria.
        self.weighted_best_guess_formulas: list of dicts,
                          for ith self.ion, ith weighted_best_guess_formulas is the suggested formulas
                          with a weighting based on how many times it was picked.
        self.weighted_best_guess_compounds: list of dicts,
                          for ith self.ion, ith weighted_best_guess_compounds is the suggested ChemSpider objects
        self.weighted_best_guess_names: list of dicts,
                          for ith self.ion, ith weighted_best_guess_compounds is the suggested common names
        self.error_on_assignment: float list, where ith new formulas are suggested, the ith error on assignment
                          is calculated and stored here.
        """

        self.reagent_ion = reagent_ion
        self.kendrick_bases = kendrick_bases
        self.kendrick_base_mass = [self._Mass_calculator(mass) for mass in self.kendrick_bases]
        self.ion = None
        self.sumFormula = None
        self.x0 = None
        self.exact_mass_minus_reagent_ion = None
        self.integer_mass = None
        self.integer_mass_minus_reagent_ion = None
        self.mass_defect = None
        self._KEM = {}
        self._KMD = {}
        self._KMD_error = {}
        self.KMD_matches = {}
        self.weighted_best_guess_formulas = None
        self.weighted_best_guess_compounds = None
        self.weighted_best_guess_names = None
        self.error_on_assignment = None


    def _Remove_punctuation(self, string):

        """Remove punctuation from string."""

        return string.replace("-","").replace("'","").replace("_","")


    def _Count_elements(self, moiety):

        """Returns dict with keys of elements and values of number of elements."""

#         ans = collections.OrderedDict({})
#         data = self._Get_molecule_params(moiety)
#         contents = data['parts'][0]['ea']
#         for block in contents:
#             ans[block['element']] = block['number']
#         return ans

        ret = []
        ans = re.findall(r'([A-Z][a-z]*)(\d*)', moiety)
        for i, value in enumerate(ans):
            val1 = value[1]
            if len(val1) == 0: val1 = 1
            ret.append((unicode(value[0]), int(val1)))

        return collections.OrderedDict(ret)


    def _Counted_elements_to_formula(self, ordered_dict):

        """Takes an ordered dict of counted elements from formula and returns a string."""

        string = ""
        for key in ordered_dict:
            string += key
            string += str(ordered_dict[key])

        return string

    def _Remove_reagent_ion(self, formula):

        """Takes a chemical formula string and removes the reagent ion."""

        split_formula = self._Count_elements(formula)
        for element in split_formula:
            if element == self.reagent_ion:
                del split_formula[element]

        return split_formula


    def _Add_reagent_ion(self, formula, mass):

        """Takes a str formula and adds reagent ion on difference between calc'd mass and given mass."""

        formula_mass = self._Mass_calculator(self._Remove_punctuation(formula))
        diff = abs(int(np.round(formula_mass)) - int(np.round(mass)))
        reagent_ion_mass = int(np.round(self._Mass_calculator(self.reagent_ion)))
        n_reagent_ions = diff / reagent_ion_mass

        if n_reagent_ions == 0:
            formula_with_reagent_ion = formula
        else:
            split_formula = self._Count_elements(formula)
            split_formula.update({self.reagent_ion : n_reagent_ions})

            updated_formula = collections.OrderedDict()
            for k, v in sorted(split_formula.items()):
                updated_formula[k] = v

            formula_with_reagent_ion = self._Counted_elements_to_formula(updated_formula)

        return formula_with_reagent_ion


#     def _Get_molecule_params(self, molecule):

#         """Hit chemcalc API to get molecule parameters."""

#         import urllib, urllib2, json

#         params = urllib.urlencode({'mf': molecule})
#         response = urllib2.urlopen('http://chemcalc.org/chemcalc/mf',params)
#         data = json.loads(response.read())

#         return data


    def _Mass_calculator(self, moiety):

        """Calculates exact mass for a given formula (string)."""

        return pyteomics.mass.calculate_mass(formula = moiety)


    def _Is_hydrocarbon(self, moiety):

        """Returns true if moiety is a hydrocarbon."""

        moiety = moiety.replace("C","")
        moiety = moiety.replace("H","")
        try:
            int(moiety)
            ishydrocarbon = True
        except ValueError:
            ishydrocarbon = False

        return ishydrocarbon


    def _Kendrick_bases_apart(self, amu1, amu2, kendrick_base):

        """How many kendrick_base amu1 and amu2 are different from each other."""

        ans = (int(np.round(amu1)) - int(np.round(amu2))) / np.round(kendrick_base)
        ans = ans if ans % 1 == 0.0 else None

        return ans


    def _F_amu_apart(self, amu1, amu2, kendrick_base_mass):

        """
        Returns True if difference between amu1 and
        amu2 when divided by moiety is equal to 0.
        """

        return abs(int(amu1) - int(amu2)) % int(kendrick_base_mass) == 0


    def _F_within_error(self, KMD1, KMD1e, KMD2, KMD2e):

        """Returns True if KMD1 +/- KMD1e is within KMD2 +/- KMD2e."""

        KMD1plus = KMD1 + KMD1e
        KMD2plus = KMD2 + KMD2e
        KMD1minus = KMD1 - KMD1e
        KMD2minus = KMD2 - KMD2e

        return (KMD2minus < KMD1minus < KMD2plus) or (KMD2minus < KMD1plus < KMD2plus)


    def Load_peak_list(self, directory, fname, **kwargs):

        """Load peaklist into init variables."""

        self.ion = []
        self.sumFormula = []
        self.x0 = []
        self.KMD_matches = {}

        df = pd.read_csv(directory+fname, header=0, **kwargs)
        self.ion = df["ion"].values.tolist()
        self.sumFormula = df['sumFormula'].values.tolist()
        self.x0 = df['x0'].values.tolist()

        del df


    def _Calc_mass_defect(self, reagent_ion_exact_mass):

        """ Populates init variables."""

        print "Calculating mass defects."

        reagent_ion_integer_mass = np.round(reagent_ion_exact_mass).astype(int)
        self.integer_mass = [int(np.round(n)) for n in self.x0]
        self.integer_mass_minus_reagent_ion = [x - reagent_ion_integer_mass if
                                               x > reagent_ion_integer_mass else
                                               x for x in self.integer_mass]
        self.exact_mass_minus_reagent_ion = [x - reagent_ion_integer_mass if
                                             x > reagent_ion_integer_mass else
                                             x for x in self.x0]
        self.mass_defect = [x - y for x, y in zip(self.exact_mass_minus_reagent_ion,
                                                  self.integer_mass_minus_reagent_ion)]


    def _Calc_kendrick_mass_defect(self, kendrick_base, kendrick_base_mass):

        """
        Calculate the kendrick mass defect for the kendrick base argument.
        Requries that _Calc_mass_defect has been run first.
        """

        assert len(self.mass_defect) > 0, '''No data in self.mass_defect.
                                             Did _Calc_mass_defect run correctly?'''

        print "calculating %s kendrick mass defects" % kendrick_base

        error = (20/1e6)
        self._KMD_error[kendrick_base] = [x * error for x in self.exact_mass_minus_reagent_ion]

        normalise = (np.round(kendrick_base_mass) / kendrick_base_mass)
        kem = [x * normalise for x in self.exact_mass_minus_reagent_ion]
        kmd = [x - y for x, y in zip(self.integer_mass_minus_reagent_ion, kem)]
        self._KEM[kendrick_base] = kem
        self._KMD[kendrick_base] = kmd


    def _Match_peaks_on_kmd(self, kendrick_base):

        """
        For each peak in self.ion, the conditions for peaks that match for
        the passed kendrick base are evaluated. Where the conditions are true,
        extract the names of those peaks and append (as a list) into self.KMD_matches

        Must run self._Calc_kendrick_mass_defect first to ensure
        the relevent columns are there:
            + ion
            + integerMass
            + KMD_<kendrick_base>
            + KMD_<kendrick_base>_error
        """

        assert len(self._KMD[kendrick_base]) > 0, '''No data in
        self._KMD[%s]. Did _Calc_kendrick_mass_defect run correctly?''' % kendrick_base

        print "matching %s kendrick mass defects" % kendrick_base

        self.KMD_matches[kendrick_base] = ["NaN" for x in self.ion]
        kendrick_base_integer_mass = int(np.round(self._Mass_calculator(kendrick_base)))

        # make bool arrays where criteria for matching is met
        amuApart1 = []
        for y in self.integer_mass:
            for x in self.integer_mass:
                amuApart1.append(self._F_amu_apart(x, y, kendrick_base_integer_mass))

        withinError1 = []
        for y, z in zip(self._KMD[kendrick_base], self._KMD_error[kendrick_base]):
            for x, w in zip(self._KMD[kendrick_base], self._KMD_error[kendrick_base]):
                withinError1.append(self._F_within_error(x, w, y, z))


        isMatch = [a and b for a, b in zip(amuApart1, withinError1)]  # both criteria are met

        for i, name in enumerate(self.ion):

            '''loop step makes i the starting point of bool array
            for each row then mask peaknames where isMatch is true'''
            mask = isMatch[i * len(self.ion) : (i * len(self.ion)) + len(self.ion)]
            matches = [x for j, x in enumerate(self.ion) if mask[j]]

            if len(matches) > 0:
                toAppend = str(list(matches)).replace('[',"").replace(']',"").replace("'","")
                self.KMD_matches[kendrick_base][i] = toAppend



    def _Calc_unknown_formula(self, known_formula, unknown_mass, kendrick_base):

        """
        Returns estimated formula for unknown_mass based on known_formula and kendrick_base.
        known_formula : string
        unknown_mass  : float
        kendrick_base : string
        """

        # get mass and elements for known arguments and collate them.
        formula_mass = self._Mass_calculator(self._Remove_punctuation(known_formula))
        kendrick_base_mass = self._Mass_calculator(self._Remove_punctuation(kendrick_base))

        formula_elements = self._Count_elements(self._Remove_punctuation(known_formula))
        kendrick_base_elements = self._Count_elements(kendrick_base)
        all_elements = formula_elements.copy()
        all_elements.update(kendrick_base_elements)

        kmd_units = self._Kendrick_bases_apart(formula_mass, unknown_mass, kendrick_base_mass)
        if not kmd_units:
            estimated = "Not integer kmd_units away."
        else:

            # new dictionary that multiplies the kendrick_bases
            # elements by how many repeating kmd bases there are
            kendrick_base_update = collections.Counter()
            for el in kendrick_base_elements:
                kendrick_base_update[el] = int(kendrick_base_elements[el] * -kmd_units)

            # update unknown_formula_elements to contain suggested formula
            unknown_formula_elements = collections.Counter()
            unknown_formula_elements.update(kendrick_base_update)

            unknown_formula_elements.update(formula_elements)

            for k,v in unknown_formula_elements.items():
                if v == 0:
                    del unknown_formula_elements[k]    # get rid of any 0 values

            # if the suggested formula have negative subscript it cant be real
            if sum(1 for number in unknown_formula_elements.values() if number < 0) > 0:
                estimated = "Can't have negative subscript in formula."
            else:
                estimated = ''.join("%s%r" % (key,val) for (key,val) in unknown_formula_elements.iteritems())

        return str(estimated)



    def _Likely_formulas(self, suggested_formulas, unknown_mass):

        """
        Takes a list of estimated formulas and returns the most likely candidates.
        The results are collated and returned as a dictionary with unique formula keys
        and frequency that key is suggested as values.

        This is the function that decides if what the solver has returned is rubbish or not.
        """

        # remove reagent ion
        suggested_formulas = [self._Counted_elements_to_formula(
                                  self._Remove_reagent_ion(suggested_formula)) for
                              suggested_formula in suggested_formulas]

        # if chemspider has no result, get rid
        potential_formulas = [x for x in suggested_formulas if len(cs.simple_search_by_formula(x)) > 0]
        #if only a hydrocarbon, get rid
        potential_formulas = [y for y in potential_formulas if not self._Is_hydrocarbon(y)]

        # count how many times a formula is guessed and make a dictionary, then apend to list
        assoc_formula = []
        for potential_formula in potential_formulas:
            assoc_formula.append({potential_formula: potential_formulas.count(potential_formula)})

        # Collect where formulas are repeated
        likely_formulas = collections.defaultdict(int)
        for formula in assoc_formula:
            key = self._Add_reagent_ion(formula.keys()[0], unknown_mass)
            likely_formulas[key] += formula.values()[0]


        return dict(likely_formulas)


    def _Weighted_best_guesses(self, best_guesses):

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


    def _Known_in_any_kmd_match(self, dict_of_array_vals):

        """
        Returns index mask array where the - character (indicating a known species)
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


    def _Get_KMD_match_suggestions(self):

        """
        Works out what the matching peak could be based on the kendrick_base
        and returns a list of the unknown name, unknown mass and chemspider objects
        """

        print "Getting suggested structures. This may take some time."

        self.weighted_best_guess_formulas = ["NaN" for x in self.ion]
        self.weighted_best_guess_compounds = ["NaN" for x in self.ion]
        self.weighted_best_guess_names = ["NaN" for x in self.ion]

        # no point iterating over every row. we are only interested where the
        # ion is unknown but at least one of the matches IS known.
        unknown_ions = [i for i, x in enumerate(self.ion) if "unknown" in x]
        known_matches = self._Known_in_any_kmd_match(self.KMD_matches)
        common = sorted(list(set(unknown_ions) & set(known_matches)))

        unknown_names = [self.ion[x] for x in common]
        unknown_masses = [self.x0[x] for x in common]

        results = []
        # for each unknown peak
        for i, elem in enumerate(common):

            print "matching %d of %d unknowns" % (i, len(common))

            unknown_mass = unknown_masses[i]
            unknown_name = unknown_names[i]
            kmd_matches = [self.KMD_matches[key][elem].split(', ') for key in self.kendrick_bases]
            kmd_matches = map(lambda x:
                              [match for match in kmd_matches[x] if "-" in match],
                              range(len(kmd_matches))
                             )

            # for each known match in that KMD_XXX_MATCH column
            suggested_formulas = []

            # iterate over the matches and their bases
            for pair in zip(self.kendrick_bases, kmd_matches):

                kendrick_base = pair[0]
                kmd_match = pair[1]
                suggested = []
                #calculate the suggested formulas for each kmd match
                map(lambda x:
                        suggested.append(self._Calc_unknown_formula(str(x),
                            float(unknown_mass), kendrick_base)),
                    kmd_match)
                suggested_formulas.append(suggested)

            # get the suggested formulas based on logic applied in self._Likely_formulas
            likely_formulas = map( lambda x:
                                   self._Likely_formulas(suggested_formulas[x], unknown_mass),
                               range(len(suggested_formulas)))

            # Collapse best guesses list into dict that counts unique guesses and their frequency
            weighted_best_guess = self._Weighted_best_guesses(likely_formulas)
            if len(weighted_best_guess) > 0:
                self.weighted_best_guess_formulas[elem] = weighted_best_guess

#             print "unknown_mass: ", unknown_mass
#             print "unknown_name: ", unknown_name
#             print "kmd_matches: ", kmd_matches
#             print "suggested_formulas:", suggested_formulas
#             print "likely_formulas:", likely_formulas
#             print "weighted_best_guess:", weighted_best_guess
#             print "\n"

            # Get common names of weighted best_guesses
            compounds = []
            common_names = []
            for formula in weighted_best_guess.keys():

                formula = self._Counted_elements_to_formula(self._Remove_reagent_ion(formula))

                CScompounds = cs.simple_search_by_formula(formula)
                compounds.append(self._Get_visible_compounds(CScompounds))
                CS_names = self._Get_visible_compounds(CScompounds, names=True)
                common_names.append(CS_names)

            self.weighted_best_guess_compounds[elem] = compounds
            self.weighted_best_guess_names[elem] = common_names

            results.append((unknown_name, unknown_mass, weighted_best_guess, compounds))

        return results


    def _Get_visible_compounds(self, list_of_csCompounds, names=False):

        """Takes a list of cs.Compounds and returns their
        names if they are visible by CIMS."""

        ls = []
        for compound in list_of_csCompounds:
            try:
                condition_met = self._Condition_for_visible(compound)
                if condition_met:
                    if names:
                        ls.append(compound.common_name)
                    else:
                        ls.append(compound)
                    if len(ls) > 4:
                        break
                else:
                    ls.append("Structure not visible by %s CIMS" % (self.reagent_ion))
            except KeyError:
                ls.append("No Common Name")

        return ls


    def _Condition_for_visible(self, compound):

        """
        Condition for whether the suggested molecule is visible by CIMS.
        + Can't be dimer.
        + Can't have overall + or - charge.
        + Cant have partial charge.
        + No wierd character in name that is odd encoding of a greek letter
            (indicates non typical oxidation state)
        """

        neg_count = 0
        pos_count = 0
        for char in compound.smiles:
            if char == "-": neg_count += 1
            if char == "+": pos_count += 1

        return ("." not in compound.smiles) and (neg_count == pos_count) and ("$" not in compound.common_name) and ("{" not in compound.common_name) and ("^" not in compound.common_name)


    def _Error_on_assignment(self):

        """Provides ppm error on assignment of unknown peak with sugggested formula."""


        self.error_on_assignment = ["NaN" for x in self.ion]
        for i, item in enumerate(self.weighted_best_guess_formulas):
            try:
                exact_x0 = self._Mass_calculator(item.keys()[0])
                self.error_on_assignment = 1e6 * ((self.x0[i] - exact_x0) / exact_x0)
            except ValueError:
                pass
            except AttributeError:
                pass



    def New_peaklist(self):

        """Return new peaklist with updated assignments."""

        new_peaklist = pd.DataFrame()

        new_peaklist['ion'] = bonfire.ion
        new_peaklist['sumFormula'] = bonfire.sumFormula
        new_peaklist['x0'] = bonfire.x0

        for i, item in enumerate(bonfire.weighted_best_guess_formulas):
            try:
                new_peaklist.loc[i, ['ion']] = item.keys()[0]
                #exact_x0 = bonfire._Mass_calculator(item.keys()[0])
                #new_peaklist.loc[i, ['x0']] = exact_x0

            except ValueError:
                pass
            except AttributeError:
                pass

        return new_peaklist



    def Run(self):


        # Calculate mass defect
        self._Calc_mass_defect(self._Mass_calculator(self.reagent_ion))

        # # Calculate kendrick mass defects
        [self._Calc_kendrick_mass_defect(base, self.kendrick_base_mass[i]) for i, base in enumerate(self.kendrick_bases)]

        # Stage 1. Match kendrick mass defects - still unknown
        [self._Match_peaks_on_kmd(base) for base in self.kendrick_bases]

        # Stage2. Get suggested formulas and compounds for matches
        self._Get_KMD_match_suggestions()
        self._Error_on_assignment()






if __name__ == '__main__':


    start_time = time.time()

    i = 0

    while i < 6:

        #initialise peak list object
        plo = Peaklist("I", ["CH2","O","CHO","OH","CO","CO2","NO2","H2","NO3"])

        # Load data into dataframe
        if i > 0:
            fname = "Bonfire14PeakList_kmd%i.txt" % i
        else:
            fname = "Bonfire14PeakList.txt"

        print "loading %s" % fname

        plo.Load_peak_list("/home/mbexkmp3/Documents/Bonfire14/", fname, sep="\t")

        # Calculates mass defect and kendrick mass defects for
        # passed kendrick bases. Matches kendrick mass defects
        # and gets suggested formulas and error on assignment.
        bonfire.Run()

        i += 1

        # Generate new peaklist based on suggestions
        plo.New_peaklist().to_csv(
            "/home/mbexkmp3/Documents/Bonfire14/Bonfire14PeakList_kmd%i.txt" % i ,
            sep="\t", index=False)


    print "--- Time taken: %s s ---" % (round(time.time() - start_time, 3))
