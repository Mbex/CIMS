#!/usr/bin/python

from chemspipy import ChemSpider
from pyteomics import mass
from itertools import cycle
import re
import collections
import os
import time
import array
import pyteomics
import sys
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

class Peaklist(object):

    def __init__(self, reagent_ion, kendrick_bases):

        """
        self.cs: str, ChemSpider instance
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
        self.mass_defect_minus_reagent_ion : float list, mass defect of peaks (with reagent ion removed)
        self._KEM: dict, keys=kendrick base, vals=float list, normalised mass to new kendrick base
        self._KMD: dict, keys=kendrick base, vals=float list, mass defect normalised to new kendrick base
        self._KMD_error: dict, keys=kendrick base, vals=float list, mass defect error to new kendrick base
        self.KMD_matches: dict, keys=kendrick base, vals=str list,
                          for ith self.ion, ith KMD_matches[kendrick_base] are peaks that meet matching c
                          criteria.
        self.suggested_formulas: list of dicts,
                          for ith self.ion, ith suggested_formulas is the suggested formulas
                          with a weighting based on how many times it was picked.
        self.suggested_names: list of dicts,
                          for ith self.ion, ith suggested_compounds is the suggested common names
        self.suggested_errors: list of dicts,
                          for ith self.ion, ith suggested_compounds is the suggested error of the assignment
        """

        self.cs = ChemSpider(os.environ["CHEM_SECURITY_TOKEN"])
        self.directory = ""
        self.fname = ""
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
        self.mass_defect_minus_reagent_ion = None
        self._KEM = {}
        self._KMD = {}
        self._KMD_error = {}
        self.KMD_matches = {}
        self.suggested_formulas = {}
        self.suggested_names = {}
        self.suggested_errors = {}


    def _Flatten(self, list_of_lists):

        """flattens a list of lists into a list"""

        return [item for sublist in list_of_lists for item in sublist]


    def _Remove_punctuation(self, string):

        """Remove punctuation from string."""

        return string.replace("-","").replace("'","").replace("_","")


    def _Count_elements(self, moiety):

        """Returns dict with keys of elements and values of number of elements."""

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
        amu2 when divided by moiety mass is equal to 0.
        """

        return abs(int(amu1) - int(amu2)) % int(kendrick_base_mass) == 0


    def _F_within_error(self, KMD1, KMD1e, KMD2, KMD2e):

        """Returns True if KMD1 +/- KMD1e is within KMD2 +/- KMD2e."""

        KMD1_ulim = KMD1 + KMD1e
        KMD2_ulim = KMD2 + KMD2e
        KMD1_llim = KMD1 - KMD1e
        KMD2_llim = KMD2 - KMD2e

        type1 = ((KMD1_llim >= KMD2_llim) and (KMD1_llim <= KMD2_ulim))
        type4 = ((KMD1_llim <= KMD2_llim) and (KMD1_ulim >= KMD2_ulim))
        type2 = ((KMD1_ulim >= KMD2_llim) and (KMD1_ulim <= KMD2_ulim))
        type3 = ((KMD1_ulim <= KMD2_ulim) and (KMD1_llim >= KMD2_llim))

        return type1 or type2 or type3 or type4

    def Load_peak_list(self, directory, fname, **kwargs):

        """Load peaklist into init variables."""

        self.ion = []
        self.sumFormula = []
        self.x0 = []
        self.KMD_matches = {}
        self.directory = directory
        self.fname = fname

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
                                               x >= reagent_ion_integer_mass else
                                               x for x in self.integer_mass]
        self.exact_mass_minus_reagent_ion = [round(x - reagent_ion_exact_mass, 6) if
                                             x >= reagent_ion_exact_mass else
                                             x for x in self.x0]

        self.mass_defect_minus_reagent_ion = [round(y - x, 6) for x, y in zip(self.integer_mass_minus_reagent_ion, self.exact_mass_minus_reagent_ion)]
        self.mass_defect = [round(y - x, 6) for x, y in zip(self.integer_mass, self.x0)]


    def _Calc_kendrick_mass_defect(self, kendrick_base, kendrick_base_mass):

        """
        Calculate the kendrick mass defect for the kendrick base argument.
        Requries that _Calc_mass_defect has been run first.
        """

        assert len(self.mass_defect_minus_reagent_ion) > 0, '''No data in self.mass_defect. Did _Calc_mass_defect run correctly?'''

        print "calculating %s kendrick mass defects" % kendrick_base

        error = (20/1e6)
        self._KMD_error[kendrick_base] = [round(x * error, 6) for x in self.exact_mass_minus_reagent_ion]

        normalise = (np.round(kendrick_base_mass) / kendrick_base_mass)
        kem = [round(x * normalise, 6) for x in self.exact_mass_minus_reagent_ion]
        kmd = [round(y - x, 6) for x, y in zip(self.integer_mass_minus_reagent_ion, kem)]
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
        amuApart1 = [self._F_amu_apart(x, y, kendrick_base_integer_mass) for x, y in zip(self.integer_mass, self.integer_mass)]
        withinError1 = [self._F_within_error(y[0], y[1], x[0], x[1]) for x, y in zip((zip(self._KMD[kendrick_base], self._KMD_error[kendrick_base])), (zip(self._KMD[kendrick_base], self._KMD_error[kendrick_base])))]
        isMatch = [a and b for a, b in zip(amuApart1, withinError1)]  # both criteria are met

        for i, name in enumerate(self.ion):

            '''loop step makes i the starting point of bool array
            for each row then mask peaknames where isMatch is true'''

            mask = isMatch[i * len(self.ion) : (i * len(self.ion)) + len(self.ion)]
            matches = [x for j, x in enumerate(self.ion) if mask[j]]
            self.KMD_matches[kendrick_base][i] = [x for x in matches if (name != x) and (len(matches) > 0)]



    def _Calc_unknown_formula(self, known_formula, unknown_mass, kendrick_base):

        """
        Returns estimated formula for unknown_mass based on known_formula and kendrick_base.
        known_formula : string
        unknown_mass  : int
        kendrick_base : string
        """

        unknown_mass = int(round(unknown_mass))

        if not isinstance(unknown_mass, int):
            raise Exception("unknown_mass should be an integer")

        # get mass and elements for known arguments and collate them.
        formula_mass = self._Mass_calculator(self._Remove_punctuation(known_formula))
        kendrick_base_mass = self._Mass_calculator(self._Remove_punctuation(kendrick_base))

        formula_elements = self._Count_elements(self._Remove_punctuation(known_formula))
        kendrick_base_elements = self._Count_elements(kendrick_base)
        all_elements = formula_elements.copy()
        all_elements.update(kendrick_base_elements)

        kmd_units = self._Kendrick_bases_apart(formula_mass, unknown_mass, kendrick_base_mass)
        if not kmd_units:
            estimated = "error - Not integer kmd_units away."
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
                estimated = "error - Can't have negative subscript in formula."
            else:
                estimated = ''.join("%s%r" % (key,val) for (key,val) in unknown_formula_elements.iteritems())

        return str(estimated)


    def _Is_formula_realistic(self, suggested_formula):

        """
        Takes a suggested formula and returns true if it passes the conditions posed here.
        Takes into account realistic structure and visiblity by CIMS.
        returns Bool.
        This is the function that decides if what the solver has returned is rubbish or not.
        """

        suggested_formula = self._Remove_reagent_ion(suggested_formula)
        suggested_formula = self._Counted_elements_to_formula(suggested_formula)
        list_of_compounds = self.cs.simple_search_by_formula(suggested_formula)

        return_val = 1
        if len(list_of_compounds) < 1:
            return_val = 0
        if self._Is_hydrocarbon(suggested_formula):
            return_val = 0

        return return_val


    def _Matched(self, kendrick_base, known=True):

        """
        Returns list of KMD matches for the passed kendrick base
        containin; known matches i.e. if a '-' is present in the name
        indicating it has been assigned a formula; or unknown matches
        if the '-' character is not present
        """

        matched = []
        for matches in self.KMD_matches[kendrick_base]:
            if known:
                matched.append([x for x in matches if "-" in x])
            else:
                matched.append([x for x in matches if not "-" in x])

        matched = list(set(self._Flatten(matched)))

        return matched


    def _Get_KMD_match_suggestions(self):

        """
        Works out what the matching peak could be based on the kendrick_base
        and returns a list of the unknown name, unknown mass and chemspider objects
        """

        print "Identifying compounds. This may take some time."

        # for each kendrick base
        for kendrick_base in self.kendrick_bases:

            self.suggested_formulas[kendrick_base] = [[] for x in self.ion]
            self.suggested_errors[kendrick_base] = [[] for x in self.ion]
            self.suggested_names[kendrick_base] = [[] for x in self.ion]

            print "Looking at matches for %s" % kendrick_base

            # find those KMD matches that are unknown
            matched_unknowns = self._Matched(kendrick_base, known=False)

            # for each of those uknowns
            for i, unknown in enumerate(matched_unknowns):

                print "matching %d of %d unknowns" % (i+1, len(matched_unknowns))

                # find the KMD matches that are known
                matched_knowns = self._Matched(kendrick_base, known=True)

                # for each of those knowns
                guesses_formulas = []
                guesses_names = []
                guesses_errors = []

                for known in matched_knowns:

                    unknown_index = self.ion.index(unknown)
                    unknown_mass = self.x0[unknown_index]

                    # get the suggested structure of the unknown that the known points to
                    suggested_formula = self._Calc_unknown_formula(known, unknown_mass, kendrick_base)

                    # igor if error is returned
                    if not suggested_formula[0:5] == "error":

                        # return true if structure is realistic
                        realistic = self._Is_formula_realistic(suggested_formula)

                        if not realistic:
                            pass

                        # Get formula without reagent ion
                        suggested_formula_noRI = self._Counted_elements_to_formula(self._Remove_reagent_ion(suggested_formula))

                        CScompounds = self.cs.simple_search_by_formula(suggested_formula_noRI)
                        suggested_names = self._Get_visible_compounds(CScompounds, names=True)

                        if suggested_names == []:
                            pass
                        else:

                            error = self._Error_on_assignment(suggested_formula, unknown_mass)

                            # put known and suggested pairs into list
                            guesses_formulas.append((known, suggested_formula))
                            guesses_errors.append((suggested_formula, error))
                            guesses_names.append((known, suggested_names))

                    else:
                        pass

                # put lists of tuples of the known and suggested unknown
                # into an object variable at the same index as the unknown
                self.suggested_formulas[kendrick_base][unknown_index] = guesses_formulas
                self.suggested_errors[kendrick_base][unknown_index] = guesses_errors
                self.suggested_names[kendrick_base][unknown_index] = guesses_names

        print "finished matching"

    def _Error_on_assignment(self, suggested_formula, unknown_mass):

        """Provides ppm error on assignment of unknown peak with sugggested formula."""

        exact_x0 = self._Mass_calculator(suggested_formula)
        error = 1e6 * ((unknown_mass - exact_x0) / exact_x0)

        return error


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


    def _Weighted_guesses(self, list_of_tuple_pairs):

        """
        Takes a list of tuple pairs where 0th element of the tuple
        is the suggesting known compound and the 1st element of the
        tuple is the unknown compound suggestion. The elements in the
        passed list represent how many times a suggestion is made.

        returns a dictionary of suggestion keys and frequency vals.
        """

        keys = set([x[1] for x in list_of_tuple_pairs])
        string_list = [x[1] for x in list_of_tuple_pairs]

        dat = {}
        for key in keys:
            dat[key] = string_list.count(key)

        return dat


    def New_peaklist(self):

        """Return new peaklist with updated assignments."""

        new_peaklist = pd.DataFrame()

        new_peaklist['ion'] = self.ion
        new_peaklist['sumFormula'] = self.sumFormula
        new_peaklist['x0'] = self.x0

        for kendrick_base in self.suggested_formulas:

            suggested_formulas = self.suggested_formulas[kendrick_base]

            for i, entry in enumerate(suggested_formulas):

                # find highest frequency of a suggestion
                weighted_guess = self._Weighted_guesses(entry)
                try:
                    most_weight = max(weighted_guess, key=weighted_guess.get)
                    new_x0 = self._Mass_calculator(most_weight)

                    new_peaklist.loc[i, ['ion']] = most_weight+"-"
                    new_peaklist.loc[i, ['x0']] = new_x0

                except ValueError:
                    pass # weighted_guess returns {}

        return new_peaklist


    def Assignment_info(self):

        """Collect estimated assignment info into dataframe"""

        assignment_info = pd.DataFrame(columns = ["ion", "kendrick base", "suggested by", "suggested formula", "error / ppm", "names"])

        i = 0
        for kendrick_base in self.kendrick_bases:

            matched_unknowns = sorted(self._Matched(kendrick_base, known=False))

            for matched_unknown in matched_unknowns:

                unknown_index = self.ion.index(matched_unknown)
                formulas = self.suggested_formulas[kendrick_base][unknown_index]
                errors = self.suggested_errors[kendrick_base][unknown_index]
                names = self.suggested_names[kendrick_base][unknown_index]

                for formula, error, name in zip(formulas, errors, names):

                    assignment_info.loc[i] = [matched_unknown, kendrick_base, formula[0].encode('utf-8'), formula[1].encode('utf-8'), round(error[1],4), str([x.encode('utf-8') for x in name[1]]).replace("'","").replace("[","").replace("]","")]
                    i +=1

        assignment_info = assignment_info.sort_values(by=["ion","kendrick base"])

        return assignment_info


    def Run(self):

        # Calculate mass defect
        self._Calc_mass_defect(self._Mass_calculator(self.reagent_ion))

        # # Calculate kendrick mass defects
        [self._Calc_kendrick_mass_defect(base, self.kendrick_base_mass[i]) for i, base in enumerate(self.kendrick_bases)]

        # Stage 1. Match kendrick mass defects - still unknown
        [self._Match_peaks_on_kmd(base) for base in self.kendrick_bases]

        # Stage2. Get suggested formulas and compounds for matches
        self._Get_KMD_match_suggestions()
