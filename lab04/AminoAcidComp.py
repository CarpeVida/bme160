#!/usr/bin/env python3
# Name: Brittney Wick (bwick)
# Group Members: Xander Pearce, Marisa Levario, Daniel Schmelter

'''
Program reads in a protein sequence, calculates the physical-chemical properties of a protein
sequence and prints out the following: total number of amino acids, total molecular weight,
molar extinction coefficient, mass extinction coefficient, theoretical isoelectric point (pI),
and amino acid composition.

Input: Protein sequence

Output: Prints out the number of amino acids, total molecular weight,
molar extinction coefficient, mass extinction coefficient,
theoretical isoelectric point (pI), and amino acid composition
'''

class ProteinParam:
    '''Class calculates the physical-chemical properties of a protein sequence.'''

    aa2mw = {                                                                                   # Table for calculating molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
            'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
            'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,               # As written, these are accessed as class attributes, for example:
            'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,               # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O
            'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
            }

    mwH2O = 18.015

# the __init__ method requires a protein string to be provided, either as a
# string or list of strings that will be concatenated

    def __init__ (self, protein):
        '''Builds a dictionary to store amino acid compostion.'''
        l = ''.join(protein).split()
        self.protString = ''.join(l).upper()
        self.aa_Comp = dict.fromkeys(list(ProteinParam.aa2mw.keys()),0)                         # Creates a new dictionary for the Amino Acid composition with starting values of 0.

        clean_string = ''                                                                       # Cleans protString to eliminate characters not included in aa_Comp dictionary.
        for aa in self.protString:
            if aa in ProteinParam.aa2mw:
                self.aa_Comp[aa] += 1
                clean_string += aa

        self.protString = clean_string                                                          # Renames clean_string as protString for easy use.

    def aaCount (self):
        '''Counts number of amino acids in string.'''
        return len(self.protString)                                                             # Returns the length of the new cleaned input amino acid string.

    def aaComposition (self):
        '''Describes composition of amino acid string.'''
        return self.aa_Comp                                                                                     # Returns amino acid compostion dictionary.

    def molecularWeight (self):
        '''Calculates the molecular weight of the amino acid string.'''
        if self.aaCount() > 0:
            temp_weight = 0
            for aa in self.protString:                                                                  # For-loop that iterates through the protein string and adds
                temp_weight += ProteinParam.aa2mw[aa]                                                   # together the weight of each amino acid in the string.
            return temp_weight
        #     mol_weight = (temp_weight) - (ProteinParam.mwH2O)*(self.aaCount() - 1)                      # Equation for the total molecular weight of protein string.
        #     return mol_weight
        # else:
        #     mol_weight = 0                                                                              # In the case that there are no amino acids present in the string, return 0.
        #     return mol_weight


# Please do not modify any of the following. This will produce a standard output that can be parsed

import sys

for inString in sys.stdin :
    myParamMaker = ProteinParam(inString)
    myAAnumber = myParamMaker.aaCount()
    print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))                             # Prints out each of the physical-chemical properties of the protein sequence.
    print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
    # print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
    print ("AA Molar Percentage:   MolarMass:    MassOfAAPercentage:       MassPercentage")
    myAAcomposition = myParamMaker.aaComposition()
    keys = list(myAAcomposition.keys())
    keys.sort()
    if myAAnumber == 0 : myAAnumber = 1                                                             # Handles the case where no AA are present
    aa_mass = 0
    total_mass = 129.216
    mass_percent = 0
    mass_percent_sum = 0
    for key in keys:
        percent_comp = myAAcomposition[key] / myAAnumber
        aa_mass = percent_comp * ProteinParam.aa2mw[key]
        #total_mass += percent_comp * ProteinParam.aa2mw[key]
        mass_percent = aa_mass/total_mass
        mass_percent_sum += mass_percent
        print ("{} = {:.2%} \t\t\t\t{} \t\t\t {:.3}\t\t\t\t\t\t{:.3}".format(key,  percent_comp, ProteinParam.aa2mw[key], aa_mass, mass_percent))
    print("Weighted average mass of all amino acids is {}, sum of mass percentages is {:.3}".format(total_mass, mass_percent_sum))
