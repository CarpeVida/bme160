#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: Nadja DiMartino (ndimarti)
 

"""
ProteinParam calculates statistics about a given protein sequence.
Input: Input a sequence of one letter amino acid codes like VLSPADKTNVKAAW
Output: Number of amino acids, molecular weight, molar extinction, mass extinction,
theoretical pI and amino acid composition.
"""


class ProteinParam :

    """
    This class contains many protein sequence analysis functions.
    """
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

# the __init__ method requires a protein string to be provided, either as a
# string or list of strings that will be concatenated

    def __init__ (self, protein):
        """
        Init parses the sequence and accepts only good amino acids.
        It also builds a few lists for later use.
        """
        upperProt = str.upper(protein)
        splitAminos = []
        allowedAminos = self.aa2mw.keys()
        for char in upperProt:
            if char in allowedAminos:
                splitAminos.append(char)
                
        l = ''.join(splitAminos).split()
        #joins strings with no separator adding: ''
        #splits joined 'protein' into a char array, assigns to l
        self.protString = ''.join(l).upper()
        #initialize protString to joined l and makes uppercase

        self.AAsum = 0
        self.compDict = {}
        #creates an empty dictionary object to represent # of each amino Acid
        for aminoAcid in self.aa2mw.keys():
        #initializes aminoAcid, loops through aa2mw.keys
            self.compDict[aminoAcid] = self.protString.count(aminoAcid)
            #assigns aminoAcid to the key, counts self.protString for the amino acid letter

    def aaCount (self):
        """
        aaCount takes in the joined, arrayed, and uppercased protString from input.
        It reutrns the number of characters in the given protString.
        """
        #return (len(self.protString))
        for count in self.compDict.keys():
            self.AAsum += self.compDict[count]
        return self.AAsum
            
        #returns len of input self's protString object

    def aaComposition (self) :
        """
        aaComposition takes the modified input, protString, and returns a dictionary.
        The one letter amino acid code is the key, scraped from aa2mw using for.
        The value is the number of that amino acid character in the string using count. 
        """
        return self.compDict


    def pI (self):
        """
        This function references _charge_ to find the pH of the protein with lowest charge.
        It takes no input, except indirectly the input sequence, and outputs pI.
        """
        if self.AAsum == 0:
            return 0.0
        else:
            temppH = 0.0
            bestpH = 42.0
            bestCharge = 999
            while( temppH <= 14):   #while pH is within normal range, loops through all 1400 numbers
                tempCharge = self._charge_(temppH)  #assigns temp charge to newly assigned pH's charge
                if abs(tempCharge) <= abs(bestCharge):    # and tempCharge >= 0:
                    bestpH = temppH
                    bestCharge = tempCharge
                temppH += 0.01
            return bestpH


    def _charge_ (self, pH):
        """
        _charge_ function takes in a float pH as an argument.
        Returns a float charge of the protein input at a specified pH.
        """
        """
        posCharge = 0.0
        negCharge = 0.0
        top = 0
        bottom = 0
        for aminoAcid in aa2chargePos.keys() and aaNterm:
            top = compDict[aminoAcid] * 10**aa2chargePos[aminoAcid]
            bottom = 10**aa2chargePos[aminoAcid]+(10 ** pH)
            posCharge += top/bottom
        for aminoAcid in aa2chargeNeg.keys() and aaCterm:
            top = compDict[aminoAcid] * 10**aa2chargeNeg[aminoAcid]
            bottom = 10**aa2chargeNeg[aminoAcid]+(10 ** pH)
            negCharge += top/bottom
        return posCharge-negCharge
        """
        posCharge = 0.0
        negCharge = 0.0
        top = 0.0
        bottom = 0.0
        for aminoAcid in self.protString:   #cycle through string, each loop initiallizing one letter to aminoAcid
            if aminoAcid in self.aa2chargePos.keys():   #if one letter code is in keys of positive AA dictionary
                top = 10 ** self.aa2chargePos[aminoAcid]  #initialize top to 10^pKa
                bottom = 10 ** self.aa2chargePos[aminoAcid] + (10 ** pH)  #initialize bottom to 10^pKa + pH
                posCharge += top/bottom
                
            elif aminoAcid in self.aa2chargeNeg.keys():
                top = 10**pH
                bottom = 10**self.aa2chargeNeg[aminoAcid]+ (10 ** pH)
                negCharge += top/bottom
            else:
                pass
        posCharge += (10**self.aaNterm)/(10**self.aaNterm + 10**pH)
        negCharge += (10**pH)/(10**self.aaCterm+10**pH)
        return posCharge - negCharge

    def molarExtinction (self):
        """
        molarExtinction calculates the light absorbance at 280nm using Y,W,C content
        """
        return float(self.protString.count('Y')*self.aa2abs280['Y']
        + self.protString.count('W')*self.aa2abs280['W']
        + self.protString.count('C')*self.aa2abs280['C'])



    def massExtinction (self):
        """
        massExtiction does not cause the die off of many species as the name suggests.
        It calculates light absorbance divided by molecular weight of the input protein
        """
        myMW =  self.molecularWeight()
        if myMW != 0:
            return self.molarExtinction() / myMW
        else:
            return 0.0

    def molecularWeight (self):
        """
        molecularWeight returns the protein's calculated molecular weight in g/mol
        This is calculated by summing individual AA weights and subtracting waters from hydrolysis
        """
        if self.AAsum == 0:
            return 0.0
        else:
            total = 0
            waterLoss = (len(self.protString) - 1)
            for newAmino in self.protString:
                total += self.aa2mw.get(newAmino)
            return total-(waterLoss*self.mwH2O)

# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys      #imports a module which directs traffic flows; stdin, stdout, stderr
for inString in sys.stdin :     #loops over sys.stdin and initializes to inString
    myParamMaker = ProteinParam(inString)   #cast inString as ProteinParam class so we can use contained funtions
    myAAnumber = myParamMaker.aaCount()     #initializes myAAnumber by calling .aaCount funtion
    print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
    print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
    print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
    print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
    print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
    print ("Amino acid composition:")
    myAAcomposition = myParamMaker.aaComposition()  #aaComp dictionary copied to myAAcomposition
    keys = list(myAAcomposition.keys()) #makes list from keys in dict
    keys.sort()
    if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
    for key in keys :
        print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))

