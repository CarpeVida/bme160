#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: Nadja DiMartino (ndimarti)

import sys


class ProteinParam:
    """
    ProteinParam calculates statistics about a given protein sequence.
    Input: Input a sequence of one letter amino acid codes like VLSPADKTNVKAAW
    Output: Number of amino acids, molecular weight, molar extinction, mass extinction,
    theoretical pI and amino acid composition.
    """
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein):
        """
        Init takes in the protein sequence and accepts only good amino acids.
        It also builds a few lists for later use.
        """
        upperProt = str.upper(protein)
        splitAminos = []
        allowedAminos = self.aa2mw.keys()
        for char in upperProt:
            if char in allowedAminos:
                splitAminos.append(char)

        list = ''.join(splitAminos).split()
        # joins strings with no separator adding: ''
        # splits joined 'protein' into a char array, assigns to l
        self.protString = ''.join(list).upper()
        # initialize protString to joined l and makes uppercase

        self.compDict = {}
        # creates an empty dictionary object to represent # of each amino Acid
        for aminoAcid in self.aa2mw.keys():
            # initializes aminoAcid, loops through aa2mw.keys
            self.compDict[aminoAcid] = self.protString.count(aminoAcid)
            # assigns aminoAcid to the key, counts self.protString for the amino acid letter

    def aaCount(self):
        """
        aaCount takes in the joined, arrayed, and uppercased protString from input.
        It reutrns the number of characters in the given protString.
        """
        #print(len(self.protString))
        return (len(self.protString))
        # returns len of input self's protString object

    def aaComposition(self):
        """
        aaComposition takes the modified input, protString, and returns a dictionary.
        The one letter amino acid code is the key, scraped from aa2mw using for.
        The value is the number of that amino acid character in the string using count.
        """
        #print("gets to aaComp")
        return self.compDict

    def molarExtinction(self):
        """
        molarExtinction calculates the light absorbance at 280nm using Y,W,C content
        """
        #print("gets to molarExtinction")
        return float(self.protString.count('Y') * self.aa2abs280['Y']
                     + self.protString.count('W') * self.aa2abs280['W']
                     + self.protString.count('C') * self.aa2abs280['C'])

    def massExtinction(self):
        """
        massExtiction does not cause the die off of many species as the name suggests.
        It calculates light absorbance divided by molecular weight of the input protein
        """
        myMW = self.molecularWeight()
        if myMW != 0:
            return self.molarExtinction() / myMW
        else:
            return 0.0

    def molecularWeight(self):
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
            return total - (waterLoss * self.mwH2O)
    def pI(self):
        """
        This function references _charge_ to find the pH of the protein with lowest charge.
        It takes no input, except indirectly the input sequence, and outputs pI.
        """
        #print(len(self.protString))

        if len(self.protString) == 0:
            return 0.0
        else:
            temppH = 0.0
            bestpH = 99.0
            bestCharge = 999
            while (temppH <= 14):  # while pH is within normal range, loops through all 1400 numbers
                tempCharge = self._charge_(temppH)  # assigns temp charge to newly assigned pH's charge
                if abs(tempCharge) <= abs(bestCharge):  # and tempCharge >= 0:
                    bestpH = temppH
                    bestCharge = tempCharge
                temppH += 0.01
            print(bestpH)
            return bestpH

    def _charge_(self, pH):
        """
        _charge_ function takes in a float pH as an argument.
        Returns a float charge of the protein input at a specified pH.
        """
        posCharge = 0.0
        negCharge = 0.0
        top = 0.0
        bottom = 0.0
        for aminoAcid in self.protString:  # cycle through string, each loop initiallizing one letter to aminoAcid
            if aminoAcid in self.aa2chargePos.keys():  # if one letter code is in keys of positive AA dictionary
                top = 10 ** self.aa2chargePos[aminoAcid]  # initialize top to 10^pKa
                bottom = 10 ** self.aa2chargePos[aminoAcid] + (10 ** pH)  # initialize bottom to 10^pKa + pH
                posCharge += top / bottom

            elif aminoAcid in self.aa2chargeNeg.keys():
                top = 10 ** pH
                bottom = 10 ** self.aa2chargeNeg[aminoAcid] + (10 ** pH)
                negCharge += top / bottom
            else:
                pass
        posCharge += (10 ** self.aaNterm) / (10 ** self.aaNterm + 10 ** pH)
        negCharge += (10 ** pH) / (10 ** self.aaCterm + 10 ** pH)
        #print(posCharge - negCharge)
        return posCharge - negCharge


class NucParams:
    """
    This class contains a series of analysis functions.
    """
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}
    validBases = {"A", "T", "C", "G", "U", "N"}
    validAminos= {'A', 'G', 'M', 'S', 'C', 'H', 'N', 'T', 'D', 'I', 'P', 'V', 'E', 'K', 'Q', 'W', 'F', 'L', 'R', 'Y'}
    aaComp = {}
    aminoList = []
    aminoAcidString = ''
    nucComp = {}
    codonComp = {}

    def __init__(self, nuc_seq):
        """
        Initializes a dictionary of one letter keys with 0 values
        Initializes only a valid list of nucleic acids.
        Takes input sequence of upper or lowercase ACGTUN,


        self.aaComp = {v:k for k,v in self.rnaCodonTable.items()}
        self.aaComp = {x:0 for x in self.aaComp}

        self.codonComp = {key:0 for key in self.rnaCodonTable.keys()}
        self.nucComp = {key:0 for key in self.validBases}
        self.addSequence(nuc_seq)
        """
        self.nuc_seq = self.stripSequence(nuc_seq)

    def initializeObjects(self):
        """
        Links all the methods to the object
        """
        self.toAminoAcid()
        self.aaComposition()
        self.nucComposition()
        self.codonComposition()

    def addSequence (self, this_sequence):
        """
        addSeq takes in a string sequence, runs it through uppercaser
        makes a nucleotide base count dictionary
        Parses this_sequence for codons and fills codon count dictionary codonComp

        this_sequence = self.stripSequence(this_sequence)

        #codon parser that fills codonComp count dictionary
        newAA =  [this_sequence[i:i+3] for i in range(0, len(this_sequence), 3)]
        for codons in newAA:
            codons = codons.replace('T','U')
            if codons in self.codonComp:
                self.codonComp[codons] += 1

                aa = self.rnaCodonTable[codons]
                self.aaComp[aa] += 1
        """
        #links the sequence to class, validates, removes spaces, uppercases
        self.nuc_seq += self.stripSequence(this_sequence)


    def stripSequence(self, sequence):
        """
        Uppercases, validates, and removes spaces in the sequence
        Returns as string
        """
        splitSeq = []
        for char in sequence:
            charUp = char.upper()
            if charUp in self.validBases:
                splitSeq.append(charUp)

        return ''.join(sequence.split()).upper()

    def aaComposition(self):
        """
        Returns a dictionary of counts of each of the 20 amino acids

        Goes through valid aminos list,
        Initialize dict with each amino as key,
        value initialized as count of amino acids in list
        Must have aminoList initialized first
        """
        for amino in self.validAminos:
            self.aaComp[amino] = self.aminoList.count(amino)

    def toAminoAcid(self):
        """
        toAminoAcid converts the object's nucleotide string into amino acids
        """
        aminoList = []
        #codon parser that fills newAA with count dictionary
        codonList = [self.nuc_seq[i:i+3] for i in range(0, len(self.nuc_seq), 3)]
        self.codonList = codonList
        isRna = -1
        for codon in codonList:
            if codon in self.rnaCodonTable:
                aminoList.append(self.rnaCodonTable[codon])
            elif codon in self.dnaCodonTable:
                aminoList.append(self.dnaCodonTable[codon])
        self.aminoList = aminoList
        self.aminoAcidString = ''.join(self.aminoList).upper()


    def nucComposition(self):
        """
        Goes through valid bases list
        counts each one's contents in self.nuc_seq and fills nucComp dictionary
        """
        for base in self.validBases:
            self.nucComp[base] = self.nuc_seq.count(base)

    def codonComposition(self):
        """
        returns rna codon count dictionary
        Discards codons with N bases in the codons
        """

        for codon in self.rnaCodonTable.keys():
            if 'N' in codon: continue
            self.codonComp[codon] = self.codonList.count(codon)

    def nucCount(self):
        """
        nucCount returns the total length on the input seq
        """
        return len(self.nuc_seq)

class FastAreader:
    """
    Fasta reader takes in a file, initialize filename to fname
    then returns header and seq objects
    """
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        Fsequence = ''

        with self.doOpen() as fileH:

            header = ''
            Fsequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, Fsequence
                    header = line[1:].rstrip()
                    Fsequence = ''
                else:
                    Fsequence += ''.join(line.rstrip().split()).upper()
        yield header, Fsequence

