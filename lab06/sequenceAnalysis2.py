#!/usr/bin/env python3
# Name: Louis Chatta (lchatta)
# Group Members: (lchatta, cjavery)

class ProteinParam :
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
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}

    # the __init__ method requires a protein string to be provided, either as a
    # string or list of strings that will be concatenated

    """
    Here I loop through all the characters in the protein and and check if
    they are allowed input, if so add to the splitAminos array. Then create
    the protein string
    """

    def __init__ (self, protein):
        splitAminos = []
        allowed_aminos = self.aa2mw.keys()
        for char in protein:
            if char in allowed_aminos:
                splitAminos.append(char)

        l = ''.join(splitAminos).split()
        self.protString = ''.join(l).upper()



#Here the Count of the protein is returned
    def aaCount (self):
        return len(self.protString)

#Here the PI is calculated by finding the ph at which the net charge is closest to 0 and non-negative
    def pI (self):
        ph = 0.00
        lowestCharge = (99999,0.0)
        while(ph <= 14.0):
            currentCharge = self._charge_(ph)
            if currentCharge < lowestCharge[0] and self._charge_(ph) >=0:
                lowestCharge = (currentCharge,ph)
            ph += 0.01

        return lowestCharge[1]

    """
    Here the dictionary of amino acids to their count is created and returned
    """
    def aaComposition (self) :
        compDict = {}
        for amino in self.aa2mw.keys():
            compDict[amino] = self.protString.count(amino)
        return compDict

    """
    Here the charge is calculated via the formula provided. Moreover, the acutual work is
    done in the dotPh
    """
    def _charge_ (self, pH):
        firstCharge = self.dotPh(['R','K','H'],pH) + (10**self.aaNterm)/(10**self.aaNterm+10**pH)
        secondCharge = self.dotPh(['D','E','C','Y'],pH) + (10**pH)/(10**self.aaCterm+10**pH)
        return firstCharge - secondCharge
    """
    Here the inner sum of each of the aminos is calculated
    """
    def dotPh(self,allowed_aminos,pH):
        totalFirstCharge = 0
        for amino in self.protString:
            if amino in allowed_aminos:
                if amino in self.aa2chargePos.keys():
                    t = 10**self.aa2chargePos[amino]
                    b = 10**self.aa2chargePos[amino] + 10**pH
                    totalFirstCharge += t/b
                elif amino in self.aa2chargeNeg.keys():

                    t = 10**pH
                    b = 10**self.aa2chargeNeg[amino] + 10**pH
                    totalFirstCharge += t/b
                else:
                    t = 10**0
                    b = 10**0 + 10**pH
                    totalFirstCharge += t/b
        return totalFirstCharge

    """
    Here the associated count of the aminos is multiplied by the associated constant
    according to the formula provided.
    """
    def molarExtinction (self):
       return float(self.protString.count('Y')*self.aa2abs280['Y']
        + self.protString.count('W')*self.aa2abs280['W']
        + self.protString.count('C')*self.aa2abs280['C'])

    def massExtinction (self):
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    """
    Here the total molecularWeight is calculated by looping
    through the sum of the molecularWeights associated to the amino acid
    and then subtracting the waters for the peptide bonds.
    As can be seen in the fomula.
    """
    def molecularWeight (self):
        firstMolecWeight = 0.0
        length = sum(len(s) for s in self.protString)

        for amino in self.protString:
            firstMolecWeight += self.aa2mw[amino]
        return (firstMolecWeight - self.mwH2O * (length-1))

class NucParams:
    """
    This class handles reading and processing data about a given genome sequence
    methods:
    __init__ (self)
    addSequence (self, thisSequence)
    aaComposition(self)
    nucComposition(self)
    codonComposition(self)
    nucCount(self)

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
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    allowedBases = {'A','C','G','T','U','N','"','P'}
    aa = {'A', 'G', 'M', 'S', 'C', 'H', 'N', 'T', 'D', 'I', 'P', 'V', 'E', 'K', 'Q', 'W', 'F', 'L', 'R', 'Y'}
    aaComp = {}
    aminoAcid = []
    aminoAcidString = ''
    nucComp = {}
    codonComp = {}

    def __init__ (self,seq):
        self.seq = self.stripSequence(seq)

    #builds all neccessary properties
    def buildNucleotide(self):
        self._toAminoAcid()
        self._aaComposition()
        self._nucComposition()
        self._codonComposition()

    #adds a genome sequence to the NucParams
    def addSequence (self, thisSequence):
        self.seq += self.stripSequence(thisSequence)

    #cleans a string of unwanted bases, returns uppercase genome string
    def stripSequence(self, protein):
        splitSeq = []
        for char in protein:
            charUp = char.upper()
            # if charUp in self.allowedBases:
            if charUp == '-': continue
            splitSeq.append(charUp)

        return ''.join(''.join(splitSeq).split()).upper()

    #returns a dictionary with counts of the 20 amino acids
    def _aaComposition(self):
        for amino in ProteinParam.aa2mw.keys():
            self.aaComp[amino] = self.aminoAcid.count(amino)

    #converts sequence to amino acid
    def _toAminoAcid(self):
        aminoAcid = []
        newAA =  [self.seq[i:i+3] for i in range(0, len(self.seq), 3)]
        self.codons = newAA
        isRna = -1
        for codon in newAA:
            if codon in self.rnaCodonTable:
                aminoAcid.append(self.rnaCodonTable[codon])
            elif codon in self.dnaCodonTable:
                aminoAcid.append(self.dnaCodonTable[codon])

        self.aminoAcid = aminoAcid
        self.aminoAcidString = ''.join(aminoAcid ).upper()

    #returns a dictionary with counts of valid nucleotides
    def _nucComposition(self):
        for base in self.allowedBases:
            self.nucComp[base] = self.seq.count(base)

    #returns a dictionary with counts of valid codons
    def _codonComposition(self):
        for amino in self.rnaCodonTable.keys():
            if 'N' in amino: continue
            self.codonComp[amino] = self.codons.count(amino)

    #returns a dictionary with counts of the 20 amino acids
    def aaComposition(self):
        return self.aaComp

    #converts sequence to amino acid
    def toAminoAcid(self):
        return self.aminoAcidString

    #returns a dictionary with counts of valid nucleotides
    def nucComposition(self):
        return self.nucComp

    #returns a dictionary with counts of valid codons
    def codonComposition(self):
        return self.codonComp

    #returns integer count of sum of valid nucleotides
    def nucCount(self):
        return len(self.seq)



class FastAreader :
    '''
    Class to provide reading of a file containing one or more FASTA
    formatted sequences:
    object instantiation:
    FastAreader(<file name>):

    object attributes:
    fname: the initial file name

    methods:
    readFasta() : returns header and sequence as strings.
    Author: David Bernick
    Date: April 19, 2013
    '''
    def __init__ (self, fname):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def readFasta (self):
        '''
        using filename given in init, returns each included FastA record
        as 2 strings - header and sequence.
        whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header.
        '''
        header = ''
        sequence = ''

        with open(self.fname) as fileH:
            # initialize return containers
            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence

    def readFastaStdIn (self,inp):
        '''
        using filename given in init, returns each included FastA record
        as 2 strings - header and sequence.
        whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header.
        '''
        header = ''
        sequence = ''

        with inp as fileH:
            # initialize return containers
            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence
