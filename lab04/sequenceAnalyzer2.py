
#!/usr/bin/env python3

# Name: Nadja DiMartino (ndimarti)

# Group Members: Dan Schmeltee (dschmelt)


class ProteinParam:

    """
    This program calculates the physical-chemical properties of a protein sequence.
    Input: String of Amino Acids. Created in help by Dan Schmeltee (dschmelt).
    Output: Amino Acid Composition, Theoretical isoelectric point (pI), the Molar extinction coefficient,
    and Mass extinction coefficient.


    EX: VLSPADKTNVKAAW
    """

# These tables are for calculating:

#   molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
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

        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189}
    mwH2O = 18.015

    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}



    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}

    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}

    aaNterm = 9.69

    aaCterm = 2.34

    myaaCount = 0

# the __init__ method requires a protein string to be provided, either as a
# string or list of strings that will be concatenated

    def __init__ (self, protein):
        """
        This is where the amino acid dictionary is stored. If the inputted sequences contains
        one of the amino acids, it adds a tally.
        """
        l = ''.join(protein).split()

        self.protString = ''.join(l).upper()

        self.aaComp = {'A' : 0, 'R' : 0, 'N' : 0, 'D' : 0, 'C' : 0, 'E' : 0, 'Q'
                   : 0, 'G' : 0, 'H' : 0, 'I' : 0, 'L' : 0, 'K' : 0, 'M' : 0,
                   'F' : 0, 'P' : 0, 'S' : 0, 'T' : 0, 'W' : 0, 'Y' : 0, 'V' : 0,
                    '-' : 0}

        for a in protein:
            if a in self.aaComp:
                self.aaComp[a] += 1



        self.myaaCount = sum(self.aaComp.values())

    def aaCount (self):
        """
        Here is where the count of the protein sequence is returned.
        """

        return self.myaaCount



    def pI (self):
        """The Theoretical Isoelectric point is calculated by finiding the pH at which the net
        charge is closet to 0 and non-negative.
        """
        pH = 0.00
        chargeClosestTo0 = (99999,0.0)
        while (pH <= 14.00):
            self.currentCharge = abs(self._charge_(pH))
            if self.currentCharge < chargeClosestTo0[0]:
                chargeClosestTo0 = (self.currentCharge, pH)
            pH += 0.01

        return chargeClosestTo0[1]


    def aaComposition (self) :
        """
        Since the dictionary is already in __init__, it just returns the dictionary.
        """
        return self.aaComp



    def _charge_ (self, pH):
        """
        The charge is calculated using the formula in the for loop.
        """
        self.posCharge = 0
        self.negCharge = 0


        for aa in self.aa2chargePos.keys():

                self.posCharge += (self.aaComp[aa] * 10**self.aa2chargePos[aa]) / (10**self.aa2chargePos[aa] + 10**pH)

        for aa in self.aa2chargeNeg.keys():
                self.negCharge += 10**pH * self.aaComp[aa] / (10**self.aa2chargeNeg[aa] + 10**pH)



        sumPosCharge = self.posCharge + (10**self.aaNterm)/(10**self.aaNterm+10**pH)
        sumNegCharge = self.negCharge + (10**pH)/(10**self.aaCterm+10**pH)

        self.totalCharge = sumPosCharge - sumNegCharge
        return self.totalCharge



    def molarExtinction (self):
        """
        We create a formula using amino acid absorbance at 280nm to find the molar extinction coefficient,
         which is used for the mass extinction coefficient
         """
        extinction = float(0)
        for aa in self.aa2abs280.keys():
            extinction += self.aaComp[aa] * self.aa2abs280[aa]

        return extinction


    def massExtinction (self):
        myMW =  self.molecularWeight()

        return self.molarExtinction() / myMW if myMW else 0.0



    def molecularWeight (self):
        """
        Total molecular weight is calculated by adding the cumulative molecular weight and subtracting the water weight.

        """
        sumMolecularWeight = 0.0
        for aa in self.protString:
            if aa in self.aa2mw:
                sumMolecularWeight += self.aa2mw[aa]

        molecularWeight = sumMolecularWeight - ((self.mwH2O)*(self.myaaCount - 1))
       #(myAAnumber * self.aa2mw) - ( self.mwH2O['a']*(myAAnumber - 1))

        return molecularWeight


class NucParams:
    """
    This class handles processing data from a given sequence
    Input: DNA and RNA sequence

    Methods:
    __init__(self, seq)
    addSequence(self, protein)
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

    allowedBases = {'A' ,'C','G','T', 'U','N'}



    def __init__ (self,seq):
        self.aaComp = {'A' : 0, 'R' : 0, 'N' : 0, 'D' : 0, 'C' : 0, 'E' : 0, 'Q'
                   : 0, 'G' : 0, 'H' : 0, 'I' : 0, 'L' : 0, 'K' : 0, 'M' : 0,
                   'F' : 0, 'P' : 0, 'S' : 0, 'T' : 0, 'W' : 0, 'Y' : 0, 'V' : 0,
                    '-' : 0}

        self.codonComp = {key:0 for key in self.rnaCodonTable.keys()}
        self.nucleotideComp = {key:0 for key in self.allowedBases}
        self.addSequence(seq)

    def addSequence (self, sequence):
        """This cleans a sequence of unwanted bases and returns string in uppercase.
        It also adds it to the NucParams"""
        sequence = self.seqCleaner(sequence)

        for base in self.allowedBases:
            self.nucleotideComp[base] += sequence.count(base)

        aminoAcid = []
        newAA =  [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        for codons in newAA:
            codons = codons.replace('T','U')
            if codons in self.codonComp:
                self.codonComp[codons] += 1

                aa = self.rnaCodonTable[codons]
                self.aaComp[aa] += 1

    def seqCleaner(self, sequence):
        return ''.join(sequence.split()).upper()


    def aaComposition(self):
        """
        returns dictionary with count of all amino acids
        """
        return self.aaComp

    def nucComposition(self):
        """
        returns a dictionary with counts of valid bases
        """

        return self.nucleotideComp


    def codonComposition(self): ####This is the error site
        """
        returns dictionary with count of codons
        """

        return self.codonComp



    def nucCount(self):
        """
        returns integer count of sum of valid nucleotides
        """
        return sum(self.nucleotideComp.values())



import sys

class FastAreader :

    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)


    def readFasta (self):

        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header,sequence

# presumed object instantiation and example usage

# myReader = FastAreader ('testTiny.fa');

# for head, seq in myReader.readFasta() :

#     print (head,seq)
