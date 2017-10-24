#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: None


class sequenceAnalysis:
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
    validAminos = {value in rnaCodonTable.items()}
    aminoComp = {}
    aminoAcid = []
    aminoAcidString = ''
    nucComp = {}
    codonComp = {}

    def __init__(self):
        """
        Takes in one string and formats it
        """

    self.seq = self.stripSequnce(seq)


def addSequence(self, thisSequence):
    """
    Adds a sequence to be read by NucParams
    """


self.seq += self.stripSequnce(thisSequence)


def stripSequence(self, sequence)
    """
    Iterates through input sequence, makes list all upper case,
    Adds valid aminos to a list iff base is valid
    Returns the formatted, verified sequence
    """


splitSeq = []
for char in sequence:
    charUpper = char.upper()
    if charUpper in self.validBases:
        splitSeq.append(charUpper)

return ''.join(''.join(splitSeq).split()).upper()


def aaComposition(self):
    """
    Goes through valid amino acids list, 
    counts each one's contents in aminoAcid and fills aminoComp dictionary
    """


for amino in self.validAminos():
    self.aminoComp[amino] = self.aminoAcid.count(aminos)


def nucComposition(self):
    """
    Goes through valid bases list 
    counts each one's contents in self.seq and fills nucComp dictionary
    """


for nucleotide in self.validBases:
    self.nucComp
    {nucleotide} = self.seq.count(nucleotide)


def codonComposition(self):
    """
    returns dictionary codonComp with counts of valid codons
    """


self.codonComp = self.rnaCodonTable
self.codonComp = {x: 0 for x in self.codonComp}

"""
for temp in self.codonComp:
    self.codonComp.item[] = 0

for amino in self.rnaCodonTable.keys():
    if "N" in amino: continue
    self.codonComp[amino] = self.codons.count(amino)
"""
multiplier = 1
for self.seq:
    if "N" not in self.seq[0 * multiplier: 2 * multiplier]
        self.codonComp[0 * multiplier: 2 * multiplier].item += 1
    multiplier += 1


def nucCount(self):
    """

    """


return len(self.seq)

import sys


class FastAreader:
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
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


        # presumed object instantiation and example usage
        # myReader = FastAreader ('testTiny.fa');
        # for head, seq in myReader.readFasta() :
        #     print (head,seq)