#!/usr/bin/env python3
# Name: Nadja DiMartino (ndimarti)
# Group Members: Dan Schmeltee (dschmelt)

"""
This program outputs the sequence length, GC content, and relative codon frequency of a sequence.


"""

from sequenceAnalyzer2 import NucParams, FastAreader, ProteinParam


class genomeAnalyzer:
    def __init__(self, filename='testGenome.fa'):


            self.genomeAnalysis(filename)


    def genomeAnalysis(self, filename):

        """
        The 'Main' function. Does all printing and calculates codon usage.
        """
        print('Reading {}....'.format(filename))
        self.myReader = FastAreader(filename)

        length, nucParams = self.sequenceLength(self.myReader)
        print('sequence length = {:.2f}Mb'.format(length),
              "\n")

        self.gc = self.gcContent(nucParams)
        print('GC content = {:.1f}%'.format(self.gc),
            " \n")



        #Calculate codon usage
        codonCount = nucParams.codonComposition()
        for codons in sorted(codonCount): #Makes sure all codons are alphabetically sorted
            total = codonCount[codons]
            aa = nucParams.rnaCodonTable[codons]
            aaCount = nucParams.aaComposition()[aa]

            if aaCount != 0:
                finalTotal = (total/aaCount) * 100
            else:
                finalTotal = (total/1) * 100

            print('{} : {} {:5.1f}% ({:6d})'.format(codons, aa, finalTotal, total))



    ##determining the GC Content
    def gcContent(self, nucParams):

        nucComp = nucParams.nucComposition()
        gc = nucComp['G'] + nucComp['C']
        gc = gc/nucParams.nucCount() * 100

        return gc



    ##determining sequence length
    def sequenceLength(self, myReader):
        nucParams = NucParams('')
        for head, seq in myReader.readFasta():
            nucParams.addSequence(seq)

        length = nucParams.nucCount()/1000000
        return length , nucParams



genomeAnalyzer()
