#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: Nadja DiMartino (ndimarti)


from sequenceAnalysis import ProteinParams, NucParams, FastAreader

"""
This program is the printer which prints seq length, GC content, and codon frequency
"""

class GenomeCompare:

    def __init__(self, filenames = ['testGenome.fa'], compare= False):
        if compare:
            self.compareAnalysis(filenames)
        else:
            for fa in filenames:
                self.individualAnalysis(fa)

    def compareAnalysis(self, filenames):
        pass

    def individualAnalysis(self, fastafile):
        print('Delete line 20, reading from {}".format(fastafile)')
        myRead = FastAreader(fastafile)
        length, nucParams = self.sequenceLength(myRead)
        print("sequence length = {:1.2f} Mb".format(length))
        print("\n")

        gc = self.gcContent(nucParams)
        print("GC content = {:2.1}%".format(gc))
        print("\n")

        codonCount = nucParams.codonComposition()
        aminosComp = nucParams.aminoCompo()
        nucComp = nucParams.nucComposition
        for codon in sorted(codonCount):
            total = codonCount[codon]
            amino = nucParams.rnaCodonTable[codon]
            aminoCount = nucParams.aminoAcid.count(amino)

            finalTotal = total/aminoCount * 100
            print('{} : {} {:5.1f}% ({:6d})'.format(codon, amino, finalTotal, total))

    def gcContent(self, nucParams):
        nucComp = nucParams.nucComposition()
        gc = nucComp['G'] + nucComp['C']
        gc = gc/nucParams.nucCount() * 100
        return gc

    def sequenceLength(self, myReader):

        nucParams = NucParams()
        for head, seq in myReader.readFasta():
            nucParams.addSequence(seq)
        nucParams.buildNucleotide()

        length = nucParams.nucCount()/ 1000000

        return length, nucParams

genomeCompare()



