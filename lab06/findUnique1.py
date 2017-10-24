#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: None

"""
Find uniques must take fasta files in from command line
parses each of the tRna sequences into sets and subsets
finds the union of each subset in other tRNA files
if found, removes subset from the set

It parses through the sequence data
    Makes allHeadsList of all heads,
    Makes seqSubseqDict of all sequences and their subsequences
    Makes allElseSubseq list of all but seq
    Find unique by doing subsequence - allElseSubseq
    Save seqUniqueDict for each sequence
    Remove unneccesary superstrings
    Print with header, seq, and unique minimums with dot spacers
"""

from collections import OrderedDict
import sys
allSubSequencesList = []
UniEssList = []

class ExternalMethods:
    def get_substrings(string):
        """
        This function splits and generates substrings
        """
        length = len(string)
        for firstPos in range(length):
            for endPos in range(firstPos+1, length +1):
                yield(string[firstPos: endPos])


    def stripSequence(sequence):
        """
        Uppercases, validates, and removes spaces in the sequence
        Returns as string
        """
        splitSeq = []
        for char in sequence:
            charUp = char.upper()
            #if charUp in self.validBases:
            if charUp == '-': continue
            splitSeq.append(charUp)

        #returns string of joined list elements
        return ''.join(''.join(splitSeq).split())

    def readFastaStdIn(input):
        """
        Takes in file from init
        Removes > from header and all whitespace
        Yields header and sequence strings
        """

        header = ''
        sequence = ''

        with input as fileH:
            # declare return objects
            header = ''
            sequence = ''

            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            #header is saved, get rest of seq
            #until next header is found
            #yield results and await next call
            #next call resumes at yield point
            #which is where we have next header

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header, sequence



class FindUnique:
    """
    This class takes in bases, head, and list of subsequences
    It parses through the sequence data
    Makes allHeadsList of all heads,
    Makes seqSubseqDict of all sequences and their subsequences
    Makes allElseSubseq list of all but seq
    Find unique by doing subsequence - allElseSubseq
    Save seqUniqueDict for each sequence
    init obj of Class, run superRemover
    Remove unneccesary superstrings
    Print with header, seq, and unique minimums with dot spacers
    """
    def __init__(self, head, tRNAseq, fullUniquesubSeq):
        """
        Gets things initialized
        """
        #self.orderedSubSequences = {}  #iterater
        self.header = head #saves the file's header
        self.tRNAseq = tRNAseq
        self.subSequences = fullUniquesubSeq
        self.uniqueEssentialSubstrings = set()
        self.uniqueEssentialSubstrings=self.removeAllSuperStrings(list(fullUniquesubSeq))
        #print(str(self.uniqueEssentialSubstrings))
        #print(self.uniqueNonSuperSubstrings)
        self.sortableHeader = head #.split('|')[1]

    def get_substrings(self, string):
        """
        This is the function that splits and generates substrings
        """
        length = len(string)
        for firstPos in range(length):
            for endPos in range(firstPos+1, length +1):
                yield(string[firstPos: endPos])

    def get_biggs(self, uniqueSubSeqList, originalSeq):
        """
        get_biggs appends list with extensions forward of each subseq
        then appends list extensions in rev of each subseq
        it returns a list of all
        """
        allBigs = set()


        #count = 0
        for eachUniqueSub in uniqueSubSeqList:
            bigs = set()
            startPos = originalSeq.find(eachUniqueSub)
            stopPos = startPos + len(eachUniqueSub)
            for i in range(0,startPos+1):
                for j in range(stopPos, len(originalSeq)+1):
                    bigs.add(originalSeq[i:j])
            #print(len(allBigs))
            #count += 1
            #print(count)
            try:
                bigs.remove(str(eachUniqueSub))
            except KeyError:
                pass
            allBigs = allBigs | bigs
            #print(len(bigs))
            # AT AAAATTTT

        return allBigs


    def removeAllSuperStrings(self, uniqueSubseqList):
        """
        This code takes input subSeqList
        Finds each element in rest of list
        Searches one by one and removes non-minimal substrings
        """
        bigsSet = set(self.get_biggs(uniqueSubseqList, self.tRNAseq))
        #print(len(uniqueSubseqList))
        #print("Sep\n\n")
        #print(len(bigsSet))
        uniEssList = set(uniqueSubseqList) - bigsSet
        #print("\n\n\n" + str(uniEssList) + "\n\n\n")
        return uniEssList

    def printSubstrings(self):
        """
        This function prints substrings to stdout
        """
        sortDic = OrderedDict()
        print(self.header)
        print(self.tRNAseq)
        #print(self.uniqueEssentialSubstrings)
        #print(len(self.subSequences))
        #print((self.uniqueEssentialSubstrings))
        for essentialUnique in self.uniqueEssentialSubstrings:
            pos = self.tRNAseq.find(essentialUnique)
            sortDic[pos] = essentialUnique

        finsortDic = OrderedDict(sorted(sortDic.items()))
        for it in finsortDic:
            print("."*it + finsortDic[it])

    """def get_subStrings(self, string):
        """"""
        This function takes in a string sequence
        Returns a list of subsequences
        """"""
        subStringList = []

        for i in range(len(string)):
            for j in range(len(string)):
                #if len(string[i:j]) > 7: continue
                self.orderedSubSequences[string[i:j +1]] = i

                subStringList.append(string[i:j+1])
                #j += 1
            #i += 1
        return subStringList
    """

def main():
    """
    Main function takes in file
    Parses arguments and runs substring sorter

    """

    InputControl()
    """
    #checks to see if file is entered, if not, print error and exit
    if len(sys.argv) != 2:
        print("Please enter a file after the program")
        sys.exit()
    #if file is entered, initialize it to infile
    if len(sys.argv) == 2:
        infile = open(sys.argv[1])
    """
def InputControl():
    for closedfile in sys.stdin:
        infile = open(str.rstrip(closedfile, "\n"))
        #print(orfReader.readFastaStdIn(infile))
        #create lists for future storage
        allHeadsList= []
        subSequence = []
        seqSubseqDict = OrderedDict()
        seqUniqueSubseqDict = {}
        """
        takes infile, cleans, makes subseq
        saves cleanseq and subseq to dictionary
        makes list of all heads
        makes list of all subseqs
        """
        for head, seq in ExternalMethods.readFastaStdIn(infile):
            cleanSeq = ExternalMethods.stripSequence(seq)
            subSequence = list(ExternalMethods.get_substrings(cleanSeq))
            seqSubseqDict[cleanSeq] = subSequence
            allHeadsList.append(head)
            allSubSequencesList.append((subSequence))

        #print("seqSubSeqDict: " + str(seqSubseqDict.keys()))
        i = 0
        tRNAList = []

        #works, seqSubseqDict is filled with 22 sets
        for thisSeq in seqSubseqDict:
            #print("seq: " + thisSeq)
            allElseSubseq = []
            for otherSeq in seqSubseqDict:
                if otherSeq != thisSeq:
                    allElseSubseq.append(set(seqSubseqDict[otherSeq]))
            #print((seqSubseqDict[thisSeq]))
            #print(len(allElseSubseq))
            uniqueSubList = set(seqSubseqDict[thisSeq]) - set.union(*allElseSubseq)
            seqUniqueSubseqDict[thisSeq] = uniqueSubList
            #print(str(i) + ": " + str(len(uniqueSubList)))
            #print(thisSeq + str(len(seqUniqueSubseqDict[thisSeq])))
            tRNAList.append(FindUnique(allHeadsList[i].lstrip(), thisSeq, uniqueSubList))
            i += 1

        #print("SeqSub len: " + str(len((seqSubseqDict))))
        print("Length of tRNAlist: " + str(len(tRNAList)))
        for trna in sorted(tRNAList, key=lambda x:x.sortableHeader):
            trna.printSubstrings()

main()
"""
Potential Pseudocode

class UniqueFinder():
    class Unique:
        pass
def main():
    uniques = UniqueFinder()
    for head,seq in seqs:
        uniques.add(head,seq)

    for unique in uniques.getUniques() :
        unique.print()
"""
