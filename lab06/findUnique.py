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

    def get_substringList(file, seqSubseqDict):
        allHeadsList = []
        for head, seq in ExternalMethods.readFastaStdIn(file):
            cleanSeq = ExternalMethods.stripSequence(seq)
            subSequence = list(ExternalMethods.get_substrings(cleanSeq))
            seqSubseqDict[cleanSeq] = subSequence
            allHeadsList.append(head)
            allSubSequencesList.append((subSequence))
        return seqSubseqDict, allHeadsList, allSubSequencesList

    def stripSequence(sequence):
        """
        Uppercases, validates, and removes spaces in the sequence
        Returns as string
        """
        splitSeq = []
        for char in sequence:
            charUp = char.upper()
            if charUp == '-': continue
            splitSeq.append(charUp)

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

    def Printer(tRNAList):
        for trna in sorted(tRNAList, key=lambda x:x.sortableHeader):
            trna.printSubstrings()

    def Intersector(seqSubseqDict,seqUniqueSubseqDict, allHeadsList):
        #make seqsubseqDict and compute intersection
        i = 0
        tRNAList = []
        for thisSeq in seqSubseqDict:
            allElseSubseq = []
            for otherSeq in seqSubseqDict:
                if otherSeq != thisSeq:
                    allElseSubseq.append(set(seqSubseqDict[otherSeq]))
            uniqueSubList = set(seqSubseqDict[thisSeq]) - set.union(*allElseSubseq)
            seqUniqueSubseqDict[thisSeq] = uniqueSubList
            tRNAList.append(FindUnique(allHeadsList[i].lstrip(), thisSeq, uniqueSubList))
            i += 1
        return tRNAList

def InputStdin():
    """
    This function takes input, cleans sequences and headers,makes subseq list and prints
    """
    for closedfile in sys.stdin:
        infile = open(str.rstrip(closedfile, "\n"))
        subSequence = []
        seqSubseqDict = OrderedDict()
        seqUniqueSubseqDict = {}
        break
    return infile, seqSubseqDict, seqUniqueSubseqDict


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
        self.header = head #saves the file's header
        self.tRNAseq = tRNAseq
        self.subSequences = fullUniquesubSeq
        self.uniqueEssentialSubstrings = set()
        self.uniqueEssentialSubstrings=self.removeAllSuperStrings(list(fullUniquesubSeq))
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

        for eachUniqueSub in uniqueSubSeqList:
            bigs = set()
            startPos = originalSeq.find(eachUniqueSub)
            stopPos = startPos + len(eachUniqueSub)
            for i in range(0,startPos+1):
                for j in range(stopPos, len(originalSeq)+1):
                    bigs.add(originalSeq[i:j])
            try:
                bigs.remove(str(eachUniqueSub))
            except KeyError:
                pass
            allBigs = allBigs | bigs
        return allBigs

    def removeAllSuperStrings(self, uniqueSubseqList):
        """
        This code takes input subSeqList
        Finds each element in rest of list
        Searches one by one and removes non-minimal substrings
        """
        bigsSet = set(self.get_biggs(uniqueSubseqList, self.tRNAseq))

        uniEssList = set(uniqueSubseqList) - bigsSet
        return uniEssList

    def printSubstrings(self):
        """
        This function prints substrings to stdout
        """
        sortDic = OrderedDict()
        print(self.header)
        print(self.tRNAseq)
        for essentialUnique in self.uniqueEssentialSubstrings:
            pos = self.tRNAseq.find(essentialUnique)
            sortDic[pos] = essentialUnique

        finsortDic = OrderedDict(sorted(sortDic.items()))
        for it in finsortDic:
            print("."*it + finsortDic[it])


def main():
    """
    Main function takes in file
    Parses arguments and runs substring sorter

    """
    infile, seqSubseqDict, seqUniqueSubseqDict = InputStdin()
    seqSubseqDict, allHeadsList, allSubSequencesList = \
            ExternalMethods.get_substringList(infile, seqSubseqDict)

    tRNAList = ExternalMethods.Intersector(seqSubseqDict, seqUniqueSubseqDict, allHeadsList)
    ExternalMethods.Printer(tRNAList)
main()
