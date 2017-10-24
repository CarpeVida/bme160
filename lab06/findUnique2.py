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

import sys
allSubSequencesList = []


class SubstringGenerator:
    """
    This class is an external substring generator
    """
    def get_substrings(self, string):
        """
        This function splits and generates substrings
        """
        length = len(string)
        for firstPos in range(length):
            for endPos in range(firstPos+1, length +1):
                yield(string[firstPos: endPos])


class FindUnique:
    """
    This class takes in bases, head, and list of subsequences
    It parses through the sequence data
    Makes allHeadsList of all heads,
    Makes seqSubseqDict of all sequences and their subsequences
    Makes allElseSubseq list of all but seq
    Find unique by doing subsequence - allElseSubseq
    Save seqUniqueDict for each sequence
    Remove unneccesary superstrings
    Print with header, seq, and unique minimums with dot spacers
    """
    def __init__(self, head, bases, fullsubSequences):
        """
        Gets things initialized
        """
        #self.orderedSubSequences = {}  #iterater
        self.header = head #saves the file's header
        self.tRNAseq = bases
        self.subSequences = fullsubSequences
        self.uniqueNonSuperSubstrings =\
            self.removeAllSuperStrings(list(fullsubSequences))
        #print(self.uniqueNonSuperSubstrings)
        #self.sortableHeader = head.split('|')[0]

    def get_substrings(self, string):
        """
        This is the function that splits and generates substrings
        """
        length = len(string)
        for firstPos in range(length):
            for endPos in range(firstPos+1, length +1):
                yield(string[firstPos: endPos])

    def removeAllSuperStrings(self, subSeqList):
        """
        This code takes input subSeqList
        Finds each element in rest of list
        Searches one by one and removes non-minimal substrings
        """
        all_strings = subSeqList[:]
        index = 0
        for eachSubstring in subSeqList:
            for iterator in range(1, len(subSeqList)):
                oneSubseq = subSeqList[iterator]
                if eachSubstring != oneSubseq and eachSubstring in oneSubseq:
                    if oneSubseq in all_strings:
                        all_strings.remove(oneSubseq)
                        index +=1
        #print("Index: " + str(index))
        return all_strings

    def printSubstrings(self):
        """
        This function prints substrings to stdout
        """
        print(self.header)
        print(self.tRNAseq)
        #print(self.uniqueNonSuperSubstrings)
        #print(len(self.subSequences))
        for subseq in self.subSequences:
            if subseq in self.uniqueNonSuperSubstrings:
                filler =''
                for i in range(0, len(subseq)):
                    filler += ('.')
                #print('{} {} \n'.format(fmtString, subseq[0]))

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
from sequenceAnalysis import NucParams, FastAreader
from collections import defaultdict

def main():
    """
    Main function takes in file
    Parses arguments and runs substring sorter

    """

    #checks to see if file is entered, if not, print error and exit
    if len(sys.argv) != 2:
        print("Please enter a file after the program")
        sys.exit()
    #if file is entered, initialize it to infile
    if len(sys.argv) == 2:
        infile = open(sys.argv[1])


        #print(orfReader.readFastaStdIn(infile))
        #create lists for future storage
        allHeadsList= []
        subSequence = []
        seqSubseqDict = {}
        seqUniqueSubseqDict = {}

        substringMaker = SubstringGenerator()
        nucParam = NucParams('')
        orfReader = FastAreader()

        #takes infile, saves cleanseq and subseq to dictionary
        #makes list of all subseqs
        #makes list of all heads
        for head, seq in orfReader.readFastaStdIn(infile):
            cleanSeq = str(nucParam.stripSequence(seq))
            subSequence = list(substringMaker.get_substrings(cleanSeq))
            seqSubseqDict.update({cleanSeq:subSequence})
            allHeadsList.append(head)
            allSubSequencesList.append((subSequence))

        #make all else subseq list
        #subseq - allElse = uniques
        #save them to seqUnique dic
        i = 0
        #print(seqSubseqDict)
        for seq in seqSubseqDict:
            allElseSubseq = []
            for otherSeq in seqSubseqDict:
                if otherSeq != seq:
                    allElseSubseq.append((seqSubseqDict[otherSeq]))

            tempListSubseq= seqSubseqDict[seq]
            print(type(tempListSubseq))
            setTempListSubseq = set(tempListSubseq)
            setAllElse = set()
            for item in allElseSubseq:
                setAllElse.add(item)
            uniqueSubList = setTempListSubseq - setAllElse
            seqUniqueSubseqDict = {seq: uniqueSubList}
            FindUnique(seq, allHeadsList[i], seqSubseqDict[seq])
            i += 1
            FindUnique(allHeadsList[0], seq, seqSubseqDict[seq])

        tRNAList = []

        headLen = len(allHeadsList)
        baseIndex = 0
        thisHeadPos=0
        thisBasesPos=1
        #print("SeqSub len: " + str(len((seqSubseqDict))))
        for thisBase in seqSubseqDict:
            #tempSeq = allSubSequencesList.copy()
            #deletes the current sequence
            #del tempSeq[baseIndex]

            #tempUnique = list()

            #joins all the sequences into tempUnique
            #tempUnique.extend((tempSeq))
            tRNAList.append(FindUnique(allHeadsList[baseIndex].lstrip(), thisBase, seqSubseqDict[thisBase]))
            baseIndex += 1
            if thisBasesPos+1 < headLen and thisHeadPos < headLen:
                        thisHeadPos +=1
                        thisBasesPos +=1
        #print(len(tRNAList))
        for trna in tRNAList:
            trna.printSubstrings()

main()
