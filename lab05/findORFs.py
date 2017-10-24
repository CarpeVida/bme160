#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: Nadja DiMartino (ndimarti)
"""
findORFs class takes in a Fasta file, parses out codons,
counts distances between start and stop codons and returns a list of ORFs
pseudocode:
    loop through seq
        loop through frames (-3,-2,-1,0,1,2)
            parse codon
            pos += 1
            if codon = start
                add it to list
            if codon = stop and no stop codons between it and last start
                save orf start and stop positions
        sort orfs by length
        save sorted list to output file

"""

outFile = None
import sequenceAnalysis
import sys


class CommandLine():
    """
    Handles the command line, usage and help requests.

    CommandLine uses argparse, to implement a standard command line argument parser
    with various argument options,a standard usage and help,
    and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    available within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    """

    def __init__(self, inOpts=None):
        """
        CommandLine constructor.
        Implements argparse to interpret the command line argv string using argparse.
        """

        import argparse
        self.parser = argparse.ArgumentParser(
             description = 'Program prolog - a brief description of what this thing does',
             epilog = 'Program epilog - '
                      'some other stuff you feel compelled to say',
             add_help = True, #default is True
             prefix_chars = '-',
             usage = '%(prog)s [options] -option1[default] <input >output'
             )
        self.parser.add_argument('inFile', action='store', help='input file name')
        self.parser.add_argument('outFile', action='store', help='output file name')
        self.parser.add_argument(
            '-lG', '--longestGene', action='store', nargs='?', const=True,
            default=False, help='longest Gene in an ORF')
        self.parser.add_argument(
            '-mG', '--minGene', type=int, choices= range(0, 2000),
            action = 'store', help='minimum Gene length')
        self.parser.add_argument(
            '-s', '--start', action = 'append', nargs='?', help='start Codon')
        self.parser.add_argument(
            '-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


class OrfFinder():
    """
    OrfFinder iterates through fasta files
    Gives orf positions in a dictionary.
    """
    frameList =[]

    def __init__(self, seq, header, minGene=0, start=["ATG"],
                 longestGene=False, revSeq=False):
        """
        init assigns all the options to self objects and sets defaults
        """
        self.seq = seq
        self.revSeq = revSeq
        self.header = header
        self.minGene = minGene
        if start:
            self.start = start
        else:
            self.start = "ATG"
        if minGene:
            self.minGene = minGene
        self.longestGene = longestGene
        self.frameList = []



    def findOrfs (self, file, reverse=False):
        """
        find Orfs on top strand and return list of Orfs
        remember to handle the dangling start and stop cases
        """

        startList = []      #temp list of start positions
        parsedSeq = None
        if reverse:         #checks to see which strand is being analyzed
            parsedSeq = self.revSeq
        else:
            parsedSeq = self.seq

        for frame in (0,1,2):       #loops through 0 1 and 2, tracks frame
            startFound = False
            for position in range(frame, len(parsedSeq), 3):
                # sets p to equal start position of each codon char
                codon = ''.join(parsedSeq[position: position+3])
                # initalize codon to = three characters
                if codon in self.start:
                    #checks if it's a start codon
                    startList.append(position)
                    #add start position to list
                    startFound = True
                    #marks start as fond before a stop
                if codon in ('TAG','TAA','TGA') and startFound:
                    if reverse:
                        self.saveOrf(len(parsedSeq) - (position+3) + 1,
                                        len(parsedSeq) - (startList[0]+1) + 1,
                                        position-(startList[0])+3,
                                     -(frame+1),file)
                    else:
                        self.saveOrf(startList[0]+1, #nucleotide startposition
                                        position+3,     #nuc stop position
                                        position-(startList[0])+3, #length
                                        frame+1,file)       #frame, file
                    startList = []  #clears start list
                    startFound = False

    def findRevOrfs (self, file):
        """
        find Orfs on the bottom strand and return that list of Orfs
        remember to fixup the orfList so that it refers to top strand
        coordinates and the rev frames
        """
        self.findOrfs(file, reverse=True)

    def clearDict(self):
        """
        clears list of saved frames
        """
        self.frameList = []

    def saveOrf (self, start, stop, length, frame, file):
        """
        Takes in all the orf data and adds it to the output list
        """
        self.frameList = []
        if length >= self.minGene:
            self.frameList.append(dict(frame=frame,start=start,stop=stop,
                                       length=length))

    def writeFramesToFile(self, file):
        """
        Writes formated output to a file argument
        """
        file.write('{} \n'.format(self.header))
        for frame in sorted(self.frameList,
                            key=lambda k: k['length'], reverse=True):
            formatString = '{:+d} {:>5d}..{:>5d} {:>5d}\n'.format\
                (frame['frame'],frame['start'],frame['stop'],frame['length'])
            file.write(formatString)

########################################################################
# Main
# Here is the main program
########################################################################

def main(myCommandLine=None):
    """
    Implements the Usage exception handler that can be raised from anywhere in process.
    """
    # myCommandLine.args.inFile has the input file name
    # myCommandLine.args.outFile has the output file name
    # myCommandLine.args.longestGene is True if only the longest Gene is desired

    if myCommandLine is None:
        myCommandLine = CommandLine([ 'tass2.fa',
                                      'tass2ORFdata-ATG-100.txt',
                                      '--longestGene',
                                      '--start=ATG',
                                      '--minGene=100'])
    else :

        myCommandLine = CommandLine(myCommandLine)
    # myCommandLine.args.startList is a list of start codons
    # myCommandLine.args.minGene is the minimum Gene length to include

    # Parse the options
    myCommandLine.args.inFile
    outFile = myCommandLine.args.outFile
    myCommandLine.args.longestGene
    myCommandLine.args.start
    myCommandLine.args.minGene

    # Parse inFile and initialize it to be read
    orfReader = sequenceAnalysis.FastAreader(myCommandLine.args.inFile)
    open(myCommandLine.args.outFile, 'w').close()
    f = open(outFile, 'a')

    # loop through sequences in the file
    # find ORFs and write them to file
    for head, seq in orfReader.readFasta():
        nucParams = sequenceAnalysis.NucParams('')
        nucParams.addSequence(seq)
        nucParams.initializeObjects()
        bases = list((''.join(nucParams.codonList)))
        reverseBases = getReverseComplement(bases)

        finder = OrfFinder(bases,head,minGene=myCommandLine.args.minGene, longestGene=myCommandLine.args.longestGene, start=myCommandLine.args.start,revSeq=reverseBases)
        finder.findOrfs(f)
        finder.findRevOrfs(f)
        finder.writeFramesToFile(f)

    # myReader = FastAreader(myCommandLine.args.inFile)
    # for head, seq in FastAreader(myCommandLine.args.inFile).readFasta() :
    #    print (head,seq)
        #call find orfs in this seq iterater

    # Clear the file if created previously,
    # and open it


def getReverseComplement(forwardSeq):
    '''
    getReverseStrand returns the reverse complementsequence
    '''
    complement = []
    #complement the sequence
    for base in forwardSeq:
        if base == 'A':
            complement.append('T')
        if base == 'T':
            complement.append('A')
        if base == 'C':
            complement.append('G')
        if base == 'G':
            complement.append('C')
    # return reversed
    return complement[::-1]
if __name__ =="__main__":
    # send all the command line objects to be parsed except for prog name
    main(sys.argv[1:])
