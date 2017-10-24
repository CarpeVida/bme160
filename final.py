# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 12:37:25 2017

@author: Carmelle
"""


"""
1. Pass args from command line to program
2. 
"""

import csv
#import matplotlib.pyplot as plt
#import numpy as np

class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
     '''
    def __init__(self, inOpts=None) :
            '''
            CommandLine constructor.
            Implements a parser to interpret the command line argv string using argparse.
            '''
            
            import argparse
            self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                                 epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                                 add_help = True, #default is True 
                                                 prefix_chars = '-', 
                                                 usage = '%(prog)s [options] -option1[default] <input >output'
                                                 )
            self.parser.add_argument('-p', '--pValue', type=float, action = 'store', help='desired p value', default=0.05)
            self.parser.add_argument('-i', '--input', action = 'store', help = 'input file')
            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class tsvReader:
    def __init__(self, myCommandline):
        self.pVal = myCommandline.args.pValue
        self.inputFile = myCommandline.args.input
        self.buildDict(self.inputFile)
        
    def buildDict(self, some_input):
        self.psiDictionary = {}
        with open(some_input) as file:
            reader = csv.reader(file, dialect = 'excel-tab')
            next(reader, None)
            for row in reader:
                try:
                    if float(row[510]) < self.pVal:
                        self.psiDictionary[row[2]] = float(row[510])
                except IndexError:
                    pass
        #self.plotResults(self.psiDictionary)
        print (self.psiDictionary)
    def plotResults(self, psi_val):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
        allData = self.psiDictionary.values()
        axes[0].violinplot(allData, showmeans = True, showmedians = True)
        axes[0].set_title('KRAS')
        plt.legend(self.psiDictionary.keys())
        plt.show()
        
        
        
        
        


def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.  

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = tsvReader(myCommandLine)


main()