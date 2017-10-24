#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: none
"""
This program asks for a DNA sequence input and returns string length, numberof As, Cs, Gs, and Ts in the DNA string.
It does not filter input.
"""
class dnaString (str):  #define class of objects dnaString which has all the string functions
    def __new__(self,s):    #creates itself
        return str.__new__(self,s.upper())      
    def length (self):      #defines function within class
        return (len(self))  #returns char length of argument's string input
    def getA (self):        #defines function within class
        num_A = self.count("A") #initiallize num_A object with the counts of "A" in the input argument
        return ((num_A))    #returns num_A
    def getT (self):    #defines function within class
        num_T = self.count("T") #initiallize num_T object with the counts of "T" in the input argument
        return ((num_T))    #returns num_T
    def getC (self):    #defines function within class
        num_C = self.count("C") #initiallize num_C object with the counts of "C" in the input argument
        return ((num_C))    #returns num_C
    def getG (self):    #defines function within class
        num_G = self.count("G") #initiallize num_G object with the counts of "G" in the input argument
        return ((num_G))    #returns num_G
dna = input("Enter a dna sequence: ")   #prints string, takes 1 input
coolString = dnaString(dna) #assigns object coolString to dnaString class
#prints length funtion 
print ("Your sequence is {0} nucotides long with the following breakdown of bases:" .format(coolString.length()))
print ("number As = {0}" .format(coolString.getA()))    #prints function
print ("number Cs = {0}" .format(coolString.getC()))    #prints function
print ("number Gs = {0}" .format(coolString.getG()))    #prints function
print ("number Ts = {0}" .format(coolString.getT()))    #prints function
