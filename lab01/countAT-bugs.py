#!/usr/bin/env python3
# Name: David Bernick(dbernick)
# Group Members: none
 
class dnaString (str):
    def __new__(self,s):
        return str.__new__(self,s.upper())      
    def length (self):
        return (length(self))
    def getAT (self):
        num_A = self.count(A)
        num_T = self.count("T")
        return ((num_A + num_T)/ self.length() )
 
dna = input("Enter a dna sequence: ")
coolString = dnaString(dna)
 
print ("AT content = {0:0.1f}".format(coolString.getAT()) )
