#!/usr/bin/env python3
# Name: Daniel Schmelter(dschmelt)
# Group Members: none

"""
This program takes in a string and determines the percentage of A and T characters.
There were three bugs in the original program.
On line 9, the count number of characters function len(self) was replaced by length(self).
On line 11, self.count(A) should have been self.count("A") to count A characters instead of an object named A
On line 18, the print formatting had only one sig fig because it said {0:0.1f} and needed to be {0:0.3f}
"""

"""
defines new class object dnaString which has same function as strings
"""
class dnaString (str):
    """
    defines new funtion which takes string
    """
    def __new__(self,s):
        """
        returns object resulting from
        """
        return str.__new__(self,s.upper())      
    def length (self):
        return (len(self))
        #funtion length is invalid and should have been the function len 
    def getAT (self):
        num_A = self.count("A")
        #argument A was not in quotes and thus referred to an unassigned object
        num_T = self.count("T")
        return ((num_A + num_T)/ self.length() )

 
dna = input("Enter a dna sequence: ") #print string and take input
coolString = dnaString(dna)   #initialize cool string to do the dna function
 
print ("AT content = {0:0.3f}".format(coolString.getAT()) )
#The format number before f was set to 1, resulting in only one sig fig output instead of 3
