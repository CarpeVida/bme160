#!/usr/bin/env python3
# Name: Daniel Schmelter(dschmelt)
# Group Members: None
 
"""
This is my first python program which prints out a few sentences of text.
""" 
 
class announcer:    #create class announcer with argument "str"
    """
    Prints Name
    """
    def printName (self):    #within, define function printMe with argument self
        print ('My name is Daniel Schmelter.')
        #when printName is called, print out argument in self
    """
    Prints User
    """
    def printUser(self):
        print('My username is dschmelt.')
    """
    Prints Roles
    """
    def printRole(self):
        print('I am an undergraduate student.')
    """
    Prints Major
    """
    def printMajor(self):
        print('My major is Biomolecular Engineering and Bioinformatics.')
    """
    Prints WhyTake
    """
    def printWhyTake(self):
        print("I'm taking this class because I want to understand how to manipulate large datasets and cross reference reads against template DNA.")
    """
    Prints Interests
    """
    def printInterests(self):
            print("I'm interested in gene therapy, transgenic organisms, understanding aging, and curing diseases.")
    """
    Prints Experience
    """
    def printExperience(self):
        print("I have prior programming experience using Java, HTML, R, and C. (I did't like C either).")
        
#Creates announcer type of announcer with the
student = announcer()
student.printName()
student.printUser()
student.printRole()
student.printMajor()
student.printWhyTake()
student.printInterests()
student.printExperience()
