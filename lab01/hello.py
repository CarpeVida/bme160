#!/usr/bin/env python3
# Name: Daniel Schmelter(dschmelt)
# Group Members: None
 
"""
This is my first python program which prints out the string Hello DanielSchmelter.
""" 
 
class announcer (str):    #create class announcer with argument "str"
    def printMe (self):    #within, define function printMe with argument self
        print (self)       #when printMe is called, print out argument in self
 
student = announcer ('Hello World!')
#Creates announcer type of announcer with the
#self is the instance of the object that's doing the thing
student.printMe()
