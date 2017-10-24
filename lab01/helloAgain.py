#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: none
 
 
"""
This is a demonstration program that shows off some basic interactive features
of a python program. The program asks the user to enter a name, builds an object
with that name, and then has the new object introduce itself.
 
input: a string of arbitrary length, which is used to name the new person object
output: greeting printed to screen
"""
 
class person:
#define new class object person
    def __init__(self,name,pet):
    #defintes initial function of person which takes a name and a pet argument
        self.myName = name
        #defintes attribute myName in person as name object
        self.myPet = pet
        #defintes attribute myPet in person as pet object
    def introduce (self):
    #creates function within person called introduce, which prints below.
        print ("Hi there, I am {0} and I like {1}s!".format(self.myName, self.myPet))
        #prints out text with variables {0} and {1} being filled by format parameters 
 
name = input( "What is my name? : " )
#print out text, take input, assign input to name
pet = input( "What is my favorite kind of pet? : ")
#same as above
newPerson = person (name, pet)
#create a new person with arguments name and pet
newPerson.introduce()
#run person's introduce funtion
