#!/usr/bin/env python3
# Name: Daniel Schmelter (dschmelt)
# Group Members: Britney Wick

import random

class ChorismateSimulation:

        def __init__(self, quantityAroH, quantityTrpE, quantityPabA):
            self.chorismateMolecules = 100
            self.enzyme =["quantityEnzyme", "Km", "AmountProduct"]
            self.AroH = [quantityAroH, 85, 0]  #Prephenate, goes to Y and F
            self.TrpE = [quantityTrpE, 7, 0]   #anthranilate goes to Tryptophan
            self.PabA = [quantityPabA, 30, 0]  #folate or acetaminophen
            self.fourABH = [1, 20, 0]
            self.nhoA = [1, 0.5, 0]

        def simulation(self):

            prephenate = 0
            anthranilate = 0
            paraAminoBenzoate = 0
            AroHrate = 1/self.AroH[1]
            TrpErate = 1/self.TrpE[1]
            PabArate = 1/self.PabA[1]
            totalProductCount = 0

            while totalProductCount < self.chorismateMolecules:
                num = random.random()
                if num < AroHrate:
                    if (totalProductCount + AroHQuantity) <= self.chorismateMolecules:
                        self.AroH[2] += AroHQuantity
                        totalProductCount += AroHQuantity
                    else:
                        self.AroH[2] += int(random.random() * AroHQuantity)
                if AroHrate < num < (AroHrate + TrpErate):
                    if (totalProductCount + TrpEQuantity) <= self.chorismateMolecules:
                        self.TrpE[2] += TrpEQuantity
                        totalProductCount += TrpEQuantity
                    else:
                        self.TrpE[2] += int(random.random() * TrpEQuantity)
                if (AroHrate + TrpErate) < num < (AroHrate + TrpErate + PabArate):
                    if (totalProductCount + PabAQuantity) <= self.chorismateMolecules:
                        self.PabA[2] += PabAQuantity
                        totalProductCount += PabAQuantity
                    else:
                        self.PabA[2] += int(random.random() * PabAQuantity)

        def results (self):
            print("\t\t\t\t\t\t\t\t\t\t\t" + str(self.enzyme))
            print ("The quantity of your enzymes are the following: {} , {} , {}".format(self.AroH, self.TrpE, self.PabA))


AroHQuantity = int(input('Please input the quantity of your AroH: '))
TrpEQuantity = int(input('Please input the quantity of your TrpE: '))
PabAQuantity = int(input('Please input the quantity of your PabA: '))

newsimulation = ChorismateSimulation(AroHQuantity, TrpEQuantity, PabAQuantity)
print("Starting quantities:")
newsimulation.results()
newsimulation.simulation()
print("You ran the simulation.")
newsimulation.results()
