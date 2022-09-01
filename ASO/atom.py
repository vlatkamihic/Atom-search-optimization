from turtle import position
import numpy
import math
import sys


class Atom(object):
    
    def __init__(self, minP, maxP, d, fitValue):

        self.dimention = d
        self.fitValue = fitValue
        self.mass = float(0.0000001)
        self.mi = float(0.001)
        self.position = list(numpy.random.uniform(minP, maxP, size=(d)))
        self.velocity = [0] * d

    def updateAtom(self, new):
        self.fitValue = new.fitValue
        self.mass = new.mass
        self.mi = new.mi
        #print("Old position: %s, new position: %s" % (str(self.fitValue), str(new.fitValue)))
        for i in range(self.dimention):

            self.position[i] = new.position[i]
            self.velocity[i] = new.velocity[i]
    
    def updateParams(self, position, velocity):
        #print("Old position: %s, new position: %s, fit: %s" % (str(self.position), str(position), str(self.fitValue)))
        for i in range(self.dimention):

            self.position[i] = position[i]
            self.velocity[i] = velocity[i]
    
    def calculateMi(self, fitWorse, fitBest):
        #print("FitValue %s, fitBest %s, fitWorst%s" % (str(self.fitValue), str(fitBest), str(fitWorse)))
        
        try:
            if(fitWorse == fitBest):
                mi = math.exp(-(self.fitValue-fitBest)/(1))
            #mi = 100
            else:
                mi = math.exp(-(self.fitValue-fitBest)/(fitWorse-fitBest))
        except OverflowError:
            mi = sys.float_info.max
                #mi = 100
        #print("MI: %s" % (str(mi)))
        self.mi = mi
        return mi
    
    def calculateMass(self, mi_sum, fitBest, fitWorse):
        mi = self.calculateMi(fitWorse, fitBest)
        #print("Best: %s Worst: %s MI: %s MI_SUM: %s" % (str(fitBest), str(fitWorse), str(mi), str(mi_sum)))
        self.mass = mi / mi_sum
        return self.mass
