from turtle import position
import numpy
import math
import sys


class Atom(object):
    
    def __init__(self, minP, maxP, minV, maxV, d, fitValue, depth, multiplierWeight):

        self.dimention = d
        self.depthWeight = depth
        self.multiplierWeight = multiplierWeight
        self.fitValue = fitValue
        self.mass = 0.0000001
        self.position = list(numpy.random.uniform(minP, maxP, size=(d)))
        #print(self.position)
        #self.velocity = list(numpy.random.uniform(minV, maxV, size=(d)))
        self.velocity = [0] * d

    def getPosition(self):
        return self.position
    #def __repr__(self):
     #   return str((self.dimention, self.depthWeight, self.multiplierWeight, self.fitValue, self.mass, self.position, self.velocity))

    def updateAtom(self, new):
        self.fitValue = new.fitValue
        self.mass = new.mass
        print("Old position: %s, new position: %s" % (str(self.fitValue), str(new.fitValue)))
        for i in range(self.dimention):

            self.position[i] = new.position[i]
            self.velocity[i] = new.velocity[i]
    
    def updateParams(self, position, velocity):
        for i in range(self.dimention):

            self.position[i] = position[i]
            self.velocity[i] = velocity[i]

    def updateFitValue(self, fitValue):
        self.fitValue = fitValue
    
    def calculateMi(self, fitWorse, fitBest):
        #print("FitValue %s, fitBest %s, fitWorst%s" % (str(self.fitValue), str(fitBest), str(fitWorse)))
        if(fitBest > 99999 or fitWorse == fitBest):
            mi = sys.float_info.max
        else:
            try:
                mi = math.exp(-(self.fitValue-fitBest)/(fitWorse-fitBest))
            except OverflowError:
                mi = sys.float_info.max
        #print("MI: %s" % (str(mi)))
        return mi
    
    def calculateMass(self, mass_sum, fitBest, fitWorse):
        mi = self.calculateMi(fitWorse, fitBest)
        #print("MI: %s" % (str(mi)))
        self.mass = mi / mass_sum
        return self.mass

    def calculateK(self, cur_iter, num_iter):
        return self.dimention - (self.dimention - 2) * (math.sqrt(cur_iter / num_iter))

