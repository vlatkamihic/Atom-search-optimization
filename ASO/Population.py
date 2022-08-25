from hashlib import new
from telnetlib import GA
import numpy
from atom import Atom
import math
from scipy.spatial import distance
import random
import sys

def calcSumDis(pDim, kBest, i):
    sum = 0
    for atom in kBest:
        sum += atom.position[i] - pDim
    return sum

class Population(object):

    atoms = []
    def __init__(self, pop_size):
        self.pop_size = pop_size
        self.fit_best = float(100)
        self.fit_worst = 0
        self.mass_sum = float(0.1)
        self.g_0 = 1.1
        self.u = 1.24
        

    # using property decorator
    # a getter function
    def __getattribute__(self,name):
        return object.__getattribute__(self, name)
    

    def generatePopulation(self, minP, maxP, minV, maxV, d, depth, multiplierWeight):
        for i in range(self.pop_size):
            fitValue = sys.float_info.max
            self.atoms.append(Atom(minP, maxP, minV, maxV, d, fitValue, depth, multiplierWeight))
        #self.calculateMassSum()


    def calculateMassSum(self):
        mass_sum = 0
        for atom in self.atoms:
            mass_sum += atom.mass
        #print(self.fit_worst)
        #print("Mass sum: %s" % (str(mass_sum)))
        self.mass_sum = mass_sum
        return mass_sum
    
    def calculateKbest(self, curr_atom, cur_iter, num_iter):
        newlist = self.atoms
        for i in range(len(newlist)):
            for j in range(len(newlist)-1):
                if(newlist[i].fitValue > newlist[j].fitValue):
                    temp = newlist[i]
                    newlist[i] = newlist[j]
                    newlist[j] = temp
        K = curr_atom.calculateK(cur_iter, num_iter)
        kBest = []
        for i in range(min(round(K), self.pop_size)):
            kBest.append(newlist[i])
      
        self.fit_best = newlist[0].fitValue
        self.fit_worst = newlist[self.pop_size - 1].fitValue
        return kBest

    def calculateDriftFactor(self, cur_iter, num_iter):
        return 0.1 * math.sin((math.pi / 2) * (cur_iter / num_iter))

    def calculateAcceleration(self, alpha, beta, curr_atom, cur_iter, num_iter):
        
        depth = alpha * (math.pow((1-(cur_iter-1)/num_iter), 3)) * math.exp(-(20*cur_iter)/(num_iter))
        lagrang = beta * math.exp(-(20*cur_iter)/num_iter)
        
        h_min = self.g_0 + self.calculateDriftFactor(cur_iter, num_iter)
        h_max = self.u
        
        kBest = self.calculateKbest(curr_atom, cur_iter, num_iter)
        
        f_i = numpy.zeros(curr_atom.dimention)
        
        pCurAtom = curr_atom.position
        
        for atom in kBest:
            i = 0
            for pDim in atom.position:
                dist = pDim - pCurAtom[i] # r_ij
                sum_dis = calcSumDis(pDim, kBest, i)
                
                sigma = math.sqrt(math.pow(dist, 2) + math.pow((sum_dis/curr_atom.calculateK(cur_iter, num_iter)), 2))
                r_sigma = dist / sigma

                h_ij = 0
                if(r_sigma < h_min):
                    h_ij = h_min
                elif(r_sigma > h_max):
                    h_ij = h_max
                else:
                    h_ij = r_sigma
                
                f_ij = -depth * (2 * math.pow(h_ij, 13) - math.pow(h_ij, 7))
                #print(f_ij)
                f_i[i] += random.uniform(0,1) * f_ij
                i = i + 1

        g_i = numpy.zeros(curr_atom.dimention)
        pBest = kBest[0].position
        acceleration = numpy.zeros(curr_atom.dimention)
        i = 0
        m_i = curr_atom.calculateMass(self.calculateMassSum(), self.fit_best, self.fit_worst)

        for position in pCurAtom:
            g_i[i] = lagrang * (pBest[i] - position)

            acceleration[i] = (g_i[i] / m_i) + ( f_i[i] / m_i )

            i = i + 1
        #print(acceleration)
        return acceleration

    def calculateVelocityAndPosition(self, alpha, beta, curr_atom, cur_iter, num_iter):
        acceleration = self.calculateAcceleration(alpha, beta, curr_atom, cur_iter, num_iter)
        velocity = curr_atom.velocity

        d = curr_atom.dimention
        for i in range(d):
            velocity[i] = random.uniform(0, 1) * velocity[i]  + acceleration[i]
        
        position = curr_atom.position

        for i in range(d):
            position[i] = position[i] + velocity[i]

        curr_atom.updateParams(position, velocity)

    def updateAtom(self, index, fitValue):
        self.atoms[index].fitValue = fitValue

    def getFitValues(self):
        gFits = []
        for atom in self.atoms:
            gFits.append(atom.fitValue)

        return gFits

    def getAtom(self, index):
        list = Population.atoms
        return list[index]