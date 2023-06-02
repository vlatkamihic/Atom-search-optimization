from hashlib import new
from telnetlib import GA
import numpy as np
from atom import Atom
import math
from scipy.spatial import distance
import random
import sys

def copyAtoms(atoms):
    copy = []
    for atom in atoms:
        copy.append(atom)
    return copy

def calcSumDis(pDim, kBest, i):
    sum = 0
    for atom in kBest:
        sum += atom.position[i] - pDim
    return sum


def fromRange(value, min, max):
    newValue = value
    if value < min:
        newValue = min
    elif value > max:
        newValue = max
    else:
        newValue = value
    return newValue


class Population():

    atoms = []
    def __init__(self, pop_size, f_min, f_max):
        self.pop_size = pop_size
        self.f_min = f_min
        self.f_max = f_max
        #self.fit_best = sys.float_info.max
        self.fit_best = float(sys.float_info.max)
        self.fit_worst = float(0)
        self.mi_sum = float(0)
        self.g_0 = 1.1
        self.u = 2.4 #1.24
        
    

    def generatePopulation(self, d):
        for i in range(self.pop_size):
            fitValue = sys.float_info.max
            self.atoms.append(Atom(self.f_min, self.f_max, d, fitValue))

    def evaluatePopulation(self, costFunction):

        for i, atom in enumerate(self.atoms):
            pos = atom.position
            #print(pos)
            self.updateAtom(i, costFunction(pos))
            #print("%s. iteration fitValue: %s" % (str(i), str(costFunction(pos))))

    def calculateK(self, cur_iter, num_iter):
        return int(self.pop_size - (self.pop_size - 2) * (math.sqrt(cur_iter / num_iter))) + 1

    def calculateMiSum(self):
        mi_sum = 0
        for atom in self.atoms:
            mi_sum += atom.mi
        self.mi_sum = mi_sum
        return mi_sum
    
    def calculateKbest(self, curr_atom, cur_iter, num_iter):
        K = self.calculateK(cur_iter, num_iter)
        kTemp = copyAtoms(self.atoms)
        kTemp = sorted(kTemp, key=lambda atom: atom.mass, reverse=True)

        thisIndex = kTemp.index(curr_atom)
        n = len(kTemp) - thisIndex
        kBest = copyAtoms(kTemp)[thisIndex:min(K, n)]
   
        
        return kBest

    def calculateDriftFactor(self, cur_iter, num_iter):
        return 0.1 * math.sin((math.pi / 2) * (cur_iter / num_iter))

    def calculateConstraintForce(self, curr_atom, best_atom, beta, dimention):
        
        g_i = beta * (best_atom.position[dimention] - curr_atom.position[dimention])
        
        return g_i



    def calculateHij(self, curr_atom, atom, d, cur_iter, num_iter, r_ij):
        h_min = self.g_0 + self.calculateDriftFactor(cur_iter, num_iter)
        h_max = self.u

        kBest = self.calculateKbest(curr_atom, cur_iter, num_iter)

        mk_average = np.mean([atom.position[d] for atom in kBest])
        length_scale = np.linalg.norm(curr_atom.position[d] - mk_average)

        if(length_scale == 0):
            length_scale = 1

        h_ij = 0
        if((r_ij/length_scale) < h_min):
            h_ij = h_min
        elif((r_ij/length_scale) > h_max):
            h_ij = h_max
        else:
            h_ij = r_ij/length_scale
        #print("hStuff: h_ij "+ str(h_ij) + " h_min " + str(h_min) + " h_max " + str(h_max))
        return h_ij

    

    def calculateForce(self, curr_atom, atom, d, cur_iter, num_iter):
        c = (1 - cur_iter / num_iter) ** 3
        eps = 2 ** (-52)
        radius = np.linalg.norm(curr_atom.position[d] - atom.position[d])
        h_ij = self.calculateHij(curr_atom, atom, d, cur_iter, num_iter, radius)
        
        potential = c * (12 * (-h_ij) ** (-13) - 6 * (-h_ij) ** (-7))
        f_ij = np.random.uniform(0,1) * potential * ((atom.position[d] - curr_atom.position[d]) / (radius + eps))
        
        return f_ij


    def calculateInteractionForce(self, alfa, curr_atom, cur_iter, num_iter, dimention):

        
        kBest = self.calculateKbest(curr_atom, cur_iter, num_iter)
        
        f_temp = 0
        for atom in kBest:
            f_temp += self.calculateForce(curr_atom, atom, dimention, cur_iter, num_iter)
        f_i = alfa * f_temp
        #print("hStuff: curr_atom" + str(cur_iter))
        #print("Force(f_i): "+str(f_i))
        return f_i
            

    def calculateAcceleration(self, alpha, beta, curr_atom, best_atom, cur_iter, num_iter):
        
        acceleration = []
        lagrangian = np.exp(-20.0 * (cur_iter + 1) / num_iter)
        #m_i = curr_atom.calculateMass(self.mi_sum, self.fit_best, self.fit_worst)
        m_i = curr_atom.mass

        for d in range(curr_atom.dimention):
            
            #print("Dimention--"+str(d)+"--")
            f_i = self.calculateInteractionForce(alpha, curr_atom, cur_iter, num_iter, d)
            g_i = self.calculateConstraintForce(curr_atom, best_atom, beta, d)

            a = lagrangian * (f_i + g_i) / m_i
            #print("Acceleration: "+str(a))
            #print("-------------------------------------")
            acceleration.append(a)
        #print("forceStuff: Total force on " + str(atom_num) + " is "  + str(total_force))
        #print("F: %s, G: %s, a: %s" % (str(f_i), str(g_i), str(acceleration)))

        return acceleration

    def calculateVelocityAndPosition(self, alpha, beta, curr_atom, best_atom, cur_iter, num_iter):

        acceleration = self.calculateAcceleration(alpha, beta, curr_atom, best_atom, cur_iter, num_iter)
        #print("Acceleration: "+str(acceleration))
        velocity = []
        position = []
        w = 0.9 * 0.5 * (cur_iter / num_iter)
        c1 = -10 * ((cur_iter / num_iter) ** 2)
        c2 = 1 -(-10 * math.pow(cur_iter / num_iter, 2))
        for d in range(curr_atom.dimention):
            #v = random.uniform(0,1) * curr_atom.velocity[d] + acceleration[d]
            rand = random.uniform(0, 1)
            v = w*rand*curr_atom.velocity[d] + c1*rand*acceleration[d] + c2*rand*(best_atom.position[d]-curr_atom.position[d])
            
            velocity.append(fromRange(v, 0.1*self.f_min, 0.1*self.f_max))
            
            position.append(curr_atom.position[d] + velocity[d])
        
        #curr_atom.updateParams(position, velocity)


        return position, velocity


    def calculateMassForEachAtom(self):
        for atom in self.atoms:
            atom.calculateMass(self.mi_sum, self.fit_best, self.fit_worst)

    def updateAtom(self, index, fitValue):
        self.atoms[index].fitValue = fitValue

    def getFitValues(self):
        gFits = []
        for atom in self.atoms:
            #print("GetFITS fit values: "+str(atom.fitValue))
            gFits.append(atom.fitValue)
        #print("Generation fit values "+ str(gFits))
        return gFits

    def getFitWorst(self):

        for atom in self.atoms:
            if(atom.fitValue > self.fit_worst):
                self.fit_worst = atom.fitValue
        #print("FIT WORST: "+str(self.fit_worst))
        return self.fit_worst
    
    def getFitBest(self):

        for atom in self.atoms:
            if(atom.fitValue < self.fit_best):
                self.fit_best = atom.fitValue
        #print("FIT BEST: "+str(self.fit_best))
        return self.fit_best
    
    def getBestInd(self, d):
        best_ind = Atom(-5, 5, d, sys.float_info.max)
        for atom in self.atoms:
            if(atom.fitValue == self.fit_best):
                best_ind.updateAtom(atom)
        #print("FIT BEST: "+str(self.fit_best))
        return best_ind


    def updatePopulation(self):
        self.calculateMiSum()
        self.calculateMassForEachAtom()
        self.getFitWorst()
        self.getFitBest()

    def __del__(self):
  # body of destructor
        self.atoms.clear()





'''
def calculateLengthScale(self, x_i, K, kBest, j):
        
        x_ij = x_i[j]

        
        sum = 0
        for k in kBest:
            x_ij_temp = k.position[j]
            sum += x_ij_temp
        
        secElement = sum / K
        lengthScale = float(math.sqrt(math.pow(x_ij, 2)+math.pow(secElement, 2)))
        if(lengthScale == 0):
            lengthScale = 1
        #print("Length scale(sigma): "+str(lengthScale))
        return lengthScale
'''
#kBest
'''
        for atom in self.atoms: newlist.append(atom)
        for i in range(len(newlist)):
            for j in range(len(newlist)):
                if(i == j):
                    continue
                if(newlist[i].fitValue < newlist[j].fitValue):
                    temp = newlist[j]
                    newlist[j] = newlist[i]
                    newlist[i] = temp
        #print("K: "+str(K))
        kBest = []
        for i in range(round(K)):
            kBest.append(newlist[i])
      
        #self.fit_best = newlist[0].fitValue
        #self.fit_worst = newlist[self.pop_size - 1].fitValue
        '''