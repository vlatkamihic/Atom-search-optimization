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

def calculateNorm(x_i, x_j):
    sum_abs = 0
    for i in x_i:
        sum_abs += math.pow(i, 2)
    for j in x_j:
        sum_abs += math.pow(j, 2)

    norm = math.sqrt(sum_abs)
    return norm

def euclideanDistance(x_i, x_j):
    sum_abs = 0
    n = len(x_i)
    for i in range(n):
        sum_abs += math.pow((x_i[i]-x_j[i]), 2)

    dist = math.sqrt(sum_abs)
    return dist

def fromRange(value, min, max):
    newValue = value
    if value < min:
        newValue = min
    elif value > max:
        newValue = max
    else:
        newValue = value
    return newValue


class Population(object):

    atoms = []
    def __init__(self, pop_size, f_min, f_max):
        self.pop_size = pop_size
        self.f_min = f_min
        self.f_max = f_max
        #self.fit_best = sys.float_info.max
        self.fit_best = float(sys.float_info.max)
        self.fit_worst = float(0)
        self.mi_sum = 1
        self.g_0 = 1.1
        self.u = 1.24
        
    

    def generatePopulation(self, d, depth, multiplierWeight):
        for i in range(self.pop_size):
            fitValue = sys.float_info.max
            self.atoms.append(Atom(self.f_min, self.f_max, d, fitValue, depth, multiplierWeight))

    def evaluatePopulation(self, costFunction):

        for i, atom in enumerate(self.atoms):
            pos = atom.position
            #print(pos)
            self.updateAtom(i, costFunction(pos))
            #print("%s. iteration fitValue: %s" % (str(i), str(costFunction(pos))))

    def calculateK(self, cur_iter, num_iter):
        return self.pop_size - (self.pop_size - 2) * (math.sqrt(cur_iter / num_iter))

    def calculateMiSum(self):
        mi_sum = 0
        for atom in self.atoms:
            mi_sum += atom.mi
        self.mi_sum = mi_sum
        return mi_sum
    
    def calculateKbest(self, cur_iter, num_iter):
        newlist = []
        for atom in self.atoms: newlist.append(atom)
        for i in range(len(newlist)):
            for j in range(len(newlist)-1):
                if(newlist[i].fitValue > newlist[j].fitValue):
                    temp = newlist[i]
                    newlist[i] = newlist[j]
                    newlist[j] = temp
        K = self.calculateK(cur_iter, num_iter)
        #print("K: "+str(K))
        kBest = []
        for i in range(round(K)):
            kBest.append(newlist[i])
      
        #self.fit_best = newlist[0].fitValue
        #self.fit_worst = newlist[self.pop_size - 1].fitValue
        return kBest

    def calculateDriftFactor(self, cur_iter, num_iter):
        return 0.1 * math.sin((math.pi / 2) * (cur_iter / num_iter))

    def calculateConstraintForce(self, curr_atom, best_atom, cur_iter, num_iter, beta, dimention):
        
        lagrangian = beta * math.exp(-(20*cur_iter) / num_iter)
        g_i = lagrangian * (best_atom.position[dimention] - curr_atom.position[dimention])
        
        #print("Lagrangian: "+str(lagrangian))
        #print("Constraint: "+str(g_i))
        return g_i

    def calculateDepthFactor(self, alfa, cur_iter, num_iter):
        base = 1 - (cur_iter-1)/(num_iter)
        eta = alfa * (math.pow(base, 3)) * math.exp(-(20*cur_iter)/num_iter)
        return eta

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


    def calculateHij(self, curr_atom, atom, d, cur_iter, num_iter):
        h_min = self.g_0 + self.calculateDriftFactor(cur_iter, num_iter)
        h_max = self.u

        r_ij = euclideanDistance(curr_atom.position, atom.position)
        K = self.calculateK(cur_iter, num_iter)
        kBest = self.calculateKbest(cur_iter, num_iter)

        length_scale = self.calculateLengthScale(curr_atom.position, K, kBest, d)

        h_ij = 0
        if((r_ij/length_scale) < h_min):
            h_ij = h_min
        elif((r_ij/length_scale) > h_max):
            h_ij = h_max
        else:
            h_ij = r_ij/length_scale
        #print("h_ij: "+str(h_ij))
        return h_ij

    

    def calculateForce(self, curr_atom, atom, d, cur_iter, num_iter):
        
        h_ij = self.calculateHij(curr_atom, atom, d, cur_iter, num_iter)
        #m_i = curr_atom.calculateMass(self.mi_sum, self.fit_best, self.fit_worst)
        m_i = curr_atom.mass
        x_i = curr_atom.position
        x_j = atom.position
        euclidean_norm = calculateNorm(x_i, x_j)
        norm = (x_j[d]-x_i[d])/(euclidean_norm)

        f_ij = random.uniform(0,1) * ((2*math.pow(h_ij, 13)-math.pow(h_ij, 7))/m_i) * norm
        #f_ij = random.choice([0,1]) * ((2*math.pow(h_ij, 13)-math.pow(h_ij, 7))/m_i) * norm
        #f_ij = random.choice([0,1]) * ((2*math.pow(h_ij, 13)-math.pow(h_ij, 7))) * norm
        
        return f_ij


    def calculateInteractionForce(self, alfa, curr_atom, i, cur_iter, num_iter, dimention):

        depth = self.calculateDepthFactor(alfa, cur_iter, num_iter)
        kBest = self.calculateKbest(cur_iter, num_iter)
        f_temp = 0
        j = 0
        for atom in kBest:
            f_temp += self.calculateForce(curr_atom, atom, dimention, cur_iter, num_iter)
            j += 1
        f_i = -depth * f_temp
        #print("Force(f_i): "+str(f_i))
        return f_i
            

    def calculateAcceleration(self, alpha, beta, curr_atom, i, best_atom, cur_iter, num_iter):
        
        acceleration = []

        #m_i = curr_atom.calculateMass(self.mi_sum, self.fit_best, self.fit_worst)
        m_i = curr_atom.mass

        for d in range(curr_atom.dimention):
            
            #print("Dimention--"+str(d)+"--")
            f_i = self.calculateInteractionForce(alpha, curr_atom, i, cur_iter, num_iter, d)
            g_i = self.calculateConstraintForce(curr_atom, best_atom, cur_iter, num_iter, beta, d)
            
            constraction = g_i / m_i
            a = (f_i + constraction)
            #print("Acceleration: "+str(a))
            #print("-------------------------------------")
            acceleration.append(fromRange(a, 0.1*self.f_min, 0.1*self.f_max))
        
        #print("F: %s, G: %s, a: %s" % (str(f_i), str(g_i), str(acceleration)))

        return acceleration

    def calculateVelocityAndPosition(self, alpha, beta, curr_atom, index,  best_atom, cur_iter, num_iter):

        acceleration = self.calculateAcceleration(alpha, beta, curr_atom, index, best_atom, cur_iter, num_iter)
        #print("Acceleration: "+str(acceleration))
        velocity = []
        position = []

        for d in range(curr_atom.dimention):
            v = random.uniform(0,1) * curr_atom.velocity[d] + acceleration[d]
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
        print(gFits)
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
    
    def getBestInd(self):
        best_ind = self.atoms[0]
        for atom in self.atoms:
            if(atom.fitValue == self.fit_best):
                best_ind = atom
        #print("FIT BEST: "+str(self.fit_best))
        return best_ind


    def updatePopulation(self):
        self.calculateMiSum()
        self.calculateMassForEachAtom()
        self.getFitWorst()
        self.getFitBest()