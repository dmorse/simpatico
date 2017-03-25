#/usr/bin/python

from math import *
import random

class Configuration():

    def __init__(self, a=10.455159, b=10.455159, c=10.455159, nMolecules=80,atomradius=2, seed=24234444):
        self.a = a
        self.b = b
        self.c = c
        self.nMolecules=nMolecules
        self.molecules=[]
        self.nmolinstant=0
        self.atomradius = atomradius;
        self.seed = seed
        random.seed(self.seed)

    def createAtom(self):
        xo = random.random() * self.a
        yo = random.random() * self.b
        zo = random.random() * self.c
        origin = (xo, yo, zo)
        return origin

    def atoms(self):

        availability = 0
        atom = self.createAtom()
        self.molecules.append(atom)

        while not availability:
            atom = self.createAtom()
            availability = self.checkAvailability(atom, self.molecules) 

            if availability:
                self.molecules.append(atom)
                availability = 0

            if len(self.molecules) == self.nMolecules:
                break

    def checkAvailability(self,atom,molecules):
        for each in molecules:
            dx = atom[0] - each[0]
            dy = atom[1] - each[1]
            dz = atom[2] - each[2]
            rsq = dx*dx + dy*dy + dz*dz;
            if (rsq < self.atomradius*self.atomradius):
                return 0
        return 1;

    def writeOutput(self, filename):
        id = 0
        with open(filename, "w") as f:

            f.write("BOUNDARY" + "\n")
            f.write("\n")
            f.write("orthorhombic" + "      " + format(self.a, '.10e') + "        " + format(self.b, '.10e')+"        " + format(self.c, '.10e') + "\n")
            f.write("\n")
            f.write("MOLECULES"+'\n')
            f.write('\n')

            f.write('species'+'     '+ '0'+'\n')
            f.write('nMolecule' + '      ' +'40'+'\n'+'\n') 
            for each in self.molecules[:40]:
                f.write("molecule" + '        ' + str(id) + '\n')
                f.write('    ' + format(each[0], '.10e') + '    ' + format(each[1], '.10e') + '    ' + format(each[2], '.10e')+'\n'+'\n')
                id += 1

            id = 0
            f.write('species'+'        '+ '1'+'\n')
            f.write('nMolecule' + '        ' +'40'+'\n'+'\n') 
            for each in self.molecules[40:80]:
                f.write("molecule" + '        ' + str(id) + '\n')
                f.write('    ' + format(each[0], '.10e') + '    ' + format(each[1], '.10e') + '    ' + format(each[2], '.10e')+'\n'+'\n')
                id += 1

molconfig = Configuration()
molconfig.atoms()
molconfig.writeOutput('config')

 


