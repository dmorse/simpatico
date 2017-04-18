#/usr/bin/env python

import random

class shuffle():

    def __init__(self, filename="config",residue=200):
        self.filename = filename + "_" + str(residue)
        self.residue  = residue

    def shuffler(self):

        old_mol = 799 
        new_mol = []

        with open("config") as f:
            with open(self.filename,"w") as g:

                while 1:
                    l = f.readline()
                    l_list = list(l)

                    try:
                        if l_list[0] == 'm':
                            l = f.readline()
                            new_mol.append(l)
                    except:
                        break
                mol = 0
                while mol < self.residue/2:
                    molid = random.randint(0,old_mol)
                    g.write("molecule" + "          " + str(mol) + "\n")
                    g.write(new_mol[molid] + "\n")
                    mol += 1
                    new_mol.pop(molid)
                    old_mol -= 1

                mol = 0
                while mol < self.residue/2:
                    molid = random.randint(0,old_mol)
                    g.write("molecule" + "          " + str(mol) + "\n")
                    g.write(new_mol[molid] + "\n")
                    mol += 1
                    new_mol.pop(molid)
                    old_mol -= 1




s = shuffle(residue=230)
s.shuffler()
