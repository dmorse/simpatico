#!/usr/bin/env python

from io import *
import string, math

class ConFile(object):
    '''
    An ConFile object contains the data in a *config output
    file produced by C++ program simpatico or molmcd. 

    Derived from Dave's version for PSCF output file.
    '''
   
    def __init__(self,filename):
        '''
        PURPOSE
          Read and parse config file filename. 
          Create an ConFile object.
        ARGUMENT
          filename - string 
        COMMENT
          File is opened and closed within body of method
        '''
        self.file = open(filename,'r')
        self.io   = IO()
        self.numberSpecies = 1
        self.nMolecule = []

        # Dictionary flags indicates which sections are present
        # Keys are flag string values ('BOUNDARY' etc.), values are all 1
        self.flags = {}

        # Read configuration file
        next = 1
        while next:
            line = self.file.readline()
            if len(line) <= 0:  # eof
                next = 0
                continue

            flag = line.strip()
            if flag == '':
                continue

            # Set key in self.flags dictionary
            self.flags[flag] = 1

            if flag == 'BOUNDARY':
                self.input_boundary()
            elif flag == 'MOLECULES':
                self.input_molecules()
            else:
                next = 0

        self.file.close()
	self.file = None

        # Make sorted list of attribute names
        self.att_names = self.__dict__.keys()
        self.att_names.sort()


    def write(self, file, major=1, minor=0):
        self.file = file

        if self.flags.has_key('BOUNDARY'):
            file.write("%-15s\n\n" % 'BOUNDARY')
            self.output_boundary()

        if self.flags.has_key('MOLECULES'):
            file.write("\n%-15s\n" % 'MOLECULES')
            self.output_molecules()

	file.close()
	self.file = None


    def input_boundary(self):
        next = 1
        while next:
            line = self.file.readline().strip()
            if line == '':
                continue

            line = line.split()
            self.boundary_type = line[0]
            self.boundary = []
            if self.boundary_type == 'orthorhombic':
                for i in range(3):
                    self.boundary.append(float(line[i+1]))
            elif self.boundary_type == 'cubic':
                length = float(line[1])
                for i in range(3):
                    self.boundary.append(length)
            next = 0


    def input_molecules(self):
        self.R = []
        self.nSpecies = 0
        next = 1
        while next:
            line = self.file.readline().strip()
            if line == '':
                continue

            line = line.split()
            if line[0] == 'species':
                self.nSpecies += 1
                Rspec = self.input_species()
                self.R.append(Rspec)

            if self.nSpecies == self.numberSpecies:
                next = 0
            
    def input_species(self):
        Rspec = []
        next = 1
        while next:
            line = self.file.readline().strip()
            if line == '':
                continue

            line = line.split()
            if line[0] == 'nMolecule':
                nMol = int(line[1])

            for i in range(nMol):
                molNext = 1
                while molNext:
                    moline = self.file.readline().strip()
                    if moline == '':
                        continue
                    else:
                        molNext = 0
                Rmol = self.input_position()
                Rspec.append(Rmol)

            next = 0

        return Rspec


    def input_position(self):
        Rmol = []
        next = 1
        while next:
            line = self.file.readline()
            if len(line) == 0:
                next = 0
                continue
            line = line.strip()
            if line == '':
                next = 0
            else:
                line = line.split()
                R = [ float(x) for x in line ]
                Rmol.append(R)
        return Rmol


    def output_boundary(self):
        if self.boundary_type == 'orthorhombic':
            self.output_vec('real', 'orthorhombic', 3)
	    self.io.output_vec(self.file, 'real', self.boundary, 3, 'orthorhombic', 'R', 'B')
        elif self.boundary_type == 'cubic':
            self.io.output_var(self.file, 'real', self.boundary[0], 'cubic', 'B')


    def output_molecules(self):
        self.file.write("\n")

        for i in range(self.nSpecies):
            self.io.output_var(self.file, 'int', i, 'species', 'B')
            self.io.output_var(self.file, 'int', len(self.R[i]), 'nMolecule', 'B')
            for j in range(len(self.R[i])):
                self.file.write("\n")
      	        self.io.output_var(self.file, 'int', j, 'molecule', 'B')
                for k in range(len(self.R[i][j])):
	            self.io.output_vec(self.file, 'real', self.R[i][j][k], 3, 'R', f='N')


    # Output pbc configuration for one molecule 
    def write_pbcR(self, mfile, i, j):
        r = self.pbcR[i][j]
        for k in range(len(r)):
            self.io.output_vec(mfile, 'real', r[k], 3, 'R', f='N')


    # Configuration manipulation methods
    def build_pbcR(self):
        self.pbcR = []
        pbcR = self.pbcR
        for i in range(self.nSpecies):
            pbcR.append([])
            for j in range(len(self.R[i])):
                pbcR[i].append([])
                r = self.R[i][j]
                for k in range(len(r)):
                    pbcR[i][j].append([])
                    if k == 0:
                        for n in range(3):
                            pbcR[i][j][k].append(r[k][n])
                    else:
                        dr = []
                        for n in range(3):
                            dr.append(r[k][n])
                            dr[n] -= r[k-1][n]
                            if dr[n] < -0.5*self.boundary[n]:
                                dr[n] += self.boundary[n]
                            if dr[n] > 0.5*self.boundary[n]:
                                dr[n] -= self.boundary[n]
                            pbcR[i][j][k].append(pbcR[i][j][k-1][n])
                            pbcR[i][j][k][n] += dr[n]


    # Input methods (wrapper for self.io.input_... methods of IO)
    def input_var(self, type, comment = None, f='A'):
        return self.io.input_var(self.file, type, comment, f)

    def input_vec(self, type, n=None, comment=None, s='R',f='A'):
        return self.io.input_vec(self.file, type, n, comment, s, f)

    def input_mat(self, type, m, n=None, comment=None, s='L', f='A'):
        return self.io.input_mat(self.file, type, m, n, comment, s, f)

    # Output methods (output by name)
    def output_var(self, type, name, f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_var(self.file, type, data, name, f)

    def output_vec(self, type, name, n=None, s='R', f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_vec(self.file, type, data, n, name, s, f)

    def output_mat(self, type, name, m, n=None, s='L', f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_mat(self.file, type, data, m, n, name, s, f)

    def find_limit(self, dim):
        limit = [self.boundary[dim], 0.0]

        for i in range(self.nSpecies):
            for j in range(len(self.R[i])):
                for k in range(len(self.R[i][j])):
                    if self.R[i][j][k][dim] > limit[1]:
                        limit[1] = self.R[i][j][k][dim]
                    if self.R[i][j][k][dim] < limit[0]:
                        limit[0] = self.R[i][j][k][dim]
        return limit

    def __getitem__(self,key):
        return self.__dict__[key]

    def __str__(self):
        s = []
        for key in self.att_names:
            s.append( key +  ' : ' + str( self[key] ) )
        return string.join(s,'\n')

    def eval(self,expr1):
        '''
	Returns the value of a python expression calculated
	by using the key names of attributes of an Outfile
	as variable names. 
	'''
        for key in self.__dict__.keys():
            exec( key + '= self.' + key )
        return eval(expr1)
         
