#!/usr/bin/env python
import sys 
import math
from confile import *

# ---------------
def write_molR():
    '''
    PURPOSE
      Read and parse config file. Output molecule configurations in pbc scheme.
    COMMENT
      File is opened and closed within body of method
    '''
    if (len(sys.argv) > 3):
        x = ConFile(sys.argv[2])
        x.build_pbcR()
        for i in range(x.nSpecies):
            for j in range(len(x.R[i])):
                f = open(sys.argv[3] + '_mol_' + str(i+1) + '_' + str(j+1) + '.dat', 'w')
                x.write_pbcR(f, i, j)
                f.close()
    else:
        print "config file or output prefix not provided"


# ---------------
def deform_box():
    '''
    PURPOSE
      Read and parse config file. Affinely scale atom positions and box
      dimensions.
    COMMENT
      File is opened and closed within body of method
    '''
    # set rescaling factors
    s = []
    if len(sys.argv) > 2:
        if sys.argv[2] == 'triaxis':
            if (len(sys.argv) > 5):
                for i in range(3):
                    s.append(float(sys.argv[3+i]))
            else:
                sys.stdout.write('triaxis parameters missing\n')
                return
        elif sys.argv[2] == 'uniaxis':
            if (len(sys.argv) > 5):
                s2 = 1.0
                for i in range(3):
                    s.append(float(sys.argv[3+i]))
                    if s[i] > 0.00001:
                        s2 *= s[i]
                s2 = math.sqrt(1.0/s2)
                for i in range(3):
                    if s[i] < 0.00001:
                        s[i] = s2
                    print 's[' + str(i) + '] = ' + str(s[i])
            else:
                sys.stdout.write('uniaxis parameters missing\n')
                return
        else:
            sys.stdout.write('unknow deformation parameters\n')
            return
    else:
        sys.stdout.write('deform option or config file not provided\n')
        return

    # load config file, rescale positions and boundary
    if len(sys.argv) > 6:
        x = ConFile(sys.argv[6])
        for i in range(3):
            x.boundary[i] *= s[i]
        x.boundary_type = 'orthorhombic' 
        for i in range(x.nSpecies):
            for j in range(len(x.R[i])):
                for k in range(len(x.R[i][j])):
                    for d in range(3):
                        x.R[i][j][k][d] *= s[d]
    else:
        sys.stdout.write('no config file provided\n')
        return

    # echo rescaled positions
    if len(sys.argv) > 7:
        f = open(sys.argv[7], 'w')
        x.write(f)
        f.close()
        sys.stdout.write('new config written to ' + sys.argv[7] + '\n')
    else:
        x.write(sys.stdout)


# -----------------
def parse_option():
    '''
    PURPOSE
      Parse options and invoke appropriate subroutines.
    '''
    if (len(sys.argv) > 1):
        option = sys.argv[1]
        if option == 'molR':
            write_molR()
        elif option == 'deform':
            deform_box()
        else:
            sys.stdout.write("option not implemented\n")
    else:
        sys.stdout.write("nothing done\n")


# ------------------------
if __name__ == '__main__':
   '''
   USAGE
     Implemented commands:
     1) molR: output the periodic configurations of each molecule
        - spbc molR conf_filename output_prefix
     2) deform the molecular conformation and output to directed locations
        - spbc deform uniaxis s1 s2 s3 conf_filename output_filename
          only one number in s1 and s2 are nonzero the other two are calculated
          by the constraint of volume conservation
        - spbc deform triaxis s1 s2 s3 conf_filename output_filename
   '''
   parse_option()
