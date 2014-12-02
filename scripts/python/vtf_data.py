#!/usr/bin/python
# vtf_data.py

# read data from input = conf
f = open('config','r')
lines = f.readlines()
f.close()
		
# open output file = chain.vtf		
fo = open('chain.vtf','w')

chlength=16
atomname=8*['A']+8*['B']

lengths=lines[2].split()
Lx = lengths[1]
Ly = Lx
Lz = Lx
'''Lz = lengths[3]
Lz = lengths[3]'''
nMolecule = lines[7].split()
nrchains = int(nMolecule[1])
nrmon = nrchains*chlength

for chain in range(0,nrchains):
  mon0=chain*chlength
  nameid=0
  for mon in range(mon0,mon0+chlength-1):
        s = 'atom '+str(mon)+' radius 0.8 name '+atomname[nameid]+' res '+str(chain)+' resid '+str(chlength)+'\n'	
        nameid+=1
        fo.write(s)
	s = 'bond '+str(mon)+':'+str(mon+1)+'\n'
        fo.write(s)
  mon+=1	
  s = 'atom '+str(mon)+' radius 0.8 name '+atomname[nameid]+' res '+str(chain)+' resid '+str(chlength)+'\n'	
  fo.write(s)


fo.write('\n')
fo.write('timestep\n')	
fo.write('pbc '+Lx+' '+Ly+' '+Lz+'\n')
begin = 10
for chain in range(0,nrchains):
  for count in range(begin,begin+chlength):
     coords=lines[count].split()
     x=coords[0]
     y=coords[1]
     z=coords[2]
     s = x+' '+y+' '+z+'\n'
     fo.write(s)
  begin=begin+chlength+2
 
fo.close()

