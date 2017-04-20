#! /usr/bin/env python
template = open('paramTemplate', 'r').readlines()
from subprocess import Popen, PIPE
data = []
for alpha in ['0.800']:
#for alpha in ['{0:0>5.3f}'.format(0.2 + 0.2*x) for x in xrange(20)]:
    kenergy=[]
    renergy=[]
    selfenergy=[]
    totalenergy=[]

    #for rcutoff in ['{0:0>6.3f}'.format(1.3 + 0.1*x) for x in xrange(40)]:
    for rcutoff in ['3.000','4.000']:

       for kcutoff in ['{0:0>6.3f}'.format(5.0 + x) for x in xrange(10) ]:
       #for kcutoff in ['12.0000']:
            print alpha, rcutoff, kcutoff
            fn = open('param', 'w')
            gn = open('results','w')
            for line in template:
                fn.write(line.replace('ALPHA', alpha).replace('RSPACE', rcutoff).replace('KSPACE', kcutoff))
            fn.close()

            data.append([float(alpha),float(rcutoff),float(kcutoff)])
            p = Popen('mdSim_c -e -p param -c commands', stdout=PIPE, stderr=PIPE, shell=True)
            (out, err), id = p.communicate(), p.pid
            for line in out.split('\n'):
                if line[17:20] == '@ES':
                    k = float(line.split()[4])
                    kenergy.append(k) 
                if line[:6] == 'rSpace':
                    r = float(line.split()[3])
                    renergy.append(r)
                if line[:5] == 'self ' :
                    s = float(line.split()[3])
                    selfenergy.append(s)
                if line[:5] == 'total' : totalenergy.append(k+r-s)
            kspaceenergy=sum(kenergy)/len(kenergy)
            rspaceenergy=sum(renergy)/len(renergy)
            self=sum(selfenergy)/len(selfenergy)
            total=sum(totalenergy)/len(totalenergy)
            data[-1].append(kspaceenergy)
            data[-1].append(rspaceenergy)
            data[-1].append(self)
            data[-1].append(total)
                
print
#gn.write("alpha".format(10)+"rcutoff".format(10)+"kcutoff".format(10)+"kenergy".format(10)+
#        "renergy".format(10)+"total".format(10)+"\n")
gn.write(reduce(lambda x, y : x + y, ['{0:10}'.format(x) for x in ['alpha', 'rcutoff', 'kcutoff', 'kenergy']]) + '\n')

for x in data:
    for each in x:
        #gn.write(str(each).format(10))
        gn.write('{0:15.8f}'.format(each))
    gn.write('\n')
    print x
gn.close()
