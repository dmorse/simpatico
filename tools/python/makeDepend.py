import os
import os.path
from string import *
from TextWrapper import *

def editDepend(infile, outfile, extraDependencies=''):
   file  = open(infile, 'r')
   lines = file.readlines()
   file.close()

   # Extract target from first line
   groups   = lines[0].split(":")
   target   = groups[0]
   lines[0] = groups[1]

   # Replace target by absolute path
   # (The gcc -M option yields only the target base name)
   dirname = os.path.dirname(infile)
   if (dirname != '.' and dirname != ''):
      dirname = os.path.abspath(dirname)
   else:
      dirname = os.getcwd()
   target = os.path.normpath(os.path.join(dirname, target))

   # Replace dependencies by absolute paths
   text = TextWrapper('\\\n', 8)
   for i in range(len(lines)):
       line = lines[i]
       if line[-1] == '\n':
           line = line[:-1]
       if line[-1] == '\\':
           line = line[:-1]
       lines[i] = line
       if i == 0:
           text.append(target + ': ')
       deps = line.split()
       for dep in deps:
           path = os.path.abspath(dep)
           text.append(path)

   # Process extraDependencies (if any)
   if extraDependencies:
       deps = extraDependencies.split()
       for dep in deps:
          path = os.path.abspath(dep)
          text.append(path)
      
   

   file  = open(outfile, 'w')
   file.write(str(text))
   file.close()


def makeDepend(cpath, options, extraDependencies=''):

   # Separate source file name into components
   dir    = os.path.dirname(cpath)
   base   = os.path.basename(cpath)
   groups = base.split('.')
   base   = groups[0]
   if dir != '':
      base  = dir + os.sep + base

   # Create command
   command  = 'g++ '
   command += options
   command += ' -MM -MF '
   command += base + '.p '
   command += cpath
   #print command
   os.system(command)

   #Edit dependency file
   editDepend(base + '.p', base + '.d', extraDependencies)
   os.remove(base + '.p')
