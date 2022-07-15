
# Simpatico 

Simpatico - Simulation Package for Polymeric and Molecular Liquids

Copyright 2010 - 2020 The Regents of the University of Minnesota
Distributed under the terms of the GNU General Public License.

## Overview

Simpatico is a C++ package for Molecular Dynamics (MD), Monte Carlo 
(MC), and hybrid MC simulations of classical mechanical models of 
polymeric and molecular liquids. It has thus far been used primarily 
for simulating course-grained models of polymer liquids. 

The simpatico package contains:

- A program for parallel MD simulations, named ddSim.

- Programs for Monte Carlo (MC) and molecular dynamics (MD) 
  simulations on a single processor, named mcSim and mdSim.

The ddSim parallel MD program uses a spatial domain decomposition 
strategy to allow efficient simulation of very large systems. 
The mcSim and mdSim simulation are convenient for simulations on
a personal computer. The single-processor simulations also provide 
convenient ways to submit embarassingly parallel simulations of 
many similar systems with one system per CPU core, which can be
used to improve statistics in MC sampling or to explore a parameter
space. 

All three simulation programs (ddSim, mdSim and mcSim) provide 
very flexible, extensible facilities for on-the-fly output and 
analysis of selected physical variables during a simulation. 
The single-processor mcSim and mdSim programs can also be used 
for postprocess analysis of trajectories, by applying the same 
analysis algorithms to a sequence of configuration snapshots 
read from trajectory file.
   
Simpatico is distributed only in source code form, and so must 
be compiled from source.

## Getting the Source Code

The simpatico source code is maintained in the github repository
                                                                             
   <https://github.com/dmorse/simpatico>.

It may be obtained by using a git version control system client to 
clone the repository. To do so, enter the command:

   git clone --recursive https://github.com/dmorse/simpatico.git

Note the use of the --recursive option to the git clone command: 
This is necessary to clone some git submodules that are maintained 
in separate repositories. 

## Documentation

A recent copy of the web manual for simpatico is available online, at

   <http://dmorse.github.com/simpatico/doc/html/index.html>.

This manual provides both user and developer documentation in an
integrated form.

The web manual for simpatico is prepared using the Doxygen documentation
utility (www.doxygen.org), and can be regenerated by users. Running doxygen 
creates html web pages that are deposited in the directory named doc/html/ 
within a users copy of the simpatico working tree. These web pages are 
created from both comments extracted from the C++ source code and from a 
set of text files with the file extension *.dox. Files containing text
for the main sections of the web manual are located in the directory 
doc/manual/, and are readable in any text editor. 

After cloning the source code, you can use doxygen to generate a local 
copy of the html documentation within the simpatico/doc/html directory 
of the source code tree.  This requires that doxygen be installed on 
your computer and that the executable be in a directory in your PATH. 
The basic instructions (after doxygen is installed) are:

  - cd to the simpatico/ root directory

  - Enter "make html"

This should create many html files in the simpatico/doc/html directory. 
To begin reading the documentation, point a browser at the file 
simpatico/doc/html/index.html, which is the main page of the manual.

## Compiling

The simpatico source code is ansi standard C++, and must be compiled 
from source. The single-processor versions of the mcSim and mdSim MC 
and MD simulation program do not depend on any external libraries. 
Multi-processor programs (ddSim and multi-processor versions of mcSim 
and mdSim) require an MPI library. The build system requires the gnu 
version of the unix make utility (gmake) and a python interpreter.

Complete instructions for compiling various versions of simpatico are 
given in Sec. 2 of the online manual. Short instructions for compiling 
the default single-processor versions of mcSim and mdSim (for the 
impatient) are given below:

- Add the simpatico/bin directory to your linux command search PATH 
  environment variable.

- Add the simpatico/scripts/python directory to your PYTHONPATH 
  environment variable.

- cd to the simpatico/ root directory

- Enter "./setup" from the root directory to run a setup script
  (you only need to do this once, before compiling the first time)

- Enter "make mcMd" from the same directory

The resulting executables, named "mcSim" and "mdSim", will be installed in 
the simpatico bin/ directory. Please see the html documentation mentioned 
above for further instructions for compiling multi-processor programs that 
depend on MPI.

Instructions for enabling or disabling a variety of optional compile-time 
features are also given in the html documentation. 
