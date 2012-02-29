#ifndef DDSIM_CPP
#define DDSIM_CPP

#include <ddMd/system/System.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/interaction/Interaction.h>
#include <util/random/Random.h>
#include <util/format/Dbl.h>
#include <util/util/initStatic.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace DdMd;

int main()
{

   MPI::Init();
   Util::initStatic();
   Util::IntVector::commitMpiType();
   Util::Vector::commitMpiType();

   System system;
   system.readParam(std::cin); 

   int myRank = system.domain().gridRank();

   std::string filename("config");
   system.readConfig(filename);
   system.exchanger().exchangeAtoms();
   system.exchanger().exchangeGhosts();

   system.interaction().findNeighbors();

   double temperature = 1.0;
   system.setBoltzmannVelocities(temperature);

   // Calculate energies before integration
   double kinetic = system.kineticEnergy();
   double potential = system.pairPotentialEnergy();
   if (myRank == 0) {
      std::cout << Dbl(kinetic) << Dbl(potential) 
                << Dbl(kinetic + potential) << std::endl;
   }

   for (int i = 0; i < 5; ++i ) {

      system.integrate(1000);

      // Calculate energies after integration
      kinetic   = system.kineticEnergy();
      potential = system.pairPotentialEnergy();
      if (myRank == 0) {
         std::cout << Dbl(kinetic) << Dbl(potential) 
                   << Dbl(kinetic + potential) << std::endl;
      }
      //system.isValid();
   }

   MPI::Finalize();
}
#endif
