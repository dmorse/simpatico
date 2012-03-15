#ifndef DDSIM_CPP
#define DDSIM_CPP

#include <ddMd/system/System.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/potentials/bond/BondPotential.h>
#include <util/random/Random.h>
#include <util/format/Dbl.h>
#include <util/util/initStatic.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace DdMd;

int main()
{

   System system;
   system.readParam(std::cin); 

   system.readCommands();

   #if 0
   MPI::Init();
   Util::initStatic();
   Util::IntVector::commitMpiType();
   Util::Vector::commitMpiType();

   int myRank = system.domain().gridRank();

   std::string filename("config");
   system.readConfig(filename);
   system.exchanger().exchange();

   system.pairPotential().findNeighbors();

   double temperature = 1.0;
   system.setBoltzmannVelocities(temperature);

   // Calculate energies before integration
   double kinetic = system.kineticEnergy();
   double pair = system.pairPotentialEnergy();
   double bond = system.bondPotentialEnergy();
   if (myRank == 0) {
      std::cout << Dbl(kinetic) << Dbl(pair) << Dbl(bond)
                << Dbl(kinetic + pair + bond) << std::endl;
   }

   for (int i = 0; i < 5; ++i ) {

      system.integrate(1000);

      // Calculate energies after integration
      kinetic = system.kineticEnergy();
      pair = system.pairPotentialEnergy();
      bond = system.bondPotentialEnergy();
      if (myRank == 0) {
         std::cout << Dbl(kinetic) << Dbl(pair) << Dbl(bond)
                   << Dbl(kinetic + pair + bond) << std::endl;
      }
      //system.isValid();
   }
   #endif

   //MPI::Finalize();
}
#endif
