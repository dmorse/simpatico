#ifndef DDMD_OUTPUT_ENERGY_CPP
#define DDMD_OUTPUT_ENERGY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputEnergy.h"
#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/potentials/bond/BondPotential.h>
#ifdef INTER_ANGLE
#include <ddMd/potentials/angle/AnglePotential.h>
#endif
#ifdef INTER_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotential.h>
#endif
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   OutputEnergy::OutputEnergy(Simulation& simulation) 
    : Diagnostic(simulation),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("OutputEnergy"); }

   /*
   * Read interval and outputFileName. 
   */
   void OutputEnergy::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void OutputEnergy::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void OutputEnergy::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
   }

  
   /*
   * Reset nSample.
   */
   void OutputEnergy::clear() 
   {  nSample_ = 0;  }

   /*
   * Dump configuration to file
   */
   void OutputEnergy::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeKineticEnergy();
         sys.computePotentialEnergies();
         if (sys.domain().isMaster()) {
            double kinetic   = sys.kineticEnergy();
            double pair      = sys.pairPotential().energy();
            double potential = pair;
            double bond      = sys.bondPotential().energy();
            potential += bond;
            Log::file() << Int(iStep, 10)
                        << Dbl(kinetic, 15)
                        << Dbl(pair, 15)
                        << Dbl(bond, 15);
            #ifdef INTER_ANGLE
            if (sys.nAngleType()) {
               double angle = sys.anglePotential().energy();
               potential += angle;
               Log::file() << Dbl(angle, 15);
            }
            #endif
            #ifdef INTER_DIHEDRAL
            if (sys.nDihedralType()) {
               double dihedral  = sys.dihedralPotential().energy();
               potential += dihedral;
               Log::file() << Dbl(dihedral, 15);
            }
            #endif
            Log::file() << Dbl(kinetic + potential, 20)
                        << std::endl;
         }
         ++nSample_;
      }
   }

}
#endif 
