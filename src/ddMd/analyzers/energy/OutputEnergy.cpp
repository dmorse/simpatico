/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputEnergy.h"
#include <ddMd/potentials/pair/PairPotential.h>
#ifdef SIMP_BOND
#include <ddMd/potentials/bond/BondPotential.h>
#endif
#ifdef SIMP_ANGLE
#include <ddMd/potentials/angle/AnglePotential.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef SIMP_EXTERNAL
#include <ddMd/potentials/external/ExternalPotential.h>
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
    : Analyzer(simulation),
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
      #if 0
      std::string filename = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      #endif
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
      #if 0
      std::string filename = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      #endif
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
   * Open outputfile
   */ 
   void OutputEnergy::setup()
   {
      if (simulation().domain().isMaster()) {
         std::string filename;
         filename  = outputFileName();
         simulation().fileMaster().openOutputFile(filename, outputFile_);
      }
   }

   /*
   * Output energy to file
   */
   void OutputEnergy::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sim = simulation();
         sim.computeKineticEnergy();
         sim.computePotentialEnergies();
         if (sim.domain().isMaster()) {
            double kinetic   = sim.kineticEnergy();
            outputFile_ << Int(iStep, 10)
                        << Dbl(kinetic, 15);
            double potential = 0.0;
            double pair = sim.pairPotential().energy();
            potential += pair;
            outputFile_ << Dbl(pair, 15);
            #ifdef SIMP_BOND
            if (sim.nBondType()) {
               double bond = sim.bondPotential().energy();
               potential += bond;
               outputFile_ << Dbl(bond, 15);
            }
            #endif
            #ifdef SIMP_ANGLE
            if (sim.nAngleType()) {
               double angle = sim.anglePotential().energy();
               potential += angle;
               outputFile_ << Dbl(angle, 15);
            }
            #endif
            #ifdef SIMP_DIHEDRAL
            if (sim.nDihedralType()) {
               double dihedral  = sim.dihedralPotential().energy();
               potential += dihedral;
               outputFile_ << Dbl(dihedral, 15);
            }
            #endif
            #ifdef SIMP_EXTERNAL
            if (sim.hasExternal()) {
               double external = sim.externalPotential().energy();
               potential += external;
               outputFile_ << Dbl(external, 15);
            }
            #endif
            outputFile_ << Dbl(kinetic + potential, 20)
                        << std::endl;
         }
         ++nSample_;
      }
   }

}
