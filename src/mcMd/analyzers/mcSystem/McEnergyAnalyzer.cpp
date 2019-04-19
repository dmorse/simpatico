/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McEnergyAnalyzer.h"
#include <mcMd/analyzers/base/AverageListAnalyzer.tpp>
#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>

#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#ifdef SIMP_BOND
#include <mcMd/potentials/bond/BondPotential.h>
#endif
#ifdef SIMP_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef SIMP_EXTERNAL
#include <mcMd/potentials/external/ExternalPotential.h>
#endif

#include <util/accumulators/Average.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McEnergyAnalyzer::McEnergyAnalyzer(McSystem& system) 
    : AverageListAnalyzer<McSystem>(system),
      #ifndef SIMP_NOPAIR
      pairId_(-1),
      #endif
      #ifdef SIMP_BOND
      bondId_(-1),
      #endif
      #ifdef SIMP_ANGLE
      angleId_(-1),
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralId_(-1),
      #endif
      #ifdef SIMP_EXTERNAL
      externalId_(-1),
      #endif
      totalId_(-1)

   {  setClassName("McEnergyAnalyzer"); }

   /*
   * Read interval and outputFileName. 
   */
   void McEnergyAnalyzer::readParameters(std::istream& in) 
   {
      AverageListAnalyzer<McSystem>::readParameters(in);

      #if 0
      //#ifdef SIMP_COULOMB
      //if (sys.hasCoulombPotential()) {
      //   coulombComponents_ = false;
      //   readOptional<bool>(in, "coulombComponents", coulombComponents_);
      //}
      #endif

      McSystem& sys = system();

      // Count number of values
      int id = 0;
      #ifndef SIMP_NOPAIR
      pairId_ = id;
      ++id;
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         bondId_ = id;
         ++id;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         angleId_ = id;
         ++id;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         dihedralId_ = id;
         ++id;
      }
      #endif
      #if 0
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         coulombId_ = id;
         ++id;
      }
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         externalId_ = id;
         ++id;
      }
      #endif
      totalId_ = id;
      ++id; 
      initializeAccumulators(id);

      // Set names
      #ifndef SIMP_NOPAIR
      setName(pairId_, "pair");
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         setName(bondId_, "bond");
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         setName(angleId_, "angle");
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         setName(dihedralId_, "dihedral");
      }
      #endif
      #if 0
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         setName(coulombId_, "coulomb");
      }
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         setName(externalId_, "external");
      }
      #endif
      setName(totalId_, "total");
   }



   /*
   * Load parameters from archive when restarting. 
   */
   void McEnergyAnalyzer::loadParameters(Serializable::IArchive& ar) 
   {
      AverageListAnalyzer<McSystem>::loadParameters(ar);

      #ifndef SIMP_NOPAIR
      ar >> pairId_;
      #endif
      #ifdef SIMP_BOND
      ar >> bondId_;
      #endif
      #ifdef SIMP_ANGLE
      ar >> angleId_;
      #endif
      #ifdef SIMP_DIHEDRAL
      ar >> dihedralId_;
      #endif
      #if 0
      #ifdef SIMP_COULOMB
      ar >> coulombId_;
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      ar >> externalId_;
      #endif
      ar >> totalId_;
   }

   /*
   * Save internal state to archive.
   */
   void McEnergyAnalyzer::save(Serializable::OArchive& ar) 
   {
      AverageListAnalyzer<McSystem>::save(ar);

      #ifndef SIMP_NOPAIR
      ar << pairId_;
      #endif
      #ifdef SIMP_BOND
      ar << bondId_;
      #endif
      #ifdef SIMP_ANGLE
      ar << angleId_;
      #endif
      #ifdef SIMP_DIHEDRAL
      ar << dihedralId_;
      #endif
      #if 0
      #ifdef SIMP_COULOMB
      ar << coulombId_;
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      ar << externalId_;
      #endif
      ar << totalId_;
   }

   /*
   * Output energy to file
   */
   void McEnergyAnalyzer::compute() 
   {
      McSystem& sys = system();

      double total = 0.0;
      #ifndef SIMP_NOPAIR
      double pair = sys.pairPotential().energy();
      //outputFile_ << Dbl(pair, 15);
      setValue(pairId_, pair);
      total += pair;
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         double bond = sys.bondPotential().energy();
         //outputFile_ << Dbl(bond, 15);
         setValue(bondId_, bond);
         total += bond;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         double angle = sys.anglePotential().energy();
         //outputFile_ << Dbl(angle, 15);
         setValue(angleId_, angle);
         total += angle;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         double dihedral  = sys.dihedralPotential().energy();
         // outputFile_ << Dbl(dihedral, 15);
         setValue(dihedralId_, dihedral);
         total += dihedral;
      }
      #endif
      #if 0
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         //if (coulombComponents_) {
         //   // R-space contribution
         //   double coulombRSpace = sys.coulombPotential().rSpaceEnergy();
         //   // K-space contribution
         //   double coulombKSpace = sys.coulombPotential().kSpaceEnergy();
         //}
         double coulomb = sys.coulombPotential().energy();
         // outputFile_ << Dbl(coulomb, 15);
         setValue(coulombId_, coulomb);
         total += coulomb;
      }
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         double external = sys.externalPotential().energy();
         // outputFile_ << Dbl(external, 15);
         setValue(externalId_, external);
         total += external;
      }
      #endif

      // outputFile_ << Dbl(total, 20)
      setValue(totalId_, total);

   }

}
