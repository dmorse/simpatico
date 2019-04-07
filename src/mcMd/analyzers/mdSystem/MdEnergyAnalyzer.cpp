/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdEnergyAnalyzer.h"
#include <mcMd/analyzers/base/AverageListAnalyzer.tpp>
#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdSystem.h>

#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/MdPairPotential.h>
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
#ifdef SIMP_COULOMB
#include <mcMd/potentials/coulomb/MdCoulombPotential.h>
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
   MdEnergyAnalyzer::MdEnergyAnalyzer(MdSystem& system) 
    : AverageListAnalyzer<MdSystem>(system),
      //#ifdef SIMP_COULOMB
      //coulombComponents_(false),
      //#endif
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
      #if 0
      #ifdef SIMP_COULOMB
      coulombId_(-1),
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      externalId_(-1),
      #endif
      potentialId_(-1),
      kineticId_(-1),
      totalId_(-1)
   {  setClassName("MdEnergyAnalyzer"); }

   /*
   * Read interval and outputFileName. 
   */
   void MdEnergyAnalyzer::readParameters(std::istream& in) 
   {
      AverageListAnalyzer<MdSystem>::readParameters(in);

      #if 0
      //#ifdef SIMP_COULOMB
      //if (sys.hasCoulombPotential()) {
      //   coulombComponents_ = false;
      //   readOptional<bool>(in, "coulombComponents", coulombComponents_);
      //}
      #endif

      MdSystem& sys = system();

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
      potentialId_ = id;
      ++id; 
      kineticId_ = id;
      ++id;
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
      setName(potentialId_, "potential");
      setName(kineticId_, "kinetic");
      setName(totalId_, "total");
   }

   /*
   * Load parameters from archive when restarting. 
   */
   void MdEnergyAnalyzer::loadParameters(Serializable::IArchive& ar) 
   {
      AverageListAnalyzer<MdSystem>::loadParameters(ar);

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
      ar >> potentialId_;
      ar >> kineticId_;
      ar >> totalId_;
   }

   /*
   * Save internal state to archive.
   */
   void MdEnergyAnalyzer::save(Serializable::OArchive& ar) 
   {
      AverageListAnalyzer<MdSystem>::save(ar);

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
      ar << potentialId_;
      ar << kineticId_;
      ar << totalId_;
   }

   /*
   * Output energy to file
   */
   void MdEnergyAnalyzer::compute() 
   {
      MdSystem& sys = system();

      double potential = 0.0;
      #ifndef SIMP_NOPAIR
      double pair = sys.pairPotential().energy();
      //outputFile_ << Dbl(pair, 15);
      setValue(pairId_, pair);
      potential += pair;
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         double bond = sys.bondPotential().energy();
         //outputFile_ << Dbl(bond, 15);
         setValue(bondId_, bond);
         potential += bond;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         double angle = sys.anglePotential().energy();
         //outputFile_ << Dbl(angle, 15);
         setValue(angleId_, angle);
         potential += angle;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         double dihedral  = sys.dihedralPotential().energy();
         // outputFile_ << Dbl(dihedral, 15);
         setValue(dihedralId_, dihedral);
         potential += dihedral;
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
         potential += coulomb;
      }
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         double external = sys.externalPotential().energy();
         // outputFile_ << Dbl(external, 15);
         setValue(externalId_, external);
         potential += external;
      }
      #endif

      // outputFile_ << Dbl(potential, 20)
      setValue(potentialId_, potential);

      double kinetic = sys.kineticEnergy();
      // outputFile_ << Dbl(kinetic, 20)
      setValue(kineticId_, kinetic);

      double total = kinetic + potential;
      // outputFile_ << Dbl(total, 20) << std::endl;
      setValue(totalId_, total);

   }

}
