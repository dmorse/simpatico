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
    : AverageListAnalyzer<MdSystem>(system)
      //#ifdef SIMP_COULOMB
      //coulombComponents_(false),
      //#endif
   {  setClassName("MdEnergyAnalyzer"); }

   /*
   * Read interval and outputFileName. 
   */
   void MdEnergyAnalyzer::readParameters(std::istream& in) 
   {
      AverageListAnalyzer<MdSystem>::readParameters(in);

      #if 0
      readInterval(in);
      readOutputFileName(in);
      readNSamplePerBlock(in);
      //#ifdef SIMP_COULOMB
      //if (sys.hasCoulombPotential()) {
      //   coulombComponents_ = false;
      //   readOptional<bool>(in, "coulombComponents", coulombComponents_);
      //}
      //#endif
      #endif

      MdSystem& sys = system();

      // Count number of values
      int id = 0;
      #ifndef SIMP_NOPAIR
      ++id;
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) ++id;
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) ++id;
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) ++id;
      #endif
      #if 0
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) ++id;
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) ++id;
      #endif
      ++id; // potential
      ++id; // kinetic
      ++id; // total
      initializeAccumulators(id);

      // Set names
      id = 0;
      #ifndef SIMP_NOPAIR
      setName(id, "pair");
      ++id;
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         setName(id, "bond");
         ++id;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         setName(id, "angle");
         ++id;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         setName(id, "dihedral");
         ++id;
      }
      #endif
      #if 0
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         setName(id, "coulomb");
         ++id;
      }
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         setName(id, "external");
         ++id;
      }
      #endif
      setName(id,"potential");
      ++id;
      setName(id,"kinetic");
      ++id;
      setName(id,"total");
      ++id;

   }

   /*
   * Output energy to file
   */
   void MdEnergyAnalyzer::compute() 
   {
      MdSystem& sys = system();

      int id = 0;
      double potential = 0.0;
      #ifndef SIMP_NOPAIR
      double pair = sys.pairPotential().energy();
      //outputFile_ << Dbl(pair, 15);
      setValue(id, pair);
      ++id;
      potential += pair;
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         double bond = sys.bondPotential().energy();
         //outputFile_ << Dbl(bond, 15);
         potential += bond;
         setValue(id, bond);
         ++id;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         double angle = sys.anglePotential().energy();
         //outputFile_ << Dbl(angle, 15);
         potential += angle;
         setValue(id, angle);
         ++id;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         double dihedral  = sys.dihedralPotential().energy();
         // outputFile_ << Dbl(dihedral, 15);
         setValue(id, dihedral);
         ++id;
         potential += dihedral;
      }
      #endif
      #if 0
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         if (coulombComponents_) {
            // R-space contribution
            double coulombRSpace = sys.coulombPotential().rSpaceEnergy();
            // K-space contribution
            double coulombKSpace = sys.coulombPotential().kSpaceEnergy();
         }
         // Total
         double coulomb = sys.coulombPotential().energy();
         // outputFile_ << Dbl(coulomb, 15);
         setValue(id, coulomb);
         ++id;
         potential += coulomb;
         // outputFile_ << Dbl(dihedral, 15);
      }
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         double external = sys.externalPotential().energy();
         potential += external;
         setValue(id, external);
         ++id;
         // outputFile_ << Dbl(external, 15);
      }
      #endif

      // outputFile_ << Dbl(potential, 20)
      setValue(id, potential);
      ++id;

      double kinetic = sys.kineticEnergy();
      // outputFile_ << Dbl(kinetic, 20)
      setValue(id, kinetic);
      ++id;

      double total = kinetic + potential;
      // outputFile_ << Dbl(total, 20)
      //             << std::endl;
      setValue(id, total);
      ++id;

   }

}
