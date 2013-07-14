#ifndef MCMD_MC_CHEMICAL_POTENTIAL_H
#define MCMD_MC_CHEMICAL_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/mcSimulation/McSystem.h>         // base template parameter
#include <util/accumulators/Average.h>          // member
#include <util/random/Random.h>                 // member
#include <util/ensembles/EnergyEnsemble.h>      // inline function
#include <util/archives/Serializable.h>         // typedef used in interface

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   class Atom;

   /**
   * McEnergyAverage averages of total potential energy.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class McChemicalPotential : public SystemDiagnostic<McSystem>
   {
   
   public:

      /**   
      * Constructor.
      */
      McChemicalPotential(McSystem& system);

      /**
      * Read parameters and initialize.
      */
      virtual void readParameters(std::istream& in);

      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /* 
      * Evaluate energy per particle, and add to ensemble. 
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

      /**
      * Save state to binary file archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from a binary file archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void load(Serializable::IArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// Get parent McSystem.
      McSystem& system();

      /// Get EnergyEnsemble object of parent McSystem.
      Simulation& simulation();

      /// Get Boundary object of parent McSystem.
      Boundary& boundary();

      /// Get EnergyEnsemble object of parent McSystem.
      EnergyEnsemble& energyEnsemble();

      /// Get Random number generator of parent Simulation.
      Random& random();

      /**
      * Boltzmann weight associated with an energy difference.
      *
      * \param energy energy or energy difference.
      * \return Boltzmann weight exp(-beta*energy)
      */
      double boltzmann(double energy);

      /**
      * Configuration bias algorithm for adding an atom to a chain end.
      *
      * This function generates and computes Rosenbluth factors for nTrial 
      * trial positions, chooses one, updates the atomic position. It does
      * not add the end atom to the system cell list.
      *
      * Upon return:
      *  
      *   - rosenbluth is the nonbonded Rosenblush factor for the added atom,
      *     i.e., the sum of Boltzmann factors from nonbonded pair interactions
      *     for all nTrial_ trial positions.
      *
      *+   - energy is the total energy (bonded + nonbonded) of the new end
      *     atom in its chosen position.
      *
      * \param endPtr     ptr to new end atom, which we attempt to add
      * \param pvtPtr     end atom of current chain, next to end
      * \param bondType   type of bond connecting pvt and end atoms
      * \param rosenbluth Rosenbluth factor of added atom (out)
      * \param energy     potential energy of deleted atom (out)
      */
      void addEndAtom(Atom* endPtr, Atom* pvtPtr, int bondType,
                         double &rosenbluth, double &energy);

      /// Maximum allowed number of trial positions for a regrown atom.
      static const int MaxTrial_ = 200;

      /// Pointer to parent McSystem object.
      McSystem  *systemPtr_;

      /// Pointer to Simulation of parent System.
      Simulation  *simulationPtr_;

      /// Pointer to Boundary of parent System.
      Boundary  *boundaryPtr_;

      /// Pointer to EnergyEnsemble of parent System.
      EnergyEnsemble  *energyEnsemblePtr_;

      /// Pointer to Random of parent System.
      Random  *randomPtr_;
   
      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Actual number of trial positions for each regrown atom.
      int  nTrial_; 

      /// Actual number of trial positions for each regrown atom.
      int  speciesId_; 

      /// Actual number of trial positions for each regrown atom.
      int  nMoleculeTrial_; 
   
      /// Number of samples per block average output.
      int nSamplePerBlock_;
   
      /// Has readParam been called?
      bool isInitialized_;

   };

   // Inline methods

   /*
   * Get parent McSystem.
   */
   inline McSystem& McChemicalPotential::system()
   {  return *systemPtr_; }

   /*
   * Get Simulation object of parent McSystem.
   */
   inline Simulation& McChemicalPotential::simulation()
   {  return *simulationPtr_; }

   /*
   * Get Boundary object of parent McSystem.
   */
   inline Boundary& McChemicalPotential::boundary()
   {  return *boundaryPtr_; }

   /*
   * Get EnergyEnsemble object of parent McSystem.
   */
   inline EnergyEnsemble& McChemicalPotential::energyEnsemble()
   {  return *energyEnsemblePtr_; }

   /*
   * Get random object of parent McSystem.
   */
   inline Random& McChemicalPotential::random()
   {  return *randomPtr_; }

   /*
   * Boltzmann weight associated with an energy difference.
   */
   inline double McChemicalPotential::boltzmann(double energy)
   {  return exp(-energyEnsemblePtr_->beta()*energy); }

   // Function template

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void McChemicalPotential::serialize(Archive& ar, const unsigned int version)
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }
      ar & accumulator_;
      ar & nSamplePerBlock_;
      ar & nTrial_;
      ar & nMoleculeTrial_;
      ar & speciesId_;
   }

}
#endif 
