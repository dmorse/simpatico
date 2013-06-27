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
#include <util/archives/Serializable.h>         // typedef used in interface

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   class Atom;
   class McSystem;

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
      virtual void readParam(std::istream& in);

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

      /**
      * CFB algorithm for deleting an end atom from a flexible chain.
      *
      * This function computes the energy of an end atom and a Rosenbluth 
      * factor for removing it. It does not remove the end atom from the 
      * system cell list.
      *
      * Upon return:
      *   
      *   - rosenbluth is the nonbonded Rosenblush factor for the deleted
      *     atom, i.e., the sum of Boltzmann factors from nonbonded pair
      *     interactions for the initial position and nTrial_ - 1 trials.
      *
      *   - energy is the total energy (bonded + nonbonded) of the end atom 
      *     before it was deleted.
      *
      * \param endPtr     ptr to end atom, which we attempt to remove
      * \param pvtPtr     ptr to atom next to end, or end after removal
      * \param bondType   type of bond connecting pvt and end atoms
      * \param rosenbluth nonbonded Rosenbluth factor of deleted atom (out)
      * \param energy     total potential energy of deleted atom (out)
      */
      void deleteEndAtom(Atom* endPtr, Atom* pvtPtr, int bondType,
                         double &rosenbluth, double &energy);

      /**
      * Rosenbluth algorithm for calculationg chemical potential.
      * 
      * This function computes the chemical potential of a system of linear
      * polymers by insertion method. Rosenbluth insertion method is used to
      * insert polymers.
      *
      * Upon return:
      *
      *   - rosenbluth factor of the whole chain is returned. 
      *
      * \param rosenbluth   rosenbluth factor of the whole grown chain
      */
      void chemicalPotential(double& rosenbluth, double& energy);

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

   protected:
   
      /// Maximum allowed number of trial positions for a regrown atom.
      static const int MaxTrial_ = 20; 
   
      /// Actual number of trial positions for each regrown atom.
      int  nTrial_; 
   
      /// Grant friend access to unit test class
      //  friend class CbEndBaseTest;

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;
   
      /// Has readParam been called?
      bool isInitialized_;
   };

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
   }

}
#endif 
