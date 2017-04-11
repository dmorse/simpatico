#ifndef MCMD_WANG_LANDAU_ALT_MOVE_H
#define MCMD_WANG_LANDAU_ALT_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>   // base class
#include <util/containers/DArray.h>
#include <util/containers/Pair.h>
#include <mcMd/chemistry/Molecule.h> //testing

namespace McMd
{

   using namespace Util;
   class SpeciesMutator;
   class GeneralpolymerSG;
   class McSystem;

   /**
   * A move that changes the type of a HomopolymerSG molecule.
   *
   * \ingroup McMd_McMove_Module
   */
   class WangLandauALTMove : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      */
      WangLandauALTMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);
    
      Molecule& randomSGMolecule(int speciesId, int typeId, int flipType);
      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();
      
      virtual void output();
      // Determine whether or not to adapt the step size 
      virtual void stepAdapt();      

  protected:
      /// Integer index for molecular species.
      int speciesId_;

      double weightSize_;
      
      int Ulimit_;

      int Llimit_;
 
      int flipSubtype_;
 
      DArray<double> weights_;

      DArray<int> stateCount_;

      /// Pointer to instance of HomopolymerSG.
      GeneralpolymerSG* speciesPtr_;
      SpeciesMutator* mutatorPtr_;

      std::string outputFileName_;
      std::string initialWeightFileName_;
      double crit_;  
      DArray<int> steps_;
      DArray<double> weightTrack_;
      int stepCount_;
      int flipType_;
      std::ofstream outputFile_;
      //Total number of semigrand molecules possible
      int capacity_;
   };

}      
#endif
