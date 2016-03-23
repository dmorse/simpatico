#ifndef MCMD_GENERALPOLYMER_ALT_SEMI_GRAND_MOVE_H
#define MCMD_GENERALPOLYMER_ALT_SEMI_GRAND_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>   // base class
#include <mcMd/chemistry/Molecule.h> //testing

namespace McMd
{

   using namespace Util;

   class GeneralpolymerSG;
   class McSystem;

   /**
   * A move that changes the type of a HomopolymerSG molecule.
   *
   * \ingroup McMd_McMove_Module
   */
   class GeneralpolymerALTSemiGrandMove : public SystemMove
   {
   
   public:
   
      /**
      * Constructor. 
      */
      GeneralpolymerALTSemiGrandMove(McSystem& system);
   
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

      /**
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();

      /**
      *   Choose a randomly selected molecule of the given subtype
      */ 
   
      Molecule& randomSGMolecule(int speciesId, int nSubType, int flipType);
   protected:
   
      /// Integer index for molecular species.
      int speciesId_;

      /// Pointer to instance of HomopolymerSG.
      GeneralpolymerSG* speciesPtr_;

      /// Limit on number of subspecies molecules that exist
      int Ulimit_;
      int Llimit_;
     
      /// Select which semigrand molecule subtype that needs flipping
      int flipSubtype_;
      
      double numoWeight_;
   };

}      
#endif
