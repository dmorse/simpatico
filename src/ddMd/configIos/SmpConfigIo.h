#ifndef DDMD_SMP_CONFIG_IO_H
#define DDMD_SMP_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/configIos/ConfigIo.h>
#include <simp/species/SpeciesGroup.h>

namespace DdMd
{

   class Simulation;

   using namespace Util;
   using namespace Simp;

   /**
   * Simpatico configuration file format.
   *
   * \ingroup DdMd_ConfigIo_Module
   */
   class SmpConfigIo  : public ConfigIo
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object
      */
      SmpConfigIo(Simulation& simulation);

      /**
      * Read configuration file.
      *
      * This routine opens and reads a file on the master, and distributes
      * atom data among the processors.
      *
      * \pre  There are no atoms, ghosts, or groups.
      * \pre  AtomStorage is set for scaled / generalized coordinates
      *
      * \post Atomic coordinates are scaled / generalized
      * \post There are no ghosts
      *
      * \param file input file stream (must be open on master)
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::ifstream& file, MaskPolicy maskPolicy);

      /**
      * Write configuration file.
      *
      * This routine opens and writes a file on the master, collecting 
      * atom data from all processors. Atomic positions are written in
      * Cartesian coordinates, independent of the coordinate system used
      * by the AtomStorage.
      *
      * \param file output file stream (must be open on master)
      */
      virtual void writeConfig(std::ofstream& file);
   
   private:

      Simulation* simulationPtr_;

      Simulation& simulation()
      {  return *simulationPtr_; }

      /**
      * Read Group<N> objects from file. 
      */
      template <int N>
      void readGroups(std::ifstream& file, 
                      const char* sectionLabel, const char* nGroupLabel,
                      GroupDistributor<N>& distributor);

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      int writeGroups(std::ofstream& file, 
                      const char* sectionLabel, const char* nGroupLabel,
                      GroupStorage<N>& storage, 
                      GroupCollector<N>& collector);
   

      /**
      * Compute and send Group<N> objects for a species.
      */
      template <int N>
      void sendSpeciesGroups(int& groupId, int& firstAtomId,
                             int  nMolecule, int  nAtom, int  nGroup,
                             const DArray<SpeciesGroup<N> >& groups,
                             GroupDistributor<N>& distributor);

      #ifdef SIMP_BOND
      /**
      * Compute all bonds from species data. 
      */
      void makeBonds();
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Compute all angles from species data. 
      */
      void makeAngles();
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Compute all dihedrals from species data. 
      */
      void makeDihedrals();
      #endif

   };

}
#endif
