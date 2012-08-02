#ifndef DDMD_CONFIG_IO_H
#define DDMD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>           // base class
#include <ddMd/communicate/AtomDistributor.h>    // member 
#include <ddMd/communicate/AtomCollector.h>      // member 
#include <ddMd/communicate/GroupDistributor.h>    // member 
#include <ddMd/communicate/GroupCollector.h>     // member 
#include <util/boundary/Boundary.h>              // typedef

#include <util/containers/DArray.h>              // member

#include <ddMd/chemistry/MaskPolicy.h>

namespace DdMd
{

   class Simulation;
   class Domain;
   class AtomStorage;
   class BondStorage;
   class AngleStorage;
   class DihedralStorage;
   class Buffer;

   using namespace Util;

   /**
   * Abstract reader/writer for configuration files.
   *
   * Each concrete subclass of ConfigIo implements a specific file format
   * by implementing the readConfig and writeConfig methods. 
   *
   * \ingroup DdMd_ConfigIo_Module
   */
   class ConfigIo  : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      ConfigIo();

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      ConfigIo(Simulation& simulation);

      /**
      * Associate with related objects.
      *
      * Required iff instantiated with default constructor.
      */
      void associate(Domain& domain, Boundary& boundary,
                     AtomStorage& atomStorage,
                     BondStorage& bondStorage,
                     #ifdef INTER_ANGLE
                     AngleStorage& angleStorage,
                     #endif
                     #ifdef INTER_DIHEDRAL
                     DihedralStorage& dihedralStorage,
                     #endif
                     Buffer& buffer);

      /**
      * Read cache size and allocate memory.
      */
      virtual void readParam(std::istream& in);

      /**
      * Read cache size and allocate memory.
      *
      * \param atomCacheCapacity size of internal atom cache. 
      * \param bondCacheCapacity size of internal bond cache. 
      */
      virtual void initialize(int atomCacheCapacity = 100,
                              int bondCacheCapacity = 100
                              #ifdef INTER_ANGLE
                              , int angleCacheCapacity = 100
                              #endif
                              #ifdef INTER_DIHEDRAL                             
                              , int dihedralCacheCapacity = 100
                              #endif
                              );

      /**
      * Read a configuration file.
      *
      * This reads a configuration file from the master processor, and 
      * distributes atoms and groups among the processors. Upon entry,
      * the file must be open for reading. Upon return:
      *
      *   - Each processor owns all atoms in its domain, and no others.
      *     Each atom is owned by one and only one processor.
      *
      *   - Each processor owns every group that contains one or more
      *     of the atoms that it owns, and no others. Each group may
      *     be owned by more than one processor.
      *
      *   - All Atom positions are expressed in generalized coordinates.
      *
      *   - All Atom Mask data is set correctly for specified maskPolicy.
      *
      *   - There are no ghost atoms on any processor.
      *
      * \param file input file stream
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::istream& file, MaskPolicy maskPolicy) = 0;

      /**
      * Write configuration file.
      *
      * This routine opens and writes a file on the master, collecting atom
      * and group data from all processors. Many file formats allow atom and
      * groups to be written in arbitrary order. Atomic positions may be 
      * written in Cartesian or generalized coordinates, depending on the
      * file format.
      *
      * \param file output file stream
      */
      virtual void writeConfig(std::ostream& file) = 0;

   protected:

      /**
      * Set Mask data on all atoms.
      */
      void setAtomMasks();

      /**
      * Get the AtomDistributor by reference.
      */
      AtomDistributor& atomDistributor();

      /**
      * Get the AtomCollector by reference.
      */
      AtomCollector& atomCollector();

      /**
      * Get the AtomDistributor by reference.
      */
      GroupDistributor<2>& bondDistributor();

      /**
      * Get the bond collector by reference.
      */
      GroupCollector<2>& bondCollector();

      #ifdef INTER_ANGLE
      /**
      * Get the angle distributor by reference.
      */
      GroupDistributor<3>& angleDistributor();

      /**
      * Get the angle collector by reference.
      */
      GroupCollector<3>& angleCollector();
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Get the dihedral distributor by reference.
      */
      GroupDistributor<4>& dihedralDistributor();

      /**
      * Get the dihedral collector by reference.
      */
      GroupCollector<4>& dihedralCollector();
      #endif

      /**
      * Get the Domain by reference.
      */
      Domain& domain();

      /**
      * Get Boundary by reference.
      */
      Boundary& boundary();
   
      /**
      * Get AtomStorage by reference.
      */
      AtomStorage& atomStorage();
   
      /**
      * Get BondStorage by reference.
      */
      BondStorage& bondStorage();
  
      #ifdef INTER_ANGLE 
      /**
      * Get AngleStorage by reference.
      */
      AngleStorage& angleStorage();
      #endif
   
      #ifdef INTER_DIHEDRAL
      /**
      * Get DihedralStorage by reference.
      */
      DihedralStorage& dihedralStorage();
      #endif
   
   private:

      AtomDistributor atomDistributor_;
      AtomCollector atomCollector_;

      GroupDistributor<2> bondDistributor_;
      GroupCollector<2> bondCollector_;

      #ifdef INTER_ANGLE
      GroupDistributor<3> angleDistributor_;
      GroupCollector<3> angleCollector_;
      #endif

      #ifdef INTER_DIHEDRAL
      GroupDistributor<4> dihedralDistributor_;
      GroupCollector<4> dihedralCollector_;
      #endif

      Domain* domainPtr_;

      Boundary* boundaryPtr_;

      AtomStorage* atomStoragePtr_;

      BondStorage* bondStoragePtr_;

      #ifdef INTER_ANGLE
      AngleStorage* angleStoragePtr_;
      #endif

      #ifdef INTER_DIHEDRAL
      DihedralStorage* dihedralStoragePtr_;
      #endif

      int  atomCacheCapacity_;

      int  bondCacheCapacity_;

      #ifdef INTER_ANGLE
      int  angleCacheCapacity_;
      #endif

      #ifdef INTER_DIHEDRAL
      int  dihedralCacheCapacity_;
      #endif

      /**
      * Read Group<N> objects from file. 
      */
      template <int N>
      int readGroups(std::istream& file, 
                     const char* sectionLabel, const char* nGroupLabel,
                     GroupDistributor<N>& distributor);

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      int writeGroups(std::ostream& file, 
                      const char* sectionLabel, const char* nGroupLabel,
                      GroupStorage<N>& storage, GroupCollector<N>& collector);

   };

   // Inline method definitions

   inline AtomDistributor& ConfigIo::atomDistributor()
   { return atomDistributor_; }

   inline AtomCollector& ConfigIo::atomCollector()
   { return atomCollector_; }

   inline GroupDistributor<2>& ConfigIo::bondDistributor()
   { return bondDistributor_; }

   inline GroupCollector<2>& ConfigIo::bondCollector()
   { return bondCollector_; }

   #ifdef INTER_ANGLE
   inline GroupDistributor<3>& ConfigIo::angleDistributor()
   { return angleDistributor_; }

   inline GroupCollector<3>& ConfigIo::angleCollector()
   { return angleCollector_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline GroupDistributor<4>& ConfigIo::dihedralDistributor()
   { return dihedralDistributor_; }

   inline GroupCollector<4>& ConfigIo::dihedralCollector()
   { return dihedralCollector_; }
   #endif

   inline Domain& ConfigIo::domain()
   { return *domainPtr_; }

   inline Boundary& ConfigIo::boundary()
   { return *boundaryPtr_; }

   inline AtomStorage& ConfigIo::atomStorage()
   { return *atomStoragePtr_; }

   inline BondStorage& ConfigIo::bondStorage()
   { return *bondStoragePtr_; }

   #ifdef INTER_ANGLE
   inline AngleStorage& ConfigIo::angleStorage()
   { return *angleStoragePtr_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline DihedralStorage& ConfigIo::dihedralStorage()
   { return *dihedralStoragePtr_; }
   #endif

}
#endif
