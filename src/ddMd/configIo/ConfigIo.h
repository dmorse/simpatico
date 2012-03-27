#ifndef DDMD_CONFIG_IO_H
#define DDMD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>           // base class
#include <ddMd/communicate/AtomDistributor.h>    // member 
#include <ddMd/communicate/GroupDistributor.h>    // member 
#include <ddMd/communicate/AtomCollector.h>      // member 
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
   class Buffer;

   using namespace Util;

   /**
   * 
   */
   class ConfigIo  : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      ConfigIo();

      /**
      * Associate with related objects.
      */
      void associate(Domain& domain, Boundary& boundary,
                     AtomStorage& atomStorage,
                     BondStorage& bondStorage,
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
                              int bondCacheCapacity = 100);

      /**
      * Read configuration file.
      *
      * This routine opens and reads a file on the master,
      * and distributes atom data among the processors.
      *
      * \param filename name of configuration file.
      */
      virtual void readConfig(std::istream& file, MaskPolicy maskPolicy);

      /**
      * Write configuration file.
      *
      * This routine opens and writes a file on the master,
      * collecting atom data from all processors.
      *
      * \param filename name of output configuration file.
      */
      virtual void writeConfig(std::ostream& file);

   protected:

      /**
      * Set masks on all atoms.
      */
      void setAtomMasks();

      /**
      * Get the AtomDistributor by reference.
      */
      AtomDistributor& atomDistributor();

      /**
      * Get the AtomDistributor by reference.
      */
      GroupDistributor<2>& bondDistributor();

      /**
      * Get the AtomCollector by reference.
      */
      AtomCollector& atomCollector();

      /**
      * Get the bond collector by reference.
      */
      GroupCollector<2>& bondCollector();

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
      * Get AtomStorage by reference.
      */
      BondStorage& bondStorage();
   
   private:

      AtomDistributor  atomDistributor_;

      GroupDistributor<2> bondDistributor_;

      AtomCollector  atomCollector_;

      GroupCollector<2>  bondCollector_;

      Domain* domainPtr_;

      Boundary* boundaryPtr_;

      AtomStorage* atomStoragePtr_;

      BondStorage* bondStoragePtr_;

      int  atomCacheCapacity_;

      int  bondCacheCapacity_;

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

   inline GroupDistributor<2>& ConfigIo::bondDistributor()
   { return bondDistributor_; }

   inline AtomCollector& ConfigIo::atomCollector()
   { return atomCollector_; }

   inline GroupCollector<2>& ConfigIo::bondCollector()
   { return bondCollector_; }

   inline Domain& ConfigIo::domain()
   { return *domainPtr_; }

   inline Boundary& ConfigIo::boundary()
   { return *boundaryPtr_; }

   inline AtomStorage& ConfigIo::atomStorage()
   { return *atomStoragePtr_; }

   inline BondStorage& ConfigIo::bondStorage()
   { return *bondStoragePtr_; }

}
#endif
