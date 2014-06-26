#ifndef DDMD_TRAJECTORY_WRITER_H
#define DDMD_TRAJECTORY_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>             // base class
//#include <ddMd/simulation/Simulation.h>
#include <ddMd/communicate/AtomCollector.h>      // member 
#include <ddMd/communicate/GroupCollector.h>     // member 
#include <util/boundary/Boundary.h>              // typedef

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
   * Based class for classes that write trajectories to a single file.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class TrajectoryWriter : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      TrajectoryWriter(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~TrajectoryWriter()
      {} 
   
      /**
      * Read parameters and initialize.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Clear nSample counter.
      */
      virtual void clear();
  
      /**
      * Dump configuration to file
      *
      * \param iStep MC step index
      */
      virtual void sample(long iStep);

      /**
      * Close ouput file.
      */
      virtual void output();

   protected:

      /**
      * Write data that should appear once, at beginning of the file. 
      *
      * Called by sample when iStep == 0.
      */
      virtual void writeHeader(std::ofstream& out, long iStep) = 0;

      /**
      * Write data that should appear in every frame.
      * 
      * Called by sample on every step.
      */
      virtual void writeFrame(std::ofstream& out, long iStep) = 0;

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
      * Get the AtomCollector by reference.
      */
      AtomCollector& atomCollector();

      #ifdef INTER_BOND
      /**
      * Get BondStorage by reference.
      */
      BondStorage& bondStorage();
  
      /**
      * Get the bond collector by reference.
      */
      GroupCollector<2>& bondCollector();
      #endif

      #ifdef INTER_ANGLE
      /**
      * Get AngleStorage by reference.
      */
      AngleStorage& angleStorage();

      /**
      * Get the angle collector by reference.
      */
      GroupCollector<3>& angleCollector();
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Get DihedralStorage by reference.
      */
      DihedralStorage& dihedralStorage();

      /**
      * Get the dihedral collector by reference.
      */
      GroupCollector<4>& dihedralCollector();
      #endif

   private:
 
      // Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;
   
      /// Has readParam been called?
      long isInitialized_;
   
      // Pointers to associated objects.
      Domain* domainPtr_;
      Boundary* boundaryPtr_;

      // Storage and collectors
      AtomCollector atomCollector_;
      AtomStorage* atomStoragePtr_;
      int atomCacheCapacity_;
      #ifdef INTER_BOND
      GroupCollector<2> bondCollector_;
      BondStorage* bondStoragePtr_;
      int bondCacheCapacity_;
      #endif
      #ifdef INTER_ANGLE
      GroupCollector<3> angleCollector_;
      AngleStorage* angleStoragePtr_;
      int angleCacheCapacity_;
      #endif
      #ifdef INTER_DIHEDRAL
      GroupCollector<4> dihedralCollector_;
      DihedralStorage* dihedralStoragePtr_;
      int dihedralCacheCapacity_;
      #endif

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      int writeGroups(std::ofstream& file, 
                      const char* sectionLabel, const char* nGroupLabel,
                      GroupStorage<N>& storage, GroupCollector<N>& collector);

   };

   // Inline method definitions

   inline Domain& TrajectoryWriter::domain()
   {  return *domainPtr_; }

   inline Boundary& TrajectoryWriter::boundary()
   {  return *boundaryPtr_; }

   inline AtomStorage& TrajectoryWriter::atomStorage()
   {  return *atomStoragePtr_; }

   inline AtomCollector& TrajectoryWriter::atomCollector()
   {  return atomCollector_; }

   #ifdef INTER_BOND
   inline BondStorage& TrajectoryWriter::bondStorage()
   {  return *bondStoragePtr_; }

   inline GroupCollector<2>& TrajectoryWriter::bondCollector()
   {  return bondCollector_; }
   #endif

   #ifdef INTER_ANGLE
   inline AngleStorage& TrajectoryWriter::angleStorage()
   {  return *angleStoragePtr_; }

   inline GroupCollector<3>& TrajectoryWriter::angleCollector()
   {  return angleCollector_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline DihedralStorage& TrajectoryWriter::dihedralStorage()
   {  return *dihedralStoragePtr_; }

   inline GroupCollector<4>& TrajectoryWriter::dihedralCollector()
   {  return dihedralCollector_; }
   #endif

}
#endif 
