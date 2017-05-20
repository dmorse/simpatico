#ifndef DDMD_TRAJECTORY_WRITER_H
#define DDMD_TRAJECTORY_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>         // base class
#include <ddMd/storage/AtomStorage.h>       
#ifdef SIMP_BOND
#include <ddMd/storage/BondStorage.h>               
#endif
#ifdef SIMP_ANGLE
#include <ddMd/storage/AngleStorage.h>               
#endif
#ifdef SIMP_DIHEDRAL
#include <ddMd/storage/DihedralStorage.h>               
#endif
#include <util/boundary/Boundary.h>         // typedef

namespace DdMd
{

   class Simulation;
   class Domain;
   class AtomCollector;
   template <int N> class GroupCollector;

   using namespace Util;

   /**
   * Base class to write a trajectory to a single file.
   *
   * \ingroup DdMd_Analyzer_Trajectory_Module
   */
   class TrajectoryWriter : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object
      * \param isBinary Is the trajectory file a binary file?
      */
      TrajectoryWriter(Simulation& simulation, bool isBinary = false);
   
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
      * Open the trajectory file. 
      */
      virtual void setup();

      /**
      * Dump configuration to file
      *
      * \param iStep MC step index
      */
      virtual void sample(long iStep);

      /**
      * Clear nSample counter.
      */
      virtual void clear();
  
      /**
      * Close ouput file.
      */
      virtual void output();

   protected:

      /**
      * Is the file format binary (true) or text (false)?
      */
      bool isBinary() const;

      /**
      * Write data that should appear once, at beginning of the file. 
      *
      * Called by sample on first invocation. Default implementation is empty.
      *
      * \param out output file stream
      */
      virtual void writeHeader(std::ofstream& out)
      {};

      /**
      * Write data that should appear in every frame.
      * 
      * \param out output file stream
      * \param iStep MD time step index
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

      #ifdef SIMP_BOND
      /**
      * Get BondStorage by reference.
      */
      BondStorage& bondStorage();
  
      /**
      * Get the bond collector by reference.
      */
      GroupCollector<2>& bondCollector();
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Get AngleStorage by reference.
      */
      AngleStorage& angleStorage();

      /**
      * Get the angle collector by reference.
      */
      GroupCollector<3>& angleCollector();
      #endif

      #ifdef SIMP_DIHEDRAL
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

      /// Has readParam been called?
      long isInitialized_;
  
      /// Is the trajectory file a binary file? 
      bool isBinary_;

      // Pointers to associated Domain.
      Domain* domainPtr_;

      // Pointers to associated Boundary.
      Boundary* boundaryPtr_;

      // Pointers to storage objects.
      AtomStorage* atomStoragePtr_;
      #ifdef SIMP_BOND
      BondStorage* bondStoragePtr_;
      #endif
      #ifdef SIMP_ANGLE
      AngleStorage* angleStoragePtr_;
      #endif
      #ifdef SIMP_DIHEDRAL
      DihedralStorage* dihedralStoragePtr_;
      #endif

   };

   // Inline method definitions

   /**
   * Is the file format binary (true) or text (false)?
   */
   inline bool TrajectoryWriter::isBinary() const
   {  return isBinary_; }

   inline Domain& TrajectoryWriter::domain()
   {  return *domainPtr_; }

   inline Boundary& TrajectoryWriter::boundary()
   {  return *boundaryPtr_; }

   inline AtomStorage& TrajectoryWriter::atomStorage()
   {  return *atomStoragePtr_; }

   inline AtomCollector& TrajectoryWriter::atomCollector()
   {  return atomStoragePtr_->collector(); }

   #ifdef SIMP_BOND
   inline BondStorage& TrajectoryWriter::bondStorage()
   {  return *bondStoragePtr_; }

   inline GroupCollector<2>& TrajectoryWriter::bondCollector()
   {  return bondStoragePtr_->collector(); }
   #endif

   #ifdef SIMP_ANGLE
   inline AngleStorage& TrajectoryWriter::angleStorage()
   {  return *angleStoragePtr_; }

   inline GroupCollector<3>& TrajectoryWriter::angleCollector()
   {  return angleStoragePtr_->collector(); }
   #endif

   #ifdef SIMP_DIHEDRAL
   inline DihedralStorage& TrajectoryWriter::dihedralStorage()
   {  return *dihedralStoragePtr_; }

   inline GroupCollector<4>& TrajectoryWriter::dihedralCollector()
   {  return dihedralStoragePtr_->collector(); }
   #endif

}
#endif 
