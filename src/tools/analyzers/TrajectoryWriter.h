#ifndef TOOLS_TRAJECTORY_WRITER_H
#define TOOLS_TRAJECTORY_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <tools/analyzers/Analyzer.h>        // base class
#include <util/boundary/Boundary.h>          // typedef

namespace Tools
{

   class Processor;
   class Configuration;
   class AtomStorage;
   template <int N> class GroupStorage;

   using namespace Util;

   /**
   * Based class for classes that write trajectories to a single file.
   *
   * \ingroup Tools_Analyzer_Module
   */
   class TrajectoryWriter : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param processor parent Processor object
      * \param isBinary Is the trajectory file a binary file?
      */
      TrajectoryWriter(Processor& processor, bool isBinary);
   
      /**
      * Constructor.
      *
      * \param configuration parent Configuration object
      * \param fileMaster associated FileMaster object
      * \param isBinary Is the trajectory file a binary file?
      */
      TrajectoryWriter(Configuration& configuration, FileMaster& fileMaster,
                       bool isBinary);
   
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
  
      #if 0 
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
      #endif

      /**
      * Open the trajectory file. 
      */
      virtual void setup();

      /**
      * Write frame to file
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
      * Called by sample on every step.
      *
      * \param out output file stream
      * \param iStep MD time step index
      */
      virtual void writeFrame(std::ofstream& out, long iStep) = 0;

      /**
      * Get Boundary by reference.
      */
      Boundary& boundary();
   
      /**
      * Get AtomStorage by reference.
      */
      AtomStorage& atoms();
   
      #ifdef SIMP_BOND
      /**
      * Get bond storage by reference.
      */
      GroupStorage<2>& bonds();
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Get angle storage by reference.
      */
      GroupStorage<3>& angles();
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Get dihedral storage by reference.
      */
      GroupStorage<4>& dihedrals();
      #endif

   private:
 
      // Output file stream
      std::ofstream outputFile_;

      /// Number of frames written thus far (first is zero).
      long nSample_;
   
      /// Has readParam been called?
      long isInitialized_;
  
      /// Is the trajectory file a binary file? 
      bool isBinary_;

      // Pointers to associated Boundary.
      Boundary* boundaryPtr_;

      // Pointers to storage objects.
      AtomStorage* atomStoragePtr_;
      #ifdef SIMP_BOND
      GroupStorage<2>* bondStoragePtr_;
      #endif
      #ifdef SIMP_ANGLE
      GroupStorage<3>* angleStoragePtr_;
      #endif
      #ifdef SIMP_DIHEDRAL
      GroupStorage<4>* dihedralStoragePtr_;
      #endif

   };

   // Inline method definitions

   /**
   * Is the file format binary (true) or text (false)?
   */
   inline bool TrajectoryWriter::isBinary() const
   {  return isBinary_; }

   inline Boundary& TrajectoryWriter::boundary()
   {  return *boundaryPtr_; }

   inline AtomStorage& TrajectoryWriter::atoms()
   {  return *atomStoragePtr_; }

   #ifdef SIMP_BOND
   inline GroupStorage<2>& TrajectoryWriter::bonds()
   {  return *bondStoragePtr_; }
   #endif

   #ifdef SIMP_ANGLE
   inline GroupStorage<3>& TrajectoryWriter::angles()
   {  return *angleStoragePtr_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   inline GroupStorage<4>& TrajectoryWriter::dihedrals()
   {  return *dihedralStoragePtr_; }
   #endif

}
#endif 
