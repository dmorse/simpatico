#ifndef MDPP_PROCESSOR_H
#define MDPP_PROCESSOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>        // base class
#include <util/boundary/Boundary.h>           // member 
#include <util/containers/DSArray.h>          // member (template)
#include <mdPp/chemistry/Atom.h>              // member (template argument)
#include <mdPp/chemistry/Group.h>             // member (template argument)
#include <mdPp/configIos/ConfigIoFactory.h>   // member 
#include <mdPp/analyzers/AnalyzerManager.h>   // member 

namespace MdPp 
{

   class ConfigIo;

   using namespace Util;

   /**
   * A post-processor for analyzing outputs of MD simulations.
   */
   class Processor : public ParamComposite 
   {

   public:

      using ParamComposite::readParam;

      /**
      * Constructor
      */
      Processor();

      /**
      * Destructor
      */
      ~Processor();

      /**
      * Read parameters and analyze.
      */
      void readParameters(std::istream& in);

      /**
      * Analyze a sequence of dump files.
      */
      void analyzeDumps(std::string& filename);

      /**
      * Analyze a trajectory file.
      */
      void analyzeTrajectory(std::string& filename);

      // Mutators for use in ConfigIo classes.

      /**
      * Return pointer to location for new atom, and add to container.
      *
      * \param  global id for new atom
      * \return pointer to location of new atom
      */
      Atom* addAtom(int id);

      /**
      * Add bond to list of those read from file.
      *
      * \return pointer to location of new atom
      */
      Group<2>* addBond();

      // etc. for angles and dihedrals
  
      // Accessors, for use in Analyzer and ConfigIo classes.

      /**
      * Get the Boundary by non-const reference
      */
      Boundary& boundary();

      /**
      * Get number of atoms.
      */ 
      int nAtom() const;

      /**
      * Get an atom reference by global id.
      */
      Atom& atom(int id);

      /**
      * Initialize an iterator for atoms.
      */
      void initAtomIterator(ArrayIterator<Atom>& iter);

      /**
      * Get number of bonds.
      */ 
      int nBond() const;

      /**
      * Access a bond by id.
      */
      Group<2>& bond(int i);
   
      /**
      * Initialize a bond iterator.
      */
      void initBondIterator(ArrayIterator< Group<2> >& iter);

   private:
     
      /// Array of atom objects, added in order read from file.
      DSArray<Atom> atoms_;

      /// Array of bond objects, added in order read from file.
      DSArray< Group<2> > bonds_;

      /// Pointers to atoms indexed by ids. Missing atoms are null pointers.
      DArray<Atom*> atomPtrs_;

      /// Boundary object defines periodic boundary conditions.
      Boundary boundary_;

      /// Pointer to current ConfigIo object.
      ConfigIo* configIoPtr_;

      /// Factory for generating ConfigIo at run time.
      ConfigIoFactory configIoFactory_;

      /// Manager for analyzers
      AnalyzerManager analyzerManager_;

      /// Maximum allowed atom id + 1 (used to allocate arrays).
      int atomCapacity_;

      /// Maximum allowed bond id + 1 (used to allocate arrays).
      int bondCapacity_;

      /// etc. for angles and dihedrals

      /// String identifier for ConfigIo class name
      std::string configIoName_;

      /// Name of configuration or trajectory input file
      std::string configFileName_;

   };

   /*
   * Return pointer to location for new atom, and add to container.
   *
   * \param  global id for new atom
   * \return pointer to location of new atom
   */
   inline Atom* Processor::addAtom(int id)
   {
      int size = atoms_.size();
      atoms_.resize(size + 1);
      Atom* ptr = &atoms_[size];
      atomPtrs_[id] = ptr;
      return ptr;
   }

   /*
   * Return pointer to location for new bond, and add to container.
   *
   * \return pointer to location of new bond
   */
   inline Group<2>* Processor::addBond()
   {
      int size = bonds_.size();
      bonds_.resize(size + 1);
      return &bonds_[size];
   }

   inline Boundary& Processor::boundary() 
   {  return boundary_; }
  
   inline int Processor::nAtom() const
   {  return atoms_.size(); }
  
   inline Atom& Processor::atom(int id)
   {  return *atomPtrs_[id]; }
  
   inline int Processor::nBond() const
   {  return bonds_.size(); }
   
}
#endif
