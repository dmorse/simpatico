namespace MdPp 
{

   using namespace Util;

   class Processor : public ParamComposite 
   {

      /**
      * Read capacities used to allocate arrays and other
      * info that is not in the input file. 
      */
      readParam(std::istream& in);

      /**
      * Analyze a sequence of dump files.
      */
      analyzeDumps(std::string& filename);

      /**
      * Analyze a trajectory file.
      */
      analyzeTrajectory(std::string& filename);

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
      Bond* addBond();

      /// Access a bond by id.
      Group<2>& bond(int i);
   
      etc.
  
      // Accessors, for use in Analyzer and ConfigIo classes.

      /**
      * Get number of atoms.
      */ 
      int nAtom();

      /**
      * Access to atom by global id.
      */
      const Atom& atom(int id) const;

      /**
      * Initialize an iterator for atoms.
      */
      void initAtomIterator(ArrayIterator<Atom>& iter);

      /**
      * Initialize a bond iterator.
      */
      void initBondIterator(ArrayIterator<Bond>& iter);

   private:
     
      // Array of atoms, added in order read from file.
      DSArray<Atom> atoms_;

      // Array of bonds, added in order read from file.
      DSArray<Group<2>> bonds_;

      //etc.

      /// Pointers to atom indexed by id. Missing atoms are null pointers.
      DArray<Atom*> atomPtrs_;

      // Boundary object defines periodic boundary conditions.
      Util::Boundary boundary_;

      ConfigIo* configIoPtr_;

      ConfigIoFactory configIoFactoryPtr_;

      AnalyzerManager analyzerManager_;

      /// Maximum allowed atom id + 1 (used to allocate arrays).
      int atomCapacity_;

      /// Maximum allowed bond id + 1 (used to allocate arrays).
      int bondCapacity_;

   };

   /*
   * Return pointer to location for new atom, and add to container.
   *
   * \param  global id for new atom
   * \return pointer to location of new atom
   */
   inline Atom* addAtom(int id)
   {
      int size = atoms_.size();
      atoms_.resize(size + 1);
      atomPtrs_[id] = size;
      return &atoms_[size];
   }

   inline int nAtom() const
   {  return atoms_.size(); }
   
   inline int nBond() const
   {  return bonds_.size(); }
   
}

