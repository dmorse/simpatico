#ifndef MDPP_PROCESSOR_CPP
#define MDPP_PROCESSOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Processor.h"

namespace MdPp 
{

   /*
   * Constructor.
   */
   Processor::Processor()
    : configIoPtr_(0),
      configIoFactory_(*this),
      analyzerManager_(*this),
      newAtomPtr_(0)
   {  setClassName("Processor"); }

   /*
   * Destructor.
   */
   Processor::~Processor()
   {
      if (configIoPtr_) {
         delete configIoPtr_;
      }
   }

   /*
   * Read parameters from file.
   */
   void Processor::readParameters(std::istream& in)
   {
      read<int>(in, "atomCapacity", atomCapacity_); 
      read<int>(in, "bondCapacity", bondCapacity_); 
      // etc. for angles dihedrals

      atoms_.allocate(atomCapacity_);
      atomPtrs_.allocate(atomCapacity_);
      for (int i = 0; i < atomCapacity_; ++i) {
         atomPtrs_[i] = 0;
      }
      bonds_.allocate(bondCapacity_);
      // etc. for angles dihedrals

      readParamComposite(in, analyzerManager_);

      read<std::string>(in, "configIoName", configIoName_);
      configIoPtr_ = configIoFactory_.factory(configIoName_);
      if (configIoPtr_ == 0) {
         UTIL_THROW("Unrecognized ConfigIo subclass name");
      }

      read<std::string>(in, "configFileName", configFileName_);
   }

   /*
   * Return pointer to location for new atom.
   */
   Atom* Processor::newAtomPtr()
   {
      if (newAtomPtr_) {
         UTIL_THROW("Error: an new atom is still active");
      }
      int size = atoms_.size() + 1;
      atoms_.resize(size);
      newAtomPtr_ = &atoms_[size - 1];
      return newAtomPtr_;
   }

   /*
   * Finalize addition of new atom.
   */
   void Processor::addAtom()
   {
      if (!newAtomPtr_) {
         UTIL_THROW("Error: No active new atom");
      }
      int id = newAtomPtr_->id;
      atomPtrs_[id] = newAtomPtr_;
      newAtomPtr_ = 0;
   }

   /*
   * Return pointer to location for new bond, and add to container.
   */
   Group<2>* Processor::newBondPtr()
   {
      int size = bonds_.size();
      bonds_.resize(size + 1);
      return &bonds_[size];
   }
   
}
#endif
