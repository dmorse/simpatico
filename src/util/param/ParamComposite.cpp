#ifndef UTIL_PARAM_COMPOSITE_CPP
#define UTIL_PARAM_COMPOSITE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <cstdio>
#include <cstring>

#include "ParamComposite.h"   // class header
#include "Begin.h"
#include "End.h"
#include "Blank.h"

namespace Util
{

   /* 
   * Default Constructor.
   */
   ParamComposite::ParamComposite() 
    : ParamComponent(),
      list_(),
      isLeaf_(),
      size_(0),
      className_("ParamComposite")
   {}

   /*
   * Constructor.
   */
   ParamComposite::ParamComposite(int capacity)
    : ParamComponent(),
      list_(),
      isLeaf_(),
      size_(0),
      className_("ParamComposite")
   {
      if (capacity <= 0 ) {
         UTIL_THROW("Attempt to reserve capacity <= 0");
      }
      list_.reserve(capacity);
      isLeaf_.reserve(capacity);
   }
   
   /* 
   * Copy constructor.
   */
   ParamComposite::ParamComposite(const ParamComposite& other) 
    : ParamComponent(other),
      list_(),
      isLeaf_(),
      size_(0),
      className_(other.className_)
   {}

   /* 
   * Destructor.
   */
   ParamComposite::~ParamComposite()
   {
      if (size_ > 0) { 
         for (int i=0; i < size_; ++i) {
            if (isLeaf_[i]) {
               delete list_[i];
            }
            /* Only delete Param, Begin, End & Blank leaf objects.
            * Do NOT delete node objects here. These are instances of
            * of subclasses of ParamComposite, which must be deleted
            * by their parent object, i.e., by the object of which they 
            * are a member or, if created on the heap, the object that
            * has sole responsibility for deleting them.
            */
         }
      }
   }

   /*
   * Read parameter block, including begin and end.
   */
   void ParamComposite::readParam(std::istream &in)
   {
      readBegin(in, className().c_str());
      readParameters(in);
      readEnd(in);
   }
   
   /* 
   * Default writeParam implementation.
   */
   void ParamComposite::writeParam(std::ostream &out) 
   {
      for (int i=0; i < size_; ++i) {
         list_[i]->writeParam(out);
      }
   }
   
   /* 
   * Reset list to empty state.
   */
   void ParamComposite::resetParam() 
   {
      for (int i=0; i < size_; ++i) {
         if (isLeaf_[i]) {
            delete list_[i];
         } else {
            list_[i]->resetParam();
         }
      }
      size_ = 0;
   }
   
   // ParamComposite object
   
   /* 
   * Add a ParamComposite node to the tree.
   */
   void ParamComposite::addParamComposite(ParamComposite &child, bool next) 
   {
      child.setIndent(*this, next);
      list_.push_back(&child);
      isLeaf_.push_back(false);
      ++size_;
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         child.setParamCommunicator(paramCommunicator());
      }
      #endif
   }
   
   /* 
   * Add a ParamComposite Node, and read the contents of that ParamComposite.
   */
   void 
   ParamComposite::readParamComposite(std::istream &in, ParamComposite &child, bool next) 
   {
      addParamComposite(child, next);
      list_.back()->readParam(in);
   }
 
   // Begin
   
   /* 
   * Add a new Begin object.
   */
   Begin& ParamComposite::addBegin(const char *label) 
   {
      Begin* ptr = new Begin(label);
      list_.push_back(ptr);
      isLeaf_.push_back(true);
      ptr->setIndent(*this, false);
      ++size_;
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         ptr->setParamCommunicator(paramCommunicator());
      }
      #endif
      return *ptr;
   }
   
   /* 
   * Read the opening line of a ParamComposite.
   */
   Begin& ParamComposite::readBegin(std::istream &in, const char *label) 
   {
      Begin* ptr = &addBegin(label);
      ptr->readParam(in);
      return *ptr;
   }
   
   // End
   
   /* 
   * Add a new End object.
   */
   End& ParamComposite::addEnd() 
   {
      End* ptr = new End();
      list_.push_back(ptr);
      isLeaf_.push_back(true);
      ptr->setIndent(*this, false);
      ++size_;
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         ptr->setParamCommunicator(paramCommunicator());
      }
      #endif
      return *ptr;
   }
   
   /* 
   * Read the closing bracket of a ParamComposite.
   */
   End& ParamComposite::readEnd(std::istream &in) 
   {
      End* ptr = &addEnd();
      ptr->readParam(in);
      return *ptr;
   }

   // Blank
   
   /* 
   * Add a Blank object (a blank line).
   */
   Blank& ParamComposite::addBlank() 
   {
      Blank* ptr = new Blank();
      list_.push_back(ptr);
      isLeaf_.push_back(true);
      ++size_;
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         ptr->setParamCommunicator(paramCommunicator());
      }
      #endif
      return *ptr;
   }
   
   /* 
   * Read a blank line.
   */
   Blank& ParamComposite::readBlank(std::istream &in) 
   {
      Blank* ptr = &addBlank();
      ptr->readParam(in);
      return *ptr;
   }

   /*
   * Set class name string.
   */
   void ParamComposite::setClassName(const char * className)
   {  className_ = className; }
      

} 
#endif
