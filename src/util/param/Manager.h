#ifndef UTIL_MANAGER_H
#define UTIL_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class
#include <util/param/Factory.h>         // pointer member, used in implementation
#include <util/global.h>

#include <string>
#include <vector>

namespace Util
{

   /**
   * Template container for pointers to objects with a common base class.
   *
   * A Manager<Data> has an array of Data* pointers to Data objects, an array
   * of corresponding subclass names, and a pointer to a Factory<Data> object.
   * The default implementation of the Manager<Data>::readParam() method uses
   * an associated Factory<Data> object to recognize the class name string
   * that begins a polymorphic block in a parameter file (which must refer to
   * a known subclass of Data) and to instantiate an object of the specified
   * subclass.
   *
   * Subclasses of Manager<Data> are used to manage arrays of Species, McMove,
   * and Diagnostic objects.
   *
   * \ingroup Manager_Module
   */
   template <typename Data>
   class Manager : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Manager();

      /**
      * Destructor.
      */
      virtual ~Manager();

      /**
      * Set a SubFactory for this Manager.
      *
      * \param subfactory Referernce to a sub-Factory to be added.
      */
      void addSubfactory(Factory<Data>& subfactory);

      /**
      * Associate a Factory with this Manager.
      *
      * \param factory Reference to a Factory object
      */
      void setFactory(Factory<Data>& factory);

      /**
      * Associated a Factory with this Manager (pass by pointer).
      *
      * \param factoryPtr pointer to a Factory<Data> object.
      */
      void setFactory(Factory<Data>* factoryPtr);

      /**
      * Read and create a set of objects.
      *
      * The default implementation of this method reads a sequence of blocks
      * for different subclasses of Data, terminated by a closing bracket.
      * For each block it:
      *
      *  - reads a className string for a subclass of Data,
      *
      *  - uses factory object to create a new instance of className.
      *
      *  - invokes the readParam() method of the new object.
      *
      * The implementation of the factory must recognize all valid
      * className string values, and invoke the appropriate constructor for
      * each.  The loop over blocks terminates when it encounters a closing
      * bracket '}' surrounded by white space.
      *
      * \param in input stream
      */
      virtual void readParam(std::istream &in);

      /**
      * Load a set of objects to an output archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save a set of objects to an output archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Append a Data object to the end of the sequence.
      *
      * \param data Data object to be appended
      * \param name subclass name
      */
      void append(Data& data, const std::string& name);

      /**
      * Get logical size.
      *
      * \return logical size of this array.
      */
      int size() const;

      /**
      * Get the subclass name for object number i.
      *
      * \param   i integer index of object
      * \return  class name of managed object
      */
      std::string className(int i) const;

      /**
      * Return a reference to the factory.
      */
      Factory<Data>& factory();

      /**
      * Return true if this Manager has a Factory, false otherwise.
      */
      bool hasFactory() const;

      /**
      * Mimic C array subscripting.
      *
      * \param  i array index
      * \return reference to element i
      */
      Data& operator[] (int i) const;

   protected:

      /**
      * Read opening line: "ManagerName{"
      */
      void beginReadManager(std::istream& in);

      /**
      * Read child blocks, return when closing bracket encountered.
      */ 
      void readChildren(std::istream &in);

      /**
      * Add closing bracket to output format.
      */
      void endReadManager();

      /**
      * Create factory if necessary.
      */
      void initFactory();

      /**
      * Create an instance of the default Factory<Data> class.
      *
      * \return a pointer to a new Factory<Data> object.
      */
      virtual Factory<Data>* newDefaultFactory() const = 0;

   private:

      /// Array of pointers to Data objects.
      std::vector<Data*> ptrs_;

      /// Array of subclass names for Data objects.
      std::vector<std::string> names_;

      /// Pointer to an associated Factory<Data> object
      Factory<Data>* factoryPtr_;

      /// Allocated size of ptrs_ array.
      int  capacity_;

      /// Logical size (number of elements with initialized data).
      int  size_;

      /// True if this manager knows its own classname
      bool  hasName_;

      /// True if this manager created the object *factoryPtr_.
      bool createdFactory_;

   };

   /*
   * Constructor.
   */
   template <typename Data>
   Manager<Data>::Manager()
    : ptrs_(),
      names_(),
      factoryPtr_(0),
      capacity_(0),
      size_(0),
      hasName_(false),
      createdFactory_(false)
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   Manager<Data>::~Manager()
   {
      // Delete objects referenced by non-null elements of ptrs_[].
      for (int i=0; i < size_; ++i) {
         if (!(ptrs_[i] == 0)) {
            delete ptrs_[i];
         }
      }

      // Delete Factory if it was created by this Manager.
      if (factoryPtr_ && createdFactory_) {
         delete factoryPtr_;
      }
   }

   /*
   * Return a reference to the factory.
   */
   template <typename Data>
   Factory<Data>& Manager<Data>::factory()
   {
      initFactory();
      return *factoryPtr_;
   }

   /*
   * Add a subfactory to the current (or default) Factory.
   */
   template <typename Data>
   void Manager<Data>::addSubfactory(Factory<Data>& subfactory)
   {
      factory().addSubfactory(subfactory);
   }

   /*
   * Associate a Factory with the Manager, pass by reference.
   */
   template <typename Data>
   void Manager<Data>::setFactory(Factory<Data>& factory)
   {
      if (factoryPtr_ && createdFactory_) {
         delete factoryPtr_;
         createdFactory_ = false;
      }
      factoryPtr_ = &factory;
   }

   /*
   * Associate a factory with the Manager, passed by pointer.
   */
   template <typename Data>
   void Manager<Data>::setFactory(Factory<Data>* factoryPtr)
   {
      if (factoryPtr_ && createdFactory_) {
         delete factoryPtr_;
         createdFactory_ = false;
      }
      factoryPtr_ = factoryPtr;
   }

   /*
   * Read instructions for creating objects from file.
   */
   template <typename Data>
   void Manager<Data>::readParam(std::istream &in)
   {
      beginReadManager(in);
      readChildren(in);
      endReadManager();
   }

   /*
   * Read instructions for creating objects from file.
   */
   template <typename Data>
   void Manager<Data>::beginReadManager(std::istream &in)
   {
      std::string managerName = ParamComposite::className();
      readBegin(in, managerName.c_str());
   }

   /*
   * Read instructions for creating objects from file.
   */
   template <typename Data>
   void Manager<Data>::readChildren(std::istream &in)
   {
      // Check if a Factory exists, create one if necessary.
      initFactory();

      // Loop over managed objects
      std::string name;
      Data* typePtr;
      bool  isEnd = false;
      while (!isEnd) {

         // Read a blank line before each object
         readBlank(in);

         // Read and instantiate a new object *typePtr
         typePtr = factoryPtr_->readObject(in, *this, name, isEnd);

         if (!isEnd) {
            if (typePtr) {
               append(*typePtr, name);
            } else {
               std::string msg("Unknown subclass name: ");
               msg += name;
               UTIL_THROW(msg.c_str());
            }
         }

      }
   }
   
   /*
   * Add closing bracket.
   */
   template <typename Data>
   void Manager<Data>::endReadManager()
   {
      End* endPtr = &addEnd();
      if (ParamComponent::echo() && isParamIoProcessor()) {
         endPtr->writeParam(Log::file());
      }
   }

   /*
   * Read instructions for creating objects from file.
   */
   template <typename Data>
   void Manager<Data>::loadParameters(Serializable::IArchive &ar)
   {
      int size;
      Data* typePtr;
      std::string name;

      // Check if a Factory exists, create one if necessary.
      initFactory();

      ar >> size;
      for (int i = 0; i < size; ++i) {
         addBlank();
         name = "unknown";
         typePtr = factoryPtr_->loadObject(ar, *this, name);
         if (typePtr) {
            append(*typePtr, name);
         } else {
            std::string msg("Unknown subclass name: ");
            msg += name;
            UTIL_THROW(msg.c_str());
         }
      }
      addBlank();
   }

   /*
   * Read instructions for creating objects from file.
   */
   template <typename Data>
   void Manager<Data>::save(Serializable::OArchive &ar)
   {
      int size = ptrs_.size();
      ar << size;
      for (int i = 0; i < size; ++i) {
         ar << names_[i];
         (ptrs_[i])->save(ar);
      }
   }

   /*
   * Append an object to the Manager.
   */
   template <typename Data>
   void Manager<Data>::append(Data& data, const std::string& name)
   {
      ptrs_.push_back(&data);
      names_.push_back(name);
      ++size_;
   }

   /*
   * Return subclass name for object number i.
   */
   template <typename Data>
   std::string Manager<Data>::className(int i) const
   {
      if (i >= size_) {
         UTIL_THROW("Invalid array index\n");
      }
      return names_[i];
   }

   /*
   * Does this Manager have an associated Factory?
   */
   template <typename Data>
   bool Manager<Data>::hasFactory() const
   {
      return (factoryPtr_ != 0);
   }

   /*
   * Return logical size (number of objects appended thus far).
   */
   template <typename Data>
   inline int Manager<Data>::size() const
   {  return size_; }

   /*
   * Array subscript operator - return a reference.
   *
   * \param  i array index
   * \return reference to element i
   */
   template <typename Data>
   inline Data& Manager<Data>::operator[] (int i) const
   {
      assert(i >= 0);
      assert(i < size_);
      return *ptrs_[i];
   }

   // Protected methods

   /*
   * Initialize Factory, create if necessary.
   */
   template <typename Data>
   void Manager<Data>::initFactory()
   {
      if (!hasFactory()) {
         factoryPtr_ = newDefaultFactory();
         createdFactory_ = true;
      }
      assert(factoryPtr_);
   }

}
#endif
