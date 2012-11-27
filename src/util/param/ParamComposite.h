#ifndef UTIL_PARAM_COMPOSITE_H
#define UTIL_PARAM_COMPOSITE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>     // base class
#include <util/param/ScalarParam.h>        // template class
#include <util/param/CArrayParam.h>        // template class
#include <util/param/DArrayParam.h>        // template class
#include <util/param/FArrayParam.h>        // template class
#include <util/param/CArray2DParam.h>      // template class
#include <util/param/DMatrixParam.h>       // template class
#include <util/archives/Serializable_includes.h>   
#include <util/global.h>            

#include <vector>

namespace Util
{

   class Begin;
   class End;
   class Blank;
   template <class Data> class Factory;

   /**
   * An object that can read multiple parameters from file.
   *
   * A ParamComposite has a private array of pointers to ParamComponent 
   * objects, stored in the order in which they are read from file by the
   * readParam() method. We will refer to this as the format array. Each 
   * element of the format array may point to a Parameter object (which 
   * represents a single parameter), a Begin or End object (which represents 
   * a line containing the opening or closing bracket for a parameter block), 
   * a Blank object (i.e., a blank line), or another ParamComposite object. 
   *
   * Any class that reads a block of parameters from the input parameter file 
   * should be derived from ParamComposite. The readParam() method must read 
   * the associated parameter file block. The virtual readParameters() method,
   * if re-implemented, should read only the body of the parameter file
   * block, without opening and closing lines. The default implementation of
   * ParamComposite::readParam() reads the opening line of the block, calls 
   * readParameters() to read the body of the  parameter file block, and then
   * reads the closing line. Almost all subclasses of ParamComposite should 
   * re-implement the readParameters() method, and rely on the default 
   * implementation of ParamComposite::readParam() to add the begin and end 
   * lines. 
   *
   * The setClassName() and className() methods are used to set and get a 
   * string representing the a subclass name. The setClassName() method 
   * should be called in the constructor of each subclass of ParamComposite. 
   * The name string set in the constructor of a subclass will replace any
   * name set by a base class, because of the order in which constructors
   * are called. The class name string is used by the default implementation 
   * of ParamComposite::readParam() in order to check that the opening line 
   * of a parameter block contains the correct class name.
   * 
   * The readParameters() method of each subclass should be implemented 
   * using protected methods provided by ParamComposite to read individual 
   * parameters and "child" ParamComposite objects.  The implementation of 
   * readParameters() normally uses read< T >, which reads an individual 
   * parameter, readParamComposite, which reads a nested subblock, readBlank,
   * which reads a blank line, and more specialized methods to read arrays
   * and matrices of parameters.  Each of these methods creates a new 
   * ParamComponent of a specified type, adds a pointer to the new object 
   * to the format array, and invokes the readParam method of the new object 
   * in order to read the associated line or block of the parameter file. 
   *
   * The ParamComposite::writeParam() method uses the format array to write
   * data to a file in the same format in which it was read by a previous
   * call to readParam(). 
   *
   * \ingroup Param_Module
   */
   class ParamComposite : public ParamComponent
   {
   
   public:
  
      /** 
      * Constructor
      */
      ParamComposite();
  
      /** 
      * Copy constructor
      */
      ParamComposite(const ParamComposite& other);
  
      /** 
      * Constructor.
      *
      * Reserve space for capacity elements in the format array.
      *
      * \param capacity maximum length of parameter list
      */
      ParamComposite(int capacity);
  
      /** 
      * Virtual destructor.
      */
      virtual ~ParamComposite();
   
      /// \name Initialization methods
      //@{

      /** 
      * Read a parameter file block. 
      *
      * Inherited from ParamComponent. This method reads the entire parameter 
      * block, including the opening line "ClassName{" and the closing bracket 
      * "}". The default implementation calls the virtual readParameters 
      * method to read the body of the block, and adds Begin and End objects.
      *
      * \param in input stream for reading
      */
      virtual void readParam(std::istream &in);
   
      /** 
      * Read body of parameter block, without opening and closing lines.
      *
      * Most subclasses of ParamComposite should overload readParameters.
      * The default implementation is empty. All subclasses must overload
      * either readParameters or readParam, but not both. 
      *
      * \param in input stream for reading
      */
      virtual void readParameters(std::istream &in)
      {};

      /** 
      * Write all parameters to an output stream.
      *
      * The default implementation calls the readParam method of each
      * child in the format array.
      *
      * \param out output stream for reading
      */
      virtual void writeParam(std::ostream &out);
   
      /** 
      * Load all parameters from an archive.
      *
      * This method is inherited from Serializable. The default
      * implementation for a ParamComposite calls loadParameters, and 
      * adds Begin and End lines. All subclasses of ParamComposite 
      * should overload loadParameters.
      *
      * \param ar input/loading archive.
      */
      virtual void load(Serializable::IArchive &ar);
   
      /** 
      * Load state from archive, without adding Begin and End lines.
      *
      * This default implementation is empty, and should be overloaded
      * by all subclasses. The overloaded method should load the entire
      * internal state from the archive, including members that do not
      * appear in the parameter file format.
      *
      * \param ar input/loading archive.
      */
      virtual void loadParameters(Serializable::IArchive &ar)
      {};

      /** 
      * Saves all parameters to an archive.
      *
      * This default implementation calls the save method for all
      * items in the parameter file format array. This is not sufficient
      * for classes that contain non-volatile members that do not
      * appear in the parameter file format. 
      *
      * If a class also defines a serialize method template, which allows 
      * instances to be serialized to any type of archive, then the save 
      * method can be implemented as follows:
      * \code
      *    void save(Serializable::OArchive& ar)
      *    { ar & *this; }
      * \endcode
      *
      * \param ar output/saving archive.
      */
      virtual void save(Serializable::OArchive &ar);
   
      /**
      * Resets ParamComposite to its empty state.
      *
      * This method deletes Parameter, Begin, End, and Blank objects in the
      * format array, recursively invokes resetParam() for any ParamComposite 
      * objects in the format array, nullifies all pointers in the list, and 
      * resets the number of items in the array to 0. 
      */
      void resetParam();
   
      #ifdef UTIL_MPI
      /**
      * Set an MPI communicator for parameter IO.
      *
      * This method recursively sets the ParamCommunicator for all children.
      */
      virtual void setParamCommunicator(MPI::Intracomm& communicator);
      #endif

      //@}
      /// \name read* methods
      /// \brief Each of these methods invokes an associated add* method to 
      /// create a new ParamComponent object, and then invoke the readParam() 
      /// method of the new object to read the associated line or block of a 
      /// parameter file.
      //@{
      
      /** 
      * Add and read a child ParamComposite.
      *
      * \param in    input stream for reading
      * \param child child ParamComposite object
      * \param next  true if the indent level is one higher than parent.
      */
      void 
      readParamComposite(std::istream &in, ParamComposite &child, 
                                           bool next = true);
   
      /**  
      * Add a new Param < Type > object, and read its value.
      *
      * \param in     input stream for reading
      * \param label  Label string 
      * \param value  reference to new ScalarParam< Type >
      */
      template <typename Type>
      ScalarParam<Type>& read(std::istream &in, const char *label, Type &value);
   
      /**  
      * Add a C array parameter, and read its elements.
      *
      * \param in     input stream for reading
      * \param label  Label string for new array
      * \param value  pointer to array
      * \param n      number of elements
      * \return reference to the new CArrayParam<Type> object
      */
      template <typename Type>
      CArrayParam<Type>& 
      readCArray(std::istream &in, const char *label, Type *value, int n);
   
      /**  
      * Add a DArray < Type > parameter, and read its elements.
      *
      * \param in     input stream for reading
      * \param label  Label string for new array
      * \param array  DArray object
      * \param n      number of elements
      * \return reference to the new DArrayParam<Type> object
      */
      template <typename Type>
      DArrayParam<Type>&
      readDArray(std::istream &in, const char *label, DArray<Type>& array, int n);
   
      /**  
      * Add and read an FArray < Type, N > fixed-size array parameter.
      *
      * \param in     input stream for reading
      * \param label  Label string for new array
      * \param array  FArray object
      * \return reference to the new FArrayParam<Type, N> object
      */
      template <typename Type, int N>
      FArrayParam<Type, N>&
      readFArray(std::istream &in, const char *label, FArray<Type, N >& array);
   
      /**  
      * Add and read a CArray2DParam < Type > C two-dimensional array parameter.
      *
      * \param in     input stream for reading
      * \param label  Label string for new array
      * \param value  pointer to array
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      * \return reference to the CArray2DParam<Type> object
      */
      template <typename Type> 
      CArray2DParam<Type>&
      readCArray2D(std::istream &in, const char *label, 
                        Type *value, int m, int n);
  
      /**  
      * Add and read a DMatrix < Type > C two-dimensional matrix parameter.
      *
      * \param in     input stream for reading
      * \param label  Label string for new array
      * \param matrix DMatrix object
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      * \return reference to the DMatrixParam<Type> object
      */
      template <typename Type> 
      DMatrixParam<Type>&
      readDMatrix(std::istream &in, const char *label, 
                       DMatrix<Type>& matrix, int m, int n);
  
      /** 
      * Add and read a class label and opening bracket. 
      *
      * \param in    input stream for reading
      * \param label class name string, without trailing bracket
      * \return reference to the new Begin object
      */
      Begin& readBegin(std::istream &in, const char* label);
   
      /** 
      * Add and read the closing bracket.
      *
      * \param in input stream for reading
      * \return reference to the new End object
      */
      End& readEnd(std::istream &in);
   
      /** 
      * Add and read a new Blank object, representing a blank line.
      *
      * \param in input stream for reading
      * \return reference to the new Blank object
      */
      Blank& readBlank(std::istream &in);

      //@}
      /// \name load* methods   
      /// \brief Each of these methods invokes an associated add* method to create
      /// a new ParamComponent object, and then invokes the load() method of the 
      /// new object to load the associated parameter value from an archive.
      //@{
     
      /** 
      * Add and load a child ParamComposite.
      *
      * \param ar    input archive for loading
      * \param child child ParamComposite object
      * \param next  true if the indent level is one higher than parent.
      */
      void 
      loadParamComposite(Serializable::IArchive &ar, ParamComposite &child, 
                         bool next = true);
   
      /**  
      * Add a new Param < Type > object, and load its value.
      *
      * \param ar     archive for loading
      * \param label  Label string 
      * \param value  reference to new ScalarParam< Type >
      */
      template <typename Type>
      ScalarParam<Type>& loadParameter(Serializable::IArchive &ar, const char *label, Type &value);
   
      /**  
      * Add a C array parameter, and load its elements.
      *
      * \param ar     archive for loading
      * \param label  Label string for new array
      * \param value  pointer to array
      * \param n      number of elements
      * \return reference to the new CArrayParam<Type> object
      */
      template <typename Type>
      CArrayParam<Type>& 
      loadCArray(Serializable::IArchive &ar, const char *label, Type *value, int n);
   
      /**  
      * Add a DArray < Type > parameter, and load its elements.
      *
      * \param ar     archive for loading
      * \param label  Label string for new array
      * \param array  DArray object
      * \param n      number of elements
      * \return reference to the new DArrayParam<Type> object
      */
      template <typename Type>
      DArrayParam<Type>&
      loadDArray(Serializable::IArchive &ar, const char *label, DArray<Type>& array, int n);
   
      /**  
      * Add and load an FArray < Type, N > fixed-size array parameter.
      *
      * \param ar     archive for loading
      * \param label  Label string for new array
      * \param array  FArray object
      * \return reference to the new FArrayParam<Type, N> object
      */
      template <typename Type, int N>
      FArrayParam<Type, N>&
      loadFArray(Serializable::IArchive &ar, const char *label, FArray<Type, N >& array);
   
      /**  
      * Add and load a CArray2DParam < Type > C two-dimensional array parameter.
      *
      * \param ar     archive for loading
      * \param label  Label string for new array
      * \param value  pointer to array
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      * \return reference to the CArray2DParam<Type> object
      */
      template <typename Type> 
      CArray2DParam<Type>&
      loadCArray2D(Serializable::IArchive &ar, const char *label, 
                   Type *value, int m, int n);
  
      /**  
      * Add and load a DMatrix < Type > C two-dimensional matrix parameter.
      *
      * \param ar     archive for loading
      * \param label  Label string for new array
      * \param matrix DMatrix object
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      * \return reference to the DMatrixParam<Type> object
      */
      template <typename Type> 
      DMatrixParam<Type>& 
      loadDMatrix(Serializable::IArchive &ar, const char *label, 
                  DMatrix<Type>& matrix, int m, int n);
  
      //@}
      /// \name add* methods
      /// \brief These methods add a ParamComponent object to the format
      /// array, but do not read any data from an input stream.
      //@{
   
      /**
      * Add a child ParamComposite object to the format array.
      *
      * \param child child ParamComposite object
      * \param next  true if the indent level is one higher than parent.
      */
      void addParamComposite(ParamComposite& child, bool next = true);
   
      /**  
      * Add a new Param < Type > object to the format array.
      *
      * \param label  Label string for new ScalarParam < Type >
      * \param value  reference to parameter value
      * \return reference to the new Param<Type> object
      */
      template <typename Type>
      ScalarParam<Type>& add(const char *label, Type &value);
   
      /**  
      * Add a C array parameter.
      *
      * \param label  Label string for new array
      * \param value  pointer to array
      * \param n      number of elements
      * \return reference to the new CArrayParam<Type> object
      */
      template <typename Type>
      CArrayParam<Type>& addCArray(const char *label, Type *value, int n);
   
      /**  
      * Add a DArray<Type> parameter.
      *
      * \param label  Label string for new array
      * \param array  DArray object.
      * \param n      number of elements
      * \return reference to the new DArrayParam<Type> object
      */
      template <typename Type>
      DArrayParam<Type>& addDArray(const char *label, DArray<Type>& array, int n);
   
      /**  
      * Add a FArray<Type> parameter.
      *
      * \param label  Label string for new array
      * \param array  FArray object.
      * \return reference to the new FArrayParam<Type, N> object
      */
      template <typename Type, int N>
      FArrayParam<Type, N>& 
      addFArray(const char *label, FArray<Type, N>& array);
   
      /**  
      * Add a C two-dimensional array parameter.
      *
      * \param label  Label string for new array
      * \param value  pointer to 2D array
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      * \return reference to the new CArray2DParam<Type> object
      */
      template <typename Type>
      CArray2DParam<Type>&
      addCArray2D(const char *label, Type *value, int m, int n);
  
      /**  
      * Add a C two-dimensional array parameter.
      *
      * \param label  Label string for new array
      * \param matrix  pointer to arr
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      * \return reference to the new DMatrixParam<Type> object
      */
      template <typename Type>
      DMatrixParam<Type>&
      addDMatrix(const char *label, DMatrix<Type>& matrix, int m, int n);

      /** 
      * Add a class label and opening bracket. 
      *
      * \param label class name string, without trailing bracket
      * \return reference to the new begin object.
      */
      Begin& addBegin(const char* label);
   
      /** 
      * Add a closing bracket.
      *
      * \return reference to the new End object.
      */
      End& addEnd();
   
      /** 
      * Add a new Blank object, representing a blank line.
      *
      * \return reference to the new Blank object
      */
      Blank& addBlank();
   
      //@}
  
      /**
      * Get class name string.
      */
      std::string className() const;
      
   protected:

      /**
      * Set class name string.
      *
      * Should be set in subclass constructor. 
      */
      void setClassName(const char* className);
      
   private:
   
      /// Array of pointers to ParamComponent objects.
      std::vector<ParamComponent*> list_;    

      /// Array of booleans, elements false for ParamComposite, true otherwise.
      std::vector<bool> isLeaf_;  

      /// Number of ParamComponent objects in format array
      int size_;     

      /// Name of subclass.
      std::string className_;

   };

   // add and read method templates for scalar parameters
 
   /* 
   * Add a ScalarParam to the format array.
   */ 
   template <typename Type>
   ScalarParam<Type>& ParamComposite::add(const char *label, Type &value) 
   {
      ScalarParam<Type>* ptr = new ScalarParam<Type>(label, value);
      list_.push_back(ptr);
      isLeaf_.push_back(true);
      ++size_;
      ptr->setIndent(*this);
      #ifdef UTIL_MPI 
      if (hasParamCommunicator()) {
         ptr->setParamCommunicator(paramCommunicator());
      }
      #endif
      return *ptr;
   }
   
   /* 
   * Add a ScalarParam to the format array, and read its contents from file. 
   */ 
   template <typename Type>
   ScalarParam<Type>& 
   ParamComposite::read(std::istream &in, const char *label, Type &value) 
   {
      ScalarParam<Type>* ptr = &add<Type>(label, value);
      ptr->readParam(in);
      return *ptr;
   }
  
   /*  
   * Add a new ScalarParam< Type > object, and load its value.
   */
   template <typename Type>
   ScalarParam<Type>& 
   ParamComposite::loadParameter(Serializable::IArchive &ar, const char *label, 
                                 Type &value)
   {
      ScalarParam<Type>* ptr = &add<Type>(label, value);
      ptr->load(ar);
      return *ptr;
   }
   
   // Templates for 1D C Arrays 

   /* 
   * Add a CArrayParam object associated with a built-in array.
   */
   template <typename Type>   
   CArrayParam<Type>& 
   ParamComposite::addCArray(const char *label, Type *value, int n) 
   {
      CArrayParam<Type>* ptr = new CArrayParam<Type>(label, value, n);
      list_.push_back(ptr);
      isLeaf_.push_back(true);
      ++size_;
      ptr->setIndent(*this);
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         ptr->setParamCommunicator(paramCommunicator());
      }
      #endif
      return *ptr;
   }
   
   /* 
   * Add a a CArrayParam associated with a built-in array, and read it.
   */
   template <typename Type>   
   CArrayParam<Type>& 
   ParamComposite::readCArray(std::istream &in, 
                   const char *label, Type *value, int n) 
   {
      CArrayParam<Type>* ptr = &addCArray<Type>(label, value, n);
      ptr->readParam(in);
      return *ptr;
   }
   
   /*  
   * Add a C array parameter, and load its elements.
   */
   template <typename Type>
   CArrayParam<Type>& 
   ParamComposite::loadCArray(Serializable::IArchive &ar, const char *label, 
                              Type *value, int n)
   {
      CArrayParam<Type>* ptr = &addCArray<Type>(label, value, n);
      ptr->load(ar);
      return *ptr;
   }

   // Templates for DArray parameters
   
   /* 
   * Add a DArrayParam < Type > object associated with a DArray<Type> container.
   */
   template <typename Type>   
   DArrayParam<Type>&
   ParamComposite::addDArray(const char *label, DArray<Type>& array, int n) 
   {
      DArrayParam<Type>* ptr = new DArrayParam<Type>(label, array, n);
      list_.push_back(ptr);
      isLeaf_.push_back(true);
      ++size_;
      ptr->setIndent(*this);
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         ptr->setParamCommunicator(paramCommunicator());
      }
      #endif
      return *ptr;
   }
   
   /* 
   * Add a a DArrayParam associated with a DArray<Type> container, and read it.
   */
   template <typename Type>   
   DArrayParam<Type>& 
   ParamComposite::readDArray(std::istream &in, 
                              const char *label, DArray<Type>& array, int n) 
   {
      DArrayParam<Type>* ptr = &addDArray<Type>(label, array, n);
      ptr->readParam(in);
      return *ptr;
   }

   /*
   * Add a DArray < Type > parameter, and load its elements.
   */
   template <typename Type>
   DArrayParam<Type>&
   ParamComposite::loadDArray(Serializable::IArchive &ar, const char *label, 
                              DArray<Type>& array, int n)
   {
      DArrayParam<Type>* ptr = &addDArray<Type>(label, array, n);
      ptr->load(ar);
      return *ptr;
   }

   // Templates for fixed size FArray array objects.

   /* 
   * Add a FArrayParam <Type, N> object.
   */
   template <typename Type, int N> 
   FArrayParam<Type, N>&
   ParamComposite::addFArray(const char *label, FArray<Type, N>& array) 
   {
      FArrayParam<Type, N>* ptr = new FArrayParam<Type, N>(label, array);
      list_.push_back(ptr);
      isLeaf_.push_back(true);
      ++size_;
      ptr->setIndent(*this);
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         ptr->setParamCommunicator(paramCommunicator());
      }
      #endif
      return *ptr;
   }
   
   /* 
   * Add and read an FArray<Type, N> fixed size array.
   */
   template <typename Type, int N>   
   FArrayParam<Type, N>&
   ParamComposite::readFArray(std::istream &in, 
                   const char *label, FArray<Type, N>& array) 
   {
      FArrayParam<Type, N>* ptr = &addFArray<Type, N>(label, array);
      ptr->readParam(in);
      return *ptr;
   }

   /*
   * Add and load an FArray < Type, N > fixed-size array parameter.
   */
   template <typename Type, int N>
   FArrayParam<Type, N>&
   ParamComposite::loadFArray(Serializable::IArchive &ar, const char *label, 
                              FArray<Type, N >& array)
   {
      FArrayParam<Type, N>* ptr = &addFArray<Type, N>(label, array);
      ptr->load(ar);
      return *ptr;
   }

   // Templates for built-in two-dimensional C Arrays

   /* 
   * Add a CArray2Param associated with a built-in 2D array.
   */
   template <typename Type>   
   CArray2DParam<Type>&
   ParamComposite::addCArray2D(const char *label, Type *value, int m, int n) 
   {
      CArray2DParam<Type>* ptr = new CArray2DParam<Type>(label, value, m, n);
      list_.push_back(ptr);
      isLeaf_.push_back(true);
      ++size_;
      ptr->setIndent(*this);
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         ptr->setParamCommunicator(paramCommunicator());
      }
      #endif
      return *ptr;
   }
   
   /* 
   * Add and read a CArray2DParam associated with a built-in 2D array.
   */
   template <typename Type>   
   CArray2DParam<Type>&
   ParamComposite::readCArray2D(std::istream &in, 
                   const char *label, Type *value, int m, int n) 
   {
      CArray2DParam<Type>* ptr = &addCArray2D<Type>(label, value, m, n);
      ptr->readParam(in);
      return *ptr;
   }
 
   /*
   * Add and load a CArray2DParam < Type > C two-dimensional array parameter.
   */
   template <typename Type> 
   CArray2DParam<Type>&
   ParamComposite::loadCArray2D(Serializable::IArchive &ar, const char *label, 
                Type *value, int m, int n)
   {
      CArray2DParam<Type>* ptr = &addCArray2D<Type>(label, value, m, n);
      ptr->load(ar);
      return *ptr;
   }
  
   // Templates for DMatrix containers

   /* 
   * Add a DMatrixParam associated with a DMatrix 2D array..
   */
   template <typename Type>   
   DMatrixParam<Type>&
   ParamComposite::addDMatrix(const char *label, DMatrix<Type>& matrix, int m, int n) 
   {
      DMatrixParam<Type>* ptr = new DMatrixParam<Type>(label, matrix, m, n);
      list_.push_back(ptr);
      isLeaf_.push_back(true);
      ++size_;
      ptr->setIndent(*this);
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         ptr->setParamCommunicator(paramCommunicator());
      }
      #endif
      return *ptr;
   }
   
   /* 
   * Add and read a DMatrixParam associated with a built-in array.
   */
   template <typename Type>   
   DMatrixParam<Type>&
   ParamComposite::readDMatrix(std::istream &in, const char *label, 
                               DMatrix<Type>& matrix, int m, int n) 
   {
      DMatrixParam<Type>* ptr = &addDMatrix<Type>(label, matrix, m, n);
      ptr->readParam(in);
      return *ptr;
   }
 
   /*
   * Add and load a DMatrix < Type > C two-dimensional matrix parameter.
   */
   template <typename Type> 
   DMatrixParam<Type>& 
   ParamComposite::loadDMatrix(Serializable::IArchive &ar, const char *label, 
                               DMatrix<Type>& matrix, int m, int n)
   {
      DMatrixParam<Type>* ptr = &addDMatrix<Type>(label, matrix, m, n);
      ptr->load(ar);
      return *ptr;
   }

   /*
   * Get class name string.
   */
   inline std::string ParamComposite::className() const
   {  return className_; }
      
} 
#endif
