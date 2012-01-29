#ifndef PARAM_COMPOSITE_H
#define PARAM_COMPOSITE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>     // base class
#include <util/param/ScalarParam.h>        // template class
#include <util/param/CArrayParam.h>        // template class
#include <util/param/DArrayParam.h>        // template class
#include <util/param/FArrayParam.h>        // template class
#include <util/param/CArray2DParam.h>      // template class
#include <util/param/DMatrixParam.h>       // template class
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
   * Any class that reads parameters from the input parameter file must
   * be derived from ParamComposite. The format of the associated block 
   * of a parameter file is defined by the implementation of the virtual
   * readParam() method. The readParam() method of each subclass is 
   * implemented using protected methods provided by ParamComposite to 
   * read any number of individual parameters, and any number of "child"
   * ParamComposite objects. The fact that a readParam() method method can 
   * invoke readParam() methods of other ParamComposite objects allows a 
   * file format to be defined as a hierarchical tree, in which each 
   * ParamComposite object is a node and each parameter is a leaf.
   * 
   * Each ParamComposite has a private array of pointers to ParamComponent 
   * objects, stored in the order in which they are read from file by the
   * readParam() method. We will refer to this as the format array. Each 
   * element of the format array may point to a Parameter object (which 
   * represents a single parameter), a Begin or End object (which represents 
   * a line containing the opening or closing bracket for a parameter block), 
   * a Blank object (i.e., a blank line), or another ParamComposite object. 
   *
   * The ParamComposite::writeParam() method writes the parameter file block 
   * by sequentially invoking the writeParam() method of each of the objects 
   * in the format array. This allows the writeParam() method to reproduce 
   * the file format that is defined by the readParam() method. 
   *
   * The readParam() method of each subclass of ParamComposite should be 
   * implemented using the methods read < T >, readBegin, readEnd, readBlank, 
   * and readParamComposite.  Each of these methods creates a new object of 
   * the specified type, adds a pointer to the new object to the format array, 
   * and then invokes the the readParam method of the new object in order to 
   * read the associated line or block of the parameter file:
   *  
   * \par
   * <b>readBegin()</b>: The readBegin() method reads the opening line for the
   * block associated with a ParamComposite, and adds a pointer to a new Begin 
   * object to the format array. The opening line should contain the class name
   * followed by an opening bracket, as in "MyClass{". Each Begin object 
   * has a class name string. 
   *
   * \par
   * <b>read < T > </b>: The read < T >() method template reads a single parameter 
   * of type T from stream and adds an associated ScalaraParam<T> object object 
   * to the format array. This template may be used with any type T for which
   * an inserter (>>) and extractor (<<) are defined. Each ScalarParam object 
   * contains a Label object and a pointer to the associated parameter variable.  
   *
   * \par
   * <b> readParamComposite():</b> The readParamComposite() method adds a
   * pointer to a child ParamComposite object to the format array, and invokes
   * the readParam() method of the child in order to read a nested subblock of
   * of the parameter file.
   *
   * \par 
   * <b> readBlank()</b>: The readBlank method add a pointer to a Blank object
   * to the format array, indicating a blank line.
   *
   * \par
   * <b>readEnd()</b>: The readEnd method reads a line containing the closing
   * bracket for the block associated with a ParamComposite, and adds a pointer 
   * to a new End object to the format array. 
   *
   * Other specialized read methods are provided for 1D and 2D C arrays.
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
      * Sets maximum length of list to listCapacity.
      *
      * \param listCapacity maximum length of parameter list
      */
      ParamComposite(int listCapacity);
  
      /** 
      * Virtual destructor.
      */
      virtual ~ParamComposite();
   
      /**
      * Resets ParamComposite to its empty state.
      *
      * This method deletes Parameter, Begin, End, and Blank objects in the
      * format list, recursively invokes resetParam() for any ParamComposite 
      * objects in the list, nullifies all pointers in the list, and resets 
      * the number of items in the list to 0. 
      */
      void resetParam();
   
      /** 
      * Read all parameters from an input stream.
      *
      * Empty default implementation.
      *
      * \param in input stream for reading
      */
      virtual void readParam(std::istream &in) 
      {}
   
      /** 
      * Write all parameters to an output stream.
      *
      * This default implementation writes all parameters to file,
      * descending children recursively. 
      *
      * \param out output stream for reading
      */
      virtual void writeParam(std::ostream &out);
   
      /// \name read* methods   
      /// \brief Each of these method invokes an associated add* method to create a 
      /// new ParamComponent object, and then invoke the readParam() method of the 
      /// new object to read the associated line or block of a file.
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
      /// \name add* methods
      /// \brief These methods add a ParamComponent object to the parameter
      /// list, but do not read any data from an input stream.
      //@{
   
      /**
      * Add a child ParamComposite object to the list.
      *
      * \param child child ParamComposite object
      * \param next  true if the indent level is one higher than parent.
      */
      void addParamComposite(ParamComposite& child, bool next = true);
   
      /**  
      * Add a new Param < Type > object to the list.
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
      
   private:
   
      /// Array of pointers to ParamComponent objects.
      std::vector<ParamComponent*> list_;    

      /// Array of booleans, elements false for ParamComposite, true otherwise.
      std::vector<bool> isLeaf_;  

      /// Number of ParamComponent objects in list
      int size_;     

   };

   // add and read method templates for scalar parameters
 
   /* 
   * Add a ScalarParam to the list.
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
      if (hasParamCommunicator()) 
         ptr->setParamCommunicator(paramCommunicator());
      #endif
      return *ptr;
   }
   
   /* 
   * Add a Param to the list, and read its contents from file. 
   */ 
   template <typename Type>
   ScalarParam<Type>& 
   ParamComposite::read(std::istream &in, const char *label, Type &value) 
   {
      ScalarParam<Type>* ptr = &add<Type>(label, value);
      ptr->readParam(in);
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
         if (hasParamCommunicator()) 
            ptr->setParamCommunicator(paramCommunicator());
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
         if (hasParamCommunicator()) 
            ptr->setParamCommunicator(paramCommunicator());
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
      if (hasParamCommunicator()) 
         ptr->setParamCommunicator(paramCommunicator());
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
      if (hasParamCommunicator()) 
         ptr->setParamCommunicator(paramCommunicator());
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
      if (hasParamCommunicator()) 
         ptr->setParamCommunicator(paramCommunicator());
      #endif
      return *ptr;
   }
   
   /* 
   * Add and read a DMatrixParam associated with a built-in array.
   */
   template <typename Type>   
   DMatrixParam<Type>&
   ParamComposite::readDMatrix(std::istream &in, 
                   const char *label, DMatrix<Type>& matrix, int m, int n) 
   {
      DMatrixParam<Type>* ptr = &addDMatrix<Type>(label, matrix, m, n);
      ptr->readParam(in);
      return *ptr;
   }
 
} 
#endif
