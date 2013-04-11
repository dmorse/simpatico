#ifndef UTIL_BROADCAST_LOADER_H
#define UTIL_BROADCAST_LOADER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/mpi/MpiFileIo.h>          // base class
#include <util/archives/Serializable.h>  // requires typedef
#include <util/global.h>


namespace Util
{

   /**
   * Class to load data from input archive and broadcast to other processors.
   *
   *  \ingroup Param_Module
   */
   class BroadcastLoader : public MpiFileIo
   {

   public:

      /*
      * Constructor
      */
      BroadcastLoader(MpiFileIo& parent)
       : MpiFileIo(parent)
      {}

      /**  
      * Load and broadcast a Type value.
      *
      * \param ar     archive for loading
      * \param value  reference to a Type
      */
      template <typename Type>
      void loadParameter(Serializable::IArchive &ar, Type &value);
   
      /**  
      * Add a C array parameter, and load its elements.
      *
      * \param ar     archive for loading
      * \param value  pointer to array
      * \param n      number of elements
      * \return reference to the new CArrayParam<Type> object
      */
      template <typename Type>
      void loadCArray(Serializable::IArchive &ar, Type *value, int n);
   
      /**  
      * Load and broadcast a DArray < Type >.
      *
      * \param ar     archive for loading
      * \param array  DArray object
      * \param n      number of elements
      */
      template <typename Type>
      void loadDArray(Serializable::IArchive &ar, DArray<Type>& array, int n);
   
      /**  
      * Load and broadcast a DArray < Type >.
      *
      * \param ar     archive for loading
      * \param array  FArray object
      */
      template <typename Type, int N>
      void loadFArray(Serializable::IArchive &ar, FArray<Type, N >& array);
   
      /**  
      * Load and broadcast a 2D CArray of Type objects.
      *
      * \param ar     archive for loading
      * \param value  pointer to array
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      */
      template <typename Type> void 
      loadCArray2D(Serializable::IArchive &ar, Type *value, int m, int n);
  
      /**  
      * Load and broadcast a DMatrix<Type> object.
      *
      * \param ar     archive for loading
      * \param matrix DMatrix object
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      */
      template <typename Type> 
      void loadDMatrix(Serializable::IArchive &ar, DMatrix<Type>& matrix, 
                       int m, int n);
 
   };
 
   /*  
   * Add a new ScalarParam< Type > object, and load its value.
   */
   template <typename Type> void 
   ParamComposite::loadScalar(Serializable::IArchive &ar, Type &value)
   {
      if (isIoProcessor()) {
         ar & value;
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<Type>(ioCommunicator(), value, 0); 
      }
      #endif
   }
   
   /*  
   * Add a C array parameter, and load its elements.
   */
   template <typename Type>
   void
   ParamComposite::loadCArray(Serializable::IArchive &ar, Type* value, int n)
   {
      if (isIoProcessor()) {
         for (int i = 0; i < n; ++i) {
            ar >> value[i];
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<Type>(ioCommunicator(), value, n, 0); 
      }
      #endif
   }

   /*
   * Add a DArray < Type > parameter, and load its elements.
   */
   template <typename Type>
   void 
   ParamComposite::loadDArray(Serializable::IArchive &ar, 
                              DArray<Type>& array, int n)
   {
      if (isIoProcessor()) {
         ar >> array;
         if (array.capacity() < n) {
            UTIL_THROW("Error: DArray capacity < n");
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<Type>(ioCommunicator(), array, n, 0); 
      }
      #endif
   }

   /*
   * Add and load an FArray < Type, N > fixed-size array parameter.
   */
   template <typename Type, int N>
   void ParamComposite::loadFArray(Serializable::IArchive &ar, 
                                   FArray<Type, N >& array)
   {
      if (isIoProcessor()) {
         for (int i = 0; i < N; ++i) {
            ar >> array[i];
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<Type>(ioCommunicator(), &(array[0]), N, 0); 
      }
      #endif
   }

   /*
   * Add and load a CArray2DParam < Type > C two-dimensional array parameter.
   */
   template <typename Type> 
   void ParamComposite::loadCArray2D(Serializable::IArchive &ar, 
                                     Type *array, int m, int n)
   {
      if (isIoProcessor()) {
         int i, j; 
         for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
               ar >> array[i*n + j];
            }
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<Type>(ioCommunicator(), &(array[0]), n*m, 0); 
      }
      #endif
   }
  
   /*
   * Add and load a DMatrix < Type > C two-dimensional matrix parameter.
   */
   template <typename Type> 
   void ParamComposite::loadDMatrix(Serializable::IArchive &ar, 
                                    DMatrix<Type>& matrix, int m, int n)
   {
      if (isIoProcessor()) {
         ar >> matrix;
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<Type>(ioCommunicator(), matrix, m, n, 0); 
      }
      #endif
   }
  
} 
#endif
