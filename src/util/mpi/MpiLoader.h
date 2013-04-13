#ifndef UTIL_MPI_LOADER_H
#define UTIL_MPI_LOADER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/mpi/MpiFileIo.h>          // base class
#include <util/containers/DArray.h>      // template argument
#include <util/containers/FArray.h>      // template argument
#include <util/containers/DMatrix.h>     // template argument
#include <util/mpi/MpiTraits.h>          // Needed for implementation
#include <util/mpi/MpiSendRecv.h>        // Needed for implementation
#include <util/global.h>


namespace Util
{

   /**
   * Methods to load data from input archive and broadcast to other processors.
   *
   * Each of the "load" method templates causes the ioProcessor to load a
   * variable from the archive ar and then broadcast it to all other 
   * processors in the communicator. 
   *
   *  \ingroup Mpi_Module
   */
   template <class IArchive>
   class MpiLoader : public MpiFileIo
   {

   public:

      /*
      * Constructor
      */
      MpiLoader(MpiFileIo& parent)
       : MpiFileIo(parent)
      {}

      /**  
      * Load and broadcast a single Type value.
      *
      * \param ar     archive for loading
      * \param value  reference to a Type
      */
      template <typename Type>
      void load(IArchive &ar, Type &value);
   
      /**  
      * Load and broadcat a C array.
      *
      * \param ar     archive for loading
      * \param value  pointer to array
      * \param n      number of elements
      * \return reference to the new CArrayParam<Type> object
      */
      template <typename Type>
      void load(IArchive &ar, Type *value, int n);
   
      /**  
      * Load and broadcast a DArray < Type > container.
      *
      * \param ar     archive for loading
      * \param array  DArray object
      * \param n      number of elements
      */
      template <typename Type>
      void load(IArchive &ar, DArray<Type>& array, int n);
   
      /**  
      * Load and broadcast an FArray < Type >.
      *
      * \param ar     archive for loading
      * \param array  FArray object to be loaded
      */
      template <typename Type, int N>
      void load(IArchive &ar, FArray<Type, N >& array);
   
      /**  
      * Load and broadcast a 2D CArray of Type objects.
      *
      * \param ar     archive for loading
      * \param value  pointer to array
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      */
      template <typename Type> void 
      load(IArchive &ar, Type *value, int m, int n);
  
      /**  
      * Load and broadcast a DMatrix<Type> object.
      *
      * \param ar     archive for loading
      * \param matrix DMatrix object
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      */
      template <typename Type> 
      void load(IArchive &ar, DMatrix<Type>& matrix, int m, int n);
 
   };
 
   /*  
   * Load and broadcast a single Type value.
   */
   template <typename IArchive> 
   template <typename Type> 
   void MpiLoader<IArchive>::load(IArchive &ar, Type &value)
   {
      if (isIoProcessor()) {
         ar >> value;
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
   template <typename IArchive> 
   template <typename Type>
   void MpiLoader<IArchive>::load(IArchive &ar, Type* value, int n)
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
   * Load a DArray < Type > container.
   */
   template <typename IArchive> 
   template <typename Type>
   void 
   MpiLoader<IArchive>::load(IArchive &ar, DArray<Type>& array, int n)
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
   * Load an FArray < Type, N > fixed-size array container.
   */
   template <typename IArchive> 
   template <typename Type, int N>
   void MpiLoader<IArchive>::load(IArchive &ar, FArray<Type, N >& array)
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
   * Load a CArray2DParam < Type > C two-dimensional array parameter.
   */
   template <typename IArchive> 
   template <typename Type> 
   void 
   MpiLoader<IArchive>::load(IArchive &ar, Type *array, int m, int n)
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
   template <typename IArchive> 
   template <typename Type> 
   void MpiLoader<IArchive>::load(IArchive &ar, DMatrix<Type>& matrix, 
                                  int m, int n)
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
