
/* 
* Example of code added to a class *.h file to define an MPI type.
* Replace string CPPTYPE by the name of the relevant cpp class.
*/

// Header include -----------------------------------------------------

#include <util/global.h>
#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>
#endif

   // Within public section of class declaration ----------------------

      #ifdef UTIL_MPI

         /**
         * Commit associated MPI DataType.
         */
         static void commitMpiType();

      #endif

   // After class declaration -----------------------------------------
 
   #ifdef UTIL_MPI

      /**
      * Explicit specialization MpiTraits<CPPTYPE>.
      */
      template <>
      class MpiTraits<CPPTYPE>
      {  public:  static MPI::Datatype type; };

   #endif

