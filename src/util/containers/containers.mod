namespace Util
{

   /**
   * \defgroup Container_Module Container Templates
   * \ingroup Util_NS_Module
   *
   * A set of container and iterator class templates.
   *
   * This module contains a set of simple container templates, some of which
   * are similar to containers provided by the C++ standard library. Bounds 
   * checking of indices for all array containers can be turned on (for 
   * safety) or off (for speed) by defining or not defining the UTIL_DEBUG 
   * preprocessor macro.
   *
   * Containers templates whose name contains the string 'Array' are one 
   * dimensional array containers, much like C arrays. All such containers
   * overload the subscript [] operator so as to return an object by 
   * reference, using the same syntax as a C array or a std::vector: If
   * A is an array, then A[i] is a reference to the ith element of A.
   * 
   * Container templates whose name contains the string 'Matrix' are two 
   * dimensional arrays. These overload the (int, int) operator to access
   * elements: If M is a Matrix, then M(i, j) is a reference to the element 
   * in column j of row i of A.
   *
   * Containers templates whose name begins with a letter 'D' (such as
   * DArray, DSArray, DPArray, and DMatrix) use dynamically allocated memory. 
   * The declaration "DArray<int> A" declares a dynamically allocated array
   * of integers.  Memory must be explicitly allocated for these containers 
   * by calling the "allocate" method after the array is instantiated. 
   * Dynamically allocated containers can only be allocated once, and (unlike 
   * std::vector) are not resizable. Attempting to allocate a container more 
   * than once is as an error. 
   * 
   * Containers templates whose name begins with a letter 'D' (such as
   * FArray, FSArray, FPArray, and FMatrix) are statistically allocated
   * fixed capacity container. The capacity of each such container is 
   * determined at compile time by a template parameter or parameters. Thus, 
   * for example, "FArray<int, 4> A;" declares a fixed size array of four 
   * integers, much like declaration "int V[4]" of a fixed size C array.
   *
   * Container templates DSArray and FSArray are "sized" arrays in which all 
   * elements are stored contiguously, and which store a logical size that is 
   * less than or equal to it physical capacity.  The capacity of an array is 
   * the number of elements for which memory has been allocated, or the maximum 
   * allowable size. The valid elements of a DSArray or an FSArray are stored 
   * contiguously from index 0 to size - 1. The size of such an array is initially 
   * set to zero, and elements may be added to the end of the array by the 
   * append() method, which increments the size. 
   *
   * The container templates DPArray and FPArray are dynamically and statically
   * allocated "pointer" arrays. A pointer array is a container that stores 
   * pointers to objects that exist outside of the array, rather than actual
   * objects. The pointer arrays are similar to the "sized" arrays in that 
   * they have a logical size that is less than or equal to there capacity,
   * and that elements must be added to the end of the array by an "append"
   * method. A pointer array DPArray< T > is different from an array of pointers
   * DArray<T*>, however, because DPArray< T > overloads the [] to return a
   * reference to an object of type T, rather than the T* pointer that points
   * to that object. The pointer arrays thus use the same syntax for element
   * access as the other arrays. 
   */

   /**
   * \defgroup Array_Module Object Arrays
   *
   * One-dimensional array containers that store objects by value, and associated 
   * iterators.
   *
   * The DArray and FArray containers are simple wrappers for dynamically
   * allocated and fixed-size C arrays, respectively. Both DArray < T > and
   * FArray <T, N> overload the [] subscript operator to provide a reference
   * to a specific element of an underlying C array of T objects.
   *
   * The DSArray and FSArray containers are dynamically and statically
   * allocated arrays, respectively, that have both a fixed capacity but
   * a variable logical size, with contiguous elements.
   *
   * An RArray < T > is an Array that is intended to be used as an alias
   * for, or a shallow copy of, a target DArray, FArray or C array. An RArray
   * contains a copy of the array address and capacity of the target array,
   * where are copied by the RArray::associate() method. Like other array
   * containers, an RArray overloads the [] operator to provide access to 
   * elements as references. The destructor of an RArray does not delete 
   * the associated C array, while the destructor of a DArray does.
   *
   * A RingBuffer is a cylic buffer array for which the append() method
   * adds elements to the end of a sequence if the buffer is not full, and
   * overwrites the oldest element if it is.
   *
   * \ingroup Container_Module
   */

   /**
   * \defgroup Pointer_Array_Module Pointer Arrays
   *
   * One-dimensional array containers that store pointers, and associated 
   * iterators.
   *
   * The DPArray and SPArray class templates are dynamically and statically
   * allocated pointer arrays, respectively. Each DPArray < T > and
   * SPArray <T, N> containers contains a C array of T* pointers, rather
   * than an actual array of T objects. Both DPArray and SPArray overload
   * the [] operator so as to return a T& reference to an associated object,
   * rather than a T* pointer. An element can be added to the end of a
   * DPArray or a SPArray by the append(T&) method.  The destructor for
   * a DPArray < T > deletes the underlying array of T* pointers, but
   * not the objects to which they point.
   *
   * An ArrayStack < T > container is a finite capacity stack that is
   * implemented as a dynamically allocated array of T* pointers. Objects
   * can be pushed onto or popped off the top of the stack using the 
   * push(T&) and pop() methods. An ArrayStack can be allocated only 
   * once, and cannot be resized.
   *
   * An SSet< T > is a container that holds pointers to an unordered
   * set of T objects. It provides fast addition and removal of single
   * elements, and iteration through all elements.
   *
   * An ArraySet < T > is a container that holds pointers to a subset of
   * the elements of an associated array of T objects. The indexing of
   * elements within the container is arbitrary.
   *
   * \ingroup Container_Module
   */

   /**
   * \defgroup Matrix_Module Matrix Containers
   *
   * Two-dimensional array containers that store by objects value.
   *
   * \ingroup Container_Module
   */

   /**
   * \defgroup List_Module Linked List
   *
   * A simple linked list implementation and associated iterator.
   *
   * \ingroup Container_Module
   */

   /**
   * \defgroup Iterator_Module Iterators
   *
   * Iterators for use with particular containers.
   *
   * \ingroup Container_Module
   */

}
