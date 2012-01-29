namespace Util
{

   /**
   * \defgroup Serialize_Module Serialization
   * \ingroup Util_NS_Module
   *
   * Serialization of C++ objects to/from file or memory.
   *
   * The classes, functions and function templates in this module provide a 
   * system for serializing sequences of C++ objects, to a file or random
   * access memory, in a form that stores the full internal state and allows
   * the object to be reconstructed.  The design is loosely based on that 
   * of the Boost serialization library,
   * http://www.boost.org/doc/libs/1_48_0/libs/serialization/doc/index.html
   * but is much simpler (and less powerful) than the Boost library. 
   */

   /**
   * \defgroup Archive_Module Archives
   * \ingroup  Serialize_Module
   *
   * The definition of an archive used here is very similar 
   * to that used in the Boost serialization library. An archive class may 
   * model either a saving / output archive, to which data is saved, or a 
   * loading / input archive, from which data is loaded.  By convention, the 
   * names of saving/output archive classes end with the string OArchive and 
   * the names of loading/input archive classes end with the string IArchive. 
   * Different archive classes store serialized objects in different forms. 
   * For example, TextFileOArchive and TextFileIArchive are saving and loading 
   * archive classes, respectively, that are wrappers for ofstream or ifstream 
   * file stream objects in which data is stored in a character representation.
   * BinaryFileOArchive and BinaryFileIArchive are saving/output and 
   * loading / input archives that store data in a binary format. 
   * MemoryOArchive and MemoryIArchive are saving and loading archives that 
   * stored data in binary form in a block of random-access memory. 
   *
   * Objects may be saved to a saving archive or loaded from a loading 
   * archive using overloaded operators, using the same syntax as that of the
   * Boost library.  Each saving archive class must define method templates
   * that overload the << (insertion) and & operators. These overloaded
   * operators must be equivalent, and must save an object to the archive.
   * If ar is an instance of a saving archive, such as BinaryFileOArchive, 
   * the expressions
   * \code
   *    ar << data; 
   *    ar &  data;
   * \endcode
   * are thus equivalent, and both save the state of variable data into
   * archive ar.  Each loading archive class must instead define template 
   * methods to overload the >> (extractor) and & operator, which must be 
   * equivalent, and which must load an object from the archive. If ar is 
   * an instance of a loading archive, such as BinaryFileIArchive, then 
   * the expressions
   * \code
   *    ar >> data;
   *    ar &  data;
   * \endcode
   * are equivalent, and both load the state of variable data from archive 
   * ar.
   *
   * Objects of type T can be saved to or loaded from an archive that is
   * is an instance of class Archive if and only if the compiler can find 
   * a function named serialize with the signature
   * \code
   *     void serialize(Archive& ar, T& data, unsigned int version)
   * \endcode
   * Here, "version" is an integer index that indicates the version of the
   * archive, which is normally given by an integer member of the archive
   * class. The operator & for a class Archive is normally implemented by 
   * a method template 
   * \endcode
   *
   *  template <typename T>
   *  void Archive::operator & (T& data);
   *  { serialize(*this, data, version_); }
   *
   * \endcode
   * that simply calls the appropiate serialize method. Here, version_ is
   * an integer member of the Archive class that stores a version id. 
   * Similar templates must be provided for the << or >> operator.
   * Definitions of the serialize function for saving archive types must 
   * save (write) data, and those for loading archive types must load 
   * (read) data. Each archive class provides serialize functions for all
   * of the built-in C++ types, as well as few other common data types 
   * such as std::string. 
   * 
   * Instances of user-defined classes may also be serialized if an appropriate 
   * serialize function is provided. Serialization of instances of a class T 
   * may be enabled by defining either:
   * 
   * - A global serialize function template with a signature
   * \code
   *
   * template <class Archive>
   * inline void serialize(Archive& ar, T& data, const unsigned int version);
   * 
   * \endcode
   *
   * - A serialize method template in class T, with a signature
   * \code
   *
   * template <class Archive>
   * void T::serialize(Archive& ar, const unsigned int version);
   * 
   * \endcode
   * Note that, in either case, the archive type is left as a template 
   * parameter. If a serialize class method is defined, it is accessed
   * by the following default template serialize function:
   * \code
   *
   * template <class Archive, typename T>
   * inline void serialize(Archive& ar, T& data, const unsigned int version)
   * {  data.serialize(ar, version); }
   *
   * \endcode
   * This template, which is defined in the file serialize.h, 
   * simply calls the serialize method of class T, if one exists.
   * When the C++ compiler needs a serialize method for a particular 
   * archive type Archive and data type T, it will look first for a
   * function serialize(Archive&, T&, unsigned int) with exactly the 
   * required signature, and then for an appropriate template. If a
   * global serialize function template is defined for class T with
   * the signature described above, in which the archive type is a 
   * template parameter but the data type T is explicit, this will 
   * be chosen in preference to the default template, in which both 
   * both the archive and data types are template parameters. 
   * Serialization of built-in C++ types always uses the explicit 
   * specializations that must be defined for these types for each
   * archive class.
   *
   * The use of a single operator & to represent both output (when applied to a 
   * saving archive) and input (when applied to a loading archive), makes it 
   * possible to write a single serialize function template for each class that
   * specifies both how to save and how to load instances of that class.  For 
   * example, consider the following definition of a simple complex number 
   * class:
   * \code 
   *
   *   class Complex  {
   *   public:
   *
   *      A(double real, double imag) : real_(real), imag_(imag) {}
   * 
   *      template <class Archive>
   *      void serialize(Archive& ar, unsigned int version)
   *      { 
   *         ar & real_;
   *         ar & imag_:
   *      }
   *
   *   private:
   *
   *      double real_;
   *      double imag_;   
   * 
   *   } 
   *
   * \endcode
   * The serialize method template provides instructions for the order in which
   * to either save the two floating point members of the class to a saving 
   * archive or load them from a loading archive. The use of a template in 
   * which the archive type is a parameter allows a single serialize method 
   * to be used with any type of saving or loading archive.
   *
   * The most serious disadvantage of this system is that definition of a 
   * serialize method as a template implies that this method cannot be 
   * virtual. As a result, the serialize method template for a class cannot 
   * be accessed polymorphically, via a pointer or reference to a base class. 
   * This limitation becomes a problem in designs in which some objects are 
   * accessed only via base class pointers.  The Serializable abstract base
   * class partially solves this problem, by replacing the serialize method 
   * template by a pair of virtual save() and load() methods.
   */

   /**
   * \defgroup Serializable_Module Serializable
   *
   * Serializable is an abstract base class that provides an interface for
   * serializing objects that uses virtual functions rather than templates.  
   * Each subclass of Serializable must define virtual save() and load() 
   * methods with the following signatures:
   * \code
   * virtual void save(Serializable::OArchiveType& ar);
   * virtual void load(Serializable::IArchiveType& ar);
   * \endcode
   * The typenames Serializable::OArchiveType and Serializable::IArchiveType
   * are typedefs that define a pair of archive classes to be used for 
   * serialization. 
   *
   * The advantage of using virtual functions is that it allows these methods 
   * to be accessed polymorphically, via base class pointers or references.
   * The disadvantage is that it requires the hard-coding (via a pair of
   * typedefs) of a single type of saving and loading archive. A serialize
   * method template for a class instead allows instances of that class to
   * be serialized to or from any type of archive, but cannot be virtual 
   * because C++ does not allow a method template to be virtual. In practice, 
   * a serialize method or function template should be defined for relatively 
   * simple, non-polymorphic classes, but more complicated polymorphic types 
   * should be derived from Serializable. 
   *
   * \ingroup Serialize_Module
   */

}
