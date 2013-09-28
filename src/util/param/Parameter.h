#ifndef UTIL_PARAMETER_H
#define UTIL_PARAMETER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>
#include <util/archives/Serializable_includes.h>
#include <util/param/Label.h>
#include <util/global.h>

namespace Util
{

   /**
   * A single variable in a parameter file.
   *
   * Parameter is a base class for objects that read and write the 
   * value of a single C++ variable from or to file. The file format
   * for a parameter contains a string label followed by a value. 
   * 
   * A Parameter may be required element in a parameter file or optional,
   * depending on the value of the bool isRequired parameter of the 
   * constructor. An optional element is "active" only if an appropriate 
   * entry with the correct label has read from a parameter file, or if 
   * a value has been loaded from an archive. A required Parameter is
   * always active.  
   *
   * \ingroup Param_Module
   */
   class Parameter : public ParamComponent
   {

   public:

      /// Width of output field for a scalar variable.
      static const int Width = 20;

      /// Precision for io of floating point data field.
      static const int Precision = 12;

      /**
      * Constructor.
      *
      * \param label       label string preceding value in file format
      * \param isRequired  Is this a required parameter?
      */
      Parameter(const char *label, bool isRequired = true);

      /**
      * Destructor.
      */
      virtual ~Parameter();

      /**
      * Read a label and (if the label matches) a parameter value.
      * 
      * The parameter file format for a Parameter consists of a label string
      * followed by value. The value is read if and only if the label matches 
      * the expected value for this Parameter. If this Parameter is required 
      * and the input label not match, an error message is printed to the log 
      * file and Exception is thrown. If the Parameter is not required and the 
      * input label does not match, the label string is retained in an buffer 
      * for later processing by the readParam method of other ParamComponent
      * objects.
      *
      * Upon entry to this function, a label string is read into a input buffer 
      * if and only if the buffer is empty. The buffer is a static member of the 
      * Label class, which can retain a label between invocations of the 
      * readParameter method of different ParamComponent objects.  Once a 
      * label string is read from file, it remains in the input buffer for 
      * processing until it is matched, in which case the buffer is cleared 
      * to allow processing of the next label.
      * 
      * \param in input stream from which to read
      */
      virtual void readParam(std::istream &in);
   
      /**
      * Load from an archive.
      *
      * An optional Parameter loads the value of an isActive flag, and then loads
      * a parameter value only if the isActive is true. A required Parameter simply
      * loads the parameter value. 
      *
      * \param ar input archive from which to load
      */
      virtual void load(Serializable::IArchive& ar);
   
      /**
      * Save to an archive.
      *
      * An optional Parameter saves the value of the isActive flag, and then saves 
      * a parameter value only if the isActive is true. A required Parameter simply
      * saves its value. The label string is not saved to the archive. 
      *
      * \param ar output archive to which to save
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Return label string.
      */
      std::string label() const;

      /**
      * Is this an optional parameter?
      */
      bool isRequired() const;

      /**
      * Is this parameter active?
      */
      bool isActive() const;

   protected:

      /// Label object that contains parameter label string.
      Label label_;

      /// Is this parameter active (always true if is required).
      bool isActive_;
      
      /**
      * Read parameter value from an input stream.
      * 
      * \param in input stream from which to read
      */
      virtual void readValue(std::istream& in){}

      /**
      * Load bare parameter value from an archive.
      *
      * \param ar input archive from which to load
      */
      virtual void loadValue(Serializable::IArchive& ar){}

      /**
      * Save parameter value to an archive.
      *
      * \param ar output archive to which to save
      */
      virtual void saveValue(Serializable::OArchive& ar){}

      #ifdef UTIL_MPI
      /**
      * Broadcast parameter value within the ioCommunicator.
      *
      * \param ar output archive to which to save
      */
      virtual void bcastValue(){}
      #endif

   };

}
#endif
