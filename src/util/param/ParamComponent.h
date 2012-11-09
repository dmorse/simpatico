#ifndef UTIL_PARAM_COMPONENT_H
#define UTIL_PARAM_COMPONENT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/mpi/MpiFileIo.h>   // class member
#include <util/global.h>

#include <iostream>
#include <sstream>
#include <string>

namespace Util
{

   /**
   *  Abstract base class for classes that input and ouput parameters to file.
   *
   *  The readParam  method reads a parameter or parameter list from iostream. 
   *  The writeParam method writes a parameter or parameter list to an ostream. 
   *  The same io format should be used by write and read methods. 
   *
   *  \ingroup Param_Module
   */
   class ParamComponent
   {

   public:

      /// Destructor.
      virtual ~ParamComponent();

      /**
      * Read parameter(s) from file.
      *
      * \param in input stream
      */
      virtual void readParam(std::istream& in) = 0;

      /**
      * Read parameter(s) to file.
      *
      * \param out output stream
      */
      virtual void writeParam(std::ostream& out) = 0;

      /**
      * Nontrivial implementation provided by ParamComposite subclass.
      *
      * This empty default implementation is used by all leaf nodes.
      */
      virtual void resetParam()
      {}

      /**
      * Set indent level.
      *
      * If next=true (default) set indent level one higher than
      * that of parent. If next=false, set the same as parent.
      *
      * \param parent parent ParamComponent object
      * \param next   If true, set level one higher than for parent.
      */
      void setIndent(const ParamComponent& parent, bool next = true);

      /**
      * Return indent string for this object (string of spaces).
      */
      std::string indent() const;

      #ifdef UTIL_MPI
      /**
      * Set an MPI communicator for parameter IO.
      */
      void setParamCommunicator(MPI::Intracomm& communicator);
      #endif

      /**
      * Serialize this ParamComponent as a string.
      *
      * \param ar      saving or loading archive
      * \param version version id for archive
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      // Public static methods

      /**
      * Initialize static echo member to false.
      */
      static void initStatic();

      /**
      * Set echo parameter to true enable echoing, or false to disable.
      */
      static void setEcho(bool echo = true);

      /**
      * Get echo parameter. true = echoing enabled, false = disable.
      */
      static bool echo();

   protected:

      /**
      * Constructor.
      *
      * Protected to prevent instantiation of a conceptually abstract 
      * class. 
      *
      * On return the indent string is empty. If UTIL_MPI is defined,
      * no communicator is set upon construction.
      */
      ParamComponent();

      /**
      * Copy constructor.
      */
      ParamComponent(const ParamComponent& other);

      /**
      * Should this processor do parameter file IO?
      *
      * Always returns true if compiled without UTIL_MPI defined.
      */
      bool isParamIoProcessor() const
      {  return io_.isIoProcessor(); }

      #ifdef UTIL_MPI
      /**
      * Has an MPI communicator been set for parameter input?
      */
      bool hasParamCommunicator() const
      {  return io_.hasCommunicator(); }

      /**
      * Return the parameter communicator.
      */
      MPI::Intracomm& paramCommunicator() const
      {  return io_.communicator(); }
      #endif

   private:

      /// Indentation string, a string of spaces.
      std::string indent_;

      /// Object to identify if this processor can do file Io.
      MpiFileIo   io_;

      /// Parameter to enable (true) or disable (false) echoing.
      static bool echo_;

   // friend:

      #ifdef UTIL_MPI
      template <class T> friend class Factory;

      /*
      * Note: This friend declaration allows the Factory::readObject() 
      * method to set the Factory paramCommunicator to be that of the 
      * parent ParamComposite. 
      */
      #endif

   };

   /*
   * Serialize a ParamComponent as a string.
   */
   template <class Archive>
   void ParamComponent::serialize(Archive& ar, const unsigned int version)
   {
      std::string str;
      if (Archive::is_saving()) {
         std::ostringstream buffer;
         writeParam(buffer);
         str = buffer.str();
         // std::cout << "Saving archive" << std::endl;
         // std::cout << str << std::endl;
      } 
      ar & str;
      if (Archive::is_loading()) {
         // std::cout << "Loading archive" << std::endl;
         // std::cout << str << std::endl;
         std::istringstream buffer(str);
         readParam(buffer);
      } 
   }

} 
#endif
