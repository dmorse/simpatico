/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EnergyEnsemble.h"
#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>
#endif

#include <iostream>

namespace Simp{

   using namespace Util;

   /*
   * Constructor.
   */
   EnergyEnsemble::EnergyEnsemble(Type type)
    : temperature_(1.0),
      beta_(1.0),
      type_(type)
   {  setClassName("EnergyEnsemble"); }

   /*
   * Destructor.
   */
   EnergyEnsemble::~EnergyEnsemble()
   {}

   /*
   * Set the temperature.
   */
   void EnergyEnsemble::setTemperature(double temperature)
   {
      if (!isIsothermal()) {
	 UTIL_THROW("Must be an isothermal ensemble");
      }
      temperature_ = temperature;
      beta_        = 1.0/temperature;
   }

   /*
   * Read the type and (if necessary) temperature from file.
   */
   void EnergyEnsemble::readParameters(std::istream& in)
   {
      read<Type>(in, "type", type_);
      if (isIsothermal()) {
         read<double>(in, "temperature", temperature_);
         beta_ = 1.0/temperature_;
      }
   }

   /*
   * Load internal state from an archive.
   */
   void EnergyEnsemble::loadParameters(Serializable::IArchive &ar)
   { 
      loadParameter<Type>(ar, "type", type_);
      if (isIsothermal()) {
         loadParameter<double>(ar, "temperature", temperature_);
         ar >> beta_;
      }
   }

   /*
   * Save internal state to an archive.
   */
   void EnergyEnsemble::save(Serializable::OArchive &ar)
   { 
      ar << type_;
      if (isIsothermal()) {
         ar << temperature_;
         ar << beta_;
      }
   }

   /*
   * Extract an EnergyEnsemble::Type from an istream as a string.
   */
   std::istream& operator >> (std::istream& in, EnergyEnsemble::Type &type)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "ADIABATIC" || buffer == "adiabatic") {
         type = EnergyEnsemble::ADIABATIC;
      } else
      if (buffer == "ISOTHERMAL" || buffer == "isothermal") {
         type = EnergyEnsemble::ISOTHERMAL;
      } else {
         UTIL_THROW("Invalid EnergyEnsemble::Type value input");
      }
      return in;
   }

   /*
   * Insert a EnergyEnsemble::Type to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, const EnergyEnsemble::Type &type)
   {
      if (type == EnergyEnsemble::ADIABATIC) {
         out << "adiabatic";
      } else
      if (type == EnergyEnsemble::ISOTHERMAL) {
         out << "isothermal";
      }
      return out;
   }


   #ifdef UTIL_MPI
   /*
   * Commit MPI Datatype.
   */
   void EnergyEnsemble::commitMpiType()
   {
      MpiStructBuilder builder;
      EnergyEnsemble   object;

      builder.setBase(&object);
      builder.addMember(&object.temperature_, MPI::DOUBLE);
      builder.addMember(&object.beta_, MPI::DOUBLE);
      builder.addMember(&object.type_, MPI::INT);
      builder.commit(MpiTraits<EnergyEnsemble>::type);
      MpiTraits<EnergyEnsemble>::hasType = true;
   }
   #endif

}
#ifdef UTIL_MPI
namespace Util 
{

   /*
   * Initialize Simp::EnergyEnsemble MPI Datatype.
   */
   MPI::Datatype MpiTraits<Simp::EnergyEnsemble>::type = MPI::BYTE;
   bool MpiTraits<Simp::EnergyEnsemble>::hasType = false;

   /*
   * Initialize Simp::EnergyEnsemble::Type MPI Datatype.
   */
   MPI::Datatype MpiTraits<Simp::EnergyEnsemble::Type>::type = MPI::INT;
   bool MpiTraits<Simp::EnergyEnsemble::Type>::hasType = true;

}
#endif
