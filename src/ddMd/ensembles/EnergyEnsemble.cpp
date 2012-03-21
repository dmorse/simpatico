#ifndef DDMD_ENERGY_ENSEMBLE_CPP
#define DDMD_ENERGY_ENSEMBLE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "EnergyEnsemble.h"
#include <iostream>

#ifdef UTIL_MPI

#include <util/mpi/MpiStructBuilder.h>

namespace Util{

   /**
   * Initialize EnergyEnsemble MPI Datatype.
   */
   MPI::Datatype MpiTraits<DdMd::EnergyEnsemble>::type = MPI::BYTE;
   bool MpiTraits<DdMd::EnergyEnsemble>::hasType = false;

   /**
   * Initialize EnergyEnsemble::Type MPI Datatype.
   */
   MPI::Datatype MpiTraits<DdMd::EnergyEnsemble::Type>::type = MPI::INT;
   bool MpiTraits<DdMd::EnergyEnsemble::Type>::hasType = true;

}

#endif

namespace DdMd
{

   using namespace Util;

   /**
   * Constructor.
   */
   EnergyEnsemble::EnergyEnsemble(Type type)
    : temperature_(1.0),
      beta_(1.0),
      type_(type)
   {}

   /**
   * Set the temperature.
   */
   void  EnergyEnsemble::setTemperature(double temperature)
   {
      if (!isIsothermal()) {
	 UTIL_THROW("Must be an isothermal ensemble");
      }
      temperature_ = temperature;
      beta_        = 1.0/temperature;
   }

   /**
   * Read the type and (if necessary) temperature from file.
   */
   void EnergyEnsemble::readParam(std::istream& in)
   {
      readBegin(in, "EnergyEnsemble");
      read<Type>(in, "type", type_);
      if (isIsothermal()) {
         read<double>(in, "temperature", temperature_);
         beta_ = 1.0/temperature_;
      }
      readEnd(in);
   }

   /*
   * Extract an EnergyEnsemble::Type from an istream as a string.
   */
   std::istream& operator>>(std::istream& in, EnergyEnsemble::Type &type)
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

   /**
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
      builder.commit(Util::MpiTraits<EnergyEnsemble>::type);
      Util::MpiTraits<DdMd::EnergyEnsemble>::hasType = true;
   }

   #endif

}
#endif
