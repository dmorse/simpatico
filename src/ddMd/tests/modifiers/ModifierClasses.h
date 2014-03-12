#ifndef DDMD_TESTS_MODIFIER_CLASSES_H
#define DDMD_TESTS_MODIFIER_CLASSES_H

#include <ddMd/modifiers/Modifier.h>

namespace DdMd {

   using namespace Util;

   class ModifierA : public Modifier 
   {
   public:
   
      ModifierA() 
       : Modifier()
      {
         setClassName("ModifierA"); 
         set(Modifier::Flags::PostIntegrate1); 
      }
   
      ~ModifierA() 
      {} 

      void readParameters(std::istream& in) 
      {
         readInterval(in);
      }

      void postIntegrate1(long iStep)
      {
         std::cout << "Calling ModifierA::postIntegrate1" << std::endl;
      }
   
   };

}
#endif //ifndef DDMD_MODIFIER_CLASSES_H
