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
         //setClassName("ModifierA"); 
         //set(Modifier::Flags::PostIntegrate1); 
      }
   
      ~ModifierA() 
      {} 

      #if 0 
      void readParameters(std::istream& in) 
      {
         readInterval(in);
      }
      #endif

      void postIntegrate1(long iStep)
      {}
   
   };

}
#endif //ifndef DDMD_MODIFIER_CLASSES_H
