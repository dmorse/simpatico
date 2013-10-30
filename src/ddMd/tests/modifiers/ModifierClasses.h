#ifndef DDMD_MODIFIER_CLASSES_H
#define DDMD_MODIFIER_CLASSES_H

#include <ddMd/modifiers/Modifier.h>

namespace DdMd {

   using namespace Util;

   class ModifierA : public Modifier 
   {
   public:
   
      ModifierA() :
         Modifier()
      { set(Modifier::Flags::PostIntegrate1); }
     
      void postIntegrate1(long iStep)
      {}
   
   };

}
#endif //ifndef DDMD_MODIFIER_CLASSES_H
