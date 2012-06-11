#include "ContainersTestComposite.h"
#include <util/util/initStatic.h>

int main() 
{
   Util::initStatic(); 
   ContainersTestComposite test;
   test.run();
}
