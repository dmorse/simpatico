#include "ChemistryTestComposite.h"
#include <util/misc/initStatic.h>

int main(int argc, char** argv) 
{
   Util::initStatic();
   ChemistryTestComposite runner;
   runner.run();
}
