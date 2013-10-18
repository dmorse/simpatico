#include "InteractionTestComposite.h"

int main(int argc, char* argv[])
{
   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }

   InteractionTestComposite runner;
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
    }
   runner.run();
}
