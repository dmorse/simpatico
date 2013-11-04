#!/usr/bin/env python

from FileEditor import *

editor = FileEditor()
editor.setIsTest(True)

#editor.setBlockSize(5)
#editor.setFilter(r'^ *namespace MolMcD *\n *{')
#editor.setOld(r'MolMcD *\n *{( *\n)*')
#editor.setNew(r'Util\n{\n\n')
#for directory in mcmd_directories:
#   editor.editFiles(directory, "*.h")
#   editor.editFiles(directory, "*.cpp")

#editor.setFilter(r'#define +([A-Z]+_)+H *$')
#editor.setOld(r'#define +')
#editor.setNew(r'#define ')
#editor.editFiles(".", "*.h")

#editor.setFilter(r'#ifndef +([A-Z]+_)+H *$')
#editor.setOld(r'#ifndef +')
#editor.setNew(r'#ifndef ')
#editor.editFiles(".", "*.h")

#editor.setFilter(r'MCMD_MCMD_')
#editor.setOld(r'MCMD_MCMD_')
#editor.setNew(r'MCMD_')
#editor.editFiles(".", "*.h")
#editor.editFiles(".", "*.cpp")

editor.setFilter(r'#include <mcMd/misc/FileMaster.h')
editor.setOld(r'<mcMd')
editor.setNew(r'<util')
editor.editFiles(".", "*.h")
editor.editFiles(".", "*.cpp")

#editor.setFilter(r'setParam\(')
#editor.setOld(r'setParam\(')
#editor.setNew(r'initialize(')
#editor.editFiles(".", "*.h")
#editor.editFiles(".", "*.cpp")
#editor.editFiles(".", "*.tpp")
