#-----------------------------------------------------------------------
# This makefile file defines:
#
#   - The name $(PYTHON_LIB) of the python module 
#
# This file is included by the makefile in the python subdirectory
#-----------------------------------------------------------------------

ifdef PYTHON_FLAG
python_LIB=$(SRC_DIR)/../lib/_simpatico.so
endif

#-----------------------------------------------------------------------
