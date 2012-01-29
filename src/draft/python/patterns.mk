# ---------------------------------------------------------------------
# File: src/mcMd/patterns.mk
#
# This makefile contains the pattern rules used to compile all sources
# files in the directory tree rooted at directory src/python, which
# contains the source code for the McMd namespace. It is included by
# all makefile files in this directory tree. 
#
#-----------------------------------------------------------------------
# Compilation pattern rules

# Path(s) to search for header files. 
INCLUDES= -I$(SRC_DIR)

PYTHON_DEFS=$(UTIL_DEFS) -fPIC

ifdef HOOMD_FLAG
INCLUDES+= -I$(HOOMD_INSTALL_PATH)/include
HOOMD_LIB=$(HOOMD_INSTALL_PATH)/lib/hoomd/python-module/hoomd.so
endif

ifdef PYTHON_FLAG
INCLUDES+= -I$(PYTHON_INCLUDE_PATH)
endif


# Rule to compile all class source (*.cpp) files.
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(PYTHON_DEFS) -c -o $@ $<
	$(MAKEDEP)

# Rule to compile *.cc main programs for unit tests. 
%.o:%.cc
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) $(PYTHON_DEFS) -c -o $@ $<
	$(MAKEDEP)

# Note: The main program files for unit tests must use a file suffix *.cc,
# while all other source files must use *.cpp. 

