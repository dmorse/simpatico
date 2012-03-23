#-----------------------------------------------------------------------
# This makefile file defines path variables that specify local installation
# paths
#-----------------------------------------------------------------------

# if HOOMD_FLAG is defined, the path to the HOOMD installation needs to be 
# set here
HOOMD_INSTALL_PATH=${HOME}/hoomd-install

# the path to python needs to be set here
PYTHON_INCLUDE_PATH=/sw/keeneland/python/2.6.4/centos5.4_gnu4.1.2/include/python2.6

# If HOOMD_FLAG is defined, the path to cuda needs to be set here
CUDA_INSTALL_PATH=/sw/keeneland/cuda/4.1/linux_binary/

# the path to the boost libraries needs to be set here (if not in the
# standard paths)
BOOST_INCLUDE_PATH=/sw/keeneland/boost/1.44.0/centos5.5_gnu4.4.0/include

#-----------------------------------------------------------------------
# everything below this needs not to be changed
#-----------------------------------------------------------------------

INCLUDES+= -I$(HOOMD_INSTALL_PATH)/include -I$(CUDA_INSTALL_PATH)/include -I$(PYTHON_INCLUDE_PATH)
ifdef BOOST_INCLUDE_PATH
INCLUDES+= -I$(BOOST_INCLUDE_PATH)
endif
HOOMD_LIB=$(HOOMD_INSTALL_PATH)/lib/hoomd/python-module/hoomd.so
LIBS+=$(HOOMD_LIB)
MCMD_DEFS+= -DENABLE_CUDA -DSINGLE_PRECISION

