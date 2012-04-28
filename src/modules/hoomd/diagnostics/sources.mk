modules_hoomd_diagnostics_SRCS=$(SRC_DIR)/modules/hoomd/diagnostics/GPUStructureFactorGrid.cpp \
                               $(SRC_DIR)/modules/hoomd/diagnostics/HoomdDiagnosticFactory.cpp

modules_hoomd_diagnostics_NVCC_SRCS=$(SRC_DIR)/modules/hoomd/diagnostics/GPUStructureFactorGrid.cu 

modules_hoomd_diagnostics_OBJS=$(modules_hoomd_diagnostics_SRCS:.cpp=.o)

modules_hoomd_diagnostics_NVCC_OBJS=$(modules_hoomd_diagnostics_NVCC_SRCS:.cu=.cu.o)

