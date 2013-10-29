modules_hoomd_diagnostics_=\
    modules/hoomd/diagnostics/GPUStructureFactorGrid.cpp \
    modules/hoomd/diagnostics/HoomdDiagnosticFactory.cpp

modules_hoomd_diagnostics_NVCC_=\
    modules/hoomd/diagnostics/GPUStructureFactorGrid.cu 

modules_hoomd_diagnostics_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_diagnostics_))
modules_hoomd_diagnostics_OBJS=\
    $(addprefix $(OBJ_DIR)/, $(modules_hoomd_diagnostics_:.cpp=.o))

modules_hoomd_diagnostics_NVCC_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_diagnostics_NVCC_))
modules_hoomd_diagnostics_NVCC_OBJS=\
    $(addprefix $(OBJ_DIR)/, $(modules_hoomd_diagnostics_NVCC_:.cu=.cu.o))

