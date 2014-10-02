modules_hoomd_analyzers_=\
    modules/hoomd/analyzers/GPUStructureFactorGrid.cpp \
    modules/hoomd/analyzers/HoomdAnalyzerFactory.cpp

modules_hoomd_analyzers_NVCC_=\
    modules/hoomd/analyzers/GPUStructureFactorGrid.cu 

modules_hoomd_analyzers_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_analyzers_))
modules_hoomd_analyzers_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_analyzers_:.cpp=.o))

modules_hoomd_analyzers_NVCC_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_analyzers_NVCC_))
modules_hoomd_analyzers_NVCC_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_analyzers_NVCC_:.cu=.cu.o))

