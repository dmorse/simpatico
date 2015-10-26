ddMd_analyzers_scattering_=\
     ddMd/analyzers/scattering/StructureFactor.cpp\
     ddMd/analyzers/scattering/StructureFactorGrid.cpp\
     ddMd/analyzers/scattering/VanHove.cpp

ddMd_analyzers_scattering_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_analyzers_scattering_))
ddMd_analyzers_scattering_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_analyzers_scattering_:.cpp=.o))

