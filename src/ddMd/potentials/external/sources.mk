ddMd_potentials_external_=\
    ddMd/potentials/external/ExternalPotential.cpp \
    ddMd/potentials/external/ExternalFactory.cpp 

ddMd_potentials_external_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_potentials_external_))
ddMd_potentials_external_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_potentials_external_:.cpp=.o))

