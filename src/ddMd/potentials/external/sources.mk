ddMd_potentials_external_SRCS=\
    $(SRC_DIR)/ddMd/potentials/external/ExternalPotential.cpp \
    $(SRC_DIR)/ddMd/potentials/external/ExternalFactory.cpp 

ddMd_potentials_external_OBJS=$(ddMd_potentials_external_SRCS:.cpp=.o)

