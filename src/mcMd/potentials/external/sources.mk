mcMd_potentials_external_=\
    mcMd/potentials/external/ExternalFactory.cpp 

mcMd_potentials_external_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_external_))
mcMd_potentials_external_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_external_:.cpp=.o))

