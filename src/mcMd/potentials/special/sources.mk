mcMd_potentials_special_=\
   mcMd/potentials/special/SpecialFactory.cpp \
   mcMd/potentials/special/MdExternalFacade.cpp

mcMd_potentials_special_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_special_))
mcMd_potentials_special_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_special_:.cpp=.o))

