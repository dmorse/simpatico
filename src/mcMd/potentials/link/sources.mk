mcMd_potentials_link_=mcMd/potentials/link/LinkFactory.cpp 

mcMd_potentials_link_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_link_))
mcMd_potentials_link_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_link_:.cpp=.o))

