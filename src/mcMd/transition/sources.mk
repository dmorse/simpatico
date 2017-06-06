mcMd_transition_=\
    mcMd/transition/MdParticle.cpp \
    mcMd/transition/MdSnapShot.cpp 

mcMd_transition_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_transition_))
mcMd_transition_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_transition_:.cpp=.o))

