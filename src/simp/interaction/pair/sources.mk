simp_interaction_pair_=\
    simp/interaction/pair/DpdPair.cpp \
    simp/interaction/pair/LJPair.cpp \
    simp/interaction/pair/WcaPair.cpp 

simp_interaction_pair_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_interaction_pair_))
simp_interaction_pair_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_interaction_pair_:.cpp=.o))

