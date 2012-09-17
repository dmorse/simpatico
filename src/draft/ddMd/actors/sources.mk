ddMd_actors_SRCS=\
     $(SRC_DIR)/ddMd/actors/Actor.cpp\
     $(SRC_DIR)/ddMd/actors/ActorManager.cpp\
     $(SRC_DIR)/ddMd/actors/ActorFactory.cpp

ddMd_actors_OBJS=$(ddMd_actors_SRCS:.cpp=.o)
