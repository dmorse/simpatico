ddMd_actors_=\
     ddMd/actors/Actor.cpp \
     ddMd/actors/ActorManager.cpp 
     #ddMd/actors/ActorFactory.cpp

ddMd_actors_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_actors_))
ddMd_actors_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_actors_:.cpp=.o))

