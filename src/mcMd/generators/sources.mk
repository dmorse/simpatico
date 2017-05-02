
mcMd_generators_= \
    mcMd/generators/Generator.cpp \
    mcMd/generators/PointGenerator.cpp \
    mcMd/generators/generatorFactory.cpp 

ifdef SIMP_BOND
mcMd_generators_+= \
    mcMd/generators/LinearGenerator.cpp 
endif

mcMd_generators_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_generators_))
mcMd_generators_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_generators_:.cpp=.o))

