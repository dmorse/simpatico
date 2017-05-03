mcMd_potentials_misc_=\
    mcMd/potentials/misc/EnergyCalculator.cpp \
    mcMd/potentials/misc/StressCalculator.cpp 

mcMd_potentials_misc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_misc_))
mcMd_potentials_misc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_misc_:.cpp=.o))

