
MDSIM=$(SRC_DIR)/main/mdSim
MCSIM=$(SRC_DIR)/main/mcSim
DDSIM=$(SRC_DIR)/main/ddSim

main_SRCS=$(MDSIM).cpp $(MCSIM).cpp
ifdef UTIL_MPI
main_SRCS+=$(DDSIM).cpp
endif

main_OBJS=$(main_SRCS:.cpp=.o)

