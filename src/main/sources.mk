
MDSIM=$(SRC_DIR)/main/mdSim
MCSIM=$(SRC_DIR)/main/mcSim
DDSIM=$(SRC_DIR)/main/ddSim

main_SRCS=$(SRC_DIR)/main/mdSim.cpp \
          $(SRC_DIR)/main/mcSim.cpp \
          $(SRC_DIR)/main/ddSim.cpp 

main_OBJS=$(main_SRCS:.cpp=.o)

