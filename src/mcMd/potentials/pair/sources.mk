mcMd_potentials_pair_SRCS=$(SRC_DIR)/mcMd/potentials/pair/DpdPair.cpp \
    $(SRC_DIR)/mcMd/potentials/pair/LJPair.cpp \
    $(SRC_DIR)/mcMd/potentials/pair/McPairPotential.cpp \
    $(SRC_DIR)/mcMd/potentials/pair/MdPairPotential.cpp \
    $(SRC_DIR)/mcMd/potentials/pair/PairFactory.cpp 

mcMd_potentials_pair_OBJS=$(mcMd_potentials_pair_SRCS:.cpp=.o)

