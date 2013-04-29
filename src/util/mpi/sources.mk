
util_mpi_SRCS=$(SRC_DIR)/util/mpi/MpiFileIo.cpp \
    $(SRC_DIR)/util/mpi/MpiTraits.cpp 

ifdef UTIL_MPI
util_mpi_SRCS+=\
    $(SRC_DIR)/util/mpi/MpiLogger.cpp \
    $(SRC_DIR)/util/mpi/MpiSendRecv.cpp \
    $(SRC_DIR)/util/mpi/MpiStructBuilder.cpp
endif

util_mpi_OBJS=$(util_mpi_SRCS:.cpp=.o)

