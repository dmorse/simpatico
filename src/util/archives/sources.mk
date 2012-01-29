
util_archives_SRCS=\
    $(SRC_DIR)/util/archives/MemoryOArchive.cpp \
    $(SRC_DIR)/util/archives/MemoryIArchive.cpp \
    $(SRC_DIR)/util/archives/MemoryCounter.cpp \
    $(SRC_DIR)/util/archives/BinaryFileOArchive.cpp \
    $(SRC_DIR)/util/archives/BinaryFileIArchive.cpp \
    $(SRC_DIR)/util/archives/TextFileOArchive.cpp \
    $(SRC_DIR)/util/archives/TextFileIArchive.cpp 

util_archives_OBJS=$(util_archives_SRCS:.cpp=.o)

