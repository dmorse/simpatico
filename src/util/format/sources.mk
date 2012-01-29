
util_format_SRCS=$(SRC_DIR)/util/format/Bool.cpp \
    $(SRC_DIR)/util/format/Dbl.cpp $(SRC_DIR)/util/format/Format.cpp \
    $(SRC_DIR)/util/format/Int.cpp $(SRC_DIR)/util/format/Lng.cpp \
    $(SRC_DIR)/util/format/Str.cpp $(SRC_DIR)/util/format/write.cpp 

util_format_OBJS=$(util_format_SRCS:.cpp=.o)

