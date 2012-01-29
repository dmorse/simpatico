
util_param_SRCS=$(SRC_DIR)/util/param/Begin.cpp \
    $(SRC_DIR)/util/param/Blank.cpp $(SRC_DIR)/util/param/End.cpp \
    $(SRC_DIR)/util/param/Label.cpp $(SRC_DIR)/util/param/ParamComponent.cpp \
    $(SRC_DIR)/util/param/ParamComposite.cpp \
    $(SRC_DIR)/util/param/Parameter.cpp 

util_param_OBJS=$(util_param_SRCS:.cpp=.o)

