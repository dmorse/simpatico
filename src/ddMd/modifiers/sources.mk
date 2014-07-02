ddMd_modifiers_=\
     ddMd/modifiers/Modifier.cpp \
     ddMd/modifiers/ModifierManager.cpp \
     ddMd/modifiers/ModifierFactory.cpp

ddMd_modifiers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_modifiers_))
ddMd_modifiers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_modifiers_:.cpp=.o))

