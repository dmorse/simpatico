ddMd_misc_=\
   ddMd/misc/DdTimer.cpp

ddMd_misc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_misc_))
ddMd_misc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_misc_:.cpp=.o))

