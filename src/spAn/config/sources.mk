spAn_config_=\
   spAn/config/ConfigReader.cpp \
   spAn/config/ConfigReaderFactory.cpp \
   spAn/config/DdMdConfigReader.cpp \
   spAn/config/HoomdConfigReader.cpp \
   spAn/config/TypeMap.cpp \
   spAn/config/ConfigWriter.cpp 

spAn_config_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_config_))
spAn_config_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_config_:.cpp=.o))

