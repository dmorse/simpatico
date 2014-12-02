tools_config_=\
   tools/config/ConfigReader.cpp \
   tools/config/ConfigReaderFactory.cpp \
   tools/config/DdMdConfigReader.cpp \
   tools/config/HoomdConfigReader.cpp \
   tools/config/TypeMap.cpp \
   tools/config/ConfigWriter.cpp \
   tools/config/ConfigWriterFactory.cpp \
   tools/config/DdMdConfigWriter.cpp \
   tools/config/HoomdConfigWriter.cpp \

tools_config_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_config_))
tools_config_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_config_:.cpp=.o))

