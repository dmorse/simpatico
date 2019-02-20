mdPp_config_=\
   mdPp/config/ConfigReader.cpp \
   mdPp/config/ConfigReaderFactory.cpp \
   mdPp/config/DdMdConfigReader.cpp \
   mdPp/config/SmpConfigReader.cpp \
   mdPp/config/HoomdConfigReader.cpp \
   mdPp/config/TypeMap.cpp \
   mdPp/config/ConfigWriter.cpp \
   mdPp/config/ConfigWriterFactory.cpp \
   mdPp/config/DdMdConfigWriter.cpp \
   mdPp/config/HoomdConfigWriter.cpp \

mdPp_config_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_config_))
mdPp_config_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mdPp_config_:.cpp=.o))

