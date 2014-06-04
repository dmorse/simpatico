ddMd_sp_analyzers_=\
     mdCf/analyzers/SpAnalyzer.cpp \
     mdCf/analyzers/SpAnalyzerFactory.cpp \
     mdCf/analyzers/SpAnalyzerManager.cpp \
     mdCf/analyzers/AtomMSD.cpp 

ddMd_sp_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_sp_analyzers_))
ddMd_sp_analyzers_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_sp_analyzers_:.cpp=.o))

