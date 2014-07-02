ddMd_sp_analyzers_=\
     ddMd/sp/analyzers/SpAnalyzer.cpp \
     ddMd/sp/analyzers/SpAnalyzerFactory.cpp \
     ddMd/sp/analyzers/SpAnalyzerManager.cpp \
     ddMd/sp/analyzers/AtomMSD.cpp 

ddMd_sp_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_sp_analyzers_))
ddMd_sp_analyzers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_sp_analyzers_:.cpp=.o))

