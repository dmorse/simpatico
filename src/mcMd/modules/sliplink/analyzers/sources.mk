mcMd_modules_sliplink_analyzers_= \
    mcMd/modules/sliplink/analyzers/SliplinkMdAnalyzerFactory.cpp \
    mcMd/modules/sliplink/analyzers/SliplinkMcAnalyzerFactory.cpp \
    mcMd/modules/sliplink/analyzers/Crosslinker.cpp \
    mcMd/modules/sliplink/analyzers/EndtoEnd.cpp \
    mcMd/modules/sliplink/analyzers/EndtoEndXYZ.cpp \
    mcMd/modules/sliplink/analyzers/G1MSD.cpp \
    mcMd/modules/sliplink/analyzers/InterIntraLink.cpp \
    mcMd/modules/sliplink/analyzers/LinkLTPos.cpp \
    mcMd/modules/sliplink/analyzers/LinkLengthDist.cpp \
    mcMd/modules/sliplink/analyzers/LinkLifeTime.cpp \
    mcMd/modules/sliplink/analyzers/LinkMSD.cpp \
    mcMd/modules/sliplink/analyzers/NLinkAverage.cpp \
    mcMd/modules/sliplink/analyzers/SSChainDist.cpp \
    mcMd/modules/sliplink/analyzers/VelProf.cpp 

mcMd_modules_sliplink_analyzers_SRCS=\
  $(addprefix $(SRC_DIR)/, $(mcMd_modules_sliplink_analyzers_))
mcMd_modules_sliplink_analyzers_OBJS=\
  $(addprefix $(BLD_DIR)/, $(mcMd_modules_sliplink_analyzers_:.cpp=.o))

