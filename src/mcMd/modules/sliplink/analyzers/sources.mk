modules_sliplink_analyzers_= \
    modules/sliplink/analyzers/Crosslinker.cpp \
    modules/sliplink/analyzers/EndtoEnd.cpp \
    modules/sliplink/analyzers/EndtoEndXYZ.cpp \
    modules/sliplink/analyzers/G1MSD.cpp \
    modules/sliplink/analyzers/InterIntraLink.cpp \
    modules/sliplink/analyzers/LinkLTPos.cpp \
    modules/sliplink/analyzers/LinkLengthDist.cpp \
    modules/sliplink/analyzers/LinkLifeTime.cpp \
    modules/sliplink/analyzers/LinkMSD.cpp \
    modules/sliplink/analyzers/NLinkAverage.cpp \
    modules/sliplink/analyzers/SSChainDist.cpp \
    modules/sliplink/analyzers/VelProf.cpp \
    modules/sliplink/analyzers/SliplinkMdAnalyzerFactory.cpp \
    modules/sliplink/analyzers/SliplinkMcAnalyzerFactory.cpp 

modules_sliplink_analyzers_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_sliplink_analyzers_))
modules_sliplink_analyzers_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_sliplink_analyzers_:.cpp=.o))

