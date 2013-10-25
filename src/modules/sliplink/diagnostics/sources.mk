modules_sliplink_diagnostics_= \
    modules/sliplink/diagnostics/Crosslinker.cpp \
    modules/sliplink/diagnostics/EndtoEnd.cpp \
    modules/sliplink/diagnostics/EndtoEndXYZ.cpp \
    modules/sliplink/diagnostics/G1MSD.cpp \
    modules/sliplink/diagnostics/InterIntraLink.cpp \
    modules/sliplink/diagnostics/LinkLTPos.cpp \
    modules/sliplink/diagnostics/LinkLengthDist.cpp \
    modules/sliplink/diagnostics/LinkLifeTime.cpp \
    modules/sliplink/diagnostics/LinkMSD.cpp \
    modules/sliplink/diagnostics/NLinkAverage.cpp \
    modules/sliplink/diagnostics/SSChainDist.cpp \
    modules/sliplink/diagnostics/VelProf.cpp \
    modules/sliplink/diagnostics/SliplinkMdDiagnosticFactory.cpp \
    modules/sliplink/diagnostics/SliplinkMcDiagnosticFactory.cpp 

modules_sliplink_diagnostics_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_sliplink_diagnostics_))
modules_sliplink_diagnostics_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_sliplink_diagnostics_:.cpp=.o))

