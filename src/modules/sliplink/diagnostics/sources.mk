modules_sliplink_diagnostics_SRCS= \
    $(SRC_DIR)/modules/sliplink/diagnostics/Crosslinker.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/EndtoEnd.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/EndtoEndXYZ.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/G1MSD.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/InterIntraLink.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/LinkLTPos.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/LinkLengthDist.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/LinkLifeTime.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/LinkMSD.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/NLinkAverage.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/SSChainDist.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/VelProf.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/SliplinkMdDiagnosticFactory.cpp \
    $(SRC_DIR)/modules/sliplink/diagnostics/SliplinkMcDiagnosticFactory.cpp 

modules_sliplink_diagnostics_OBJS=$(modules_sliplink_diagnostics_SRCS:.cpp=.o)

