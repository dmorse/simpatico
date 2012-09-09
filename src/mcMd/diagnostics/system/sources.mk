mcMd_diagnostics_system_SRCS=$(SRC_DIR)/mcMd/diagnostics/system/AtomMSD.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/BondLengthDist.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/BlockRadiusGyration.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/ComMSD.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/CompositionProfile.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/DumpConfig.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/IntraPairAutoCorr.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/IntraStructureFactor.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/LinearRouseAutoCorr.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/RadiusGyration.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/RDF.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/RingRouseAutoCorr.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/StructureFactor.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/StructureFactorGrid.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/StructureFactorP.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/StructureFactorPGrid.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/SystemDiagnosticFactory.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/VanHove.cpp \
    $(SRC_DIR)/mcMd/diagnostics/system/BoundaryAverage.cpp 

ifdef MCMD_PERTURB
ifdef UTIL_MPI
mcMd_diagnostics_system_SRCS+=\
    $(SRC_DIR)/mcMd/diagnostics/system/MigratingVanHove.cpp
endif
endif

mcMd_diagnostics_system_OBJS=$(mcMd_diagnostics_system_SRCS:.cpp=.o)

