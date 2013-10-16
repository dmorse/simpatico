mcMd_diagnostics_system_=\
    mcMd/diagnostics/system/AtomMSD.cpp \
    mcMd/diagnostics/system/BondLengthDist.cpp \
    mcMd/diagnostics/system/BlockRadiusGyration.cpp \
    mcMd/diagnostics/system/ComMSD.cpp \
    mcMd/diagnostics/system/CompositionProfile.cpp \
    mcMd/diagnostics/system/DumpConfig.cpp \
    mcMd/diagnostics/system/IntraPairAutoCorr.cpp \
    mcMd/diagnostics/system/IntraStructureFactor.cpp \
    mcMd/diagnostics/system/IntraStructureFactorGrid.cpp \
    mcMd/diagnostics/system/LinearRouseAutoCorr.cpp \
    mcMd/diagnostics/system/RadiusGyration.cpp \
    mcMd/diagnostics/system/RDF.cpp \
    mcMd/diagnostics/system/RingRouseAutoCorr.cpp \
    mcMd/diagnostics/system/StructureFactor.cpp \
    mcMd/diagnostics/system/StructureFactorGrid.cpp \
    mcMd/diagnostics/system/StructureFactorP.cpp \
    mcMd/diagnostics/system/StructureFactorPGrid.cpp \
    mcMd/diagnostics/system/SystemDiagnosticFactory.cpp \
    mcMd/diagnostics/system/VanHove.cpp \
    mcMd/diagnostics/system/BoundaryAverage.cpp 

ifdef MCMD_PERTURB
ifdef UTIL_MPI
mcMd_diagnostics_system_+=\
    mcMd/diagnostics/system/MigratingVanHove.cpp
endif
endif

mcMd_diagnostics_system_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_diagnostics_system_))
mcMd_diagnostics_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_diagnostics_system_:.cpp=.o))

