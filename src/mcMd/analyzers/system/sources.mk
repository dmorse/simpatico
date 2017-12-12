mcMd_analyzers_system_=\
    mcMd/analyzers/system/AtomMSD.cpp \
    mcMd/analyzers/system/ComMSD.cpp \
    mcMd/analyzers/system/CompositionProfile.cpp \
    mcMd/analyzers/system/RDF.cpp \
    mcMd/analyzers/system/StructureFactor.cpp \
    mcMd/analyzers/system/StructureFactorGrid.cpp \
    mcMd/analyzers/system/StructureFactorP.cpp \
    mcMd/analyzers/system/StructureFactorPGrid.cpp \
    mcMd/analyzers/system/SystemAnalyzerFactory.cpp \
    mcMd/analyzers/system/VanHove.cpp \
    mcMd/analyzers/system/BoundaryAverage.cpp \
    mcMd/analyzers/system/ConfigWriter.cpp 

ifdef SIMP_BOND
mcMd_analyzers_system_+=\
    mcMd/analyzers/system/BondLengthDist.cpp \
    mcMd/analyzers/system/BlockRadiusGyration.cpp \
    mcMd/analyzers/system/IntraPairAutoCorr.cpp \
    mcMd/analyzers/system/IntraStructureFactor.cpp \
    mcMd/analyzers/system/IntraStructureFactorGrid.cpp \
    mcMd/analyzers/system/LinearRouseAutoCorr.cpp \
    mcMd/analyzers/system/RadiusGyration.cpp \
    mcMd/analyzers/system/RingRouseAutoCorr.cpp \
    mcMd/analyzers/system/ClusterIdentifier.cpp \
    mcMd/analyzers/system/ClusterHistogram.cpp \
    mcMd/analyzers/system/Cluster.cpp 
endif

mcMd_analyzers_system_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_system_))
mcMd_analyzers_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_system_:.cpp=.o))

