mcMd_mdIntegrators_=\
    mcMd/mdIntegrators/MdIntegrator.cpp \
    mcMd/mdIntegrators/MdIntegratorFactory.cpp \
    mcMd/mdIntegrators/NveVvIntegrator.cpp \
    mcMd/mdIntegrators/NvtDpdVvIntegrator.cpp \
    mcMd/mdIntegrators/NvtNhIntegrator.cpp \
    mcMd/mdIntegrators/NphIntegrator.cpp

mcMd_mdIntegrators_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mdIntegrators_))
mcMd_mdIntegrators_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_mdIntegrators_:.cpp=.o))

