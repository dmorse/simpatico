mcMd_mdIntegrators_SRCS=$(SRC_DIR)/mcMd/mdIntegrators/MdIntegrator.cpp \
    $(SRC_DIR)/mcMd/mdIntegrators/MdIntegratorFactory.cpp \
    $(SRC_DIR)/mcMd/mdIntegrators/NveVvIntegrator.cpp \
    $(SRC_DIR)/mcMd/mdIntegrators/NvtDpdVvIntegrator.cpp \
    $(SRC_DIR)/mcMd/mdIntegrators/NvtNhIntegrator.cpp \
    $(SRC_DIR)/mcMd/mdIntegrators/NphIntegrator.cpp

mcMd_mdIntegrators_OBJS=$(mcMd_mdIntegrators_SRCS:.cpp=.o)

