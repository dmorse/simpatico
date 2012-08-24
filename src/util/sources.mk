include $(SRC_DIR)/util/util/sources.mk
include $(SRC_DIR)/util/format/sources.mk
include $(SRC_DIR)/util/memory/sources.mk
include $(SRC_DIR)/util/containers/sources.mk
include $(SRC_DIR)/util/mpi/sources.mk
include $(SRC_DIR)/util/signal/sources.mk
include $(SRC_DIR)/util/param/sources.mk
include $(SRC_DIR)/util/math/sources.mk
include $(SRC_DIR)/util/space/sources.mk
include $(SRC_DIR)/util/random/sources.mk
include $(SRC_DIR)/util/boundary/sources.mk
include $(SRC_DIR)/util/crystal/sources.mk
include $(SRC_DIR)/util/ensembles/sources.mk
include $(SRC_DIR)/util/accumulators/sources.mk
include $(SRC_DIR)/util/archives/sources.mk

util_SRCS=$(util_util_SRCS) $(util_format_SRCS) \
    $(util_memory_SRCS) $(util_containers_SRCS) $(util_mpi_SRCS) \
    $(util_signal_SRCS) $(util_param_SRCS) $(util_math_SRCS) \
    $(util_space_SRCS) $(util_random_SRCS) $(util_boundary_SRCS) \
    $(util_crystal_SRCS) $(util_ensembles_SRCS) \
    $(util_accumulators_SRCS) $(util_archives_SRCS)

util_OBJS=$(util_SRCS:.cpp=.o)

$(util_LIB): $(util_OBJS)
	$(AR) rcs $(util_LIB) $(util_OBJS)

