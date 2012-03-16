SRC_DIR =..

include $(SRC_DIR)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/python/defines.mk
include $(SRC_DIR)/python/patterns.mk
include $(SRC_DIR)/python/sources.mk
include $(SRC_DIR)/python/cpp/sources.mk

all: $(python_OBJS) $(python_LIB)

clean:
	rm -f $(python_OBJS) $(python_OBJS:.o=.d) $(python_LIB)

clean-deps:
	rm -f $(python_OBJS:.o=.d)

ifdef HOOMD_FLAG
PLUGIN_INSTALL_DIR=$(HOOMD_INSTALL_PATH)/lib/hoomd/python-module/hoomd_plugins/simpatico

install:
	$(shell [ -d "$(PLUGIN_INSTALL_DIR)" ] || mkdir -p $(PLUGIN_INSTALL_DIR))
	cp $(python_LIB) $(PLUGIN_INSTALL_DIR)
	cp -r pymodule/* $(PLUGIN_INSTALL_DIR)
endif

-include $(python_OBJS:.o=.d)
