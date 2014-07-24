#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include makefiles

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/inter/config.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk

#-----------------------------------------------------------------------
# Main targets

all: $(inter_OBJS) $(inter_LIB)

clean:
	rm -f $(inter_OBJS) $(inter_OBJS:.o=.d) $(inter_LIB)
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	rm -f */*.o */*.d
	rm -f lib*.a

#-----------------------------------------------------------------------
# Include dependency files

-include $(inter_OBJS:.o=.d)

