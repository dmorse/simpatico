SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk

all: $(util_OBJS) $(util_LIB)

clean:
	rm -f $(util_OBJS) $(util_OBJS:.o=.d) $(util_LIB)
	cd tests; make clean
	rm -f */*.o */*/*.o */*/*/*.o
	rm -f */*.d */*/*.d */*/*/*.d

# Removes all generated files except util/defines.mk
veryclean:
	make clean
	-rm -f lib*.a
	-rm -f boundary/Boundary.h

-include $(util_OBJS:.o=.d)

