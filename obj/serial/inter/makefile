OBJ_DIR_REL =..

include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk

all: $(inter_OBJS) $(inter_LIB)

clean:
	rm -f $(inter_OBJS) $(inter_LIB)
	rm -f $(inter_OBJS:.o=.d)
	cd tests; make clean

veryclean:
	make clean
	rm -f */*.o */*/*.o
	rm -f */*.d */*/*.d
	rm -f lib*.a

-include $(inter_OBJS:.o=.d)

