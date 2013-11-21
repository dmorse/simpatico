OBJ_DIR_REL =..

include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk

all: $(util_OBJS) $(util_LIB)

clean:
	rm -f $(util_OBJS) $(util_OBJS:.o=.d)
	rm -f $(util_LIB)
	cd tests; make clean

veryclean:
	make clean
	rm -f *.o */*.o */*/*.o 
	rm -f *.d */*.d */*/*.d 
	rm -f lib*.a

-include $(util_OBJS:.o=.d)

