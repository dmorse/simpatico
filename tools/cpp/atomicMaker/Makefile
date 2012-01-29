SRC_DIR=../../../src

include $(SRC_DIR)/makefile.all
include $(SRC_DIR)/mcMd/makefile.inc
include $(SRC_DIR)/util/makefile.inc

clean:	
	-rm AtomicMaker

# ==========================================================================

ATOMICMAKER=AtomicMaker

$(ATOMICMAKER): $(ATOMICMAKER).cpp $(mcMd_OBJS) $(util_OBJS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $(ATOMICMAKER) $(ATOMICMAKER).cpp \
                          $(mcMd_OBJS) $(util_OBJS)

all: $(ATOMICMAKER)

