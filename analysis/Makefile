SOURCES := 	./src/main.c \
		./src/utils.c \
		./src/xmem.c \
		./src/xstring.c \
		./src/parse_ini.c \
		./src/dconfObj.c \
		./src/outgal.c \
		./src/read_trees.c \
		./src/run.c \
		./src/zas_lists.c \
		./src/phys_const.c \
		./src/statistics.c \
		./src/build_index_tree_walking.c \
		./src/derive_properties.c \
		./src/selection.c \
		./src/uv.c \
		./src/operations_on_properties.c \
		./src/statistics_evolution.c \
		./src/selection_evolution.c \
		./src/sort.c \
		./src/communicate_mpi.c \
		./src/tree_correct.c \
		./src/get_grid_properties.c
		
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := analyse

USE-MPI=YES

include common.mk


.PHONY: all clean clena celan celna

all: $(SOURCES) $(EXECUTABLE)

celan celna clena:clean


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
