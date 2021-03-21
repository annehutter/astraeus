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
		./src/tree_correct.c \
		./src/get_grid_properties.c \
		./src/communicate_mpi.c \
		./src/get_further_derived_properties.c \
		./src/SNfeedback.c \
		./src/ion_emissivity.c \
		./src/smoothing_grid.c \
		./src/map_gal_to_grid.c \
		./src/find_median.c \
		./src/statistics_median.c \
		./src/selection_median.c \
		./src/statistics_analytic.c \
		./src/selection_analytic.c \
		./src/numdens_analytic.c \
		./src/statistics_evolution_analytic.c \
		./src/selection_evolution_analytic.c \
		./src/output_lists.c
		
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := analyse

USE-MPI=YES
# GCC-DEBUG=YES

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
