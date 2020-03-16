SOURCES := 	./src/main.c \
		./src/utils.c \
		./src/xmem.c \
		./src/xstring.c \
		./src/parse_ini.c \
		./src/output.c \
		./src/outgal.c \
		./src/dconfObj.c \
		./src/halo_tree.c \
		./src/gal_gtree.c \
		./src/read_trees.c \
		./src/copy_halo_gal.c \
		./src/run.c \
		./src/nion.c \
		./src/nion_gal.c \
		./src/domain.c \
		./src/sort.c \
		./src/redshift_list.c \
		./src/photion_gal.c \
		./src/cutandresort.c \
		./src/evol_gal.c \
		./src/radiative_feedback.c \
		./src/starformation_SNfeedback.c \
		./src/ion_emissivity.c \
		./src/fesc.c \
		./src/check_cifogParam.c \
		./src/cifog/confObj.c \
		./src/cifog/grid.c \
		./src/cifog/input_nion.c \
		./src/cifog/fraction_q.c \
		./src/cifog/filtering.c \
		./src/cifog/phys_const.c \
		./src/cifog/self_shielding.c \
		./src/cifog/density_distribution.c \
		./src/cifog/recombination.c \
		./src/cifog/mean_free_path.c \
		./src/cifog/convolution_fftw.c \
		./src/cifog/input_redshifts.c\
		./src/cifog/input_grid.c \
		./src/cifog/photion_background.c \
		./src/cifog/redshift_tools.c \
		./src/cifog/checks.c \
		./src/cifog/cifog_delphi.c
		
		
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := astraeus

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
