include ../../Makefile.in

program = run

#=====================================================================#
# Object Files
#=====================================================================#

include OBJECTS

# Algebra
F90FLAGS += -I../../lib/
LDFLAGS += -L../../lib/ ../../lib/libcoarrays.a -llapack

#=====================================================================#
# Targets
#=====================================================================#

all: $(program)

$(program): $(objects)
	$(F90) $(objects) -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so $(program)

#=====================================================================#
# Rule
#=====================================================================#

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<
