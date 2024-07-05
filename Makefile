###############################
###   MAKEFILE FOR VP_PIC   ###
###############################

# IMPORTANT: Please check which compiler you have before compiling the code.
# I have several options here but it might be that you have a different compiler
# so you will need to change the flags. Hopefully it will not be very hard to do.

# Simple macro for new line (there must be two blank lines inside for it to work!).

define newline


endef

# Figure out on which machine and architecture we are.
# Notice that "uname" returns "Linux" on a linux box
# and "Darwin" on a Mac.

HOST := $(shell hostname -s)
ARCH := $(shell uname)

# Directories where make will search for files.

vpath %.f90 src:src
vpath %.mod objs
vpath %.o objs

vpath VP_PIC exe

# Select the FORTRAN compiler.

# If we have MPI figure out which compiler it uses.
# The "mpi** -show" shows the actual mpi call
# explicitly, including the compiler called and the flags. 
# This is piped into a sed call which matches all letters
# and numbers in first word, but it will not match extra
# characters like "-" or "." (this is good, as it will eliminate
# version number and things like that at the end of the
# compiler name).  We then pipe this again into 'basename',

ifeq ($(shell which gfortran | sed -n "s:.*\/gfortran:gfortran:p"),gfortran)
FC := gfortran
else

ifeq ($(shell which ifort | sed -n "s:.*\/ifort:ifort:p"),ifort)
FC := ifort
else

# No Fortran compiler found.

$(error $(newline)Makefile error: No FORTRAN compiler found)

# Close all the if statements.  This is ugly, but it avoids
# using "else if" that older versions of make don't allow.

endif
endif

# Compiler and compilation flags for the gfortran compiler.
#
# -O3                  Optimization.
# -ffree-form          Fortran free form format.
# -Jobjs               Put module files in subdirectory "objs".
# -funroll-loops       Unroll loops.
# -fopenmp             Use OpenMP calls for optimization.

ifeq ($(FC),gfortran)

#FLAGS := -O3 -ffree-form -Jobjs -funroll-loops -fopenmp -w  -Wall#-fallow-argument-mismatch
#FLAGS := -O3 -Wall -I. -Wno-unused-dummy-argument -fopenmp
#FLAGS := -Ofast -ffast-math -ftree-vectorize -funroll-loops -fopenmp
#FLAGS := -O2 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace

else

# For the current version of the INTEL compiler "ifort".
#
# -O3                  Optimization.
# -free                Fortran free form format.
# -shared-intel        Dynamically look for libraries.
# -nowarn              Suppress warnings.
# -module objs         Put module files in subdirectory "objs".

ifeq ($(FC),ifort)

#FLAGS := -O3 -axCORE-AVX512,CORE-AVX2 -xAVX -free -shared-intel -module objs -qopenmp -g -trace -check all -fpe0 -warn -traceback -debug extended 
#FLAGS := -O3 -xhost -free -shared-intel -nowarn -module objs -qopenmp -llapack -lblas
FLAGS := -O3 -xhost -free -shared-intel -module objs -qopenmp -align array64byte
else

# Unknown compiler.  If we get here it means that mpi uses
# a compiler that is not considered above so the flags are not
# known.

$(error $(newline)Makefile error: Unknown FORTRAN compiler: $(FC))

# Close all the if statements.  This is ugly, but it avoids
# using "else if" that older versions of make don't allow.

endif
endif

# Fortran Compiler
FF := $(FC)


# Object files corresponding to Fortran modules. I separate
# them from the rest to be sure they are compiled first.

MODS = parameters.o arrays.o utils.o functions.o


# This line automatically looks for all f90 files in
# the "src" directory.  No need to add new files by hand.

FILES := $(patsubst %.f90,%.o,$(wildcard src/*.f90))
FILES += $(patsubst %.f90,%.o,$(wildcard src/*/*.f90))


# Object files.

OBJS := $(notdir $(FILES))

OBJS := $(filter-out $(MODS),$(sort $(OBJS)))

# Define pattern rule to compile f90 files.
#
# -I objs           Look for object files in subdirectory "objs".
# -c                Do not link, just create object file.
# -o                Name of object file follows.
#  $<               Source file (in this case %.f90)
#  $@               Target file (in this case %.o)

%.o : %.f90
	@ echo "COMPILING FILE: $(notdir $<)"
	@ $(FF) $(FLAGS) -I objs -c $< -o objs/$@
	@ echo

# By default, do nothing for targets for which no rule is specified.
# Basically, this stops make from complaining about the fact that the
# f90 files created by the perl scripts are not there to begin with.

.DEFAULT : ; @

# Main make target.  The prerequisites are done from left to right,
# and are defined below.

start :  hello dir compile link .timeend

# Hello message.

hello :
	@ touch .timestart
	@ /bin/rm -f .config; echo $(ARCH) > .config
	@ echo
	@ echo "****************************"
	@ echo "***   COMPILING VP_PIC   ***"
	@ echo "****************************"
	@ echo
	@ echo $(WW)
	@ echo "FORTRAN COMPILER =" $(FC)
	@ echo

# Create subdirectories objs and exe.

dir :
	@ mkdir -p exe; mkdir -p objs

# Targets to run perl scripts that create files for parameter
# assignment, for dealing with arrays and for system calls.
# Notice that the scripts need to be run only if the prerequisite
# files have changed since the last time we compiled.  Since the
# make system can only compare the dates on files, I use this trick:
# I create empty files ".perl*" in the "objs" directory after running
# each perl script.  The next time we compile the code, the dates on
# these files can be compared with the dates of the prerequisite files.


# Compile and link object files.

compile : $(OBJS)

link : VP_PIC

VP_PIC : $(OBJS)
	@ echo
	@ echo "LINKING ..."
	@ echo
	@ echo
	cd objs; $(FF) $(FLAGS) $(MODS) $(OBJS) -o ../exe/VP_PIC
	@ echo
	@ echo
	@ echo  "COMPILATION DONE!"
	@ echo
	@ echo
	@ touch .timeend

# Create object files.  Here we create all object files in one go
# using the rule defined above.  Notice that no prerequisites are
# needed apart from the default ones *.f90 defined by the compilation
# rule above.

$(MODS) :
$(OBJS) : $(MODS)

# Up to date message.  Here I make use of the same trick described
# above for the perl scripts.  To know if everything is up to date, I
# compare the date on the empty file ".timestart" which is touched as
# soon as make starts, with the date of the file ".timeend" which is
# touched only after creating the executable.  In this way, if the
# executable was not created (because it is more recent than the
# object files, which are in turn more recent that the f90 files,
# i.e. the executable is up to date), then the file ".timeend" will be
# older than the file ".timestart" and the message will be output.

.timeend : .timestart
	@ echo "EXECUTABLE IS UP TO DATE, NOTHING TO COMPILE."
	@ echo
	@ echo

# clean deletes subdirectory objs, timing files and executable file,
# but leaves everything else in subdirectory exe alone.

clean :
	@ /bin/rm -f -r objs .time* .config
	@ /bin/rm -f exe/VP_PIC

# cleanobj deletes subdirectory "objs" and timing files
# but leaves subdirectory "exe" alone.

cleanobj :
	@ /bin/rm -f -r objs .time*

# veryclean is the same as clean, but also deletes the subdirectory "exe"
# with ALL its contents!

veryclean :
	@ /bin/rm -f -r objs exe .time* .config
	@ /bin/rm -f src/auto/*.f90 src/auto/*.inc

# real8, real10 and rea16 change all declarations of real numbers
# in the code to the corresponding values.

# Help.

help:
	@ echo " "
	@ echo "make            - Compiles code and copies par files to subdirectory 'exe'"
	@ echo "make clean      - Deletes object files, executable, and all automatically generated Fortran files"
	@ echo "make cleanobj   - Deletes object files"
	@ echo "make veryclean  - Same as 'make clean', but also deletes subdirectory 'exe' with all its contents, BEWARE!"
	@ echo "make real8      - Changes all Fortran sources to use real8  variables"
	@ echo "make real10     - Changes all Fortran sources to use real10 variables"
	@ echo "make real16     - Changes all Fortran sources to use real16 variables"
	@ echo " "

