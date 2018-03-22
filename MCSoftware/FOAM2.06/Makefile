#----------------------------------------------------------------------
# NOTES: This is usually needed before execution
#	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}
# If you use ROOT then remember to define env. variable ROOTSYS
#----------------------------------------------------------------------
#       make Demo-run
#	make DemoPers
# without ROOT (modify ROOT_DEF.h amd make clean)
#	make DemoNR-run
#----------------------------------------------------------------------
.SUFFIXES:  .o .c .cxx .f .l  .s  .sh  .h  .a
#======================================================================
#  Root 
RPATH =$(ROOTSYS)
#RPATH = /afs/cern.ch/user/j/jadachsp/private/root
#RPATH = /afs/cern.ch/na49/library.4/ROOT/new/i386_redhat60/root
RLIB=   -L$(RPATH)/lib -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint \
	-lPostscript -lMatrix -lPhysics -lm -ldl -lpthread -rdynamic
#======================================================================
#####  GNU c++ compiler
CXXFLAGS    = -g  -ansi -Wall -fPIC -D_REENTRANT -I$(ROOTSYS)/include
CXX         =  g++
CLD         =  g++
#-----------------------------
# Compile all c++ files
.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<
#======================================================================
#####  f77 compiler
# forbid underscore completely       -fno-second-underscore
# forbid double underscore only      -fno-underscoring
FFLAGS  =  -g -O2 -fno-second-underscore -fno-automatic -ffixed-line-length-132 -fbackslash 
F77 = gfortran
#  Fortran libraries for linker
FLIB= -lgfortran -lgfortranbegin 
CLF         =  $(F77)
# Compile all f77 files
.f.o:
	$(F77) $(FFLAGS) -c $<
#--------------------------------------------------------------------
# make Demo-run		make Demo-map
#--------------------------------------------------------------------
Demo.o:		      TFDISTR.o   TFOAM.o     TFMAXWT.o
#
TFDISTR.o:	TFDISTR.h	TFOAM_INTEGRAND.h
TFOAM.o:	TFOAM.h     	TFOAM_INTEGRAND.h  TRND.o
TFMAXWT.o:	TFMAXWT.h
TRND.o:		TRND.h
TFOAM_INTEGRAND.o:	TFOAM_INTEGRAND.h
TFHST.o: 	TFHST.h
#------------------------------------------------------------------------------
HEADERS = TRND.h TRandf.h TFMAXWT.h TFOAM_INTEGRAND.h TFOAM.h TFDISTR.h
OBJECTS = TRND.o TRandf.o TFMAXWT.o TFOAM_INTEGRAND.o TFOAM.o TFDISTR.o
OBJECTF = $(OBJECTS) circe2.o
FoamDict.cxx: $(HEADERS) FoamLinkDef.h
	@echo "Generating dictionary $@..."
	rootcint -f $@ -c $^
libFoam206.so: $(OBJECTS) FoamDict.o FoamDict.cxx
	g++ -shared  -o libFoam206.so $(OBJECTS) FoamDict.o
#------------------------------------------------------------------------------
Demo.exe:                  Demo.o $(OBJECTF)  libFoam206.so
	$(CLD) -o Demo.exe Demo.o $(OBJECTF)  libFoam206.so $(RLIB)  $(FLIB) 
DemoPers.exe:                  DemoPers.o $(OBJECTF)  libFoam206.so
	$(CLD) -o DemoPers.exe DemoPers.o $(OBJECTF)  libFoam206.so  $(RLIB)  $(FLIB) 
#------------------------------------------------------------------------------
Demo-run: Demo.exe
	(rm -f output-Demo; time ./Demo.exe > output-Demo; echo "See output-Demo for results")
Demo-map: Demo-run
	(latex Demo-map.tex; dvips Demo-map.dvi -o; gv Demo-map.ps) 
DemoPers: Demo-run DemoPers.exe
	(./DemoPers.exe )
#----------------------------------------------------------------------
line_count:
	(cat *.cxx *.h | wc -l )
#==============================================================================
# No ROOT version, note that TFHST replaces THF of ROOT
# dont forget to modify ROOT_DEF.h
# make DemoNR-run -f Makefile
DemoNR.exe:                  Demo.o $(OBJECTF) TFHST.o
	$(CLD) -o DemoNR.exe Demo.o $(OBJECTF) TFHST.o $(FLIB)
DemoNR-run: DemoNR.exe
	(rm -f output-Demo; time ./DemoNR.exe > output-Demo; echo "See output-Demo for results.")
#==============================================================================
clean:
	(rm -f *.o *.a *~ *.exe core *.gz lmap *.out *.dvi *.aux *.log FoamDict* libFoam206.so )
#------------------------------------------------------------------------------
