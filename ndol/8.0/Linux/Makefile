#
# Makefile for NDOL2015 - v. 7.0.3
#
# The following variables will be set when the Makefile is invoked with
# a specific machine type.  The variables correspond to:
# CPP     => The C preprocessor on the machine
# DEBUG   => The flag used to invoke compiling in debug mode
# EXE     => The name of the executable that will be generated
# FC      => The name of the Fortran compiler
# FFLAGS  => A list of compile options
# LDFLAGS => A list of link or load options
# EXTRAS  => 'Special' routines used on this machine
#
#CPP	= cpp -P -traditional
DEBUG	= 
EXE     = ndol2015-7.0.4.x
FC	= gfortran
#FFLAGS  = -O3 -march=pentium -fbounds-check -fno-automatic
FFLAGS  = -O3 -march=native -fbounds-check -fno-automatic
#LDFLAGS	= second.o
LDFLAGS = -static 
#LDFLAGS =
MAKE	= make
EXTRAS  =
#
# This variable is set to the list of "include" files which is just SIZES.
#
INCLUDES	= ndoldim.inc 
#
# OBJS is the list of all object files that need to be linked into the
# executable.  If a new source file is added, then the corresponding
# object file name must be added to the list below.
#
OBJS	= \
app.o \
bdata.o \
cardin.o \
cimat1.o \
ciout.o \
cisym.o \
degen.o \
dipm.o \
ener.o \
excindo.o \
excite.o \
exmat.o \
filen1.o \
filen5.o\
filen6.o \
geome.o \
getgeo.o \
gettxt.o \
guessalan.o \
guessp.o \
hyf.o \
init.o \
input1.o \
input5.o \
input6.o \
intout.o \
iter.o \
main.o \
matnis.o \
movlap.o \
mumu.o \
muorb.f \
peglag.o \
pegleg.o \
peglig.o \
peglog.o \
qrdiag.o \
reada.o \
readint.o \
rmopac.o \
scfmat.o \
scfout.o \
scfq.o \
selec.o \
theoga.o \
ucase.o
#
# For Sun arcstation 2 - version 4.1.3 of OS
#all:
#	@echo "Compiling for SPARC..."
#	@$(MAKE) FC='$(SPARCFC)' FFLAGS='$(SPARCFFLAGS)' CPP='$(SPARCCPP)' \
#	        EXE='$(SPARCEXE)' EXTRAS='$(SPARCEXTRAS)' linkit
#
# For Linux
all:
	@echo "Compiling for LINUX..."
	@$(MAKE) FC='$(FC)' FFLAGS='$(FFLAGS)' \
		EXE='$(EXE)'  linkit

#
# The default method of how to make object files from Fortran files.
#
.f.o:
#	$(CPP) $< TEMP.f
#	$(FC) -c $(FFLAGS) TEMP.f
	$(FC) -c $(FFLAGS) $<
#	/bin/mv TEMP.o $*.o
#	/bin/rm TEMP.[of]
#
#
linkit: $(OBJS)
	@echo " "
	@echo "Linking $(EXE)..."
	$(FC) $(LDFLAGS) -o $(EXE) $(OBJS)
	@echo " "
	@ls -l $(EXE)
#
# This removes all object files when "make clean" is entered by the user.
#
clean:
	rm -f *.o *.il
#
# Next comes all of the source file names and their dependencies (if any)
# and specific compile options (if any).  If no specific compile options
# are listed, then the default specified above is used.  If a new source file
# is added to NDOL, then an entry should be included below
app.o:		app.f      
bdata.o: 	bdata.f	$(INCLUDES) 
cardin.o:	cardin.f
cimat1.o:	cimat1.f	 $(INCLUDES)
ciout.o:	ciout.f	$(INCLUDES)
cisym.o:	cisym.f	$(INCLUDES)
degen.o:	degen.f	$(INCLUDES)
dipm.o:   	dipm.f	$(INCLUDES)   
ener.o:   	ener.f	$(INCLUDES)   
excindo.o:	excindo.f	$(INCLUDES)
excite.o:	excite.f	$(INCLUDES)
exmat.o: 	exmat.f	$(INCLUDES) 
filen1.o:	filen1.f	$(INCLUDES)
filen5.o:	filen5.f
filen6.o:       filen6.f
geome.o: 	geome.f	$(INCLUDES) 
getgeo.o:	getgeo.f	$(INCLUDES)
gettxt.o:	gettxt.f
guessp.o:	guessp.f	$(INCLUDES)
guessalan.o:    guessalan.f     $(INCLUDES)
hyf.o:   	hyf.f		$(INCLUDES)   
init.o:  	init.f	$(INCLUDES)  
input1.o:	input1.f	$(INCLUDES)
input5.o:	input5.f	$(INCLUDES)
input6.o:	input6.f	$(INCLUDES)
intout.o:	intout.f	$(INCLUDES)
iter.o:  	iter.f	$(INCLUDES)  
main.o:  	main.f	$(INCLUDES)  
matnis.o:	matnis.f	$(INCLUDES)
movlap.o:	movlap.f	$(INCLUDES)
mumu.o:  	mumu.f	$(INCLUDES)  
muorb.o:	muorb.f	$(INCLUDES)
peglag.o:	peglag.f	$(INCLUDES)
pegleg.o:	pegleg.f	$(INCLUDES)
peglig.o:	peglig.f	$(INCLUDES)
peglog.o:	peglog.f	$(INCLUDES)
qrdiag.o: 	qrdiag.f	$(INCLUDES) 
reada.o:  	reada.f  
readint.o:	readint.f	$(INCLUDES)
rmopac.o:	rmopac.f	$(INCLUDES)
scfmat.o:	scfmat.f	$(INCLUDES)
scfout.o:	scfout.f	$(INCLUDES)
scfq.o:  	scfq.f	$(INCLUDES)  
selec.o: 	selec.f	$(INCLUDES) 
theoga.o:	theoga.f	$(INCLUDES)
ucase.o:	ucase.f
