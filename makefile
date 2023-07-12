SRC1 = invtci.F gtm.F plot.F limits.F intdde.F 

SRC2 = shader.f limits.f cspack.f drawcl.f pltpl.f \
       pltl.f pltlp.f fexp.f drwvec.f

SINGLE_OR_DOUBLE         = "double precision"
#SINGLE_OR_DOUBLE        = "real"

FC = ncargcc
CPPFLAGS = -DSINGLE_OR_DOUBLE=$(SINGLE_OR_DOUBLE)
OBJ =	$(SRC1:.F=.o)  $(SRC2:.f=.o)

LIBS = -lslatec
GLIB = -lplotf -lncarg -lncarg_gks -lncarg_loc
GLIBsmooth = -lplotf -dashsmooth

FFLAGS = -u

invtc: $(OBJ)
	ncargcc $(OBJ) $(FFLAGS) $(LIBS) $(GLIBsmooth) -o $@

.F.f :
	/lib/cpp -C -P $(CPPFLAGS) $< > $*.f

inv.f : $(SRC1) $(SRC2)
	cat $(SRC1) $(SRC2)  > temp
	/lib/cpp -C -P $(CPPFLAGS) temp $@
	rm temp

invtci.f : invtci.F

intdde.f : intdde.F

gtm.f      : gtm.F

erf.f      : erf.F

plot2.f    : plot2.F

limits.f   : limits.F

# set suffixes list:
# -----------------
.SUFFIXES :
.SUFFIXES : .o .F .f

