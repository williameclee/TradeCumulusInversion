ncargf77 -c -o invtci.o invtci.F
ncargf77 -c -o gtm.o gtm.F
ncargf77 -c -o plot.o plot.F
ncargf77 -c -o limits.o limits.F
ncargf77 -c -o intdde.o intdde.F

ncargf77 -c -o shader.o shader.f
ncargf77 -c -o limits.o limits.f
ncargf77 -c -o cspack.o cspack.f
ncargf77 -c -o drawcl.o drawcl.f
ncargf77 -c -o pltpl.o pltpl.f
ncargf77 -c -o pltl.o pltl.f
ncargf77 -c -o pltlp.o pltlp.f
ncargf77 -c -o fexp.o fexp.f
ncargf77 -c -o drwvec.o drwvec.f

ncargf77 invtci.o gtm.o plot.o limits.o intdde.o shader.o cspack.o drawcl.o pltpl.o pltl.o pltlp.o fexp.o drwvec.o -o invtci_main.o