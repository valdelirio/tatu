dipolos: Anderson.o filtros.o escolhadofiltro.o dehx.o dehy.o dev.o InvtgcDpls.o
	gfortran -o dipolos.x Anderson.o filtros.o escolhadofiltro.o dehx.o dehy.o dev.o InvtgcDpls.o

Anderson.o: Anderson.for
	gfortran -c Anderson.for

filtros.o: filtros.f90
	gfortran -c filtros.f90

escolhadofiltro.o: escolhadofiltro.f90
	gfortran -c escolhadofiltro.f90

dehx.o: dehx.f90
	gfortran -c dehx.f90

dehy.o: dehy.f90
	gfortran -c dehy.f90

dev.o: dev.f90
	gfortran -c dev.f90

InvtgcDpls.o: InvtgcDpls.f90
	gfortran -c InvtgcDpls.f90
