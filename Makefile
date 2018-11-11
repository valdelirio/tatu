mod1d: Anderson.o computationalstuff.o filtros.o escolhadofiltro.o dehx.o dehy.o dev.o dmhx.o dmhy.o dmv.o csv_file.o main1D.o
	gfortran -o mod1d.x Anderson.o computationalstuff.o filtros.o escolhadofiltro.o dehx.o dehy.o dev.o dmhx.o dmhy.o dmv.o csv_file.o main1D.o

Anderson.o: Anderson.for
	gfortran -std=legacy -c Anderson.for

computationalstuff.o: computationalstuff.f08
	gfortran -c computationalstuff.f08

filtros.o: filtros.f08
	gfortran -c filtros.f08

escolhadofiltro.o: escolhadofiltro.f08
	gfortran -c escolhadofiltro.f08

dehx.o: dehx.f08
	gfortran -c dehx.f08

dehy.o: dehy.f08
	gfortran -c dehy.f08

dev.o: dev.f08
	gfortran -c dev.f08

dmhx.o: dmhx.f08
	gfortran -c dmhx.f08

dmhy.o: dmhy.f08
	gfortran -c dmhy.f08

dmv.o: dmv.f08
	gfortran -c dmv.f08

csv_file.o: csv_file.f08
		gfortran -c csv_file.f08

main1D.o: main1D.f08
	gfortran -c main1D.f08

clean:
	rm *.o *.mod *.x
