main : test.o drift.o orbel_xv.o 
	gfortran test.o drift.o orbel_xv.o -o main

test.o : test.f90
	gfortran test.f90 -c test.o

drift.o : drift.f90
	gfortran drift.f90 -c drift.o

orbel_xv.o : orbel_xv.f90
	gfortran orbel_xv.f90 -c orbel_xv.o

clean :
	rm *.o main 

.PHONY : clean
