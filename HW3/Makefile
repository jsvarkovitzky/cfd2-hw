# $CLASSHG/codes/fortran/roots/Makefile

OBJECTS = projection_method.o

.PHONY: flow plot clean

flow: projection_method.exe
	./projection_method.exe

projection_method.exe: $(OBJECTS)
	gfortran $(OBJECTS) -o projection_method.exe

plot: 
	@echo "Assuming .txt files are correct!"
	python makeplot.py

%.o : %.f90
	gfortran -c $< 

clean:
	rm -f *.o *.exe *.mod





