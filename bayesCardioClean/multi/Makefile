CC       = gcc
OBJ      = \
Amazigh/mathMatrix.o Amazigh/mathVector.o Amazigh/openFile.o Amazigh/plot2D.o \
Amazigh/printFile.o Amazigh/showFile.o Amazigh/supportC.o Amazigh/solvers.o \
Solvers/problemDef.o Solvers/solversAQ.o Solvers/RV-GOMEA/RV-GOMEA.o \
Solvers/util/FOS.o Solvers/util/SO_optimization.o Solvers/util/Tools.o \
Solvers/CMAES/src/boundary_transformation.o Solvers/CMAES/src/cmaes.o \
Solvers/LBFGSB/driver1.o Solvers/LBFGSB/lbfgsb.o Solvers/LBFGSB/linesearch.o Solvers/LBFGSB/linpack.o \
Solvers/LBFGSB/miniCBLAS.o Solvers/LBFGSB/print.o Solvers/LBFGSB/subalgorithms.o Solvers/LBFGSB/timer.o \
Gaussian/gaussian.o Solvers/functions.o Bayes/bayesTest.o \
Cardio/cell.o Cardio/FEM.o Cardio/FEMTests.o Cardio/cardio.o Cardio/geometry.o \
home.o
LINKOBJ  = \
Amazigh/mathMatrix.o Amazigh/mathVector.o Amazigh/openFile.o Amazigh/plot2D.o \
Amazigh/printFile.o Amazigh/showFile.o Amazigh/supportC.o Amazigh/solvers.o \
Solvers/problemDef.o Solvers/solversAQ.o Solvers/RV-GOMEA/RV-GOMEA.o \
Solvers/util/FOS.o Solvers/util/SO_optimization.o Solvers/util/Tools.o \
Solvers/CMAES/src/boundary_transformation.o Solvers/CMAES/src/cmaes.o \
Solvers/LBFGSB/driver1.o Solvers/LBFGSB/lbfgsb.o Solvers/LBFGSB/linesearch.o Solvers/LBFGSB/linpack.o \
Solvers/LBFGSB/miniCBLAS.o Solvers/LBFGSB/print.o Solvers/LBFGSB/subalgorithms.o Solvers/LBFGSB/timer.o \
Gaussian/gaussian.o Solvers/functions.o Bayes/bayesTest.o \
Cardio/cell.o Cardio/FEM.o Cardio/FEMTests.o Cardio/cardio.o Cardio/geometry.o \
home.o
CFLAGS   = -g -Wall -lm -std=c99
RM       = -rm -f

all: bayesOpt

bayesOpt: $(OBJ) 
	gcc $(OBJ) -o bayesOpt -lm
#Amazigh

Amazigh/mathMatrix.o: Amazigh/mathMatrix.c
	$(CC) -c Amazigh/mathMatrix.c -o Amazigh/mathMatrix.o $(CFLAGS)

Amazigh/mathVector.o: Amazigh/mathVector.c
	$(CC) -c Amazigh/mathVector.c -o Amazigh/mathVector.o $(CFLAGS)

Amazigh/openFile.o: Amazigh/openFile.c
	$(CC) -c Amazigh/openFile.c -o Amazigh/openFile.o $(CFLAGS)

Amazigh/plot2D.o: Amazigh/plot2D.c
	$(CC) -c Amazigh/plot2D.c -o Amazigh/plot2D.o $(CFLAGS)

Amazigh/printFile.o: Amazigh/printFile.c
	$(CC) -c Amazigh/printFile.c -o Amazigh/printFile.o $(CFLAGS)

Amazigh/showFile.o: Amazigh/showFile.c
	$(CC) -c Amazigh/showFile.c -o Amazigh/showFile.o $(CFLAGS)

Amazigh/supportC.o: Amazigh/supportC.c
	$(CC) -c Amazigh/supportC.c -o Amazigh/supportC.o $(CFLAGS)

Amazigh/solvers.o: Amazigh/solvers.c
	$(CC) -c Amazigh/solvers.c -o Amazigh/solvers.o $(CFLAGS)
#CMAES

Solvers/solversAQ.o: Solvers/solversAQ.c
	$(CC) -c Solvers/solversAQ.c -o Solvers/solversAQ.o $(CFLAGS)

Solvers/CMAES/src/cmaes.o: Solvers/CMAES/src/cmaes.c
	$(CC) -c Solvers/CMAES/src/cmaes.c -o Solvers/CMAES/src/cmaes.o $(CFLAGS)
	
Solvers/CMAES/src/boundary_transformation.o: Solvers/CMAES/src/boundary_transformation.c
	$(CC) -c Solvers/CMAES/src/boundary_transformation.c -o Solvers/CMAES/src/boundary_transformation.o $(CFLAGS)
#GOMEA
	
Solvers/RV-GOMEA/RV-GOMEA.o: Solvers/RV-GOMEA/RV-GOMEA.c
	$(CC) -c Solvers/RV-GOMEA/RV-GOMEA.c -o Solvers/RV-GOMEA/RV-GOMEA.o $(CFLAGS)

Solvers/util/FOS.o: Solvers/util/FOS.c
	$(CC) -c Solvers/util/FOS.c -o Solvers/util/FOS.o $(CFLAGS)

Solvers/util/SO_optimization.o: Solvers/util/SO_optimization.c
	$(CC) -c Solvers/util/SO_optimization.c -o Solvers/util/SO_optimization.o $(CFLAGS)

Solvers/util/Tools.o: Solvers/util/Tools.c
	$(CC) -c Solvers/util/Tools.c -o Solvers/util/Tools.o $(CFLAGS)

Solvers/problemDef.o: Solvers/problemDef.c
	$(CC) -c Solvers/problemDef.c -o Solvers/problemDef.o $(CFLAGS)
#LBFGS

Solvers/LBFGSB/driver1.o: Solvers/LBFGSB/driver1.c
	$(CC) -c Solvers/LBFGSB/driver1.c -o Solvers/LBFGSB/driver1.o $(CFLAGS)

Solvers/LBFGSB/lbfgsb.o: Solvers/LBFGSB/lbfgsb.c
	$(CC) -c Solvers/LBFGSB/lbfgsb.c -o Solvers/LBFGSB/lbfgsb.o $(CFLAGS)

Solvers/LBFGSB/linesearch.o: Solvers/LBFGSB/linesearch.c
	$(CC) -c Solvers/LBFGSB/linesearch.c -o Solvers/LBFGSB/linesearch.o $(CFLAGS)

Solvers/LBFGSB/linpack.o: Solvers/LBFGSB/linpack.c
	$(CC) -c Solvers/LBFGSB/linpack.c -o Solvers/LBFGSB/linpack.o $(CFLAGS)

Solvers/LBFGSB/miniCBLAS.o: Solvers/LBFGSB/miniCBLAS.c
	$(CC) -c Solvers/LBFGSB/miniCBLAS.c -o Solvers/LBFGSB/miniCBLAS.o $(CFLAGS)

Solvers/LBFGSB/print.o: Solvers/LBFGSB/print.c
	$(CC) -c Solvers/LBFGSB/print.c -o Solvers/LBFGSB/print.o $(CFLAGS)

Solvers/LBFGSB/subalgorithms.o: Solvers/LBFGSB/subalgorithms.c
	$(CC) -c Solvers/LBFGSB/subalgorithms.c -o Solvers/LBFGSB/subalgorithms.o $(CFLAGS)

Solvers/LBFGSB/timer.o: Solvers/LBFGSB/timer.c
	$(CC) -c Solvers/LBFGSB/timer.c -o Solvers/LBFGSB/timer.o $(CFLAGS)
#Bayes

Solvers/functions.o: Solvers/functions.c
	$(CC) -c Solvers/functions.c -o Solvers/functions.o $(CFLAGS)

Gaussian/gaussian.o: Gaussian/gaussian.c
	$(CC) -c Gaussian/gaussian.c -o Gaussian/gaussian.o $(CFLAGS)

Bayes/bayesTest.o: Bayes/bayesTest.c
	$(CC) -c Bayes/bayesTest.c -o Bayes/bayesTest.o $(CFLAGS)
#Cardio

Cardio/cardio.o: Cardio/cardio.c
	$(CC) -c Cardio/cardio.c -o Cardio/cardio.o $(CFLAGS)

Cardio/FEM.o: Cardio/FEM.c
	$(CC) -c Cardio/FEM.c -o Cardio/FEM.o $(CFLAGS)

Cardio/FEMTests.o: Cardio/FEMTests.c
	$(CC) -c Cardio/FEMTests.c -o Cardio/FEMTests.o $(CFLAGS)

Cardio/cell.o: Cardio/cell.c
	$(CC) -c Cardio/cell.c -o Cardio/cell.o $(CFLAGS)
	
Cardio/geometry.o: Cardio/geometry.c
	$(CC) -c Cardio/geometry.c -o Cardio/geometry.o $(CFLAGS)
#home

home.o: home.c
	$(CC) -c home.c -o home.o $(CFLAGS)
	
clean :
	rm bayesOpt $(OBJ)
