#find the correct mpifort command
FC1 := $(shell { command -v mpifort; } )
FC2 := $(shell { command -v mpif90; } )

FC := mpifort
ifeq ($(strip $(FC1)),)
  FC := mpif90
endif

#master
objs_master := mydata.o master.o

#worker
objs_worker := mydata.o worker.o



all: build

build: master worker

master: mydata.o master.o
	$(FC) $(FFLAGS) $^ -o $@

worker: mydata.o worker.o
	$(FC) $(FFLAGS) $^ -o $@

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	@rm -f *.mod *.o worker master

run:
	@mpirun -n $(nm) master : -n $(nw) worker