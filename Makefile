FC = mpif90 -O2 ## gfotran
#FC = mpiifort -O3 -xHOST -ipo -ip ## intel
#FC = mpifrtpx -O3 -Kfast ##FX100@Nagoya

LN = #-llapack -lblas
#LN = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm

VPATH = src:object
SRC = $(shell cd src ;ls *.f90 ;cd ..)
OBJ = $(SRC:.f90=.o)
OBJ_dir = $(addprefix object/,$(OBJ))

PROG = spin-boson

$(PROG):global_variables.o CTEF_module.o $(OBJ)
	$(FC) -o $(PROG) $(OBJ_dir) $(LN)

main.o:main.f90
	$(FC) -c $< $(LN);mv $@  object 
%.o:%.f90
	$(FC) -c $< $(LN);mv $@  object 


clean:
	rm  -f  object/*.o  *.mod ${PROG}
clean_complete:
	rm  -f *~  */*~ */*/*~ object/*.o  */*.mod *.mod ${PROG} */#* *.out *.log
