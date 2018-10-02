
 #cuda compiler and flags
NVCC = nvcc
FLAG1 = -arch 'compute_35' -code 'sm_35' #-maxrregcount 32
PROF = --ptxas-options=-v  -lineinfo
FLAG2 = --use_fast_math
MAT = -ftz=true -prec-div=false
#CACHE = -Xptxas -dlcm=cg # for caching only in L2
#PTX = -Xptxas -abi=no #for removing the ABI that consent printf in the stack

# linking libraries and includes for heterogeneous cuda-fortran programming
CUDA = /usr/local/cuda-8.0
INC = -I$(CUDA)/include
LIB = -L$(CUDA)/lib64 -lcudart -lc -lstdc++ -lcurand #-lcuda

LIBD =   libacml.a
COMP  = mpif90 -132 -r8
EXEC = g++
FLAGC = -c

turb5:  decomp_2d.o sfft.o ders_n.o main_comp.o fix_var.o mom.o solver.o param.txt read_tables.o memory_copy.o spline.o mc_gpu.o param.h NarrowBand.h
	$(COMP) $(DBGc) $(OPT) main_comp.o  decomp_2d.o io.o fix_var.o mom.o solver.o sfft.o ders_n.o read_tables.o memory_copy.o spline.o mc_gpu.o $(LIBD) $(INC) $(LIB) -o volh
main_comp.o : main_comp.f param.txt
	$(COMP) $(DBGc) $(FLAGC) $(OPT) main_comp.f 
mom.o : mom.f
	$(COMP) $(DBGc) $(FLAGC) $(OPT) mom.f 
spline.o: spline.cpp
	$(EXEC) $(FLAGC) $(OPT) spline.cpp
mc_gpu.o: mc_gpu.cu NarrowBand.h param.h
	$(NVCC) $(DBG) $(OPT) $(MAT) $(PROF) $(PTX) $(FLAG2) $(CACHE)  $(FLAGC) $(OPT) mc_gpu.cu $(FLAG1)
read_tables.o: read_tables.cpp
	$(EXEC) $(DBGc)  $(FLAGC) $(OPT) read_tables.cpp
memory_copy.o: memory_copy.cu
	$(NVCC) $(DBG) $(OPT) $(MAT) $(PROF) $(PTX) $(FLAG2) $(CACHE)  $(FLAGC) $(OPT) memory_copy.cu $(FLAG1)
solver.o : solver.f param.txt
	$(COMP) $(DBGc) $(FLAGC) $(OPT) solver.f 
sfft.o : sfft.f
	$(COMP) $(DBGc) $(FLAGC) $(OPT) sfft.f
ders_n.o: ders_n.f
	$(COMP) $(DBGc) $(FLAGC) $(OPT) ders_n.f
fix_var.o : fix_var.f
	$(COMP) $(DBGc) $(FLAGC) $(OPT) fix_var.f
decomp_2d.o: decomp_2d.f90 
	$(COMP)  -DDOUBLE_PREC -DOVERWRITE -O3  -cpp  -c decomp_2d.f90 io.f90
clean:
	rm -rf *.o *.mod
