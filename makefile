CC 	   = mpic++
#FFTW_LOC = /opt/seas/pkg/gcc/fftw3/mpi/double/3.3.7
FFTW_LOC = ${HOME}/Install/fftw3
EIGEN_LOC = ${HOME}/Install/eigen
CFLAGS     = -I${FFTW_LOC}/include -I${EIGEN_LOC} -O3 -Wno-unused-result -Wno-write-strings
LIBS      = -lm -lfftw3_mpi -lfftw3 -O3 -L${FFTW_LOC}/lib

#CC 	   = icpc 
#FFTW_LOC = ${HOME}/Install/fftw3a-threaded
#CFLAGS     = -O3 -openmp -I${FFTW_LOC}/include
#LIBS      = -openmp -O3 -lfftw3_threads -lfftw3 -lpthread -L${FFTW_LOC}/lib

#############################################################################
# nothing should be changed below here

SRCS = main.cpp array_utils.cpp die.cpp  random.cpp grid_utils.cpp \
			 fftw_mpi_wrappers.cpp initialize.cpp config_utils.cpp io_utils.cpp \
			 update_positions.cpp forces.cpp integ_utils.cpp read_input.cpp \
			 bonded.cpp calc_unb.cpp  charge_forces.cpp anneal_utils.cpp \
			 mpi_utils.cpp communicate_utils.cpp angles.cpp \
			 orientation_utils.cpp ms_forces.cpp \
			 pair_style.cpp pair_style_gaussian.cpp \
       
       
			 


OBJS = ${SRCS:.cpp=.o}

.cpp.o:
	${CC} ${CFLAGS} ${DFLAGS} -c  $<

dmft_mpi:  ${OBJS}
	$(CC) ${CFLAGS} ${DFLAGS} -o $@ ${OBJS} $(LIBS)

clean:
	rm -f *.o
	rm -f dmft_mpi
	rm -f *~

