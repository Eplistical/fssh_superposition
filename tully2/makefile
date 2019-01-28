CXX = g++ 
OPT = -O3 
LIBS = -lboost_program_options -lfftw3 -lopenblas -lpthread -lgfortran 

all : exact_1d fssh_1d fssh_1d_phase_corr

exact_1d: exact_1d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_1d: fssh_1d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_1d_phase_corr: fssh_1d_phase_corr.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

