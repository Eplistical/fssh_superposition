CXX = g++ 
OPT = -O3 
LIBS = -lboost_program_options -lfftw3 -lopenblas -lpthread -lgfortran 

all : exact_1d fssh_1d ehrenfest_1d exact_1d_model2 fssh_1d_model2 ehrenfest_1d_model2 fssh_1d_model2_phase_corr

exact_1d: exact_1d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_1d: fssh_1d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_1d: ehrenfest_1d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_1d_model2: exact_1d_model2.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_1d_model2: fssh_1d_model2.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_1d_model2: ehrenfest_1d_model2.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_1d_model2_phase_corr: fssh_1d_model2_phase_corr.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)
