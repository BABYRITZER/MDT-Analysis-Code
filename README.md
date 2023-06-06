# MDT Analysis Code
 
to build analysis: 
g++ -Wall -Wextra event_dict.cpp main.cpp fitt0s.cpp radiustimefunction.cpp line_fitting.cpp gransac_implementations.cpp load_from_root_file.cpp -o analysis `root-config --cflags --libs` -lMinuit -fopenmp

to build the chamber viewer:
g++ -Wall -Wextra event_dict.cpp radiustimefunction.cpp chambers.cpp -o drawchambs `root-config --cflags --libs` -lMinuit -fopenmp