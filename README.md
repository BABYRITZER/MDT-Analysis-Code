# MDT Analysis Code
 
## To build analysis run the command: 

`` g++ event_dict.cpp main.cpp fitt0s.cpp radiustimefunction.cpp line_fitting.cpp gransac_implementations.cpp load_from_root_file.cpp efficiency_calc.cpp -o analysis `root-config --cflags --
libs` -lMinuit -fopenmp ``

## To build the chamber viewer run the command: 

`` g++ -Wall -Wextra event_dict.cpp radiustimefunction.cpp chambers.cpp -o drawchambs `root-config --cflags --libs` -lMinuit -fopenmp ``

## To build the efficiency graph maker run the command:

`` g++ eff_graphs.cpp -o effgraphs `root-config --cflags --libs` ``

## The RANSAC implementation was taken from:
https://github.com/drsrinathsridhar/GRANSAC/. 
