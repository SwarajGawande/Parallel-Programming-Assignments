g++ -std=c++11 -fopenmp dataP.cpp -o data 
mpic++ -std=c++17 -o code code.cpp -fopenmp
g++ -o cnvrt convert.cpp
