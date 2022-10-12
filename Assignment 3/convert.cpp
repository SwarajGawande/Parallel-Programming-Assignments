#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <bits/stdc++.h>

using namespace std;

int main(int argc, char* argv[]){

    std::string infile = "user_output_file.dat";
    FILE * ifp = fopen(infile.c_str(), "rb");

    string outfile=argv[1];
    ofstream file;
    file.open(outfile.c_str());
    int buf;
    int D=stoi(argv[2]);
    int count=0;
    while(fread(&buf, 4, 1, ifp) != 0)
    {
        count++;
        string s=to_string(buf);
        file << (s);
        if (count == D)
        {
            file<<"\n";
            count = 0;
        }
        else{
            file<<" ";
        }
    }
    file.close();
    fclose(ifp);

    return 0;
}

