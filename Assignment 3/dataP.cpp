#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

int main(int argc, char* argv[]){

    std::string inpath = argv[1];
    std::string outpath = argv[2];
    std::string infiles[2]={"vect.txt","index.txt"};
    std::string outfiles[2]={"output.dat","outi.dat"};
    #pragma omp parallel for shared(infiles, outfiles, inpath, outpath)
    for (int id =0; id < 2; id++)
    {
    // int id=omp_get_thread_num();
    string infile=inpath+infiles[id];
    ifstream file;
    file.open(infile.c_str());
    int c=0;
    int count=0;
    string line;
    getline(file,line);
    cout<<"reading first lin \n";
    while(c<line.length()){
        if (line[c]==' '){
            count++;
        }
	//cout<<line[c];
	c++;
    }
    int d=count+1;
    file.seekg(0);
    string word;
    cout<<d<<endl;
	FILE *of;
   string outfile=outpath+outfiles[id];
    const char* s = (outfile).c_str();
    of = fopen(s, "wb");

    int l=0;
    count=0;
    while (file >> word)
    {
        if (count==d-1){
            count=0;
            l++;
            //cout<<l<<endl;
        }
        else
        count++;
        if(id==1)
        {
	        int f = stoi(word);
	        fwrite(&f , 4, 1, of);
        }
        else
        {
            float f = stof(word);
	        fwrite(&f , 4, 1, of);
        }
    }
    file.close();
    fclose(of);
    string str = to_string(id);
    s = (str+"_mdata.dat").c_str();
    cout << "s: " << s << endl;
    of = fopen(s, "wb");
    int bufi;
    bufi=l;
    cout << "bufi: " << bufi << endl;
    cout << "of: " << of << endl;
    fwrite(&bufi , 4, 1, of);
    bufi=d;
    fwrite(&bufi , 4, 1, of);
    cout<<l<<','<<d<<endl;
    }
    #pragma omp taskwait
    return 0;
}

