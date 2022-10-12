#include <string>
#include <assert.h>
#include <bits/stdc++.h>
#include <utility>
#include <mpi.h>
#include "randomizer.hpp"
using namespace std;

bool comparator(pair<int,int> p1,pair<int,int> p2){
	return p1.second>p2.second || (p1.second==p2.second && p1.first<p2.first);
}

class exec{
	int user;
	int initial;
	int current;
	int num_steps;
	int num_walks;
	unordered_map<int,int> scores;
	
	public:
	
	exec(int userv,int initialv,int currentv,int num_stepsv,int num_walksv,unordered_map<int,int> scoresv){
		user=userv;
		initial=initialv;
		current=currentv;
		num_steps=num_stepsv;
		num_walks=num_walksv;
		scores=scoresv;
	};
	
	void setNew(int initialv){
		initial=initialv;
		current=initialv;
	};
	
	vector<pair<int,int>> getOrderedScores(){
		vector<pair<int,int>> pairs;
		scores[user]=0;
		for (auto iter : scores){
			pairs.push_back(make_pair(iter.first,iter.second));
		}
		sort(pairs.begin(),pairs.end(),comparator);
		return pairs;
	};
	
	void GenScores(vector<vector<int>> edges, Randomizer random_generator){
		//cout<<user<<" "<<initial<<" "<<num_steps<<" "<<num_walks<<"\n";
		for(int j=0;j<num_walks;j++){
			current=initial;
			for(int i=0;i<num_steps;i++){
				int num_child=edges[current].size();
				if(num_child > 0){
					//Called only once in each step of random walk using the original node id 
					//for which we are calculating the recommendations
					int next_step = random_generator.get_random_value(user);
					//Random number indicates restart
					if(next_step<0){
						current=initial;
						//std::cout << "Restart \n";
					}else{
						//Deciding next step based on the output of randomizer which was already called
						int child = next_step % num_child; //is the number of child of the current node
						current=edges[current][child];
						//std::string s =  std::to_string(child)+"\n";
						//std::cout << s;
					}
				}
				else{
					current=initial;
					//std::cout << "Restart \n";
				}
				if (scores.find(current)==scores.end()){
					scores[current]=1;
				}
				else{
					scores[current]++;
				}
			}
		}
		scores[initial]=numeric_limits<int>::min();
	};
};


int main(int argc, char* argv[]){
    assert(argc > 8);
    std::string graph_file = argv[1];
    int num_nodes = std::stoi(argv[2]);
    int num_edges = std::stoi(argv[3]);
    float restart_prob = std::stof(argv[4]);
    int num_steps = std::stoi(argv[5]);
    int num_walks = std::stoi(argv[6]);
    int num_rec = std::stoi(argv[7]);
    int seed = std::stoi(argv[8]);
    
    //Only one randomizer object should be used per MPI rank, and all should have same seed
    Randomizer random_generator(seed, num_nodes, restart_prob);
    auto start = std::chrono::steady_clock::now();
    
    vector<vector<int>> edges;
    for(int i=0;i<num_nodes;i++){
    	vector<int> v;
    	edges.push_back(v);
    }
    
    
    FILE* file_;
    unsigned char buffer[4];
    file_ = fopen(graph_file.c_str(), "rb");
    cout << "started reading"<<endl;
    int i=0;
    while (i<num_edges) // to read file
    {
        // function used to read the contents of file
        fread(buffer, 1, 4, file_);
        int from = (int)buffer[3] | (int)buffer[2]<<8 | (int)buffer[1]<<16 | (int)buffer[0]<<24;
        fread(buffer, 1, 4, file_);
        int to = (int)buffer[3] | (int)buffer[2]<<8 | (int)buffer[1]<<16 | (int)buffer[0]<<24;
        edges[from].push_back(to);
        //cout<<i<<","<<from<<","<<to<<endl;
        i++;
    }
    fclose(file_);
    
    int rank, size;

    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
//    if (rank==0){
//    	FILE* file_out_=fopen("output.dat","wb");
//    	for(i=0;i<num_nodes*(2*num_rec+1)*4;i++){
//    		char* s="NULL";
//			fwrite(s,4,1,file_out_);	
//    	}
//    fclose(file_out_);
//    }
//    
//    cout<<"empty file created\n";
    
    //string s="output"+to_string(rank)+".dat";
    MPI_Barrier( MPI_COMM_WORLD);
    
    //MPI_Barrier(MPI_COMM_WORLD);
    FILE* file_out=fopen("output.dat","wb");
    cout<<"file opened\n";
    for(int n=(rank*num_nodes)/size;n<num_nodes && n<((rank+1)*num_nodes)/size;n++){
    	vector<pair<int,int>> scores;
    	unordered_map<int,int> count;
    	exec iter=exec(n,0,0,num_steps,num_walks,count);
    	for(int j=0;j<edges[n].size();j++){
    		iter.setNew(edges[n][j]);
			iter.GenScores(edges,random_generator);
			//cout<<n<<","<<edges[n][j]<<" done\n";
		}
		scores=iter.getOrderedScores();
		//cout<<n<<" "<<edges[n].size()<<",";
		fseek (file_out, n*(2*num_rec+1)*4, SEEK_SET);
		int m=edges[n].size();
		unsigned char buffer[4];
		buffer[0]=m>>24;
		buffer[1]=m>>16;
		buffer[2]=m>>8;
		buffer[3]=m;
		fwrite(buffer,4,1,file_out);
		for(int j=0;j<num_rec;j++){
			if (j<scores.size() && scores[j].second>0){
				//cout<<scores[j].first<<" "<<scores[j].second<<",";
				m=scores[j].first;
				buffer[0]=m>>24;
				buffer[1]=m>>16;
				buffer[2]=m>>8;
				buffer[3]=m;
				fwrite(buffer,4,1,file_out);
				m=scores[j].second;
				buffer[0]=m>>24;
				buffer[1]=m>>16;
				buffer[2]=m>>8;
				buffer[3]=m;
				fwrite(buffer,4,1,file_out);
			}
			else{
				//cout<<"NULL NULL"<<",";
				unsigned char* buffer;
				char* s="NULL";
				fwrite("NULL",4,1,file_out);
				fwrite("NULL",4,1,file_out);
			}
		}
		//cout<<"\n";
		//MPI_Barrier( MPI_COMM_WORLD);
	
	//MPI_Barrier( MPI_COMM_WORLD);
	}
	fclose(file_out);
	cout<<"seperate done\n";
	MPI_Finalize();
	auto stop = std::chrono::steady_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
	cout<<duration.count()<<endl;
}
