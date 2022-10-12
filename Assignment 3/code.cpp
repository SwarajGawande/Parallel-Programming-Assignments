#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

float cosine_dist(vector<float> &a, vector<float> &b);

bool cmp(const vector<int>& a, const vector<int>& b){
    return (a[a.size() - 1] < b[b.size() - 1]);
}

void heap_trim(vector<pair<float, int>> *v, int k)
{
    while(v->size() > k)
    {
        pop_heap(v->begin(), v->end());
        v->pop_back();
    }
}

vector<pair<float, int>> SearchLayer(vector<float> &q, vector<pair<float, int>> &candidates, vector<int> &indptr, vector<int> &index, vector<int> &level_offset, int lc, unordered_map<int, bool> &visited, vector<vector<float>> &vect, int topk, int user)
{
    vector<pair<float, int>> top_k(candidates);
    while (candidates.size() > 0)
    {
        int ep = candidates[0].second;
        //cout << "ep : " << ep << endl;
        pop_heap(candidates.begin(), candidates.end());
        candidates.pop_back();
        int start = indptr[ep] + level_offset[lc];
        int end = indptr[ep] + level_offset[lc + 1];
        //cout << "start " << start << " end " << end << endl;
        
        for (int p_ = start; p_ < end; p_++)
        {   
            int px=index[p_];
            //cout << "test_px: " << px << endl;
            if (visited.find(px) != visited.end() || px == -1)
            {
                continue;
            }
            visited[px] = true;
            float _dist = cosine_dist(q, vect[px]);
            //cout << _dist <<  " " << px << endl;
            if (_dist > top_k.front().first && top_k.size() >= topk)
            {
                continue;
            }
            top_k.push_back({_dist, px});
            push_heap(top_k.begin(), top_k.end());
            heap_trim(&top_k, topk);
            candidates.push_back({_dist, px});
            //cout << "px: " << px << endl;
            push_heap(candidates.begin(), candidates.end());
        }
    }
    
    return top_k;
}

float cosine_dist(vector<float> &a, vector<float> &b)
{
    float c_dist = 0.0;
    float mag_a = 0.0;
    float mag_b = 0.0;
    int sz = a.size();

    for (int i = 0; i < sz; i++)
    {
        c_dist += a[i]*b[i];
        mag_a += a[i]*a[i];
        mag_b += b[i]*b[i];
    }

    mag_a = sqrt(mag_a);
    mag_b = sqrt(mag_b);

    return 1.0 - (c_dist / (mag_a * mag_b));
}


vector<pair<float, int>> QueryHNSW(vector<float> &q, int topk, int ep, vector<int> &indptr, vector<int> &index, vector<int> &level_offset, int max_level, vector<vector<float>> &vect, int user)
{
    vector<pair<float, int>> top_k;
    // for (int i=0;i<topk-1;i++)
    // top_k.push_back({2.0,-1});
    top_k.push_back({cosine_dist(q, vect[ep]), ep});
    make_heap(top_k.begin(), top_k.end());
    //push_heap(top_k.begin(), top_k.end());

    unordered_map<int, bool> visited;
    visited[ep] = true;
    for (int level = max_level; level >=0; level--)
    {
        //cout << level << endl;
        top_k = SearchLayer(q, top_k, indptr, index, level_offset, level, visited, vect, topk, user);
    }
    return top_k;
}


int main(int argc, char **argv)
{
    int rank, size;

    int topk = stoi(argv[2]);

    string out_path = argv[1];
    string bin_filename = out_path + "/output.dat";
    string ep_file = out_path + "/ep.txt";
    string index_file = out_path + "/outi.dat";
    string indptr_file = out_path + "/indptr.txt";
    string level_file = out_path + "/level.txt";
    string level_offset_file = out_path + "/level_offset.txt";
    string max_level_file = out_path + "/max_level.txt";
    string user_file = argv[3];
    string meta_data = out_path + "/0_mdata.dat";
    string user_output_file = "user_output_file.dat";


    const char *file = bin_filename.c_str();
    const char *file1 = ep_file.c_str();
    const char *file2 = index_file.c_str();
    const char *file3 = indptr_file.c_str();
    const char *file4 = level_file.c_str();
    const char *file5 = level_offset_file.c_str();
    const char *file6 = max_level_file.c_str();
    const char *file7 = user_file.c_str();
    const char *file8 = meta_data.c_str();
    const char *file9 = user_output_file.c_str();

    FILE *ifp;

    ifp = fopen(file8, "rb");
    int D;
    fread(&D, 4, 1, ifp);
    fread(&D, 4, 1, ifp);
    fclose(ifp);

    //cout << "0" << endl;

    ifp = fopen(file, "rb");

    vector<vector<float>> v_node;

    float buf;

    int count = 0;

    vector<float> row(D);
    while(fread(&buf, 4, 1, ifp) != 0)
    {
        row[count] = buf;
        count++;
        if (count == D)
        {
            v_node.push_back(row);
            count = 0;
        }
    }

    fclose(ifp);

    fstream epfile;
    epfile.open(file1);
    int ep;
    string eps;
    epfile >> eps;
    ep=stoi(eps);
    epfile.close();

    ifp = fopen(file2, "rb");
    int idx;
    vector<int> index;

    while(fread(&idx, 4, 1, ifp) != 0)
    {
        index.push_back(idx);
    }
   
    fclose(ifp);


    fstream ipfile;
    ipfile.open(file3);
    string ips;

    int indptr_idx;
    vector<int> indptr;

    while(ipfile >> ips)
    {
        indptr_idx = stoi(ips);
        indptr.push_back(indptr_idx);
    }

    ipfile.close();

    fstream lofile;
    lofile.open(file5);
    string los;

    int level_offset_idx;
    vector<int> level_offset;

    while(lofile >> los)
    {
        level_offset_idx = stoi(los);
        level_offset.push_back(level_offset_idx);
    }

    lofile.close();

    fstream maxfile;
    maxfile.open(file6);
    int max_level;
    string ms;
    maxfile >> ms;
    max_level = stoi(ms);
    maxfile.close();

    vector<vector<float>> user;

    float user_buf;

    int user_count = 0;
    
    int U = 0;

    vector<float> user_row(D);
    fstream user_file1;
    user_file1.open(file7);
    while(user_file1 >> user_buf)
    {
        user_row[user_count] = user_buf;
        user_count++;
        if (user_count == D)
        {
            user.push_back(user_row);
            user_count = 0;
            U++;
        }
    }

    user_file1.close();
    fstream file_user("user_output_file.dat" , ios::out|std::ios::binary);
    // Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int lb,rb;
    lb = (rank*U)/size;
    rb = ((rank+1)*U)/size;

    // lb = 0;
    // rb = U;
    //ifp = fopen(file9, "wb");
    vector<vector<pair<float,int>>> top_kv;
    vector<vector<int>> top_ks;
    //int offset=0;
    #pragma omp parallel for 
    for (int i = lb; i < rb; i++)
    {
       // cout << i << " started" << endl;
        vector<pair<float, int>> top_k = QueryHNSW(user[i], topk, ep, indptr, index, level_offset, max_level, v_node, i);
	top_k.push_back({3, i});
	sort(top_k.begin(), top_k.end());
	vector<int> tempv;
	for (int j = 0; j < top_k.size(); j++)
	{
		tempv.push_back(top_k[j].second);
	}
	
        #pragma omp critical
        {
 		top_ks.push_back(tempv);    
        }
        //cout<<i <<" completed"<<endl;
    }
    cout<<"writing now\n";
    #pragma omp taskwait
    MPI_Barrier(MPI_COMM_WORLD);
    sort(top_ks.begin(), top_ks.end(), cmp);

    //for (int i = lb; i < rb; i++)
    //{
      //  sort(top_kv[i].begin(), top_kv[i].end());
        //for (int j = 0; j < topk; j++)
        //{
            //int offset = topk*i*4 + 4*j;
           // fseek(ifp, offset, SEEK_SET);
            //fwrite(&top_kv[i-lb][j].second, 1, 4, ifp);
            //offset+=4;
            //cout<<i<<", "<<j<<" "<<top_kv[i][j].first<<","<<top_kv[i][j].second<<endl;
        //}
    //}

    file_user.seekg(4*((rank*U)/size)*topk, ios::beg);
    for (int i = 0; i < top_ks.size(); i++)
    {
        file_user.write(reinterpret_cast<const char*>(&top_ks[i][0]), topk*sizeof(int));
    	//cout<<""<<endl;
	}
	
    //fclose(ifp);
    file_user.close();
    

    MPI_Finalize();

    return 0;
}



