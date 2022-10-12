#include <iostream>
#include <bits/stdc++.h>
using namespace std;

// float calcAvg(int m,int n, vector<vector<tuple<int,int,int>>> v, int x,int y){
//     int sum=0;
//     for(int i=0;i<m;i++){
//         for (int j=0;j<n;j++){
//             int R=get<0>(v[x+i][y+j]);
//             int G=get<1>(v[x+i][y+j]);
//             int B=get<2>(v[x+i][y+j]);
//             sum+=(R+G+B);
//         }
//     }
//     float avg=sum/((float) (m*n*3));
//     return avg;
// }

__host__ void calcAvg(int m,int n, int* v, int x,int y, float* avg, int N){
    int sum=0;
    for(int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            int R=v[(x+i)*3*N+(y+j)*3];
            int G=v[(x+i)*3*N+(y+j)*3+1];
            int B=v[(x+i)*3*N+(y+j)*3+2];
            sum+=(R+G+B);
            //cout<<R<<","<<G<<","<<B<<endl;
        }
    }
    *avg=sum/((float) (m*n*3));
}


pair<int,int> rotate(int x,int y, float angle){
   double const PI = 3.14159265358979323;
   float radians = angle * (PI / 180.0f);   // convert degrees to radians
   int nx = x * cos(radians) - y * sin(radians); 
   int ny = x * sin(radians) + y * cos(radians);
   return make_pair(nx,ny);
}

void rotateimg(int m,int n, float angle, float * v){
    for(int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            pair<int,int> p = rotate(i,j,angle);
            v[i*2*n+j*2]=p.first;
            v[i*2*n+j*2+1]=p.second;
            //cout<<i<<","<<j<<endl;
        }
    }
}

__global__ void query(int M,int N, int m, int n, float avgq, float* poso, float* posp, float* posn, int* img, int* imgq, float* soln, float th1, float th2){
    int index=blockIdx.x*blockDim.x+threadIdx.x;
    int i = index/N;
    int j = index% N;
    float avgo=0.0;
    float avgp=0.0;
    float avgn=0.0;
    bool valid = false;
    //0 degree
    if (!(i+poso[2*n*(m)-2]>M-1 || j+poso[2*n*(m)-1]>N-1)){
         int sum=0;
         int xmin=i;
         int ymin=j;
         int xmax=i+m;
         int ymax=j+n;
         for(int mi=xmin;mi<xmax;mi++){
             for (int ni=ymin;ni<ymax;ni++){
                int R=img[(M-1-mi)*3*N+(ni)*3];
                int G=img[(M-1-mi)*3*N+(ni)*3+1];
                int B=img[(M-1-mi)*3*N+(ni)*3+2];
                sum+=(R+G+B);
                //cout<<R<<","<<G<<","<<B<<endl;
            }
         }
        avgo=sum/((float) (m*n*3));
        if(abs(avgo-avgq)<=th2){
            long sum=0;
            for(int mi=0;mi<m;mi++){
                for(int ni=0;ni<n;ni++){
                    int R=abs(img[(M-1-i-mi)*3*N+(j+ni)*3]-imgq[(m-1-mi)*n*3+ni*3]);
                    int G=abs(img[(M-1-i-mi)*3*N+(j+ni)*3+1]-imgq[(m-1-mi)*n*3+ni*3+1]);
                    int B=abs(img[(M-1-i-mi)*3*N+(j+ni)*3+2]-imgq[(m-1-mi)*n*3+ni*3+2]);
                    sum=sum+R*R+G*G+B*B;
                }
            }
            float rmsg=sqrt(sum/((float)m*n*3));
            if (rmsg<=th1){
                valid=true;
                soln[3*index]=rmsg;
                soln[3*index+1]=avgo;
                soln[3*index+2]=1;
                //printf("%d,%d,%f,%d\n",i,j,rmsg,index);
            }
        }
    }

    // //45 degree

    // if (!( i+posp[2*(n-1)] <0 || j+posp[2*(n-1)+1]>N-1 || i+posp[2*n*m-2]>M-1 || j+posp[2*n*m-1]>N-1 || i+posp[2*n*(m-1)]>M-1 || j+posp[2*n*(m-1)+1]<0)){
    //     int sum=0;
    //     int xmin=i;
    //     int ymin=(int) j-(m/sqrt(2.0));
    //     int xmax=(int) i+((m+n)/sqrt(2.0))+1;
    //     int ymax=(int) j+(n/sqrt(2.0))+1;
    //     if(i==985 && j==239)
    //     printf("%d,%d,%d,%d\n",xmin,xmax,ymin,ymax);
    //     for(int mi=xmin;mi<xmax;mi++){
    //          for (int ni=ymin;ni<ymax;ni++){
    //             int R=img[(M-1-mi)*3*N+(ni)*3];
    //             int G=img[(M-1-mi)*3*N+(ni)*3+1];
    //             int B=img[(M-1-mi)*3*N+(ni)*3+2];
    //             sum+=(R+G+B);
    //             //cout<<R<<","<<G<<","<<B<<endl;
    //         }
    //     }
    //     avgp=sum/((float) ((xmax-xmin+1)*(ymax-ymin+1)*3));
    //     if(i==985 && j==239)
    //     printf("%d,%d,%f\n",i,j,avgp);
    //     if(abs(avgp-avgq)<th2){
    //         long sum=0;
    //         for(int mi=0;mi<m;mi++){
    //             for(int ni=0;ni<n;ni++){
    //                 int x=(int) i+posp[mi*n*2+ni*2];
    //                 int y=(int) j+posp[mi*n*2+ni*2+1];
    //                 float R=0;
    //                 float G=0;
    //                 float B=0;
    //                 for(int k=0;k<3;k++){
    //                     int x_=x+k/2;
    //                     int y_=y+k%2;
    //                     R+=(img[(M-1-x_)*3*N+3*y_]*sqrt(pow(x_-i-posp[mi*n*2+ni*2],2)+pow(y_-j-posp[mi*n*2+ni*2+1],2)));
    //                     G+=(img[(M-1-x_)*3*N+3*y_+1]*sqrt(pow(x_-i-posp[mi*n*2+ni*2],2)+pow(y_-j-posp[mi*n*2+ni*2+1],2)));
    //                     B+=(img[(M-1-x_)*3*N+3*y_+2]*sqrt(pow(x_-i-posp[mi*n*2+ni*2],2)+pow(y_-j-posp[mi*n*2+ni*2+1],2)));
    //                 }
    //                 R=abs(R-imgq[(m-1-mi)*n*3+ni*3]);
    //                 G=abs(G-imgq[(m-1-mi)*n*3+ni*3+1]);
    //                 B=abs(B-imgq[(m-1-mi)*n*3+ni*3+2]);
    //                 sum=sum+R*R+G*G+B*B;
    //             }
    //         }
    //         float rmsg=sqrt(sum/((float)m*n*3));
    //         if (rmsg<th1 && (!valid || rmsg<soln[3*index])){
    //             valid=true;
    //             soln[index]=rmsg;
    //             printf("%d,%d,%f\n",i,j,rmsg);
    //         }
    //     }
    // }

    // //-45 degree
    // if (!( i+posn[2*(n-1)] <0 || j+posn[2*(n-1)+1]>N-1 || i+posn[2*n*m-2]>M-1 || j+posn[2*n*m-1]>N-1 || i+posn[2*n*(m-1)]>M-1 || j+posn[2*n*(m-1)+1]<0)){
        
    //     int sum=0;
    //     int xmin=(int) i-(n/sqrt(2.0));
    //     int ymin=(int) j;
    //     int xmax=(int) i+(m/sqrt(2.0))+1;
    //     int ymax=(int) j+((m+n)/sqrt(2.0))+1;
    //     for(int mi=xmin;mi<xmax;mi++){
    //          for (int ni=ymin;ni<ymax;ni++){
    //             int R=img[(M-1-mi)*3*N+(ni)*3];
    //             int G=img[(M-1-mi)*3*N+(ni)*3+1];
    //             int B=img[(M-1-mi)*3*N+(ni)*3+2];
    //             sum+=(R+G+B);
    //             //cout<<R<<","<<G<<","<<B<<endl;
    //         }
    //     }
    //     avgn=sum/((float) (m*n*3));
    //     if(abs(avgn-avgq)<th2){
    //         long sum=0;
    //         for(int mi=0;mi<m;mi++){
    //             for(int ni=0;ni<n;ni++){
    //                 int x=(int) i+posn[mi*n*2+ni*2];
    //                 int y=(int) j+posn[mi*n*2+ni*2+1];
    //                 float R=0;
    //                 float G=0;
    //                 float B=0;
    //                 for(int k=0;k<3;k++){
    //                     int x_=x+k/2;
    //                     int y_=y+k%2;
    //                     R+=(img[x_*3*N+3*y_]*sqrt(pow(x_-i+posn[mi*n*2+ni*2],2)+pow(y_-j+posn[mi*n*2+ni*2+1],2)));
    //                     G+=(img[x_*3*N+3*y_+1]*sqrt(pow(x_-i+posn[mi*n*2+ni*2],2)+pow(y_-j+posn[mi*n*2+ni*2+1],2)));
    //                     B+=(img[x_*3*N+3*y_+2]*sqrt(pow(x_-i+posn[mi*n*2+ni*2],2)+pow(y_-j+posn[mi*n*2+ni*2+1],2)));
    //                 }
    //                 R=abs(R-imgq[mi*n+ni]);
    //                 G=abs(G-imgq[mi*n+ni+1]);
    //                 B=abs(B-imgq[mi*n+ni+2]);
    //                 sum=sum+R*R+G*G+B*B;
    //             }
    //         }
    //         float rmsg=sqrt(sum/((float)m*n*3));
    //         if (rmsg<th1 && (!valid || rmsg<soln[index])){
    //             valid=true;
    //             soln[index]=rmsg;
    //         }
    //     }
    // }

    if(!valid){
        soln[3*index]=-1.0;
        soln[3*index+1]=-1.0;
        soln[3*index+2]=-1.0;
    }
}

int main(int argc, char* argv[]){
    string infile=argv[1];
    string infile_q=argv[2];
    string outfile="output.txt";
    float th1=stof(argv[3]);
    float th2=stof(argv[4]);
    int n_=stoi(argv[5]);
    ifstream inf;
    inf.open(infile.c_str(),ios::in);
    string word;
    inf>>word; 
    // cout<<word;
    int M = stoi(word);
    inf>>word;
    //cout<<word;
    int N=stoi(word);
    cout<<M<<","<<N<<endl;
    // vector<vector<tuple<int,int,int>>> img;
    int * img = new int[M*N*3];
    cout<<M*N*3<<endl;
    // for(int i=0;i<M;i++){
    //     vector<tuple<int,int,int>> v(N,make_tuple(0,0,0));
    //     img.push_back(v);
    // }
    cout<<"reading started"<<endl;
    for(int i=0;i<M;i++){
        for (int j=0;j<N;j++){
            inf>>word;
            //cout<<word<<",";
            int R= stoi(word);
            img[3*N*i+3*j]=R;
            inf>>word;
            //cout<<word<<",";
            int G= stoi(word);
            img[3*N*i+3*j+1]=G;
            inf>>word;
            //cout<<word<<endl;
            int B= stoi(word);
            img[3*N*i+3*j+2]=B;
            // img[i][j]=(make_tuple(R,G,B));
            // cout<<i<<","<<j<<endl;
        }
    }
    inf.close();
    cout<<"reading done"<<endl;
    ifstream infq;
    infq.open(infile_q.c_str(),ios::in);
    infq>>word; 
    int m = stoi(word);
    infq>>word;
    int n=stoi(word);
    cout<<m<<","<<n<<endl;
    // vector<vector<tuple<int,int,int>>> imgq;
    int * imgq = new int[m*n*3];
    // for(int i=0;i<m;i++){
    //     vector<tuple<int,int,int>> v(n,make_tuple(0,0,0));
    //     imgq.push_back(v);
    // }
    cout<<"reading started"<<endl;
    for(int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            infq>>word;
            int R= stoi(word);
            imgq[3*n*i+3*j]=R;
            infq>>word;
            int G= stoi(word);
            imgq[3*n*i+3*j+1]=G;
            infq>>word;
            int B= stoi(word);
            imgq[3*n*i+3*j+2]=B;
            //cout<<i<<","<<j<<endl;
            //cout<<R<<","<<G<<","<<B<<endl;
        }
        
        //cout<<i<<endl;
    }
    infq.close();
    cout<<"calculating avg"<<endl;
    float qavg;
    calcAvg(m,n,imgq,0,0,&qavg,n);
    cout<<"average: "<<qavg<<endl;
    cout<<"roatation started"<<endl;
    float* oimgq=new float[m*n*2];
    rotateimg(m,n,0.0,oimgq);
    float* pimgq=new float[m*n*2];
    rotateimg(m,n,45.0,oimgq);
    float* nimgq=new float[m*n*2];
    rotateimg(m,n,-45.0,oimgq);
    cout<<"rotation done"<<endl;
    int * imgc;
    int * imgcq;
    float* imgcqo;
    float* imgcqp;
    float* imgcqn;
    float *soln;
    //float* solnh=new float(M*N);
    cudaMalloc (&imgc,(size_t) M*N*3*sizeof(int));
    cudaMallocManaged(&soln,(size_t) 3*M*N*sizeof(float));
    cudaMalloc (&imgcq,(size_t) m*n*3*sizeof(int));
    cudaMalloc (&imgcqo,(size_t) m*n*3*sizeof(float));
    cudaMalloc (&imgcqp,(size_t) m*n*3*sizeof(float));
    cudaMalloc (&imgcqn,(size_t) m*n*3*sizeof(float));
    cout<<"memory allocated"<<endl;
    cudaMemcpy(imgc,img,(size_t)sizeof(int)*M*N*3,cudaMemcpyHostToDevice);
    cudaMemcpy(imgcq,imgq,(size_t)sizeof(int)*m*n*3,cudaMemcpyHostToDevice);
    cudaMemcpy(imgcqo,oimgq,sizeof(float)*m*n*3,cudaMemcpyHostToDevice);
    cudaMemcpy(imgcqp,pimgq,sizeof(float)*m*n*3,cudaMemcpyHostToDevice);
    cudaMemcpy(imgcqo,nimgq,sizeof(float)*m*n*3,cudaMemcpyHostToDevice);
    cout<<"memory copied"<<endl;
    cout<<"average: "<<qavg<<endl;
    query<<<((M*N)/1024 + 1), 1024>>>(M,N,m,n,qavg,imgcqo,imgcqp,imgcqn,imgc,imgcq,soln,th1,th2);
    cout<<"processing done"<<endl;
    //cudaMemcpy(solnh,soln,(size_t)sizeof(int)*M*N,cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    float min = 100000.0;
    int idxi=-1;
    for(int i=0;i<M;i++){
        for(int j=0;j<N;j++){
            int idx=i*N+j;
            if (soln[3*(idx)]<=min && soln[3*(idx)+2]>0){
                min=soln[3*idx];
                idxi=idx;
                //cout<<i<<","<<j<<","<<soln[3*(i*N+j)]<<","<<soln[3*(i*N+j)+1]<<","<<soln[3*(i*N+j)+2]<<endl;
            }
        }
    }
    ofstream outf;
    outf.open(outfile.c_str(),ios::out);
    outf<<idxi/N<<" ";
    outf<<idxi%N<<" ";
    outf<<soln[idxi*3+2]-1<<" ";
    outf.close();
    cudaFree(imgc);
    cudaFree(soln);
    cudaFree(imgcq);
    cudaFree(imgcqo);
    cudaFree(imgcqp);
    cudaFree(imgcqn);

}