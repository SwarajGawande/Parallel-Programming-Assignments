#include "psort.h"
#include <omp.h>

void swap(uint32_t* a,uint32_t* b){
	uint32_t temp=*b;
	*b=*a;
	*a=temp;
}

int partition(uint32_t* arr,int high,int low){
int mid = (high+low)/2;
swap(&arr[mid], &arr[high]);
uint32_t pivot = arr[high];
    int i = (low - 1);
 
    for (int j = low; j < high ; j++)
    {

        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void Quicksort(uint32_t* arr,int n,int l){
 if (l<n){
 int part=partition(arr,n,l);
 Quicksort(arr,part-1,l);
 Quicksort(arr,n,part+1);
 }
}

void seqSort(uint32_t* arr,uint32_t n){
	Quicksort(arr, n-1,0);
}


int range(uint32_t data,uint32_t* S,int p,int k){
	if (k>0){
	if (data<=S[k-1]) return false;
	}
	if (k<p-1){
	if (data>S[k]) return false;
	}
	return true;
}


void classify(uint32_t* data,uint32_t** buckets,uint32_t* counter,uint32_t ndata,int p,uint32_t* S)
{
   int numt= omp_get_num_threads();
   for(int k=0;k<p;k++){
   #pragma omp task shared(k,counter,buckets,data,S)
   {  //std::cout<<k<<"\n";
      for(uint32_t i=0; i<ndata; i++) { // Threads together share-loop through all of Data
         bool v = range(data[i],S,p,k);// For each data, find the interval of data's key,
         //std::cout<<data[i]<<","<<v<<"\n";
         if (v){
         buckets[k][counter[k]]=data[i]; // Found one key in interval v
         counter[k]++;
         }
      }
      if(counter[k] >= 2*ndata/p ) ParallelSort(buckets[k],counter[k],p);
			  else seqSort(buckets[k],counter[k]);
      //printArray(counter,p);
   }
   }
   #pragma omp taskwait
}

void ParallelSort(uint32_t *data, uint32_t n, int p)
{
    // Entry point to your sorting implementation.
    // Sorted array should be present at location pointed to by data.
    if (n<=p*p){seqSort(data,n); return;}
    
    int size=n/p;
    uint32_t* R=new uint32_t[p*p];
    uint32_t* S=new uint32_t[p-1];
    for (int i=0;i<p;i++){
    	for(int j=0;j<p;j++){
    		R[i*p+j]=data[size*i+j];
    	}
    }
    
    seqSort(R,p*p);
    
    for(int j=0;j<p-1;j++){
		S[j]=R[p*(j+1)];
	}
	uint32_t** buckets=new uint32_t*[p];
	uint32_t* counter=new uint32_t[p];
	for(int i=0;i<p;i++){
		buckets[i]=new uint32_t[n];
		counter[i]=0;
	}
	classify(data,buckets,counter,n,p,S);
	//assert_classify(buckets,counter,S,n,p);	
    uint32_t total = 0;
    for(uint32_t i=0;i<p;i++){
        for(uint32_t j=0;j<counter[i];j++){ 
        data[total++] = buckets[i][j];
        }
    }
    //std::cout<<"array:";
    //printArray(counter,p);
    //printArray(S,p-1);
}
