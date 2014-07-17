#include <cstdlib>
#include "dbtproj.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <time.h>

using namespace std;

int main(int argc, char** argv) {
    
    int nmem_blocks = atoi(argv[2]);
    clock_t start, stop;
    double t = 0.0;
    
    unsigned int nsorted_segs=0, npasses=0, nios=0,nres=0,nunique=0;
    block_t *buffer = new block_t[nmem_blocks];
    int field = atoi(argv[3]);
    
    //call Mergesort
    start = clock();
    cout<<"Running Mergesort...\n";
    MergeSort(argv[1], field, buffer, nmem_blocks,"MergesortOutput.bin", &nsorted_segs, &npasses, &nios);
    stop = clock();
    t = (double) (stop-start)/CLOCKS_PER_SEC;
    cout<<"Mergesort run:";
    cout<<"\nNumber of sorted segments:   "<<nsorted_segs;
    cout<<"\nNumber of passes:            "<<npasses;
    cout<<"\nNumber of I/Os:              "<<nios;
    cout<<"\nRuntime:                     "<<t<<"\n\n";
    
    //call EliminateDuplicates
    nios = 0;
    start = clock();
    cout<<"Running Eliminate Duplicates...\n";
    EliminateDuplicates(argv[1], field, buffer, nmem_blocks, "EliminateOutput.bin", &nunique, &nios);
    stop = clock();
    t = (double) (stop-start)/CLOCKS_PER_SEC;
    cout<<"Eliminate Duplicates run:";
    cout<<"\nNumber of unique records:     "<<nunique;
    cout<<"\nNumber of I/Os:               "<<nios;
    cout<<"\nRuntime:                      "<<t<<"\n\n";
    
    //call Mergejoin if there is a second file 
    if(argv[4] != NULL)
    {
        nios = 0;
        start = clock();
        cout<<"Running Mergejoin...\n";
        MergeJoin(argv[1], argv[4], field, buffer, nmem_blocks, "MergeJoinOutput.bin", &nres, &nios);
        stop = clock();
        t = (double) (stop-start)/CLOCKS_PER_SEC;
        cout<<"MergeJoin run:";
        cout<<"\nNumber of pairs:             "<<nres;
        cout<<"\nNumber of I/Os:              "<<nios;
        cout<<"\nRuntime:                     "<<t<<"\n\n";
    }
    
    delete[] buffer;
    return 0;
}

