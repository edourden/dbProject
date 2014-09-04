#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <fstream>
#include "dbtproj.h"

using namespace std;

//variable for the function compare
unsigned char FIELD;

/**
 * Function used for comparisons, return 0 if the field of a and b are equal, -1 if the field of b > a
 * and 1 if the field of a > b.
 */
int compare (const void * a,const void * b)
{
  if (FIELD == 1) //compare for num
  {
	  if( ((record_t*)a)->num > ((record_t*)b)->num ) return 1;
	  else if ( ((record_t*)a)->num < ((record_t*)b)->num ) return -1;
	  else return 0;

  }
  else if (FIELD == 0) //compare for recid
    {
  	  if( ((record_t*)a)->recid > ((record_t*)b)->recid ) return 1;
  	  else if ( ((record_t*)a)->recid < ((record_t*)b)->recid ) return -1;
  	  else return 0;

    }
  else if (FIELD == 2) //compare for str
      {
	  	  return strcmp(((record_t*)a)->str ,((record_t*)b)->str);
      }
  else //compare for both num and str
  {
	  if( ((record_t*)a)->num > ((record_t*)b)->num ) return 1;
	  	  	  else if ( ((record_t*)a)->num < ((record_t*)b)->num ) return -1;
	  	  	  else
	  	  	  {
	  	  		return strcmp(((record_t*)a)->str ,((record_t*)b)->str);
	  	  	  }
  }
}
 
/**
 * swaps the records buffer[blockid].entries[i] and buffer[store_block].entries[storeid]
 * @param buffer
 * @param i
 * @param storeid
 * @param blockid
 * @param store_block
 */
 void swap(block_t *buffer, int i, int storeid, int blockid, int store_block)
 {
     record_t temp;
     temp = buffer[blockid].entries[i];
     buffer[blockid].entries[i] = buffer[store_block].entries[storeid];
     buffer[store_block].entries[storeid] = temp;
 }

 void print (block_t *b, int blocks)
 {
     for (int i=0;i< blocks; i++)
     {
         for (int j=0;j<b[i].nreserved;j++)
         {
             printf("this is block id: %d, record id: %d, num: %d, str: %s\n",
	b[i].blockid, b[i].entries[j].recid, b[i].entries[j].num, b[i].entries[j].str);
         }
     }
 }
 
 /**
  * quicksort for buffer, finds the position of pivot in the buffer from start_block start_index,
  * to end_block end_index, end_index record of end_block is set as pivot
  * @param buffer
  * @param start_block
  * @param start_index
  * @param end_block
  * @param end_index
  */
 void QuicksortBuffer(block_t *buffer,int start_block, int start_index, int end_block, int end_index){
     if (start_block >= end_block && start_index >= end_index)
     {
         return;
     }
     else 
     {
         
         int blockid; //id of the block we start the sorting
         int store = start_index; //keep tracks of the position the pivot belongs to in every loop
         int store_block = start_block; //block the store pointer belongs to
         int pivot = end_index; //choose last element of the last block as pivot
         int i = start_index;
         
         //finds the right position of the element pivot in the blocks of the buffer
         //any element before the the pivot is lesser or equal than it and any element after pivot is greater than it
         for (blockid = start_block; blockid <= end_block; blockid++) {
            //for every block in the from start_block to end_block in buffer
            while ( i < buffer[blockid].nreserved ) {
                if(blockid == end_block && i == end_index)
                {//if we reached the end of blocks and indexes
                    break;
                }
                else
                {
                //for every element in the current block
                if (compare(&(buffer[blockid].entries[i]), &(buffer[end_block].entries[end_index])) <= 0) 
                {
                    //if it is lesser of or equal to pivot
                    //swap the element at i with the element at store
                    swap(buffer, i, store, blockid, store_block);
                    
                    if(store+1 == buffer[store_block].nreserved && store_block+1 <= end_block)
                    {//if store has reached the end of the block and there is another block after that
                     //make store point to the start of the next block
                        store = 0;
                        store_block++;
                    }
                    else if (store+1 < buffer[store_block].nreserved)
                    {
                     //if store has not reached the end of the block yet,increment it
                        store++;
                    }
                }
                }
                i++;
            }
            i = 0;
        }
         //store_block shows the right block and store the right index in the block for pivot
         swap(buffer, end_index, store, end_block, store_block);
         
         
         //call QuicksortBuffer for the indexes before and after pivot
         if(store == 0 && store_block > start_block)
         {//if store is the first index of the block and there are more blocks before this one 
             QuicksortBuffer(buffer, start_block, start_index, store_block-1, buffer[store_block-1].nreserved-1);
         }
         else
         {//call QB for the blocks from start to pivot
             QuicksortBuffer(buffer, start_block, start_index, store_block, store-1);
         }
         if(store == buffer[store_block].nreserved-1 && store_block < end_block)
         {//if pivot is the last of the block and there are more blocks after it
             QuicksortBuffer(buffer, store_block+1, 0 , end_block, end_index);
         }
         else
         {//call QB for the blocks from pivot to end
             QuicksortBuffer(buffer, store_block, store+1, end_block, end_index);
         }
     }
     
 } 
 
/**
 * Merges the sorted files from start_file till nmem_blocks-1, to one sorted file
 * @param buffer
 * @param start_file
 * @param nmem_blocks
 * @param outfile
 */
void MergeFiles(block_t *buffer,int start_file, int nmem_blocks, char *outfile, unsigned int *nios)
 {   
     char filename[100];
     int id = 0;
                
     //open file for output
     FILE *output = fopen(outfile, "w");
     
     //variables to keep track of the next record in every block
     int block_pointers[nmem_blocks-1];
     
     //set the variables to the first record of every block
     for(int i=0; i<nmem_blocks-1; i++)
     {
         block_pointers[i] = 0;
     }
     
     //array of the files 
     FILE *files[nmem_blocks-1];
     
     //assign each variable of the array to a file
     for(int i=0; i<nmem_blocks-1; i++)
     {
         //filenames are file1,file2,file3 etc...
         sprintf(filename, "file%d", i+start_file);
         files[i] = fopen(filename, "r");
         //read the first block of the file to buffer[i]
         fread(&buffer[i], 1, sizeof(block_t), files[i]); // read the next block
         (*nios)++;
     }
     
     //initialize variables for merge
     //min_i and min_block point to the minimum record 
     int min_i=0,min_block=0;
     //last_block_i points to the next available record in the last block of buffer
     int last_block_i = 0;
     //finished keep how many files have finished 
     int finished = 0;
     
     //merge files 
     while(finished < nmem_blocks-1)
     {//while there are unfinished files
         min_block = 0;
         min_i = block_pointers[0];
         
         for (int i = 0; i < nmem_blocks - 1; i++) 
         {//check the first element of every block to find the minimum record
            if (block_pointers[i] < buffer[i].nreserved) 
            {//if the block[i] has not reached the end of the file[i]
                if (compare(&(buffer[min_block].entries[min_i]), &(buffer[i].entries[block_pointers[i]])) > 0) 
                {//if the current min element is greater than the one compared to
                    min_block = i;
                    min_i = block_pointers[i];
                }
            }
        }
         
        //write min record to the last block of the buffer
        buffer[nmem_blocks - 1].entries[last_block_i] = buffer[min_block].entries[min_i];
        last_block_i++;
        block_pointers[min_block]++;
        
        //check if the block_pointers of the min_i block reached the end of the block
        if (block_pointers[min_block] == buffer[min_block].nreserved) 
        {//check if there is another block in the file of this block
            if (fread(&buffer[min_block], 1, sizeof (block_t), files[min_block]))
            {
                block_pointers[min_block] = 0;
            } 
            else
            {
                finished++;
                if(finished == nmem_blocks - 2)
                {//if there is only one not finished file
                        
                    //write last block to output file
                    buffer[nmem_blocks - 1].nreserved = last_block_i;
                    buffer[nmem_blocks - 1].blockid = id;
                    id++;
                    fwrite(&buffer[nmem_blocks - 1], 1, sizeof (block_t), output);
                    last_block_i = 0;
                    (*nios)++;
                        
                    //find the not finished block in the buffer
                    for(int j=0; j<nmem_blocks - 1; j++)
                    {
                        if (block_pointers[j] < buffer[j].nreserved) 
                        {//buffer[j] hasnt finished yet
                            //add every record of buffer[j] to last block and when they are finished write it to 
                            while(block_pointers[j] < buffer[j].nreserved)
                            {
                                buffer[nmem_blocks-1].entries[last_block_i] = buffer[j].entries[block_pointers[j]];
                                last_block_i++;
                                block_pointers[j]++;
                            }
                            //write last block
                            buffer[nmem_blocks - 1].nreserved = last_block_i;
                            buffer[nmem_blocks - 1].blockid = id;
                            id++;
                            fwrite(&buffer[nmem_blocks - 1], 1, sizeof (block_t), output);
                            last_block_i = 0;
                            (*nios)++;
                            
                            //write every block from file to output file
                            while (fread(&buffer[0], 1, sizeof (block_t), files[j]))
                            {
                                fwrite(&buffer[0], 1, sizeof(block_t), output);
                                (*nios) += 2;
                            } 
                        }
                    }
                    break;
                }
            }
        }
        //check if last block is full
        if (last_block_i == MAX_RECORDS_PER_BLOCK) 
        {
            //write last block to output file
            buffer[nmem_blocks - 1].nreserved = last_block_i;
            buffer[nmem_blocks - 1].blockid = id;
            id++;
            fwrite(&buffer[nmem_blocks - 1], 1, sizeof (block_t), output);
            last_block_i = 0;
            (*nios)++;
        }
     } 
     
     //close all files and delete the input files used in the merge
     fclose(output);
     for(int i=0; i<nmem_blocks-1; i++)
     {
         fclose(files[i]);
         sprintf(filename, "file%d", i+start_file);
         remove(filename);
     }
 }
 
 /**
  * returns minimum of a and b
  */
 int min(int a, int b)
 {
     if (a <= b)
         return a;
     return b;
 }


/* ----------------------------------------------------------------------------------------------------------------------
   infile: the name of the input file
   field: which field will be used for sorting: 0 is for recid, 1 is for num, 2 is for str and 3 is for both num and str
   buffer: pointer to memory buffer
   nmem_blocks: number of blocks in memory
   outfile: the name of the output file
   nsorted_segs: number of sorted segments produced (this should be set by you)
   npasses: number of passes required for sorting (this should be set by you)
   nios: number of IOs performed (this should be set by you)
   ----------------------------------------------------------------------------------------------------------------------
*/
void MergeSort (char *infile, unsigned char field, block_t *buffer, unsigned int nmem_blocks, char *outfile, unsigned int *nsorted_segs, unsigned int *npasses, unsigned int *nios)
{
    //open the infile for reading
    FILE *input,*output,*temp;
    char filename[100];
    
    //number of blocks contained in the input file
    int blocks = 0; //used to find the npasses

    //sets the global variable FIELD to be used by the function compare
    FIELD = field;
    
    //open file for input
    input = fopen(infile,"r");

    if(input != NULL) //check the file opened and exists
    {
        //number of files that we will create for every merge
        int num_of_files = 0;
        
    	while (!feof(input)) 
        {
            //fill buffer from input file
            int i = 0;
            while (i < nmem_blocks && fread(&buffer[i], 1, sizeof (block_t), input)) 
            {
                //increment i and nios
                i++;
                (*nios)++;
            }
            if (i > 0) 
            {
                //call QuicksortBuffer to sort the buffer
                QuicksortBuffer(buffer, 0, 0, i - 1, buffer[i - 1].nreserved - 1);
                
                //increment nsorted_segments
                (*nsorted_segs)++;
                
                //add i to blocks
                blocks += i;
                
                //buffer is a sorted segment,write it to a temp file
                num_of_files++;
                sprintf(filename, "file%d", num_of_files);
                temp = fopen(filename, "w");

                //write sorted blocks
                for (int j = 0; j < i; j++) {
                    (*nios)++;
                    fwrite(&buffer[j], 1, sizeof (block_t), temp); // write the block to the file
                }
                fclose(temp);
            }
            
        } 
    	fclose(input);
        
        //calculate npasses
        blocks = blocks / nmem_blocks;
        (*npasses)++;
        while(blocks > 0)
        {
            blocks = blocks / (nmem_blocks-1);
            (*npasses)++;
        }

        if (num_of_files == 1) 
        {//if only one file was created,no need for merge
            //rename it to the filename specified in outfile
            rename(filename, outfile);
        } else {

            //we have sorted the buffer and store it to num_of_files different files
            //we have to merge them until we have only one file
            int start = 1;
            int last_file = num_of_files;
            int end;
            
            while (num_of_files > 1) {
                //set end to minimum of number of files+1 or size of the buffer
                end = min(nmem_blocks, num_of_files + 1);
                //if all files can be merged set the filename of the output to outfile
                if (end == num_of_files + 1)
                    strcpy(filename, outfile);
                else
                    sprintf(filename, "file%d", last_file + 1);

                //merge files
                MergeFiles(buffer, start, end, filename, nios);
                //change variables according to how many files were merged and created
                last_file++;
                num_of_files -= end - 2;
                start += end - 1;
            }
        }
    }
    else
    {//opening the file failed
    	cout<<"Could not open the input file.\n";
    } 
}

/////////////////////////////////////////
void EliminateDuplicates (char *infile, unsigned char field, block_t *buffer, unsigned int nmem_blocks, char *outfile, unsigned int *nunique, unsigned int *nios)
{
    //call MergeSort
    unsigned int nsorted_segs=0, npasses=0;
    MergeSort(infile, field, buffer, nmem_blocks, "temp.bin", &nsorted_segs, &npasses, nios);
    
    //record keep track of the last record
    record_t previous_record;
    //set previous_record valid to 0 for the comparison with the first record
    previous_record.valid = 0;
    
    //int to keep track of the last record in last block of buffer
    int last_block_i = 0;
    //int to keep track of the last block's id everytime we write it to output file
    int id = 0;
    
    //open sorted file
    FILE *input;
    input = fopen("temp.bin", "r");
    FILE *output;
    output = fopen(outfile, "w");

    while (!feof(input)) 
    {
        //fill buffer
        int i = 0;
        while (i < nmem_blocks - 1 && fread(&buffer[i], 1, sizeof(block_t), input)) 
        {
            (*nios)++;
            i++;
        }
        
        //for every block in buffer, parse each record keeping track of previous and write to last block of buffer 
        //if previous record is different than current record
        for (int j = 0; j < i; j++) 
        {
            for (int k = 0; k < buffer[j].nreserved; k++) 
            {
                //compare returns 0 if the records have the same value in the field we specified
                if (compare(&(buffer[j].entries[k]), &(previous_record))) 
                {
                    (*nunique)++;
                    //records are different in that field
                    //write that record in the last block of buffer
                    buffer[nmem_blocks - 1].entries[last_block_i] = buffer[j].entries[k];
                    last_block_i++;

                    //check if last block is full
                    if (last_block_i == MAX_RECORDS_PER_BLOCK) 
                    {
                        //write last block to output file
                        buffer[nmem_blocks - 1].blockid = id;
                        buffer[nmem_blocks - 1].nreserved = MAX_RECORDS_PER_BLOCK;
                        fwrite(&buffer[nmem_blocks - 1], 1, sizeof (block_t), output);
                        (*nios)++;
                        id++;
                        last_block_i = 0;
                    }
                }
                //change previous record to current record
                previous_record = buffer[j].entries[k];
            }
        }
    }
    
    //if there are records in the last block not written to output file
    if(last_block_i > 0)
    {
       //write last block to output file
       buffer[nmem_blocks - 1].blockid = id;
       buffer[nmem_blocks - 1].nreserved = last_block_i;
       fwrite(&buffer[nmem_blocks - 1], 1, sizeof (block_t), output); 
       (*nios)++;
    }
    
    //close files
    fclose(input);
    fclose(output);
    remove("temp.bin");
}


///////////////////////////////////////////////////

void MergeJoin (char *infile1, char *infile2, unsigned char field, block_t *buffer, unsigned int nmem_blocks, char *outfile, unsigned int *nres, unsigned int *nios)
{
    unsigned int nsorted_segs = 0;
    unsigned int npasses = 0;

    //sort the inputs
    MergeSort(infile1, field, buffer, nmem_blocks, "temp1.bin", &nsorted_segs, &npasses, nios);
    MergeSort(infile2, field, buffer, nmem_blocks, "temp2.bin", &nsorted_segs, &npasses, nios);
    
    //open sorted file
    FILE *input1, *input2;
    input1 = fopen("temp1.bin", "r");
    input2 = fopen("temp2.bin", "r");
    //open output file
    FILE *output;
    output = fopen(outfile, "w");
    
    //fill the first positions of the buffer from input1,until there are only 2 
    //empty blocks in buffer or we reach the end of input1
    int blocks1 = 0;
    while(blocks1 < nmem_blocks-2 && fread(&buffer[blocks1], 1, sizeof (block_t), input1))
    {
        (*nios)++; 
        blocks1++;
    }
    //input1 is from 0 to blocks1,block in buffer
    
    //fill the rest of the buffer,excluding last block,from input2
    int blocks2 = blocks1;
    while(blocks2 < nmem_blocks-1 && fread(&buffer[blocks2], 1, sizeof (block_t), input2))
    {
      (*nios)++;
      blocks2++;
    }
    //input1 is from 0 to blocks1
    //input2 is from blocks1 to blocks2 
    //and the output block is from blocks2 to nmem_blocks
    
    int ignore_till_segment = 0;
    int last_block_i = 0; //keeps track of the last entry in the last block
    int id = 0; //keeps track of the id for the output blocks
    
    for(int i=0; i<blocks1; i++)
    {//i goes through the blocks from input1
        for(int r1=0; r1<buffer[i].nreserved; r1++)
        {//r1 goes through every record of the i block
            for(int j=blocks1; j<blocks2; j++)
            {//j goes through the blocks from input2
                for(int r2=0; r2<buffer[j].nreserved; r2++)
                {//r2 goes through every record of the j block
                    if (!compare(&(buffer[i].entries[r1]), &(buffer[j].entries[r2]))) 
                    {//FIELD of these record is the same
                        //write that record in the last block of buffer
                        buffer[nmem_blocks - 1].entries[last_block_i] = buffer[i].entries[r1];
                        last_block_i++;
                        //check if last block is full
                        if (last_block_i == MAX_RECORDS_PER_BLOCK) 
                        {
                            //write last block to output file
                            buffer[nmem_blocks - 1].blockid = id;
                            buffer[nmem_blocks - 1].nreserved = MAX_RECORDS_PER_BLOCK;
                            fwrite(&buffer[nmem_blocks - 1], 1, sizeof (block_t), output);
                            (*nios)++;
                            id++;
                            last_block_i = 0;
                        }
                    }
                }
            }
        }
    }
    
    
    //close the files,and delete any temp files that were created
    fclose(input1);
    fclose(input2);
    fclose(output);
    remove("temp1.bin");
    remove("temp2.bin");
}

///////////////////////////////////////////////////

int readFromInput(ifstream &input_file, block_t *buffer, unsigned int n_input_mem_blocks, unsigned int firstPosition, unsigned int *nios, unsigned int &recordsRead) {
    int i = firstPosition;//the first position of the input buffer
    int sum = 0;
    int returnValue=0;//if the input doesn't finish, the value won't change
    while (!input_file.eof() && i < firstPosition+n_input_mem_blocks) {

        input_file.read(reinterpret_cast<char*> (&buffer[i]), sizeof (block_t));
        
        sum += buffer[i].nreserved;
        
        

        if (input_file.eof()) 
        {
            if(i==firstPosition)
            {
                returnValue=1;//if the input finished earlier but was not "caught" by eof())
            }
            else
            {
                returnValue=2;//if the input finished now
            }
        }
        i++;
    }
    recordsRead = sum;
    (*nios)++;
    
    return returnValue;
}

void writeToOutput(ofstream &output_file, block_t *buffer, unsigned int firstOutputIndex, unsigned int *nios, unsigned int n_output_blocks) {
    //for every block in the output buffer
    for (int i = firstOutputIndex; i < firstOutputIndex + n_output_blocks; i++) {

        output_file.write(reinterpret_cast<const char*> (&buffer[i]), sizeof (block_t));
    }
    (*nios)++;
}

int compareRecords(record_t *record, record_t *record2, unsigned char field) {
    int returnValue;
    if (field == 0) {
        if (record->recid == record2->recid) returnValue = 0;
        else if (record->recid > record2->recid) returnValue = 1;
        else returnValue = -1;
    } else if (field == 1) {
        if (record->num == record2->num) returnValue = 0;
        else if (record->num > record2->num) returnValue = 1;
        else returnValue = -1;
    } else if (field == 2) {
        if (strcmp(record->str, record2->str) == 0) returnValue = 0;
        else if (strcmp(record->str, record2->str) > 0) returnValue = 1;
        else returnValue = -1;
    } else if (field == 3) {
        if (record->num == record2->num) {
            int strCR = strcmp(record->str, record2->str);
            if (strCR == 0) returnValue = 0;
            else if (strCR > 0) returnValue = 1;
            else returnValue = -1;
        } else if (record->num > record2->num) returnValue = 1;
        else returnValue = -1;
    }

    return returnValue;
}


void MergeJoinold(char *infile1, char *infile2, unsigned char field, block_t *buffer, unsigned int nmem_blocks, char *outfile, unsigned int *nres, unsigned int *nios)
{
    
    
    unsigned int nsorted_segs = 0;
    unsigned int npasses = 0;

    //sort the input
    MergeSort(infile1, field, buffer, nmem_blocks, "temp1.bin", &nsorted_segs, &npasses, nios);
    MergeSort(infile2, field, buffer, nmem_blocks, "temp2.bin", &nsorted_segs, &npasses, nios);
    
    
    //initialize the buffer values
    for(int i=0;i<nmem_blocks;i++)
    {
        for(int j=0;j<MAX_RECORDS_PER_BLOCK;j++)
        {
            buffer[i].entries[j].num=0;
            buffer[i].entries[j].recid=0;
            buffer[i].entries[j].valid=false;
            strcpy(buffer[i].entries[j].str,"\0");
        }
        buffer[i].nreserved=0;
    }

    //open the file streams
    ifstream input_file1;
    input_file1.open("temp1.bin", ios::binary);
    ifstream input_file2;
    input_file2.open("temp2.bin", ios::binary);
    ofstream output_file;
    output_file.open(outfile, ios::binary);


    unsigned int outputBufferSize = 1; //output buffer size in blocks
    unsigned int firstBufferSize = (nmem_blocks - outputBufferSize)/2; //the buffer for the first input in blocks
    unsigned int secondBufferSize = nmem_blocks - outputBufferSize - firstBufferSize; //the buffer for the second input in blocks

    //initialize the counters
    unsigned int firstBufferBlockCounter = 0;
    unsigned int firstBufferRecordCounter = 0;
    unsigned int secondBufferBlockCounter=firstBufferSize;
    unsigned int secondBufferRecordCounter=0;
    unsigned int outputBufferBlockCounter = firstBufferSize+secondBufferSize;
    unsigned int outputBufferRecordCounter = 0;

    //initialize the pointers
    record_t *firstBufferPointer = &(buffer[firstBufferBlockCounter].entries[firstBufferRecordCounter]);
    record_t *secondBufferPointer = &(buffer[secondBufferBlockCounter].entries[secondBufferRecordCounter]);
    record_t *outputBufferPointer = &(buffer[outputBufferBlockCounter].entries[outputBufferRecordCounter]);
    //record_t *compareRecordPointer = NULL;

    int input1Finished = 0; //shows the conditions for input1
    int input2Finished = 0; //shows the conditions for input2
    bool oneInputFinished=false;

    bool needToRead1 = true; //is true when there is the need to read more records from the first input
    bool needToRead2 = true; //is true when there is the need to read more records from the second input
    bool needToWrite = false;//is true when the output buffer is full and the records must get written to the output file

    unsigned int recordsRead1 = 0; //counts how many records were read from the first input
    unsigned int recordsRead2 = 0; //counts how many records were read from the second input
    unsigned int blocksUsed1 = 0; //counts how many blocks are being used after reading from the first input
    unsigned int blocksUsed2 = 0; //counts how many blocks are being used after reading from the second input
    unsigned int reservedInOutput=0;//counts how many records have been written in the output buffer

    
    for(int i=firstBufferSize+secondBufferSize;i<nmem_blocks;i++)
    {
        buffer[i].valid=false;
    }
    for(int i=0;i<firstBufferSize+secondBufferSize;i++)
    {
        buffer[i].valid=true;
    }
    
    unsigned int counter=0;//defines the block id of each block written to output
    
    while (!oneInputFinished) {//while none of the inputs has finished
        
        if(needToRead1)//read input1 if needed
        {
            recordsRead1 = 0;
            input1Finished = readFromInput(input_file1, buffer, firstBufferSize, 0, nios, recordsRead1);
            blocksUsed1 = (recordsRead1 / MAX_RECORDS_PER_BLOCK);
            if (recordsRead1 % MAX_RECORDS_PER_BLOCK != 0) blocksUsed1++;
            needToRead1 = false;
        }
        if(needToRead2)//read input2 if needed
        {
            recordsRead2 = 0;
            input2Finished = readFromInput(input_file2, buffer, secondBufferSize, firstBufferSize, nios, recordsRead2);
            blocksUsed2 = (recordsRead2 / MAX_RECORDS_PER_BLOCK);
            if (recordsRead2 % MAX_RECORDS_PER_BLOCK != 0) blocksUsed2++;
            needToRead2 = false;
        }
        
        if(input1Finished==1 || input2Finished==1) break;
        
        if(input1Finished==2 || input2Finished==2) oneInputFinished=true;

        bool loop = true; //becomes false when all the records of one memory buffer were parsed

        while (loop)
        {
            int check = compareRecords(firstBufferPointer,secondBufferPointer,field);
            if (check == 0) {//if the records are the same
                //copy the fields from the input part of the buffer to the output part of the buffer
                outputBufferPointer->num = firstBufferPointer->num;
                outputBufferPointer->recid = firstBufferPointer->recid;
                strcpy(outputBufferPointer->str, firstBufferPointer->str);
                outputBufferPointer->valid = firstBufferPointer->valid;
                needToWrite = true;
                reservedInOutput++;
                buffer[outputBufferBlockCounter].nreserved++;

                (*nres)++;
                outputBufferRecordCounter++;
                if (outputBufferRecordCounter >= MAX_RECORDS_PER_BLOCK) {
                    outputBufferRecordCounter = 0;
                    outputBufferBlockCounter++;

                    if (outputBufferBlockCounter >= outputBufferSize) {
                        outputBufferBlockCounter = firstBufferSize+secondBufferSize;
                        unsigned int position=firstBufferSize+secondBufferSize;
                        while(position<nmem_blocks)
                        {
                            buffer[position].blockid=counter;
                            counter++;
                            position++;
                        }

                        //write to the output file
                        writeToOutput(output_file, buffer, firstBufferSize+secondBufferSize, nios, outputBufferSize);
                        reservedInOutput==0;
                        
                        position=firstBufferSize+secondBufferSize;
                        
                        while(position<nmem_blocks)
                        {
                            buffer[position].nreserved = 0;
                            for(int j=0;j<MAX_RECORDS_PER_BLOCK;j++)
                            {
                                buffer[position].entries[j].valid=false;
                            }
                            position++;
                        }
                        
                        
                        needToWrite = false;
                    }
                }
                
                outputBufferPointer = &(buffer[outputBufferBlockCounter].entries[outputBufferRecordCounter]);
                
                firstBufferRecordCounter++;
                secondBufferRecordCounter++;
            }
            else if(check>0)//else if buffer1 record is bigger
            {
                secondBufferRecordCounter++;
            }
            else//else
            {
                firstBufferRecordCounter++;
            }
            
            if(firstBufferRecordCounter >= buffer[firstBufferBlockCounter].nreserved)
            {
                firstBufferRecordCounter=0;
                firstBufferBlockCounter++;
                if(firstBufferBlockCounter >= blocksUsed1)
                {
                    firstBufferBlockCounter = 0;
                    loop=false;
                    needToRead1=true;
                }
            }
            
            if(secondBufferRecordCounter >= buffer[secondBufferBlockCounter].nreserved)
            {
                secondBufferRecordCounter=0;
                secondBufferBlockCounter++;
                if(secondBufferBlockCounter-firstBufferSize >= blocksUsed2)
                {
                    secondBufferBlockCounter = firstBufferSize;
                    loop=false;
                    needToRead2=true;
                }
            }

            firstBufferPointer = &(buffer[firstBufferBlockCounter].entries[firstBufferRecordCounter]);
            secondBufferPointer = &(buffer[secondBufferBlockCounter].entries[secondBufferRecordCounter]);
        }
    }

    if (needToWrite) {
        unsigned int blocksUsed = (reservedInOutput / MAX_RECORDS_PER_BLOCK);
        if (reservedInOutput % MAX_RECORDS_PER_BLOCK != 0) blocksUsed++;
        writeToOutput(output_file, buffer, firstBufferSize+secondBufferSize, nios, blocksUsed);
    }
    input_file1.close();
    input_file2.close();
    output_file.close();
    
}