//#include "nctest.h"
#include "uebpgdecls.h"
//The following program creates a matrix. The inputs are matrix dimensions. It allocates a memory block of size nrows*ncols* (size of float)
//and returns an array of pointers to the allocated memory block
float*** Create3DArray(int nt, int nr, int nc)       //inputs: no. of rows and no. of colos (dimensions) of matrix
{
	float*** myMatrix = new float**[nt]; 	
	for(int i = 0;i<nt;i++)
	{
		myMatrix[i] = new float*[nr];  
		for(int j=0; j< nr; j++)
			myMatrix[i][j] = new float [nc];		
	}
	return myMatrix;

}//float** CreateMatrix

//The following program deletes a matrix passed to it (it frees up the memory block allocated to the matrix)
void DeleteMatrix(float ***A, int nt, int nr, int nc)            //input: A matrix
{
	for(int i = 0; i< nt; i++)
	{
		for(int j = 0; j< nr; j++)
			delete [] A[i][j];
		delete [] A[i];
	}
	delete [] A;

	return;
}//void DeleteMatrix

