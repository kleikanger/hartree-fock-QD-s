#ifndef NEWMATRIX_H
#define NEWMATRIX_H

//void **matrix(int row, int col, int num_bytes);
//void **tria_matrix(int n, int num_bytes);
//void free_matrix(void **matr);

#include <string.h>
#include <cstdlib>
#include <iostream>

using std::cout;
using std::endl;
using std::nothrow;

/*
   * The function                             
   * reserves dynamic memory for a two-dimensional matrix 
   * using the C++ command new . No initialization of the elements. 
   * Input data:                      
   *  int row      - number of  rows          
   *  int col      - number of columns        
   *                 element                  
   * Returns a void  **pointer to the reserved memory location.
   */

template <class T>
T **matrix(int row, int col)
{
	int num_bytes = sizeof(T);
	int i, num;
	T **pointer, *ptr;

	pointer = new(nothrow) T* [row];
	if(!pointer) 
	{
		cout << "Exception handling: Memory allocation failed";
		cout << " for "<< row << "row addresses !" << endl;
		return nullptr;
	}

	i = (row * col * num_bytes)/sizeof(char);
	pointer[0] = new(nothrow) T [i];
	if(!pointer[0]) 
	{
		cout << "Exception handling: Memory allocation failed";
		cout << " for address to " << i << " characters !" << endl;
		return nullptr;
	}
	ptr = pointer[0];
	num = col * num_bytes;
	for(i = 0; i < row; i++, ptr += num )   
	{
		pointer[i] = ptr; 
	}

return pointer;

} // end: function void **matrix()
    
/*
     * The function                         
     *      void free_matrix()              
     * releases the memory reserved by the function matrix() 
     * for the two-dimensional matrix[][] 
     * Input data:                          
     *  void far **matr - pointer to the matrix
     */

template <class T>
void free_matrix(T **matr)
{
  delete [] matr;
} // End:  function free_matrix() 

#endif
