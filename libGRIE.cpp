#ifndef LIBGRIE_CPP
#define LIBGRIE_CPP

/*
   *
   * This class is tailored to store and return the coulomb interaction energies of 
   * the quantum dots
   *
   *
   *
   *
   *
   */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "OpenFCI/classQDotInteraction.hpp"
#include "libGRIE.h"

#define RUNMINIMAL true

using std::cerr;
using std::cout;
using std::endl;
using std::nothrow;

//allocate memory for a 2d array pointer
void **mtrix(int, int, int);
/*
   *
   * Constructor
   *
   */
libGRIE::libGRIE(int* pi_mARG, int* pi_nARG, const int irARG, const int inoARG,  
		const double d_lambda,  bool b_useveff = false, bool b_energycut = true)
{/*startvimfold*/
	pi_n = pi_nARG;
	pi_m = pi_mARG;
	ir = irARG;
	ino = inoARG;
	
	ir=irARG;
	oOFCI = new QdotInteraction();
	if ((b_energycut) && (!b_useveff))
		oOFCI->setR(ir);
	else
		//must be used if effective interactions are used!? 
		oOFCI->setR(ir*2); 
	oOFCI->setLambda(d_lambda);
	if (b_useveff)
		//Eff. interactions
		oOFCI->buildEffectiveInteractionComBlocks(1); 
	else
		//Bare interactions
		oOFCI->buildInteractionComBlocks(); 
};/*endvimfold*/
/*
   *
   * Destructor
   *
   */
libGRIE::~libGRIE()
{/*startvimfold*/
	delete oOFCI;
	if (ppd_gm1) 
		delete [] ppd_gm1;
	if (ppd_gm2)
		delete [] ppd_gm2;
	if (pi_iindex)
	{
		delete [] pi_iindex;
	}
	if (pppd_cie)
	{
		delete [] **pppd_cie;
		delete [] *pppd_cie;
		delete [] pppd_cie;
	}
};/*endvimfold*/
/*
   *
   *  Declare and initialize the GRIE matrices
   *  Method only works for ino>2
   *
   */
void libGRIE::initGRIEmatrices()
{/*startvimfold*/
	int i,j;
	ppd_gm1 = (double**) mtrix(ino, 3, sizeof(double));
	ppd_gm2 = (double**) mtrix(ino, 3, sizeof(double));
	
	//Set up the GRIE martix for case 1,3
	for (i=0; i<3; i++)
	{
		ppd_gm1[0][i] = 0.0;
		ppd_gm1[1][i] = 1.0;
	}
	for (i=1; i<ino; i++)
	{
		ppd_gm1[i][0] = ppd_gm1[i-1][0] + 1.0;
		for (j=1; j<3; j++)
			ppd_gm1[i][j] = ppd_gm1[i-1][j] + ppd_gm1[i][j-1];
	}
	//Set up the GRIE martix for case 2
	for (i=0; i<3; i++)
		ppd_gm2[0][1] = 0; //never used
	
	j = 0;
	for (i=ino-1; i>=2; i--)
	{
		ppd_gm2[i][0] = j;
		j++;
	}
	ppd_gm2[2][1] = j;
	ppd_gm2[ino-1][1] = 0; //never used
	ppd_gm2[2][2] = j; 
	j--;
	for (i=3; i<ino-1; i++)
	{
		ppd_gm2[i][1] = ppd_gm2[i-1][1]+j;
		j--;
	}
	for (i=3; i<ino-1; i++)
	{
		ppd_gm2[i][2] = ppd_gm2[i-1][2] + ppd_gm2[i][1];
	}
	ppd_gm2[ino-1][2] =  ppd_gm2[ino-2][2] + ppd_gm2[ino-2][1]+1;
	
	//calc. the total nr. of unique cie's
	ino1 = ppd_gm1[ino-1][0] + ppd_gm1[ino-1][1] + ppd_gm1[ino-1][2] + 2;
	ino2 = ppd_gm2[ino-2][1] + ppd_gm2[ino-1][2] + 2;
	ino3 = ppd_gm1[ino-2][0] + ppd_gm1[ino-2][1] + 2;
}/*endvimfold*/
/*
   * 
   * Initialize arrays with the interaction elements.
   * Initializing 3D array pppd_cie[i_spinproj][n_grie][pi_iindx].
   * All enpty pppd_cie[spin_proj][n_GRIE] set to nullptr to avoid sneaky bugs.
   *
   */
void libGRIE::initInteractionElements()
{/*startvimfold*/
	int i, j, k, l, im, inumcie, igrie, istat;
	bool bex;
	double d_dir, d_exc;
	double *pd_uuuu, *pd_udud, *pd_uddu;
	double **ppd_uuuu, **ppd_udud, **ppd_uddu;

	//count the number of interaction elem.
	//NB: uuuu for i=j=k=l always 0 (very few elem).	
#if !RUNMINIMAL
	cout << endl << endl << "Counting nr. of CIE's.... ";
#endif
	inumcie = 0;
	for (j=0; j<ino; j++)
	for (k=j; k<ino; k++)
	for (l=k; l<ino; l++)
	{
		im = - pi_m[j] + pi_m[k] + pi_m[l];
		for (i=0; i<=j; i++)
			if (pi_m[i]==im)
				inumcie++;
	}
	for (j=ino-1; j>=2; j--)
	for (k=1; k<j; k++)
	for (l=k; l<ino; l++)
	{
		im = - pi_m[j] + pi_m[k] + pi_m[l];
		for (i=0; i<k; i++)
			if (pi_m[i]==im)
				inumcie++;
	}
	for (j=1; j<ino; j++)
	for (l=j; l<ino; l++)
	{
		im = pi_m[j] - pi_m[l];
		for (i=0; i<j; i++)
			if (im==0)
				inumcie++;
	}

#if !RUNMINIMAL
	cout << "finished." << endl;
	cout << "Number of CIE's:    " << inumcie << endl;
	cout << "Number of pointers: " << ino1+ino2+ino3 
		<< " Case I: " << ino1 
		<< " Case II:  " << ino2  
		<< " Case III: " << ino3 <<  endl;
#endif

	//Res. mem. for all interactionelem.
	pd_uuuu = new(nothrow) double[inumcie];
	pd_udud = new(nothrow) double[inumcie];
	pd_uddu = new(nothrow) double[inumcie];
  	if ( (!pd_uuuu) || (!pd_uddu) || (!pd_udud) ) 
	{
    	cout << "Memory allocation failed (A)";
	}
	//Res. mem. for pointers.
	ppd_uuuu = new(nothrow) double*[ino1+ino2+ino3];
	ppd_udud = new(nothrow) double*[ino1+ino2+ino3];
	ppd_uddu = new(nothrow) double*[ino1+ino2+ino3];
  	if ( (!pd_uuuu) || (!pd_uddu) || (!pd_udud) ) 
	{
    	cout << "Memory allocation failed (B)";
	}
#if !RUNMINIMAL
	cout << "Initializing CIE's.... " << endl;
#endif
	//Calculate and store all cie's
	igrie = inumcie = 0;
	for (j=0; j<ino; j++)
	for (k=j; k<ino; k++)
	for (l=k; l<ino; l++)
	{
		igrie = ppd_gm1[j][0] + ppd_gm1[k][1] + ppd_gm1[l][2];
		im = - pi_m[j] + pi_m[k] + pi_m[l];
		istat = inumcie;
		bex = 0;

		for (i=0; i<=j; i++)
			if (pi_m[i]==im)
			{
				d_dir = calcCoulIE(i, j, k, l);
				d_exc = calcCoulIE(i, j, l, k);
				pd_uuuu[inumcie] = d_dir-d_exc;
				pd_udud[inumcie] = d_dir;
				pd_uddu[inumcie] = - d_exc;
				inumcie++;
				bex = 1;
			}
		if (bex)
		{
			ppd_uuuu[igrie] = &pd_uuuu[istat];
			ppd_udud[igrie] = &pd_udud[istat];
			ppd_uddu[igrie] = &pd_uddu[istat];
		}
		else
		{
			ppd_uuuu[igrie] = nullptr;
			ppd_udud[igrie] = nullptr;
			ppd_uddu[igrie] = nullptr;
		}
	}
	for (j=ino-1; j>=2; j--)
	for (k=1; k<j; k++)
	for (l=k; l<ino; l++)
	{
		igrie = ppd_gm2[j][0] + ppd_gm2[k][1] + ppd_gm2[l][2] + ino1;
		im = - pi_m[j] + pi_m[k] + pi_m[l];
		istat = inumcie;
		bex = 0;
		
		for (i=0; i<k; i++)
			if (pi_m[i]==im)
			{	
				d_dir = calcCoulIE(i, j, k, l);
				d_exc = calcCoulIE(i, j, l, k);
				pd_uuuu[inumcie] = d_dir-d_exc;
				pd_udud[inumcie] = d_dir;
				pd_uddu[inumcie] = -d_exc;
				inumcie++;
				bex = 1;
			}
		if (bex)
		{
			ppd_uuuu[igrie] = &pd_uuuu[istat];
			ppd_udud[igrie] = &pd_udud[istat];
			ppd_uddu[igrie] = &pd_uddu[istat];
		}
		else
		{
			ppd_uuuu[igrie] = nullptr;
			ppd_udud[igrie] = nullptr;
			ppd_uddu[igrie] = nullptr;
		}
		
	}
	for (j=1; j<ino; j++)
	for (l=j; l<ino; l++)
	{
		igrie = ppd_gm1[j-1][0] + ppd_gm1[l-1][1] + ino1 + ino2;
		istat = inumcie;
		bex = 0;
		if (pi_m[j]==pi_m[l])
		{
			for (i=0; i<j; i++)
			{
				d_dir = calcCoulIE(i, j, i, l);
				d_exc = calcCoulIE(i, j, l, i);
				pd_uuuu[inumcie] = d_dir-d_exc;
				pd_udud[inumcie] = d_dir;
				pd_uddu[inumcie] = - d_exc;
				inumcie++;
				bex = 1;
			}
		}
		if (bex)
		{
			ppd_uuuu[igrie] = &pd_uuuu[istat];
			ppd_udud[igrie] = &pd_udud[istat];
			ppd_uddu[igrie] = &pd_uddu[istat];
		}
		else
		{
			ppd_uuuu[igrie] = nullptr;
			ppd_udud[igrie] = nullptr;
			ppd_uddu[igrie] = nullptr;
		}
		
	}
	//Class variable pppd_cie contains all Coul. Interaction. Elem.
	pppd_cie = new double**[3];
	pppd_cie[0] = &ppd_uuuu[0];
	pppd_cie[1] = &ppd_udud[0];
	pppd_cie[2] = &ppd_uddu[0];

#if !RUNMINIMAL
	cout << "finished." << endl;
#endif
}/*endvimfold*/
/*
   *
   * Initializes pi_iindex.
   * pi_iindex[i] contains the number of states j where j<i 
   * which has the same angular momentum as i.
   *
   */
void libGRIE::initIIndices()
{/*startvimfold*/
	pi_iindex = new int[ino];
	int i, j, m, n;

	for (m=-ir; m<=ir; m++)
	{
		n=0;
		for (i=0; i<ino; i++)
			if (pi_m[i]==m)
			{
				pi_iindex[i]=n;
				n++;
			}
	}
}/*endvimfold*/
/*
   *
   * Calc the antisymmetric coulomb element.
   * Input: 
   * 	a,b,c,d : num. of the orbs as sorted in pi_m and pi_n.
   * 	i_spinproj: 0 ={ uuuu, dddd}, 1 = {udud, dudu}, 2 = {uddu, duud}.
   *
   */
double libGRIE::calcCoulIE(const int a, const int b, 
		const int c, const int d)
{/*startvimfold*/
	//The class vars. pi_m,pi_n contains the spin and the angular mom.
	//calculate the antisymmetric matrix element.
	return oOFCI->singleElement(
				2*pi_n[a] + abs(pi_m[a]), pi_m[a],
				2*pi_n[b] + abs(pi_m[b]), pi_m[b],
				2*pi_n[c] + abs(pi_m[c]), pi_m[c],
				2*pi_n[d] + abs(pi_m[d]), pi_m[d]); 
		//coulomb(
		//		2*pi_n[a] + abs(pi_m[a]), pi_m[a],
		//		2*pi_n[b] + abs(pi_m[b]), pi_m[b],
		//		2*pi_n[c] + abs(pi_m[c]), pi_m[c],
		//		2*pi_n[d] + abs(pi_m[d]), pi_m[d]); 
}/*endvimfold*/
/*
   *
   * Calculation of the element on the fly. For testing.
   *
   */
double libGRIE::testCIE(int a, int b, int c, int d)
{/*startvimfold*/

	const int imk = pi_m[a/2];
	const int iml = pi_m[b/2];
	const int imp = pi_m[c/2];
	const int imq = pi_m[d/2];

	const int innk = 2*pi_n[a/2] + abs(imk);
	const int innl = 2*pi_n[b/2] + abs(iml);
	const int innp = 2*pi_n[c/2] + abs(imp);
	const int innq = 2*pi_n[d/2] + abs(imq);

	double dret = 0.0;
	if ((a%2==c%2) && (b%2==d%2)) //last test always true if first test is!
	{
		//Generate matrix elements on the fly
		dret += oOFCI->singleElement(innk, imk, innl, iml, innp, imp, innq, imq);
	}
	if ((a%2==d%2) && (b%2==c%2)) 
	{
		//Generate matrix elements on the fly
		dret -= oOFCI->singleElement(innk, imk, innl, iml, innq, imq, innp, imp);
	}
	return dret;
}/*endvimfold*/

//Morten H-J's code below.

/*
   * The function                             
   *      void  **matrix()                    
   * reserves dynamic memory for a two-dimensional matrix 
   * using the C++ command new . No initialization of the elements. 
   * Input data:                      
   *  int row      - number of  rows          
   *  int col      - number of columns        
   *  int num_bytes- number of bytes for each 
   *                 element                  
   * Returns a void  **pointer to the reserved memory location.                                
   */

void **mtrix(int row, int col, int num_bytes)
  {
  int      i, num;
  char     **pointer, *ptr;

  pointer = new(nothrow) char* [row];
  if(!pointer) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for "<< row << "row addresses !" << endl;
    return NULL;
  }
  i = (row * col * num_bytes)/sizeof(char);
  pointer[0] = new(nothrow) char [i];
  if(!pointer[0]) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for address to " << i << " characters !" << endl;
    return NULL;
  }
  ptr = pointer[0];
  num = col * num_bytes;
  for(i = 0; i < row; i++, ptr += num )   {
    pointer[i] = ptr; 
  }

  return  (void **)pointer;

  } // end: function void **matrix()

#endif

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=startvimfold,endvimfold
