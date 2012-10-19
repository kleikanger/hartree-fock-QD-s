#include <cmath>
#include <iostream>
#include "libGRIE.h"

class hamiltonianElem
{
	private:
		int i_numstates;
		int i_shells;
		int i_numpart;
		int *pi_n;
		int *pi_m;
		libGRIE *ogrie;

	public:
		hamiltonianElem(int i_numstatesARG, int i_shellsARG, int i_numpartARG)
		{
			i_numstates = i_numstatesARG;
			i_shells 	= i_shellsARG;
			i_numpart 	= i_numpartARG;
		
			std::cout << i_numstates << std::endl;

			//init arrays with m and n
			pi_n = new int[i_numstates];
			pi_m = new int[i_numstates];
			int i, j;
			j=-1, i=0;
			while (j<i_shells-1)
			{
			j++;
				for (int k=-j;k<=j;)
				{
					pi_n[i]=(int)(j-abs(k))/2;
					pi_m[i]=k;
					k+=2, i++;
				}
			}
			ogrie = new libGRIE(pi_m, pi_n, i_shells-1, i_numstates, 1.0, false, true);
			ogrie->initGRIEmatrices();
			ogrie->initIIndices();
			ogrie->initInteractionElements();
		};
		~hamiltonianElem()
		{
			delete [] pi_n;
			delete [] pi_m;
		};
		/*
		   * initiate the coulomb interaction elements
		   */
		double**** initCoulombElem();
		/*
		   * return pointer to array with the single particle energies
		   */
		double* singleParticletE()
		{
			double *pd_spe = new double[i_numstates];
			for (int i=0; i<i_numstates; i++)
			{
				pd_spe[i] = 1.0 + 2.0*pi_n[i/2] + fabs(pi_m[i/2]);
			}
			return pd_spe;
		};
		double cCIE(int i, int j, int k, int l)
		{
			if ((i%2+j%2==k%2+l%2) && (pi_m[i/2]+pi_m[j/2]==pi_m[k/2]+pi_m[l/2]))
			{
//				std::cout << ogrie->testCIE(9,0,1,0) << "\n";
				return ogrie->getCoulIE(i, j, k, l);
//				return ogrie->calcCoulIE(i, j, k, l);
			}
			else
				return 0.0;
		}
		/*
		   * identify all basis states with spin = 0 and anngular momentum = 0
		   * return : pointer to an (i_numstates x i_dimh) array
		   * i_dimh is calculated in the first while loop
		   */
		int **basisStates(int &i_dimh)
		{
			int i, j, is, im, k, i_temp, i_nallowed;
			int **ppi_bs;
			int pi_t[i_numpart];

			i_nallowed = 0;

			for (i=0; i<i_numpart; i++)
				pi_t[i] = i;

			while (pi_t[0]<=i_numstates-i_numpart)
			{
				
				is = im = 0;
				for (k=0; k<i_numpart; k++)
				{
					im += pi_m[pi_t[k]/2];
					is += pi_t[k]%2;
				}
				if ( (is==1) && (im==0) )
				{
					i_nallowed++;
				}

				pi_t[i_numpart-1]++;
				if (pi_t[i_numpart-1]==i_numstates)
				{
					j = i_numpart-1;
					k = 0;
					while ( (j>0) && (pi_t[j]==i_numstates-k) )
					{
						k++;
						pi_t[j-1]++;
						pi_t[j]=pi_t[j-1];
						j--;
					}
					j+=1;
					while (j<i_numpart)
					{
						pi_t[j]=pi_t[j-1]+1;
						j++;
					}
				}
			}

			ppi_bs = new int*[i_nallowed];
			for (i=0; i<i_nallowed; i++)
				ppi_bs[i] = new int[i_numpart];
			
			for (i=0; i<i_numpart; i++)
				pi_t[i] = i;
			i_nallowed = 0;
			while (pi_t[0]<=i_numstates-i_numpart)
			{
				
				is = im = 0;
				for (k=0; k<i_numpart; k++)
				{
					im += pi_m[pi_t[k]/2];
					is += pi_t[k]%2;
				}
				if ( (is==1) && (im==0) )
				{
					for (k=0; k<i_numpart; k++)
						ppi_bs[i_nallowed][k] = pi_t[k];
					i_nallowed++;
				}

				pi_t[i_numpart-1]++;
				if (pi_t[i_numpart-1]==i_numstates)
				{
					j = i_numpart-1;
					k = 0;
					while ( (j>0) && (pi_t[j]==i_numstates-k) )
					{
						k++;
						pi_t[j-1]++;
						pi_t[j]=pi_t[j-1];
						j--;
					}
					j+=1;
					while (j<i_numpart)
					{
						pi_t[j]=pi_t[j-1]+1;
						j++;
					}
				}
			}
			i_dimh = i_nallowed;
	/*	
			//PRINT	
			std::cout << "number of basis states: " << i_nallowed << std::endl;
			std::cout << "basis states: " << std::endl;
			for (i=0; i<i_nallowed; i++)
			{
				std::cout << "|";
				for (k=0; k<i_numpart; k++)
					std::cout << ppi_bs[i][k] << ",";
				std::cout << ">" << std::endl;
			}
*/
			
			return ppi_bs;
		};
};
