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
//		double**** initCoulombElem();
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
};
