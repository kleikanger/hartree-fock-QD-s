
#include "hamiltonianElem.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "dsyev.h"
#include "newmatrix.h"

using namespace std;

class hartreeFock
{
	private: 
		int i_shells, i_numpart, i_numstates;
		double d_enull;
		//single particle energies
		double *pd_spe;
		//hamiltoniam matrix
		double **ppd_h;
		//eigenvals
		double *pd_eigv;
		//eigenvectors
		double **ppd_eigvecs;
		double *pd_eigvold;
		//class handles all interaction elements
		hamiltonianElem *ohe;

	public:
		hartreeFock(int i_numpartARG, int i_shellsARG)
		{
			i_numpart = i_numpartARG;
			i_shells = i_shellsARG;
			i_numstates = i_shells*(i_shells+1);
			
			//This class stores and tabulates the coulomb interaction elements
			ohe = new hamiltonianElem(i_numstates, i_shells, i_numpart);
			//and initiates the single particle energies.
			pd_spe = ohe->singleParticletE();
	
			//reserve mem for eigvecs and eigvals
	   		pd_eigv = new double[i_numstates];
   			pd_eigvold = new double[i_numstates];
			for (int i=0; i<i_numstates; i++)
				pd_eigvold[i] = 0;
			ppd_eigvecs = matrix<double>(i_numstates, i_numstates);
			ppd_h = matrix<double>(i_numstates, i_numstates);
		};
		~hartreeFock(){};
		void runHF(int, int);
		void hfIter();
		double calcEnergy();
};

void hartreeFock::runHF(int i_numpartARG, int i_shellsARG)
{
	i_numpart = i_numpartARG;
	i_shells = i_shellsARG;
	i_numstates = i_shells*(i_shells+1);

	int i, j, k, l, m, n;

	//init the unitary matrix as the identity matrix martix	
	for (i=0; i<i_numstates; i++) 
		for (j=0; j<i_numstates; j++)
			ppd_eigvecs[i][j] = 0.0;
	for (i=0; i<i_numstates; i++)
		ppd_eigvecs[i][i] = 1.0;	
	
	int i_numiter = 0, i_converged;
	bool b_iterate = true;
	double d_ctol = 1e-4;

	while (b_iterate) 
	{
		hfIter();
		
		if (i_numiter>0)
		{	
			i_converged = 0;
			for (i = 0; i < i_numstates; i++) 
			{
				if (fabs((pd_eigv[i] - pd_eigvold[i])/pd_eigv[i]) < d_ctol) 
				{
					i_converged++;
				}
			}
			if (i_converged == i_numstates) 
			{
				b_iterate = false;
			}
			for (int i = 0; i < i_numstates; i++)
			{
				pd_eigvold[i] = pd_eigv[i]; 	
			}
		}
	
		i_numiter++;

		cout << "E:" << calcEnergy() << endl;
	}//End of while (b_iterate)

	cout << "Number of iterations = " << i_numiter <<"\n";
}
/*
   *
   *
   *
   */
double hartreeFock::calcEnergy()
{
    double d_energy = 0;
    for (int a = 0; a < i_numpart; a++) 
        for (int beta = 0; beta < i_numstates; beta++) 
		{
            d_energy += 
				ppd_eigvecs[beta][a]
				* ppd_eigvecs[beta][a]
				* pd_spe[beta];
        }

	//calc E
    for (int a = 0; a < i_numpart; a++) 
        for (int b = 0; b < i_numpart; b++)
            for (int mu = 0; mu < i_numstates; mu++) 
                for (int nu = 0; nu < i_numstates; nu++) 
                    for (int beta = 0; beta < i_numstates; beta++) 
                        for (int delta = 0; delta < i_numstates; delta++)
						{	
						d_energy  += 0.5 
							* ppd_eigvecs[mu][a] 
							* ppd_eigvecs[nu][b]
							* ppd_eigvecs[beta][a]
							* ppd_eigvecs[delta][b]
							* ohe->cCIE(mu, nu, beta, delta);
                        }
	return d_energy;
}
/*
   *
   *
   *
   */
void hartreeFock::hfIter()
{
	//First part; Coulomb energies.
	for (int alpha = 0; alpha < i_numstates; alpha++) 
		for (int gamma = alpha; gamma < i_numstates; gamma++) 
		{
			ppd_h[alpha][gamma] = 0;

			for (int a = 0; a < i_numpart; a++) 
				for (int beta = 0; beta < i_numstates; beta++) 
					for (int delta = 0; delta < i_numstates; delta++) 
					{
						ppd_h[alpha][gamma] += 
							ppd_eigvecs[beta][a] 
							* ppd_eigvecs[delta][a] 
							* ohe->cCIE(alpha, beta, gamma, delta);
					}
		}
	//Second part: Adding eigenenergies of basis states.
	for (int i=0;i<i_numstates;i++)
	{
		ppd_h[i][i] += pd_spe[i];
	}
	dsyev(ppd_h, i_numstates, pd_eigv, ppd_eigvecs);
}
/*
   *
   * MAIN: RUN HF SIMULATIONS
   *
   */
//MAIN:
int main()
{
	hartreeFock ohf(6, 10);

	cout << "2 pt" << endl;
	for (int i=3; i<11; i++)
	{
		cout << "num shells = " << i << endl;
		ohf.runHF(2, i);
	}
	
	cout << "6 pt" << endl;
	for (int i=3; i<11; i++)
	{
		cout << "num shells = " << i << endl;
		ohf.runHF(6, i);
	}
}
