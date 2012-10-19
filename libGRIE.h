#ifndef libGRIE_H
#define libGRIE_H

//ref: Simen Kvaal OpenFCI
#include "OpenFCI/classQDotInteraction.hpp"
using namespace quantumdot;

class libGRIE
{
	private:
		double **ppd_gm1 = nullptr, **ppd_gm2 = nullptr;
		double ***pppd_cie = nullptr;
		int *pi_iindex = nullptr, *pi_m, *pi_n;
		int ino, ino1, ino2, ino3, ir;
		QdotInteraction *oOFCI;
		
	public:
		libGRIE(int*, int*, const int, const int, const double, bool, bool);
		~libGRIE();
		void initGRIEmatrices();
		void initInteractionElements();
		void initIIndices();
		//double getCoulIE(int, int, int, int) const;
		double calcCoulIE(const int, const int,	const int, const int);
		double testCIE(int a, int b, int c, int d);
/*
   *
   * It is assumed that the imput orbs a, b, c, d conserves the spin and the 
   * angular momentum:
   *
   *   m(p)+m(q) = m(r)+m(s) and s(p)+s(q) = s(r)+s(s)
   *
   * Declare the file here to inline in code. 
   *
   */
		inline double getCoulIE(int a, int b, int c, int d) const
		{
			double dsign = 1;
			int i_spinproj, i_ngm, i_indx, i_temp;
			//switch s.t. a<=b and b<=c
			//check spin projection (0,1,2). 1 1/3 if tests on average.
			if (a>b)
			{
				i_temp = a;
				a = b;
				b = i_temp;
				dsign *= -1.;
			}
			if (c>d)
			{
				i_temp = c;
				c = d;
				d = i_temp;
			   dsign *= -1.;   
			}
			if (a%2==c%2)
			{
				if (b%2==c%2)
					//uuuu or dddd
					i_spinproj = 0;
				else 
					// udud or dudu
					i_spinproj = 1;
			}
			else //p==s and p!=r and symmetry means uddu
			{
				//uddu or duud
				i_spinproj = 2;
			}
			//Integer div.. Find the nr. of the orb wf. 
			a = a/2, b = b/2, c = c/2, d = d/2;
			//Find the indices of the CIE. Never more that 4 if tests
			if (a==c)  
			{
				if ( (a==b) || (c==d) )
				{
				i_indx = pi_iindex[a];
					//case I
					if (b==c)
					{
						i_ngm = ppd_gm1[b][0] + ppd_gm1[c][1] + ppd_gm1[d][2];
					}
					else
					{
						i_ngm = ppd_gm1[d][0] + ppd_gm1[c][1] + ppd_gm1[b][2];
					}
				}
				else
				{
					//case III
					i_indx = a; //all a allowed!
					if (b>d)
					{
						i_ngm = ppd_gm1[d-1][0] + ppd_gm1[b-1][1] + ino1 + ino2;
					}
					else
					{
						i_ngm = ppd_gm1[b-1][0] + ppd_gm1[d-1][1] + ino1 + ino2;
					}
				}
			}
			else if ((b<=c)||(d<=a))
			{
				//case I
				if (d<=a) //cdab
				{
					i_indx = pi_iindex[c];
					i_ngm = ppd_gm1[d][0] + ppd_gm1[a][1] + ppd_gm1[b][2];
				}
				else //abcd
				{
					i_indx = pi_iindex[a];
					i_ngm = ppd_gm1[b][0] + ppd_gm1[c][1] + ppd_gm1[d][2];
				}
			}
			else
			{
				//case II
				if (a<c) //abcd
				{
					i_indx = pi_iindex[a];
					i_ngm = ppd_gm2[b][0] + ppd_gm2[c][1] + ppd_gm2[d][2] + ino1;
				}
				else //cdab
				{
					i_indx = pi_iindex[c];
					i_ngm = ppd_gm2[d][0] + ppd_gm2[a][1] + ppd_gm2[b][2] + ino1;
				}
			}
			return dsign * pppd_cie[i_spinproj][i_ngm][i_indx];
		}
};

#endif
