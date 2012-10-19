#include "classQDotInteraction.hpp"
#include <iostream>

using std::cout;
using std::endl;
using namespace quantumdot;

int main()
{
	int im_max = 5;
	double d_lambda = 1; //l(omega)  1/sqrt(omega) ??? XXX XXX

	quantumdot::QdotInteraction *oOFCI;
	oOFCI = new QdotInteraction();
	oOFCI->setR(im_max);
	oOFCI->setLambda(d_lambda);
    oOFCI->buildInteractionComBlocks();
	//oOFCI->buildEffectiveInteractionComBlocks(1); //must be one according to Frank.

   double dElement1 = 0.0, dElement2 = 0.0;

   int iN1 = 0;
   int iM1 = 0;
   int iS1 = 1;
   int iN2 = 0;
   int iM2 = 0;
   int iS2 = 0;
   int iN3 = 1;
   int iM3 = 0;
   int iS3 = 0;
   int iN4 = 1;
   int iM4 = 0;
   int iS4 = 1;

   int iNNNN = iN1 + iN2 + iN3 + iN4;
   int iP = 1-2*(iNNNN%2);

   int iNN1 =2*iN1 + abs(iM1);
   int iNN2 =2*iN2 + abs(iM2);
   int iNN3 =2*iN3 + abs(iM3);
   int iNN4 =2*iN4 + abs(iM4);

   if ((iS1%2==iS3%2) && (iS2%2==iS4%2)) {
       dElement1 = oOFCI->singleElement(iNN1, iM1, iNN2, iM2, iNN3, iM3, iNN4, iM4);
   }
   if ((iS1%2==iS4%2) && (iS2%2==iS3%2)) {
       dElement2 = oOFCI->singleElement(iNN1, iM1, iNN2, iM2, iNN4, iM4, iNN3, iM3);
   }
   //return iP*(dElement1-dElement2);
   cout << iP*(dElement1-dElement2) << endl;
 
	cout << "tjohei" << endl;

	cout << "1 2" << endl;
    cout <<  oOFCI->singleElement(1, 1, 1, -1, 1, -1, 1, 1) << endl;
	cout << "2 1" << endl;
    cout <<  oOFCI->singleElement(1, -1, 1, 1, 1, 1, 1, -1) << endl;
	cout << "1 1" << endl;
    cout <<  oOFCI->singleElement(1, 1, 1, -1, 1, 1, 1, -1) << endl;
	cout << "2 2" << endl;
    cout <<  oOFCI->singleElement(1, -1, 1, 1, 1, -1, 1, 1) << endl;
	cout << "0 1" << endl;
    cout <<  oOFCI->singleElement(0, 0, 0, 0, 1, -1, 1, 1) << endl;
	cout << "0 2" << endl;
    cout <<  oOFCI->singleElement(0, 0, 0, 0, 1, 1, 1, -1) << endl;
	cout << "2 0" << endl;
    cout <<  oOFCI->singleElement(1, 1, 1, -1, 0, 0, 0, 0) << endl;
	cout << "1 0" << endl;
    cout <<  oOFCI->singleElement(1, -1, 1, 1, 0, 0, 0, 0) << endl;
	cout << "0 0" << endl;
    cout <<  oOFCI->singleElement(0, 0, 0, 0, 0, 0, 0, 0) << endl;
	cout << "3 3" << endl;
    cout <<  oOFCI->singleElement(2, 0, 2, 0, 2, 0, 2, 0) << endl;
	cout << "3 2" << endl;
    cout <<  oOFCI->singleElement(2, 0, 2, 0, 1, -1, 1, 1) << endl;
	cout << "3 1" << endl;
    cout <<  oOFCI->singleElement(2, 0, 2, 0, 1, 1, 1, -1) << endl;
	cout << "3 0" << endl;
    cout <<  oOFCI->singleElement(2, 0, 2, 0, 0, 0, 0, 0) << endl;
	
	cout << "tjohei" << endl;
		
	for (int i=0; i<1; i++)
	{
   		dElement1 = oOFCI->singleElement(iNN1, iM1, iNN2, iM2, iNN3, iM3, iNN4, iM4);
   		//dElement1 = oOFCI->singleElement(10, 0, 10, 0, 10, 0, 10, 0	);
   		//cout << iP*(dElement1 - dElement2) << endl;
	}
   		cout << iP*(dElement1 - dElement2) << endl;
   		cout << dElement1 << endl;
   		cout << dElement2 << endl;
}
