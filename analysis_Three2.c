#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int normalFormIDX (int* vect);
void SAP1 (int* P);
void SAP2 (int* P);
void SAP3 (int* P);
void chooseTB (int* selectedPoints,int* timeBracket, int tPrec);
int gaussianRand(int minInt, int maxInt,float center, float sigma);
int uniformRand(int minInt, int maxInt);
double doubleRand();

////////////////
// Global variables

int Nlen=5400;
int N=1000000;
int tau=10;
float centerValue = 0.5;
float sigmaValue = 0.5;

////////////////
// MAIN

int main(int argc, char *argv[])
{
  FILE *theOutputFile;
	int i,j,k,vectIDX,vectIDXFut;
	int P[3][Nlen];
	int vect[3];
	double** S;
	double*** condS;
	float prog,progStep;

  srand(time(NULL));

	// Matrix initializations

	S = (double**)malloc(sizeof(double*)*Nlen);
	for (i=0;i<Nlen;i++) {
		S[i] = (double*)malloc(sizeof(double)*64);
	}

	for (i=0;i<Nlen;i++) {
		for (j=0;j<64;j++) {
			S[i][j] = 0.0;
		}
	}

	condS = (double***)malloc(sizeof(double**)*(Nlen-tau));
	for (i=0;i<Nlen-tau;i++) {
		condS[i] = (double**)malloc(sizeof(double*)*64);
		for (j=0;j<64;j++) {
			condS[i][j] = (double*)malloc(sizeof(double)*64);
		}
	}

	for (i=0;i<Nlen-tau;i++) {
		for (j=0;j<64;j++) {
			for (k=0;k<64;k++) {
				condS[i][j][k] = 0.0;
			}
		}
	}

  // Calculating the Number Piece realizations

	progStep=0.01;
	printf("[");
	for (prog=0;prog<=1-progStep;prog+=progStep)
		printf(" ");
	printf("]\n x");fflush(stdout);
	prog=progStep;

	for (i=0;i<N;i++) {
		if ((float)i/(float)N>prog) {
			printf("x");fflush(stdout);
			prog+=progStep;
		}

		SAP1(P[0]);
		SAP2(P[1]);
		SAP3(P[2]);

		for (j=0;j<Nlen;j++) {
			for (k=0;k<3;k++) {
				vect[k]=P[k][j];
			}
			vectIDX=normalFormIDX(vect);
			S[j][vectIDX] = S[j][vectIDX] + (1/(double)N);

			if (j<Nlen-tau) {

			for (k=0;k<3;k++) {
				vect[k]=P[k][j+tau];
			}
			vectIDXFut=normalFormIDX(vect);
			condS[j][vectIDXFut][vectIDX] = condS[j][vectIDXFut][vectIDX] + (1/(double)N);
			}
		}

	}
	printf("\n");


 	// Computation of the conditional probabilities

	for (i=0;i<Nlen-tau;i++) {
		for (j=0;j<64;j++) {
			for (k=0;k<64;k++) {
				if (S[i][k]>0)
					condS[i][j][k] = condS[i][j][k]/S[i][k];
			}
		}
	}

 	// Writing output files

  theOutputFile = fopen("outputS.txt","w");
  if (theOutputFile!=NULL) {
        for (i=0;i<Nlen;i++) {
		for (j=0;j<64;j++) {
			fprintf(theOutputFile,"%.8e\n",S[i][j]);
		}
	  }
        fclose(theOutputFile);
  } else {
         printf("!! Problem opening the ouput S file... !!\n");
  }

  theOutputFile = fopen("outputCondS.txt","w");
  if (theOutputFile!=NULL) {
        for (i=0;i<Nlen-tau;i++) {
		for (j=0;j<64;j++) {
			for (k=0;k<64;k++) {
				fprintf(theOutputFile,"%.8e\n",condS[i][j][k]);
			}
		}
	  }
        fclose(theOutputFile);
  } else {
         printf("!! Problem opening the ouput condS file... !!\n");
  }
	printf("Done\n");

 	// Freeing memory

	for (i=0;i<Nlen-tau;i++) {
			for (j=0;j<64;j++)
				free(condS[i][j]);
		free(condS[i]);
	}
	free(condS);

	for (i=0;i<Nlen;i++) {
		free(S[i]);
	}
	free(S);
}


// Prime form calculation

int normalFormIDX (int* vect) {
  /*
  * vect is the vector corresponding to what the players are playing at time t
  *
  * Returns a number indicating the prime form
  *
  */
	return vect[0]+4*vect[1]+16*vect[2];
}

// Computation of the individual parts of the Number Piece

void SAP1(int* P) {
  // Computation of player 1 part
	int i,tPrec;
	int TB[9][4];
	int resTB[9][2];

	for (i=0;i<Nlen;i++) {
		P[i]=0;
	}

  // Time-brackets for this Number Piece are hard-coded in the program
	TB[0][0]=0; 		TB[0][1]=450; 		TB[0][2]=300; 		TB[0][3]=750;
	TB[1][0]=600; 		TB[1][1]=1050; 		TB[1][2]=900; 		TB[1][3]=1350;
	TB[2][0]=1200;		TB[2][1]=1650;		TB[2][2]=1500;		TB[2][3]=1950;
	TB[3][0]=1800;		TB[3][1]=2250;		TB[3][2]=2100;		TB[3][3]=2550;
	TB[4][0]=2550;		TB[4][1]=2550;		TB[4][2]=2850;		TB[4][3]=2850;
	TB[5][0]=2850;		TB[5][1]=3300;		TB[5][2]=3150;		TB[5][3]=3600;
	TB[6][0]=3450;		TB[6][1]=3900;		TB[6][2]=3750;		TB[6][3]=4200;
	TB[7][0]=4050;		TB[7][1]=4500;		TB[7][2]=4350;		TB[7][3]=4800;
	TB[8][0]=4650;		TB[8][1]=5100;		TB[8][2]=4950;		TB[8][3]=5400;


	tPrec=0;
	for (i=0;i<9;i++) {
		chooseTB(resTB[i],TB[i],tPrec);
		tPrec=resTB[i][1];
	}

	////// TB1
	for (i=resTB[0][0];i<resTB[0][1];i++)
		P[i]=3;

	////// TB2
	for (i=resTB[1][0];i<resTB[1][1];i++)
		P[i]=2;

	////// TB3
	for (i=resTB[2][0];i<resTB[2][1];i++)
		P[i]=1;

	////// TB4
	for (i=resTB[3][0];i<resTB[3][1];i++)
		P[i]=3;

	////// TB5
	for (i=resTB[4][0];i<resTB[4][1];i++)
		P[i]=3;

	////// TB6
	for (i=resTB[5][0];i<resTB[5][1];i++)
		P[i]=2;

	////// TB7
	for (i=resTB[6][0];i<resTB[6][1];i++)
		P[i]=3;

	////// TB8
	for (i=resTB[7][0];i<resTB[7][1];i++)
		P[i]=3;

	////// TB9
	for (i=resTB[8][0];i<resTB[8][1];i++)
		P[i]=1;
}

void SAP2(int* P) {
  // Computation of player 2 part
	int i,tPrec;
	int TB[9][4];
	int resTB[9][2];

	for (i=0;i<Nlen;i++) {
		P[i]=0;
	}

  // Time-brackets for this Number Piece are hard-coded in the program
	TB[0][0]=0; 		TB[0][1]=450; 		TB[0][2]=300; 		TB[0][3]=750;
	TB[1][0]=750; 		TB[1][1]=750; 		TB[1][2]=1050; 		TB[1][3]=1050;
	TB[2][0]=1050;		TB[2][1]=1500;		TB[2][2]=1350;		TB[2][3]=1800;
	TB[3][0]=1650;		TB[3][1]=2100;		TB[3][2]=1950;		TB[3][3]=2400;
	TB[4][0]=2250;		TB[4][1]=2700;		TB[4][2]=2550;		TB[4][3]=3000;
	TB[5][0]=2850;		TB[5][1]=3300;		TB[5][2]=3150;		TB[5][3]=3600;
	TB[6][0]=3450;		TB[6][1]=3900;		TB[6][2]=3750;		TB[6][3]=4200;
	TB[7][0]=4050;		TB[7][1]=4500;		TB[7][2]=4350;		TB[7][3]=4800;
	TB[8][0]=4650;		TB[8][1]=5100;		TB[8][2]=4950;		TB[8][3]=5400;


	tPrec=0;
	for (i=0;i<9;i++) {
		chooseTB(resTB[i],TB[i],tPrec);
		tPrec=resTB[i][1];
	}

	////// TB1
	for (i=resTB[0][0];i<resTB[0][1];i++)
		P[i]=1;

	////// TB2
	for (i=resTB[1][0];i<resTB[1][1];i++)
		P[i]=1;

	////// TB3
	for (i=resTB[2][0];i<resTB[2][1];i++)
		P[i]=1;

	////// TB4
	for (i=resTB[3][0];i<resTB[3][1];i++)
		P[i]=3;

	////// TB5
	for (i=resTB[4][0];i<resTB[4][1];i++)
		P[i]=2;

	////// TB6
	for (i=resTB[5][0];i<resTB[5][1];i++)
		P[i]=1;

	////// TB7
	for (i=resTB[6][0];i<resTB[6][1];i++)
		P[i]=3;

	////// TB8
	for (i=resTB[7][0];i<resTB[7][1];i++)
		P[i]=2;

	////// TB9
	for (i=resTB[8][0];i<resTB[8][1];i++)
		P[i]=2;
}

void SAP3(int* P) {
  // Computation of player 3 part
	int i,tPrec;
	int TB[9][4];
	int resTB[9][2];

	for (i=0;i<Nlen;i++) {
		P[i]=0;
	}

  // Time-brackets for this Number Piece are hard-coded in the program
	TB[0][0]=0; 		TB[0][1]=450; 		TB[0][2]=300; 		TB[0][3]=750;
	TB[1][0]=600; 		TB[1][1]=1050; 		TB[1][2]=900; 		TB[1][3]=1350;
	TB[2][0]=1200;		TB[2][1]=1650;		TB[2][2]=1500;		TB[2][3]=1950;
	TB[3][0]=1800;		TB[3][1]=2250;		TB[3][2]=2100;		TB[3][3]=2550;
	TB[4][0]=2400;		TB[4][1]=2850;		TB[4][2]=2700;		TB[4][3]=3150;
	TB[5][0]=3000;		TB[5][1]=3450;		TB[5][2]=3300;		TB[5][3]=3750;
	TB[6][0]=3600;		TB[6][1]=4050;		TB[6][2]=3900;		TB[6][3]=4350;
	TB[7][0]=4350;		TB[7][1]=4350;		TB[7][2]=4650;		TB[7][3]=4650;
	TB[8][0]=4650;		TB[8][1]=5100;		TB[8][2]=4950;		TB[8][3]=5400;


	tPrec=0;
	for (i=0;i<9;i++) {
		chooseTB(resTB[i],TB[i],tPrec);
		tPrec=resTB[i][1];
	}

	////// TB1
	for (i=resTB[0][0];i<resTB[0][1];i++)
		P[i]=2;

	////// TB2
	for (i=resTB[1][0];i<resTB[1][1];i++)
		P[i]=2;

	////// TB3
	for (i=resTB[2][0];i<resTB[2][1];i++)
		P[i]=2;

	////// TB4
	for (i=resTB[3][0];i<resTB[3][1];i++)
		P[i]=2;

	////// TB5
	for (i=resTB[4][0];i<resTB[4][1];i++)
		P[i]=2;

	////// TB6
	for (i=resTB[5][0];i<resTB[5][1];i++)
		P[i]=2;

	////// TB7
	for (i=resTB[6][0];i<resTB[6][1];i++)
		P[i]=1;

	////// TB8
	for (i=resTB[7][0];i<resTB[7][1];i++)
		P[i]=1;

	////// TB9
	for (i=resTB[8][0];i<resTB[8][1];i++)
		P[i]=3;
}


// Routines for random choice

void chooseTB(int* selectedPoints,int* timeBracket, int tPrec) {
	int ts,te;
  
  // Replace 'gaussianRand' by the appropriate routine to change the
  // probability distribution being used
	if (tPrec>timeBracket[0]) {
        ts=gaussianRand(tPrec,timeBracket[1],centerValue,sigmaValue);
	} else {
        ts=gaussianRand(timeBracket[0],timeBracket[1],centerValue,sigmaValue);
	}

	if (ts>timeBracket[2]) {
        te=gaussianRand(ts,timeBracket[3],1.0-centerValue,sigmaValue);
	} else {
        te=gaussianRand(timeBracket[2],timeBracket[3],1.0-centerValue,sigmaValue);
	}

	selectedPoints[0]=ts;
	selectedPoints[1]=te;
}

/// Gaussian distribution selection

int gaussianRand(int minInt, int maxInt,float center, float sigma) {
  /*
  * minInt, maxInt: limits of the time interval
  * center, sigma: parameters of the gaussian distribution
  *               center should be between 0.0 (minInt) and 1.0 (maxInt)
  *               sigma is indicated in proportion of the length maxInt-minInt
  *
  * Returns a number indicating the prime form
  *
  */
	double c,s,u,h,r;

	c=center*(double)minInt+(1.0-x)*(double)maxInt;
	s=sigma*((double)maxInt-(double)minInt);

	u=doubleRand()*((double)maxInt-(double)minInt)+(double)minInt;
	h=exp(-(u-c)*(u-c)/(2*s*s));
	r=doubleRand();

	while (r>h) {
		u=doubleRand()*((double)maxInt-(double)minInt)+(double)minInt;
		h=exp(-(u-c)*(u-c)/(2*s*s));
		r=doubleRand();
	}

	return (int)u;
}

int uniformRand(int minInt, int maxInt) {
	return (minInt+(int)(doubleRand()*((double)maxInt - (double)minInt)));
}


double doubleRand()
{
    return (double)rand() / (double)RAND_MAX ;
}
