#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <math.h>

int normalFormIDX (int* vect);
void SBP1 (int* P);
void SBP2 (int* P);
void SBP3 (int* P);
void SBP4 (int* P);
void chooseTB (int* selectedPoints,int* timeBracket, int tPrec);
int gaussianRand(int minInt, int maxInt,float center, float sigma);
int uniformRand(int minInt, int maxInt);
double doubleRand();

////////////////
// Global variables

int Nlen=3000;
int N=10000;
int tau=10;
float centerValue = 0.5;
float sigmaValue = 0.1;

/*
* The 49 normal forms up to T/I invariance are:
* ['0-1','1-1','2-1','2-2','2-3','2-4','2-5','2-6','3-1','3-2','3-3','3-4',
*  '3-5','3-6','3-7','3-8','3-9','3-10','3-11','3-12','4-1','4-2','4-3','4-4',
*  '4-5','4-6','4-7','4-8','4-9','4-10','4-11','4-12','4-13','4-14','4-z15',
*  '4-16','4-17','4-18','4-19','4-20','4-21','4-22','4-23','4-24','4-25',
*  '4-26','4-27','4-28','4-z29']
*
* The 70 normal forms up to T invariance are:
* ['0-1','1-1','2-1','2-2','2-3','2-4','2-5','2-6','3-1','3-2A','3-2B','3-3A',
*  '3-3B', '3-4A','3-4B','3-5A','3-5B','3-6','3-7A','3-7B','3-8A','3-8B','3-9',
*  '3-10','3-11A','3-11B','3-12','4-1','4-2A','4-2B','4-3','4-4A','4-4B','4-5A',
*  '4-5B','4-6','4-7','4-8','4-9','4-10','4-11A','4-11B','4-12A','4-12B',
*  '4-13A','4-13B','4-14A','4-14B','4-z15A','4-z15B','4-16A','4-16B','4-17',
*  '4-18A','4-18B','4-19A','4-19B','4-20','4-21','4-22A','4-22B','4-23','4-24',
*  '4-25','4-26','4-27A','4-27B','4-28','4-z29A','4-z29B']
*/

// This is the correct version of the normal forms, with the special case(0,4,7,8)
// for which the minimum is 281 (and not 401)
int data_normalForms[69] = {1,3,5,9,17,33,65,7,11,13,19,25,35,49,67,97,21, \
                            37,41,69,81,133,73,137,145,273,15,23,29,27,39,57,\
                            71,113,135,51,99,195,45,43,53,77,89,75,105,141,177,\
                            83,101,163,197,153,147,201,275,281,291,85,149,169,\
                            165,277,325,297,293,329,585,139,209};
/*
* This defines a map from the normal forms to the normal form index
* Comment one or the other depending on whether normal forms up to T/I
* invariance or normal forms up to T invariance are needed.
* Adjust the num_normalForms in consequence.
*/
// int map_normalForms[69] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20, \
//                           21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,\
//                           39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,\
//                           56,57,58,59,60,61,62,63,64,65,66,67,68};

int map_normalForms[69] = {1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12,13,14,14,15,15,\
                           16,17,18,18,19,20,21,21,22,23,23,24,24,25,26,27,28,29,\
                           30,30,31,31,32,32,33,33,34,34,35,35,36,37,37,38,38,\
                           39,40,41,41,42,43,44,45,46,46,47,48,48};

// Do not forget to add 1 for the 0-1 normal form
// (so, either 49 or 70)
int num_normalForms = 49;

//////////////////////////////////////////////////////////////
// MAIN

int main(int argc, char *argv[])
{
  FILE *theOutputFile;
	int i,j,k,vectIDX,vectIDXFut;
	int P[4][Nlen];
	int vect[4];
	float** S;
	float*** condS;
	float prog,progStep;

  srand(time(NULL));

	/* Matrix initializations
  * S[i][j] is the probability value that at time t=i, the normal form j is
  * being played.
  *
  * condS[i][j][k] is the conditional probability value that at time t=i+tau
  * the normal form j is being played, knowing that at time t=i the normal
  * form k is being played.
  */
	S = (float**)malloc(sizeof(float*)*Nlen);
	for (i=0;i<Nlen;i++) {
		S[i] = (float*)malloc(sizeof(float)*num_normalForms);
	}

	for (i=0;i<Nlen;i++) {
		for (j=0;j<num_normalForms;j++) {
			S[i][j] = 0.0;
		}
	}

	condS = (float***)malloc(sizeof(float**)*(Nlen-tau));
	for (i=0;i<Nlen-tau;i++) {
		condS[i] = (float**)malloc(sizeof(float*)*num_normalForms);
		for (j=0;j<num_normalForms;j++) {
			condS[i][j] = (float*)malloc(sizeof(float)*num_normalForms);
		}
	}

	for (i=0;i<Nlen-tau;i++) {
		for (j=0;j<num_normalForms;j++) {
			for (k=0;k<num_normalForms;k++) {
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

		SBP1(P[0]);
		SBP2(P[1]);
		SBP3(P[2]);
		SBP4(P[3]);

		for (j=0;j<Nlen;j++) {
			for (k=0;k<4;k++) {
				vect[k]=P[k][j];
			}
			vectIDX=normalFormIDX(vect);
			S[j][vectIDX] = S[j][vectIDX] + (1/(float)N);

			if (j<Nlen-tau) {

			for (k=0;k<4;k++) {
				vect[k]=P[k][j+tau];
			}
			vectIDXFut=normalFormIDX(vect);
      /* At first condS contains the join probability valuethat at time t=i+tau
      * the normal form j is being played and that at time t=i the normal
      * form k is being played. This will be later converted to a conditional
      * probability.
      */
			condS[j][vectIDXFut][vectIDX] = condS[j][vectIDXFut][vectIDX] + (1/(float)N);
			}
		}

	}
	printf("\n");


 	// Computation of the conditional probabilities

	for (i=0;i<Nlen-tau;i++) {
		for (j=0;j<num_normalForms;j++) {
			for (k=0;k<num_normalForms;k++) {
				if (S[i][k]>0)
					condS[i][j][k] = condS[i][j][k]/S[i][k];
			}
		}
	}

 	// Writing output files

  theOutputFile = fopen("outputS.txt","w");
  if (theOutputFile!=NULL) {
    for (i=0;i<Nlen;i++) {
      for (j=0;j<num_normalForms;j++) {
        fprintf(theOutputFile,"%f\n",S[i][j]);
		  }
	  }
    fclose(theOutputFile);
  } else {
    printf("!! Problem opening the ouput S file... !!\n");
  }

  theOutputFile = fopen("outputCondS.txt","w");
  if (theOutputFile!=NULL) {
    for (i=0;i<Nlen-tau;i++) {
		  for (j=0;j<num_normalForms;j++) {
			  for (k=0;k<num_normalForms;k++) {
				  fprintf(theOutputFile,"%f\n",condS[i][j][k]);
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
    for (j=0;j<num_normalForms;j++)
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
//
int normalFormIDX (int* vect) {
  /*
  * vect is the vector corresponding to what the players are playing at time t
  *     i.e. vect[i] is the pitch class value for player i at the given time t
  *
  * Returns a number indicating the prime form
  *
  */
	int i,c,IDX,min;
	unsigned short normalVector;

  IDX=0;
	if (vect[0]==-1 && vect[1]==-1 && vect[2]==-1 && vect[3]==-1) {
		return IDX;
	}

	normalVector = 0;

	for (i=0;i<4;i++) {
		if (vect[i] > -1) {
			normalVector = normalVector | (1 << vect[i]);
		}
	}

	min = 512000;
	if (normalVector<min)
		min=normalVector;

	for (i=1;i<12;i++) {
    // The following code performs a circular (left) bit-shift
    // on the first 12 bits of an unsigned short
    c = ((normalVector << 1) & 4096) == 4096;
    normalVector = ((normalVector << 1) & 4095) | c;

		if (normalVector<min)
			min=normalVector;
	}

	for (i=0;i<69;i++) {
		if (min==data_normalForms[i]) {
			return map_normalForms[i];
		}
	}

	return IDX;
}

//////////////////////////////////////////////////////////////
// Computation of the individual parts of the Number Piece

void SBP1(int* P) {
	int i,tPrec,r1,r2;
	int TB[10][4];
	int resTB[10][2];

	for (i=0;i<Nlen;i++) {
		P[i]=-1;
	}


	TB[0][0]=0; 		TB[0][1]=225; 		TB[0][2]=150; 		TB[0][3]=375;
	TB[1][0]=300; 		TB[1][1]=525; 		TB[1][2]=450; 		TB[1][3]=600+75;
	TB[2][0]=600; 		TB[2][1]=600+225; 	TB[2][2]=600+150;	TB[2][3]=600+375;
	TB[3][0]=600+300; 	TB[3][1]=600+525; 	TB[3][2]=600+450;	TB[3][3]=1200+75;
	TB[4][0]=1200+75; 	TB[4][1]=1200+75;	TB[4][2]=1200+225;	TB[4][3]=1200+225;
	TB[5][0]=1200+225; 	TB[5][1]=1200+450;	TB[5][2]=1200+375;	TB[5][3]=1800;
	TB[6][0]=1200+525; 	TB[6][1]=1800+150;	TB[6][2]=1800+75;	TB[6][3]=1800+300;
	TB[7][0]=1800+225; 	TB[7][1]=1800+450;	TB[7][2]=1800+375;	TB[7][3]=2400;
	TB[8][0]=1800+525; 	TB[8][1]=2400+150;	TB[8][2]=2400+75;	TB[8][3]=2400+300;
	TB[9][0]=2400+225; 	TB[9][1]=2400+450;	TB[9][2]=2400+375;	TB[9][3]=3000;


	tPrec=0;
	for (i=0;i<10;i++) {
		chooseTB(resTB[i],TB[i],tPrec);
		tPrec=resTB[i][1];
	}

	/////// TB1
	r1=uniformRand(resTB[0][0],resTB[0][1]);
	for (i=resTB[0][0];i<r1+1;i++)
		P[i]=10;
	for (i=r1+1;i<resTB[0][1];i++)
		P[i]=7;

	/////// TB2
	for (i=resTB[1][0];i<resTB[1][1];i++)
		P[i]=8;

	/////// TB3

	for (i=resTB[2][0];i<resTB[2][1];i++)
		P[i]=1;

	/////// TB4

	for (i=resTB[3][0];i<resTB[3][1];i++)
		P[i]=9;

	/////// TB5

	for (i=resTB[4][0];i<resTB[4][1];i++)
		P[i]=10;

	/////// TB6

	r1=uniformRand(resTB[5][0],resTB[5][1]);
	r2=uniformRand(r1,resTB[5][1]);
	for (i=resTB[5][0];i<r1;i++)
		P[i]=9;
	for (i=r2;i<resTB[5][1];i++)
		P[i]=1;

	/////// TB7

	for (i=resTB[6][0];i<resTB[6][1];i++)
		P[i]=11;

	/////// TB8

	r1=uniformRand(resTB[7][0],resTB[7][1]);
	for (i=resTB[7][0];i<r1+1;i++)
		P[i]=1;
	for (i=r1+1;i<resTB[7][1];i++)
		P[i]=3;

	/////// TB9

	r1=uniformRand(resTB[8][0],resTB[8][1]);
	r2=uniformRand(r1,resTB[8][1]);
	for (i=resTB[8][0];i<r1;i++)
		P[i]=7;
	for (i=r2;i<resTB[8][1];i++)
		P[i]=9;

	/////// TB10

	for (i=resTB[9][0];i<resTB[9][1];i++)
		P[i]=5;
}

/// DONE SBP2
void SBP2(int* P) {
	int i,tPrec,r1,r2,r3;
	int TB[10][4];
	int resTB[10][2];

	for (i=0;i<Nlen;i++) {
		P[i]=-1;
	}


	TB[0][0]=0; 		TB[0][1]=225; 		TB[0][2]=150; 		TB[0][3]=375;
	TB[1][0]=300; 		TB[1][1]=525; 		TB[1][2]=450; 		TB[1][3]=600+75;
	TB[2][0]=600; 		TB[2][1]=600+225; 	TB[2][2]=600+150;	TB[2][3]=600+375;
	TB[3][0]=600+300; 	TB[3][1]=600+525; 	TB[3][2]=600+450;	TB[3][3]=1200+75;
	TB[4][0]=1200+75; 	TB[4][1]=1200+75;	TB[4][2]=1200+225;	TB[4][3]=1200+225;
	TB[5][0]=1200+225; 	TB[5][1]=1200+450;	TB[5][2]=1200+375;	TB[5][3]=1800;
	TB[6][0]=1200+525; 	TB[6][1]=1800+150;	TB[6][2]=1800+75;	TB[6][3]=1800+300;
	TB[7][0]=1800+225; 	TB[7][1]=1800+450;	TB[7][2]=1800+375;	TB[7][3]=2400;
	TB[8][0]=1800+525; 	TB[8][1]=2400+150;	TB[8][2]=2400+75;	TB[8][3]=2400+300;
	TB[9][0]=2400+225; 	TB[9][1]=2400+450;	TB[9][2]=2400+375;	TB[9][3]=3000;


	tPrec=0;
	for (i=0;i<10;i++) {
		chooseTB(resTB[i],TB[i],tPrec);
		tPrec=resTB[i][1];
	}

	//////// TB1
	for (i=resTB[0][0];i<resTB[0][1];i++)
		P[i]=9;

	//////// TB2
	r1=uniformRand(resTB[1][0],resTB[1][1]);
	for (i=resTB[1][0];i<r1+1;i++)
		P[i]=10;
	for (i=r1+1;i<resTB[1][1];i++)
		P[i]=3;

	//////// TB3
	r1=uniformRand(resTB[2][0],resTB[2][1]);
	for (i=resTB[2][0];i<r1+1;i++)
		P[i]=3;
	for (i=r1+1;i<resTB[2][1];i++)
		P[i]=6;

	//////// TB4
	for (i=resTB[3][0];i<resTB[3][1];i++)
		P[i]=2;

	//////// TB5
	r1=uniformRand(resTB[4][0],resTB[4][1]);
	r2=uniformRand(r1,resTB[4][1]);
	r3=uniformRand(r2,resTB[4][1]);
	for (i=resTB[4][0];i<r1+1;i++)
		P[i]=1;
	for (i=r1+1;i<r2+1;i++)
		P[i]=8;
	for (i=r3;i<resTB[4][1];i++)
		P[i]=0;

	//////// TB6
	for (i=resTB[5][0];i<resTB[5][1];i++)
		P[i]=10;

	//////// TB7
	for (i=resTB[6][0];i<resTB[6][1];i++)
		P[i]=2;

	//////// TB8
	r1=uniformRand(resTB[7][0],resTB[7][1]);
	r2=uniformRand(r1,resTB[7][1]);
	for (i=resTB[7][0];i<r1+1;i++)
		P[i]=6;
	for (i=r2;i<resTB[7][1];i++)
		P[i]=2;

	//////// TB9
	r1=uniformRand(resTB[8][0],resTB[8][1]);
	r2=uniformRand(r1,resTB[8][1]);
	for (i=resTB[8][0];i<r1+1;i++)
		P[i]=6;
	for (i=r2;i<resTB[8][1];i++)
		P[i]=8;

	//////// TB10
	r1=uniformRand(resTB[9][0],resTB[9][1]);
	for (i=resTB[9][0];i<r1+1;i++)
		P[i]=8;
	for (i=r1+1;i<resTB[9][1];i++)
		P[i]=3;
}

//// DONE SBP3
void SBP3(int* P) {
	int i,tPrec,r1,r2;
	int TB[10][4];
	int resTB[10][2];

	for (i=0;i<Nlen;i++) {
		P[i]=-1;
	}


	TB[0][0]=0; 		TB[0][1]=225; 		TB[0][2]=150; 		TB[0][3]=375;
	TB[1][0]=300; 		TB[1][1]=525; 		TB[1][2]=450; 		TB[1][3]=600+75;
	TB[2][0]=600; 		TB[2][1]=600+225; 	TB[2][2]=600+150;	TB[2][3]=600+375;
	TB[3][0]=600+300; 	TB[3][1]=600+525; 	TB[3][2]=600+450;	TB[3][3]=1200+75;
	TB[4][0]=1200+75; 	TB[4][1]=1200+75;	TB[4][2]=1200+225;	TB[4][3]=1200+225;
	TB[5][0]=1200+225; 	TB[5][1]=1200+450;	TB[5][2]=1200+375;	TB[5][3]=1800;
	TB[6][0]=1200+525; 	TB[6][1]=1800+150;	TB[6][2]=1800+75;	TB[6][3]=1800+300;
	TB[7][0]=1800+225; 	TB[7][1]=1800+450;	TB[7][2]=1800+375;	TB[7][3]=2400;
	TB[8][0]=1800+525; 	TB[8][1]=2400+150;	TB[8][2]=2400+75;	TB[8][3]=2400+300;
	TB[9][0]=2400+225; 	TB[9][1]=2400+450;	TB[9][2]=2400+375;	TB[9][3]=3000;


	tPrec=0;
	for (i=0;i<10;i++) {
		chooseTB(resTB[i],TB[i],tPrec);
		tPrec=resTB[i][1];
	}

	//////// TB1
	r1=uniformRand(resTB[0][0],resTB[0][1]);
	r2=uniformRand(r1,resTB[0][1]);
	for (i=resTB[0][0];i<r1;i++)
		P[i]=9;
	for (i=r2;i<resTB[0][1];i++)
		P[i]=0;

	//////// TB2
	r1=uniformRand(resTB[1][0],resTB[1][1]);
	for (i=resTB[1][0];i<r1+1;i++)
		P[i]=7;
	for (i=r1+1;i<resTB[1][1];i++)
		P[i]=5;

	//////// TB3
	r1=uniformRand(resTB[2][0],resTB[2][1]);
	for (i=resTB[2][0];i<r1+1;i++)
		P[i]=7;
	for (i=r1+1;i<resTB[2][1];i++)
		P[i]=4;

	//////// TB4
	r1=uniformRand(resTB[3][0],resTB[3][1]);
	r2=uniformRand(r1,resTB[3][1]);
	for (i=resTB[3][0];i<r1;i++)
		P[i]=11;
	for (i=r2;i<resTB[3][1];i++)
		P[i]=11;

	//////// TB5
	for (i=resTB[4][0];i<resTB[4][1];i++)
		P[i]=5;

	//////// TB6
	for (i=resTB[5][0];i<resTB[5][1];i++)
		P[i]=11;

	//////// TB7
	for (i=resTB[6][0];i<resTB[6][1];i++)
		P[i]=3;

	//////// TB8
	r1=uniformRand(resTB[7][0],resTB[7][1]);
	for (i=resTB[7][0];i<r1+1;i++)
		P[i]=4;
	for (i=r1+1;i<resTB[7][1];i++)
		P[i]=5;

	//////// TB9
	r1=uniformRand(resTB[8][0],resTB[8][1]);
	for (i=resTB[8][0];i<r1+1;i++)
		P[i]=2;
	for (i=r1+1;i<resTB[8][1];i++)
		P[i]=10;

	//////// TB10
	for (i=resTB[9][0];i<resTB[9][1];i++)
		P[i]=4;
}

void SBP4(int* P) {
	int i,tPrec,r1,r2,r3,r4,r5;
	int TB[10][4];
	int resTB[10][2];

	for (i=0;i<Nlen;i++) {
		P[i]=-1;
	}


	TB[0][0]=0; 		TB[0][1]=225; 		TB[0][2]=150; 		TB[0][3]=375;
	TB[1][0]=300; 		TB[1][1]=525; 		TB[1][2]=450; 		TB[1][3]=600+75;
	TB[2][0]=600; 		TB[2][1]=600+225; 	TB[2][2]=600+150;	TB[2][3]=600+375;
	TB[3][0]=600+300; 	TB[3][1]=600+525; 	TB[3][2]=600+450;	TB[3][3]=1200+75;
	TB[4][0]=1200+75; 	TB[4][1]=1200+75;	TB[4][2]=1200+225;	TB[4][3]=1200+225;
	TB[5][0]=1200+225; 	TB[5][1]=1200+450;	TB[5][2]=1200+375;	TB[5][3]=1800;
	TB[6][0]=1200+525; 	TB[6][1]=1800+150;	TB[6][2]=1800+75;	TB[6][3]=1800+300;
	TB[7][0]=1800+225; 	TB[7][1]=1800+450;	TB[7][2]=1800+375;	TB[7][3]=2400;
	TB[8][0]=1800+525; 	TB[8][1]=2400+150;	TB[8][2]=2400+75;	TB[8][3]=2400+300;
	TB[9][0]=2400+225; 	TB[9][1]=2400+450;	TB[9][2]=2400+375;	TB[9][3]=3000;


	tPrec=0;
	for (i=0;i<10;i++) {
		chooseTB(resTB[i],TB[i],tPrec);
		tPrec=resTB[i][1];
	}

	///////// TB1
	r1=uniformRand(resTB[0][0],resTB[0][1]);
	r2=uniformRand(r1,resTB[0][1]);
	r3=uniformRand(r2,resTB[0][1]);
	for (i=resTB[0][0];i<r1+1;i++)
		P[i]=2;
	for (i=r1+1;i<r2+1;i++)
		P[i]=3;
	for (i=r2+1;i<r3+1;i++)
		P[i]=8;
	for (i=r3+1;i<resTB[0][1];i++)
		P[i]=1;

	///////// TB2
	r1=uniformRand(resTB[1][0],resTB[1][1]);
	r2=uniformRand(r1,resTB[1][1]);
	r3=uniformRand(r2,resTB[1][1]);
	for (i=resTB[1][0];i<r1+1;i++)
		P[i]=9;
	for (i=r2;i<r3+1;i++)
		P[i]=5;
	for (i=r3+1;i<resTB[1][1];i++)
		P[i]=6;

	///////// TB3
	r1=uniformRand(resTB[2][0],resTB[2][1]);
	r2=uniformRand(r1,resTB[2][1]);
	for (i=resTB[2][0];i<r1+1;i++)
		P[i]=1;
	for (i=r2;i<resTB[2][1];i++)
		P[i]=0;

	///////// TB4
	for (i=resTB[3][0];i<resTB[3][1];i++)
		P[i]=9;

	///////// TB5
	r1=uniformRand(resTB[4][0],resTB[4][1]);
	r2=uniformRand(r1,resTB[4][1]);
	r3=uniformRand(r2,resTB[4][1]);
	for (i=resTB[4][0];i<r1+1;i++)
		P[i]=4;
	for (i=r1+1;i<r2+1;i++)
		P[i]=5;
	for (i=r3;i<resTB[4][1];i++)
		P[i]=6;

	///////// TB6
	r1=uniformRand(resTB[5][0],resTB[5][1]);
	r2=uniformRand(r1,resTB[5][1]);
	for (i=resTB[5][0];i<r1+1;i++)
		P[i]=8;
	for (i=r2;i<resTB[5][1];i++)
		P[i]=10;

	///////// TB7
	r1=uniformRand(resTB[6][0],resTB[6][1]);
	r2=uniformRand(r1,resTB[6][1]);
	r3=uniformRand(r2,resTB[6][1]);
	r4=uniformRand(r2,resTB[6][1]);
	r5=uniformRand(r2,resTB[6][1]);
	for (i=resTB[6][0];i<r1+1;i++)
		P[i]=11;
	for (i=r1+1;i<r2+1;i++)
		P[i]=3;
	for (i=r3;i<r4+1;i++)
		P[i]=6;
	for (i=r4+1;i<r5+1;i++)
		P[i]=10;
	for (i=r5+1;i<resTB[6][1];i++)
		P[i]=0;

	///////// TB8
	r1=uniformRand(resTB[7][0],resTB[7][1]);
	r2=uniformRand(r1,resTB[7][1]);
	for (i=resTB[7][0];i<r1+1;i++)
		P[i]=11;
	for (i=r2;i<resTB[7][1];i++)
		P[i]=2;

	///////// TB9
	r1=uniformRand(resTB[8][0],resTB[8][1]);
	for (i=resTB[8][0];i<r1+1;i++)
		P[i]=0;
	for (i=r1+1;i<resTB[8][1];i++)
		P[i]=2;

	///////// TB10
	for (i=resTB[9][0];i<resTB[9][1];i++)
		P[i]=5;
}

//////////////////////////////////////////////////////////////
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

	c=center*(double)minInt+(1.0-center)*(double)maxInt;
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
