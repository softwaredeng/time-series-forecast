

//------original c++ function
//#include "stdafx.h"
#include <tchar.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>  
//#include <math.h>
#include <cmath>
using namespace std;
#define PI 3.14159265

#include "mex.h"//for matlab


double linfit(double *const x, double *const y, int n)
{
double m;//slope
double b;//intercept
int error = 0;
double sumx = 0,
sumy = 0,
sumx2 = 0,
sumxy = 0;
double dn = (double) n;
int i;
if (n <= 1) {
m = 0;
b = 0;
error = 1;
} else {
double divisor;
//error = 0;
for (i = 0; i < n; i++) {
sumx += x[i];
sumy += y[i];
sumx2 += (x[i] * x[i]);
sumxy += (x[i] * y[i]);
}
divisor = (sumx2 - ((sumx * sumx) / dn));
if (divisor != 0) {
m = (sumxy - ((sumx * sumy) / dn)) / divisor;
b = (sumy - ((m) * sumx)) / dn;
} else {
m = 0;
b = 0;
//error = 2;
}
}
return(m);
}
double entropy(double *set, int setN)// calculate the entropy of a set of class labels
{
	double * clsSet = NULL;
	double * clsSetN = NULL;// clsSet will store the class label found, initilized as a unpossible class label; clsSetN store the number of the corresponding classes, initilized as 0;
	clsSet = new double[setN]; 
	clsSetN = new double[setN];
	
	int count = 0; 

	for(int i=0;i<setN;i++)
	{
		clsSet[i]=-1000;
		clsSetN[i]=0;
	}
	
	for(int i=0;i<setN;i++)//for each class in set, firstly see if it exists; if yes, then increase the corresponding number; if no, then store. 
	{
		int found = 0;
		//search if the class exist; if exist, then found, and increase the class number onece for that class
		for(int j=0;j<count;j++)
		{
			if(clsSet[j]==set[i])
			{
				clsSetN[j]++;
				found = 1;
				break;
			}
		}
		//if not exist, then remember the class, and increase the class number 
		if(found==0)
		{
			count++;
			clsSet[count-1]=set[i];
			clsSetN[count-1]++;
		}

	}
	
	//calculate entropy
	double e=0;
	for(int i=0;i<count;i++)
	{
		//cout<<"\n"<<clsSetN[i]/setN<<"-";
		e = e - clsSetN[i]/setN*  log(double(clsSetN[i]/setN))/log(double(2));
		//cout<<clsSetN[i]/setN*  log(double(clsSetN[i]/setN))/log(double(2))<<"-";
	}
	
	
	delete [] clsSet;  delete [] clsSetN;
	clsSet = NULL; clsSetN = NULL;

	return(e);

}

void BestCut(double *v,double rangeV, double *cls, int nrow, double *cutSet, int cutSetN, double &bestE, double &bestCut,double &bestEntropy, double &bestMargin, double alpha) //v: a set of values; cls: the class corresponding to the v; nrow: row of v; cutSet: cut thresold set to be evluated;  
{

	double *set1=NULL;//two children with clsses
	double *set2=NULL;//two children will be cut
	set1 = new double[nrow];
	set2 = new double[nrow];


	double *set1V=NULL;
	double *set2V=NULL;//two children with value
	set1V = new double[nrow];
	set2V = new double[nrow];

	for(int i=0;i<cutSetN;i++)//evaluate each cut value
	{
		for(int j=0;j<nrow;j++)//initilize each division
		{
			set1[j]=-1000;set2[j]=-1000;//not possible classes numbers
			set1V[j]=-1000;set2V[j]=-1000;//not possible classes numbers
		}
		int count1 = 0; int count2=0;

		double cut = cutSet[i];//for each threold value, divide into two sets, and see the entropy gain
		double distNeibor = 1000000;//the threshold to the minimum neibor, i.e. margin
		//double distNeibor = 1000000;//the threshold to the minimum neibor, i.e. margin

		for(int j=0;j<nrow;j++)
		{
			if(v[j]<=cut)
			{
				set1[count1]=cls[j];
				set1V[count1] = v[j];		
				count1++;

				double tempDistNeibor = fabs(cut-v[j]);//the threshold to the minimum neibor, i.e. margin //reviseFeb16: remove /rangeV
				if(tempDistNeibor<distNeibor)
				{
					distNeibor = tempDistNeibor;
				}
			}
			else
			{
				set2[count2]=cls[j];
				set2V[count2] = v[j];
				count2++;

				double tempDistNeibor = fabs(cut-v[j]);//the threshold to the minimum neibor, i.e. margin
				if(tempDistNeibor<distNeibor)
				{
					distNeibor = tempDistNeibor;
				}
			}
		}
		double *set1True=NULL;double *set2True=NULL;//does not contain void values such as -1000, set1True contains the true size of the divided child
		set1True = new double[count1];set2True = new double[count2];
		//1.set value; 2, calcuate the distance to each center
		
        double centerSet0 = 0; double centerSet1 = 0; double centerSet2= 0;		
        //set0:the spliting node; set1 and set2: children nodes
		for(int j=0;j<count1;j++)
		{
			set1True[j]=set1[j];//1.set value;	
		}
		for(int j=0;j<count2;j++)
		{
			set2True[j]=set2[j];
		}
 		//double iG =-1 ;
        
        distNeibor = distNeibor;
        
        double thisBestEntropy = entropy(cls,nrow) - double(count1)/double(nrow)*entropy(set1True,count1) - double(count2)/double(nrow)*entropy(set2True,count2);
		if(thisBestEntropy>bestEntropy)//thisBestE>bestE
		{
		    bestCut = cut;
            bestE = -10;
            bestEntropy = thisBestEntropy;
            bestMargin = distNeibor;
		}               
        else if(thisBestEntropy==bestEntropy && distNeibor>bestMargin && alpha>0)//prefer larger margin
    	{
		    bestCut = cut;
            bestE = -10;
            bestEntropy = thisBestEntropy;
            bestMargin = distNeibor;
		}
		delete [] set1True;  delete [] set2True;
		set1True = NULL; set2True = NULL; // Clear a to prevent using invalid memory reference.

	}
	delete [] set1V;  delete [] set2V;
	set1V = NULL;    set2V = NULL; // Clear a to prevent using invalid memory reference.
	delete [] set1;  delete [] set2;
	set1 = NULL;    set2 = NULL; // Clear a to prevent using invalid memory reference.


}


void slideWin(double *T, double * Cls, double nrow0, double ncol0,double nMean0, double nSlope0, double nVar0, int seed, int sampleModeWSZ, int sampleModePos,double alpha, int minWin, int maxWin, double* outDouble, double *outInt)//double *x0[],double *y0[]
//void slideWin(double *T, double * Cls, int nrow, int ncol,int nMean, int nSlope)//double *x0[],double *y0[]
{
    int debug_here = 0;
	//sampleModeWSZ : 0:log2; 1:sqrt; 2:all window size;  sampleModePos: 1:sqrt; 2:all window sizes
	int nrow = int(nrow0); int nMean = int(nMean0);int nSlope = int(nSlope0);int nVar = int(nVar0);
	int ncol = int(ncol0);
	//double* outDouble = NULL; outDouble = new double[2];outDouble[0]=0;outDouble[1]=-1000; //0: bestEntropy; 1: best Cut
	//int* outInt = NULL; outInt = new double[4]; outInt[0]=-1;outInt[1]=-1;outInt[2]=-1;outInt[3]=-1;//0:bestWin; 1: bestP1; 2:bstP2; 3: bestType;


//------max value in the time series, and min Vluae
	double minTV = 100000000; double maxTV = -1000000;
		for(int i=0;i<nrow;i++)
		{
			double sum = 0;
			for(int j=0;j<=ncol;j++)
			{
				if(minTV>T[i*ncol+j])
				{
					minTV = T[i*ncol+j];
				}
				if(maxTV<T[i*ncol+j])
				{
					maxTV = T[i*ncol+j];
				}
			}
		}


		int needPosPool = 1;//if needPool = 1, then create a pStartPool and sample from it 
		int needWZPool = 1;//if needPool = 1, then create a winSizePool and sample from it 

double bestE=outDouble[0];double bestCut = outDouble[1]; 
int bestWin = int(outInt[0]); int bestP1 = int(outInt[1]); int bestP2 = int(outInt[2]); int bestType = int(outInt[3]);

double bestEntropy = -1000;
double* bestMarginV = NULL; double* bestEntropyV = NULL;
int* bestP1V = NULL;int* bestP2V = NULL;int* bestWinV = NULL;double* bestCutV = NULL;
double* bestMarginNormV = NULL;

bestMarginV = new double[3]; bestEntropyV = new double[3]; 
bestP1V = new int[3]; bestP2V = new int[3];
bestWinV = new int[3]; bestCutV = new double[3];
bestMarginNormV = new double[3];
for(int iii = 0; iii < 3; iii++)
	{
		bestMarginV[iii]=-1000;		
        bestEntropyV[iii]=-1000;
        bestP1V[iii]=-1000;	
        bestP2V[iii]=-1000;	
        bestWinV[iii]=-1000;	
        bestCutV[iii]=-1000;
        bestMarginNormV[iii]=-1000;
	}
	srand ((seed));//random seed
	//srand(time(0));
//int a = rand() % 5; //( value % 30 + 1985 ) is in the range 1985 to 2014   ( value % 100 + 1 ) is in the range 1 to 100


//generate window sizes from least:ncol without replacement, 
int least = 1;
int winSizeN;

//0:log sample; 1: sqrt window siae;  2: all window size
	if(sampleModeWSZ==0)//log sample
	{
		winSizeN = floor(log(double(ncol))/log(double(2)))+1;//log2 position, plus largest window
		needWZPool = 0;//do not need a sample, because it is just 1,2,3...
	}
	if(sampleModeWSZ==1)//sqrt random sample
	{
		winSizeN = int(ceil(sqrt(double(ncol))))+1;//sqrt + whole window size
	}
	if(sampleModeWSZ==2)//all position
	{
		winSizeN = ncol;//use all postions of number = the window size
		needWZPool = 0;//do not need a sample, because it is just 1,2,3...
	}



int* winSize = NULL;
winSize = new int[winSizeN];//this can assign length of the array dynamiclly

int* winSizePool = NULL;//winSizePool store column number from 1 to ncol
winSizePool = new int[ncol];


	if(needWZPool==1)
	{
		for(int i = 0; i < ncol; i++)
		{
			winSizePool[i]=i+1;		
		}
		int lastIndex = ncol-1; //not include the largest window 
		for(int i = 0; i < winSizeN; i++)
		{
			int randNum=0;//the number should be larger than the least number
			while(randNum<=least-1)//least-1 means the index
			{
				randNum = rand() % (ncol-i);//sample from 1 to lastIndex
			}		
			winSize[i]=winSizePool[randNum];
			//without replacement sampeling        
			std::swap(winSizePool[randNum], winSizePool[ncol-i-1]); 
			//cout<<winSize[i]<<"-";
			//cout<<winSize[i]<<"-"<<'\n';
		}
		//winSize[winSizeN-1]=ncol;//replace the last window size as the whole time series; later mask this in case the separate standardization of two classes
	}

//----for each window sizes, generate start position p1, then the end position p2 = p1 + wsz -1;

	for(int i=0;i<winSizeN;i++)
	{
		int tempWZ;//cout<<tempWZ<<"-";//<<'\n';
	
		if(sampleModeWSZ==0)//log2 and the whole window
		{
			if(i==winSizeN-1)//the last one choose the whole window
			{
				tempWZ=ncol;
			}
			else
			{
				tempWZ = int(pow(double(2),(i+1))); //
			}	
			//cout<<tempWZ<<"-";

		}
		if(sampleModeWSZ==1)
		{
			tempWZ = winSize[i]; //use sqrt of the sizes
		}
		if(sampleModeWSZ==2)
		{
			tempWZ = i+1; //use all sizes

		}
		if(tempWZ>maxWin||tempWZ<minWin)//new: constraint the window size 
          //if(tempWZ>=ncol*10/10)//new: constraint the window size     
		{
			continue;
		}



	//generate the start position of the window
	int pStartN = int(ceil(sqrt(double(ncol-tempWZ+1))));
	 //pStartN = int(((double(ncol-tempWZ+1))));
	int* pStart = NULL;	pStart = new int[pStartN];//
	int* pStartPool = NULL;	pStartPool = new int[ncol-tempWZ+1];


	if(needPosPool==1)//only sqrt sampling need this
	{
		 for(int ii = 0; ii < ncol-tempWZ+1; ii++)
		 {
			pStartPool[ii]=ii;		
		 }
		int lastIndex = ncol-1;
		//cout<<"\n"<<"-----------"<<"\n";
		for(int ii = 0; ii < pStartN; ii++)
		{
			int randNum=0;
			randNum = rand() % (ncol-tempWZ+1-ii);//from 0 -> ncol-tempWZ
			pStart[ii]=pStartPool[randNum];
			//without replacement sampeling        
			std::swap(pStartPool[randNum], pStartPool[ncol-tempWZ-ii]);
		}
	}

	//int positionStartN = ncol-tempWZ+1;//use all postion with the window size
	int positionStartN;//
	
	if(sampleModePos==2)//2: all windo sizes
	{
		positionStartN = ncol-tempWZ+1;//use all postion with the window size
	}
	if(sampleModePos==1)//1: sqrt window sizes
	{
		positionStartN = pStartN;//use sqrt of the positions
	}

	for(int j=0;j<positionStartN;j++)//
	{
		
		int p1;
		if(sampleModePos==2)
		{
			p1=j; //use all postion with the window size
		}
		if(sampleModePos==1&&needPosPool==1)
		{
			p1 = pStart[j]; //use sqrt of the positions
		}
		if(sampleModePos==1&&needPosPool==0)
		{
			p1 = j; //use all position when the window size is big enough
		}
		
		

		int p2=p1+tempWZ-1;

		double *pX = NULL; 
		pX = new double[tempWZ];//a set of x axis from p1, p1+1,...p2; this is for cacluate the slope. Since regression, so it is double
		for(int ii =0;ii<tempWZ;ii++)
		{
			pX[ii]=p1+ii+1;//pX means the true X asis, so need + 1
		}

		
		double* meanV = NULL;//winSizePool store column number from 1 to ncol; meanV is the set of mean of that segment of a time seires. 
		meanV = new double[nrow];


		double* varV = NULL;//covV is the set of covaraince of that segment of a time seires. 
		varV = new double[nrow];
		
		//----mean of each time series, i.e. mean of each row
		


		for(int ii=0;ii<nrow;ii++)
		{
			double sum = 0;
			for(int jj=p1;jj<=p2;jj++)
			{
				sum = sum+ T[ii*ncol+jj];
			}
			double mean = sum/(tempWZ);
			meanV[ii]=mean;

			//calcuate variance
			double tempVar = 0;
			for(int jj=p1;jj<=p2;jj++)
			{
				 tempVar=tempVar+(mean-T[ii*ncol+jj])*(mean-T[ii*ncol+jj]);
				 //tempVar=tempVar+fabs(mean-T[ii*ncol+jj]);
			}
            if(tempWZ>1)
            {
                varV[ii]=sqrt(tempVar/(tempWZ-1));
            }
            if(tempWZ==1)
            {
                varV[ii]=0;
            }

		}
		
        int flagEqCut=1;
        
		//get maximum and minimum mean
		double minM =meanV[0];
		double maxM =meanV[0];
		for(int ii = 0; ii<nrow;ii++)
		{
			if(meanV[ii]<minM)
			{
				minM =meanV[ii];
			}
			if(meanV[ii]>maxM)
			{
				maxM =meanV[ii];
			}
		}
		double eachM = (maxM - minM)*0.9999/(nMean+1);//put here 0.999 here, so that +eachM, can have nMean intervals
		//produce nMean number of values between minM and maxM: meanSet
		
		double* meanSet = NULL;//winSizePool store column number from 1 to ncol
		
        meanSet = new double[nMean];   
		
        for(int ii=0;ii<nMean;ii++)
		{
			//uniform sampling //int a1 = ceil((maxM-minM)*10000);int a2 = ceil(minM*10000);		//int rand1 =  rand() % a1 + a2;
			//double rand2 = double(rand1)/10000;//meanSet[ii]=rand2;
			//grid sampling
		  meanSet[ii]= minM + eachM*(ii+1);
		}	

		double thisBestE=0;double thisBestCut = -10000; double thisBestEntropy = -10000;
        double thisBestMargin = -10000;
		double rangeV = maxM - minM;//this is the difference btween the mean values of this segment
		//rangeV = maxTV - minTV;//difference between max and min value of the whole time series, thus a constant for all window
		if(flagEqCut==1){BestCut(meanV, rangeV, Cls, nrow, meanSet, nMean, thisBestE, thisBestCut, thisBestEntropy,thisBestMargin,alpha);}//v: a set of values; cls: the class corresponding to the v; nrow: row of v; cutSet: cut thresold set to be evluated;  
		//if(flagEqCut==0){BestCut(meanV, rangeV, Cls, nrow, meanV, nrow, thisBestE, thisBestCut, thisBestEntropy,alpha);}
        
        int thisType =0;
		if(thisBestEntropy>bestEntropyV[thisType])//thisBestE>bestE
		{
            bestCutV[thisType] = thisBestCut; 
            bestWinV[thisType] = tempWZ; 
            bestP1V[thisType] = p1; bestP2V[thisType] = p2; 
			//bestType = thisType;
            bestEntropyV[thisType] = thisBestEntropy; 
            bestMarginV[thisType]= thisBestMargin;
            bestMarginNormV[thisType]= thisBestMargin/rangeV;
            //debug_here=1;
		}       
        else if(thisBestEntropy==bestEntropyV[thisType] && thisBestMargin>bestMarginV[thisType] && alpha>0)//prefer larger margin
    	{
            bestCutV[thisType] = thisBestCut; 
            bestWinV[thisType] = tempWZ; 
            bestP1V[thisType] = p1; bestP2V[thisType] = p2; 
			//bestType = thisType;
            bestEntropyV[thisType] = thisBestEntropy; bestMarginV[thisType]= thisBestMargin;
            bestMarginNormV[thisType]= thisBestMargin/rangeV;
		}

		//--------------------------------use variance
		//get maximum and minimum variance
		//if(nVar>0&&tempWZ>=15)// 25 nVar must be greater than 0, and set the minimum window size = 10. so it will be less overfitting
          if(nVar>0)// 25 nVar must be greater than 0, and set the minimum window size = 10. so it will be less overfitting
		{
			double minVar =varV[0];
			double maxVar =varV[0];
			for(int ii = 0; ii<nrow;ii++)
			{
				if(varV[ii]<minVar)
				{
					minVar =varV[ii];
				}
				if(varV[ii]>maxVar)
				{
					maxVar =varV[ii];
				}
			}
			double eachVar = (maxVar - minVar)*0.9999/(nVar+1);//put here 0.999 here, so that +eachM, can have nMean intervals
			//produce nMean number of values between minM and maxM: meanSet
			
			double* varSet = NULL;//winSizePool store column number from 1 to ncol
			varSet = new double[nVar];
			for(int ii=0;ii<nVar;ii++)
			{
				//uniform sampling
				//int a1 = ceil((maxM-minM)*10000);int a2 = ceil(minM*10000);
				//int rand1 =  rand() % a1 + a2;
				//double rand2 = double(rand1)/10000;
				//meanSet[ii]=rand2;
		
				//grid sampling
				varSet[ii]= minVar + eachVar*(ii+1);			
				//cout<<"\n"<<meanSet[ii]<<"-"<<'\n';
			}
			rangeV = (maxVar - minVar);//the upper bound of variance
			//rangeV = (maxTV - minTV)*(maxTV - minTV);//the upper bound of variance. use the maxum value and minivalue of time series 
			thisBestE=-10;thisBestCut = -10000; thisBestEntropy=-10000;thisBestMargin=-1000;
			if(flagEqCut==1){BestCut(varV, rangeV, Cls, nrow, varSet, nVar, thisBestE, thisBestCut,thisBestEntropy,thisBestMargin,alpha);} //v: a set of values; cls: the class corresponding to the v; nrow: row of v; cutSet: cut thresold set to be evluated; 
            
        thisType = 1;
		if(thisBestEntropy>bestEntropyV[thisType])//thisBestE>bestE
		{
            bestCutV[thisType] = thisBestCut; 
            bestWinV[thisType] = tempWZ; 
            bestP1V[thisType] = p1; bestP2V[thisType] = p2; 
			//bestType = thisType;
            bestEntropyV[thisType] = thisBestEntropy; bestMarginV[thisType]= thisBestMargin;
            bestMarginNormV[thisType]= thisBestMargin/rangeV;
		}       
        else if(thisBestEntropy==bestEntropyV[thisType] && thisBestMargin>bestMarginV[thisType] && alpha>0)//prefer larger margin
    	{
            bestCutV[thisType] = thisBestCut; 
            bestWinV[thisType] = tempWZ; 
            bestP1V[thisType] = p1; bestP2V[thisType] = p2; 
			//bestType = thisType;
            bestEntropyV[thisType] = thisBestEntropy; bestMarginV[thisType]= thisBestMargin;
            bestMarginNormV[thisType]= thisBestMargin/rangeV;
		}           
			delete [] varSet; varSet = NULL;
		}
		//----var

int slopeFlag=1;//0:not use slope;
if(slopeFlag==1)//whether use slope
   {
		double* slopeV = NULL;
		slopeV = new double[nrow];
		
		//----slope of each time series, i.e. slope of each row
		for(int ii=0;ii<nrow;ii++)
		{
			//size_t length = sizeof x0 / sizeof x0[0];
			//double m=0;
			//double b=0;
			double * tempY = NULL; tempY = new double[tempWZ];
			for(int iii=0; iii<tempWZ; iii++)
			{
				tempY[iii] = T[ii*ncol+p1+iii];
			}
			double debug1 = (linfit(pX, tempY, tempWZ));
			double debug = atan(linfit(pX, tempY, tempWZ));
			
            if(tempWZ>1){slopeV[ii]=atan(linfit(pX, tempY, tempWZ));}
            if(tempWZ==1){slopeV[ii]=0;}
			//use arttan
			
		//	cout<<slopeV[ii]<<"-";
			delete [] tempY; tempY = NULL;
		}

	
		//get maximum and minimum slope
		double minS =slopeV[0];
		double maxS =slopeV[0];
		for(int ii = 0; ii<nrow;ii++)
		{
			if(slopeV[ii]<minS)
			{
				minS =slopeV[ii];
			}
			if(slopeV[ii]>maxS)
			{
				maxS =slopeV[ii];
			}
		}


		double eachS = (maxS - minS)*0.9999/(nSlope+1);//put here 0.999 here, so that +eachS, can have nSlope intervals

		if(maxS-minS>0.00000001)//if the max and min dif is too small, then don't use this feature
		{
			//produce nSlope number of values between minS and maxS: meanSet
			double* slopeSet = NULL;//
			slopeSet = new double[nSlope];
			for(int ii=0;ii<nSlope;ii++)
			{
				//uniform sampling
				//int a1 = ceil((maxS-minS)*10000);
				//int a2 = ceil(minS*10000);
				//int rand1 =  rand() % a1 + a2;
				//double rand2 = double(rand1)/10000;
				//slopeSet[ii]=rand2;
				
				//grid sampling
				slopeSet[ii]= minS + eachS*(ii+1);	
				//cout<<"\n"<<"----------------------"<<"\n";
				//cout<<"\n"<<slopeSet[ii]<<"-"<<'\n';
			}	

							//debug ---------
		//cout<<"--------------------------------------"<<"\n";
		//for(int ii=0;ii<nSlope;ii++)
		//{
		//		cout<<slopeSet[ii]<<"\n";
		//}

			thisBestE=-10;thisBestCut = -10000; thisBestEntropy = -10000;thisBestMargin=-1000;
			rangeV = maxS - minS;//the the difference of max and min of arctan(slope) for this segment
			//rangeV = PI;//the max: difference of the arctan(sople). Because the arctan = [-pi/2,pi/2]
            
            
			if(flagEqCut==1){BestCut(slopeV, rangeV, Cls, nrow, slopeSet, nSlope, thisBestE, thisBestCut,thisBestEntropy,thisBestMargin,alpha);} //v: a set of values; cls: the class corresponding to the v; nrow: row of v; cutSet: cut thresold set to be evluated;  

        thisType = 2;
		if(thisBestEntropy>bestEntropyV[thisType])//thisBestE>bestE
		{
            bestCutV[thisType] = thisBestCut; 
            bestWinV[thisType] = tempWZ; 
            bestP1V[thisType] = p1; bestP2V[thisType] = p2; 
			//bestType = thisType;
            bestEntropyV[thisType] = thisBestEntropy; bestMarginV[thisType]= thisBestMargin;
            bestMarginNormV[thisType]= thisBestMargin/rangeV;
		}       
        else if(thisBestEntropy==bestEntropyV[thisType] && thisBestMargin>bestMarginV[thisType] && alpha>0)//prefer larger margin
    	{
            bestCutV[thisType] = thisBestCut; 
            bestWinV[thisType] = tempWZ; 
            bestP1V[thisType] = p1; bestP2V[thisType] = p2; 
			//bestType = thisType;
            bestEntropyV[thisType] = thisBestEntropy; bestMarginV[thisType]= thisBestMargin;
            bestMarginNormV[thisType]= thisBestMargin/rangeV;
		}         

			delete [] slopeSet; slopeSet = NULL;
		}
  
			delete [] slopeV;slopeV = NULL; 
  }  
//double *v, double *cls, int nrow, double *cutSet, int cutSetN, double &bestCut, double &bestE



		delete [] pX; pX = NULL;		
		delete [] meanSet; meanSet = NULL;		
		delete [] meanV;meanV = NULL; 
		delete [] varV;varV = NULL; 
	}
		delete [] pStart;pStart = NULL; 
		delete [] pStartPool;pStartPool = NULL; 

	}

bestType=rand() % 3;         // v1 in the range 0 to 2; 
bestE = bestEntropyV[bestType];bestCut=bestCutV[bestType];
bestEntropy=bestEntropyV[bestType];bestWin=bestWinV[bestType]; 
bestP1=bestP1V[bestType]; bestP2=bestP2V[bestType]; 
double bestMarginNorm = bestMarginNormV[bestType];

for(int typeI=0; typeI<3; typeI++)
{
    if(bestEntropyV[typeI]>bestEntropy)//thisBestE>bestE
		{
            bestType=typeI;bestE = bestEntropyV[bestType];bestCut=bestCutV[bestType];
            bestEntropy=bestEntropyV[bestType];bestWin=bestWinV[bestType]; 
            bestP1=bestP1V[bestType]; bestP2=bestP2V[bestType]; 
            bestMarginNorm = bestMarginNormV[bestType];
            //debug_here=2;
		}
//       else if( bestEntropyV[typeI]==bestEntropy && bestWinV[typeI] < bestWin)        
//       //else if( bestEntropyV[typeI]==bestEntropy && bestMarginNormV[typeI] > bestMarginNorm)      
//       {
//               bestType=typeI;bestE = bestEntropyV[bestType];bestCut=bestCutV[bestType];
//               bestEntropy=bestEntropyV[bestType];bestWin=bestWinV[bestType]; 
//               bestP1=bestP1V[bestType]; bestP2=bestP2V[bestType];    
//               bestMarginNorm = bestMarginNormV[bestType];
//               //debug_here=1;
//       }    
}

outDouble[0]=bestE;  outDouble[1]=bestCut; outDouble[2]=bestEntropy;
outDouble[3]=bestMarginNormV[0];outDouble[4]=bestMarginNormV[1];outDouble[5]=bestMarginNormV[2];
outDouble[6]=bestEntropyV[0];outDouble[7]=bestEntropyV[1];outDouble[8]=bestEntropyV[2];
outDouble[9]=debug_here;outDouble[10]=bestMarginNorm;

outInt[0] = bestWin; outInt[1]=bestP1+1; outInt[2]=bestP2+1; outInt[3]= bestType;

delete [] bestMarginV;bestMarginV=NULL;		
delete [] bestMarginNormV;bestMarginNormV=NULL;	
delete [] bestEntropyV;bestEntropyV=NULL;	
delete [] bestP1V;	bestP1V=NULL;	
delete [] bestP2V;	bestP2V=NULL;	
delete [] bestWinV;	bestWinV=NULL;	
delete [] bestCutV; bestCutV=NULL;	

delete [] winSizePool;  // When done, free memory pointed to by windowSisePool
winSizePool = NULL;     // Clear a to prevent using invalid memory reference.
delete [] winSize;  // When done, free memory pointed to by a.
winSize = NULL;     // Clear a to prevent using invalid memory reference.
}

//the above is copied to Matlab to be mexed
//----above is from c++

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    
    //input
    double* T=            mxGetPr(prhs[0]);
    double* Cls=                mxGetPr(prhs[1]);
    double nrow=                mxGetScalar(prhs[2]);
    double ncol=        mxGetScalar(prhs[3]);
    double nMean=       mxGetScalar(prhs[4]);   
    double nSlope=            mxGetScalar(prhs[5]);
    double nVar=            mxGetScalar(prhs[6]);
    int seed=            mxGetScalar(prhs[7]);
    int sampleModeWSZ=            mxGetScalar(prhs[8]);
    int sampleModePos=            mxGetScalar(prhs[9]);
    double alpha=            mxGetScalar(prhs[10]);
    int minWin=            mxGetScalar(prhs[11]);
    int maxWin=            mxGetScalar(prhs[12]);
	
    
    //output
    double* outDouble;	
    double* outInt;	
    
    /* set up output arguments */
	plhs[0] = mxCreateDoubleMatrix(11,1,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(4,1,mxREAL);
	
    outDouble = mxGetPr(plhs[0]); 
    outInt = mxGetPr(plhs[1]);    

    /* call the computational routine */
    slideWin(T,Cls,nrow,ncol,nMean,nSlope,nVar,seed,sampleModeWSZ,sampleModePos,alpha,minWin,maxWin,outDouble,outInt);
    
    //clear;
}
