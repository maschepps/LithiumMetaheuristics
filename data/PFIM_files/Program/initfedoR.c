/*	Algorithme de Fedorov-Wynn Version Sylvie

Appel dans R
.C("Fedorov_R",as.char(namfich))
	namfich est un vecteur de deux noms de fichiers:
		namfile est un fichier contenant les protocoles elementaires
		namprob est une chaine contenant le nom du probleme
	ex : namfich<-c("matelem.dat","problem1")

Resultat
	ecrit un fichier resultat contenant les Mf des protocoles elementaires
		=> <namprob>.out
	et un fichier
		=> <namprob>.res
*/

#include <stdio.h> 
#include <stdlib.h>
#include <assert.h>
#include <math.h> 	/* Mathematical functions */
#include <time.h>	/* Function time used to initialise the random number generator */
#include <float.h>	/* Implementation related constants */
#include <signal.h>	/* Signal handling used to detect arithmetic errors */

#define FTOL 1e-6 /* Minimal change in criterion during optimisation step*/
#define FREQMIN 1e-8 /* Minimal frequency for a protocol to be kept */
#define EPSILON DBL_EPSILON
#define CRI_MAX 1e30 /* maximum value for criterium, used to test correct comp*/
#define NRANSI 
#define tol 0.0001
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;} 
#define TINY 1.0e-20
#define SMALL 1.0e-7
#define IA 16807  	/*definition des constantes pour ran1 */
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define NTAB 32 
#define NDIV (1+(IM-1)/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS) 
#define FREE_ARG char* 

/******************************************************************************
   function prototyping
*/
#include "classfed.h" /* Structures defined for the program*/


/******************************************************************************
   Fonctions calculatoires
*/
double matrace(int iprot,PROTOC *prot,POPROT *pop);
double lik(POPROT *pop);
void lubksb(double **a, int n, int *indx, double b[]);
int ludcmp(double **a, int n, int *indx, double *d);

/* Random number generator */
double randu(double x);
double gauss1(long *idum);
float ran1(long *idum);

/******************************************************************************
   Defining global variables (see variables_npml.txt for details)
*/
volatile sig_atomic_t ierr=0;	/* Parameter signaling a SIGFPE error */
				/* to be used with a signal handler */
long seed=-10;     	/* graine du generateur aleatoire */

/******************************************************************************
   lecture et ecriture de fichiers
*/
PROTOC *readfich_R(char *callfile,char *outfile,int nprot,
		int ndim,int ntot)
{	
	FILE *incall=fopen(callfile,"rt");
	FILE *inout=fopen(outfile,"w"); /* Empties the result file */
	int iread;
	int iprot,jtps,ij,ik,icas;
	double xread;
	PROTOC *hh;

/* reading the Splus file with the elementary protocols
Structure:
-1st line: nb_protocols nb_times nb_dimensions total_cost
-then for each protocol
	protocol number
	sampling times
	Fisher information matrix (symmetric)
*/
/*	fprintf( stderr,"Reading individual protocols from %s\n", callfile );
*/	if( incall==NULL)
	{
	fprintf( stderr,"File [%s] does not exist...\n", callfile );
	return NULL;
	}
	else
	{
/*		fprintf(stderr,"Opening %s for results \n", outfile );
*/		fprintf(inout,"RESULTS FROM FILE %s\n", outfile );
		hh = (PROTOC *)calloc(nprot, sizeof(PROTOC));	
		for (iprot=0;iprot<nprot;iprot++) {
/*		fprintf(stderr,"Protocole %d \n",iprot);*/
			fscanf(incall,"%d",&iread);
			fscanf(incall,"%d",&iread);
			PROTOC_alloc(&hh[iprot],iread,ndim);
			for (jtps=0;jtps<hh[iprot].ntps;jtps++) {
				fscanf(incall,"%lf",&xread);
				hh[iprot].tps[jtps]=xread;}
			icas=0;
			for (ij=0;ij<ndim;ij++) {
				for (ik=0;ik<=ij;ik++) {
					fscanf(incall,"%lf",&xread);
/* eco : cost comes in here Ntot/ntps
*/
	hh[iprot].fisher[icas]=xread*(double)(ntot)/(double)(hh[iprot].ntps);
					icas++;}
				for (ik=(ij+1);ik<ndim;ik++) {
					fscanf(incall,"%lf",&xread);}
					}
				}
			}
	fclose(incall);
	fclose(inout);
	return hh;
}

void ecritprot(char *outfile,POPROT *pop,int nit)
{	FILE *inout=fopen(outfile,"a"); /* Appends to file */
	int i,j;
	double xcal;

	if (inout == NULL)
	{
	fprintf(stderr,"Impossible d'ouvrir le fichier %s en lecture/écriture\n",outfile);
	return;
	}
	if(nit==(-1))
		fprintf(inout,"INITIALISATION\n");
	else
		fprintf(inout,"ITERATION  %d    \n",nit);
	fprintf(inout,"NOMBRE DE GROUPES  %d    ",pop->np);
	xcal=lik(pop);
	xcal=exp(log(pop->det)/(double)pop->ndim);
	fprintf(inout,"DET %e EFF %e\n",pop->det,xcal);
/* Ici calcul du cout dans programme initial
	fprintf(inout,"nombre total de points  %d    ",880);
	fprintf(inout,"nombre de sujets %lf \n",(double)880/4);
  ntot=800, donne comme cout maximal
  xcout=1
  cout(iprot)=ntps(iprot)/float(ntot)*xcout
  xnn=Somme_i xcout*al(i)/cout(q(i))
  ga(i) remplace ici par 2 fois frequence
*/
	for(i=0;i<pop->np;i++) 
	{
		fprintf(inout,"%d %lf %lf %d %d ",i,pop->freq[i],pop->freq[i],
			pop->num[i]+1,pop->pind[i].ntps);
		for (j=0;j<pop->pind[i].ntps;j++) 
		{
			fprintf(inout,"%lf ",pop->pind[i].tps[j]);
		}
		fprintf(inout,"\n ");
	}
	fclose(inout);
	return;
}

/******************************************************************************
   Initialisation et ajout de protocole élémentaire
*/
POPROT *initprot(PROTOC *allprot,int *protdep,double *freqdep)
{
	int i,ij,np,nmax;
	float xrand;
	POPROT *mypop;
	
/*	xrand=ran1(&seed);
*/	np=protdep[0];
/*	fprintf(stderr,"Initialisation avec les protocoles numéro ");
	for(i=1;i<np+1;i++)
		fprintf(stderr,"%d ",protdep[i]);
	fprintf(stderr,"\n");
*/
/*	if (np>allprot[0].ndim) 
		np=allprot[0].ndim;
*/	nmax=allprot[0].ndim*(allprot[0].ndim+1)/2+1;
	mypop=POPROT_alloc(np,allprot[0].ndim,nmax);
	for(i=0;i<np;i++) {
		mypop->freq[i]=freqdep[i];
		mypop->num[i]=protdep[i+1]-1;
		ij=mypop->num[i];
		PROTOC_alloc(&mypop->pind[i],allprot[ij].ntps,allprot[ij].ndim);
		PROTOC_copy(&allprot[ij],&mypop->pind[i]); /* Avec PROTOC_copy*/
		}
/*	fprintf(stderr,"\n Initialisation :\n");
	POPROT_print(mypop,0);
*/	return mypop;
}

/******************************************************************************
   Optimisation de protocole
*/

void tassement(POPROT *pop)
{	int i,idec=0,np;
	np=pop->np;
	for(i=0;i<pop->np;i++)
	{
		if(pop->freq[i]<FREQMIN)
		{
			idec++;
			np--;
		} 
		else
		{
			if (idec>0)
			{
				pop->freq[i-idec]=pop->freq[i];
				pop->num[i-idec]=pop->num[i];
/* Next: free memory allocated to pind[i-idec] and allocate memory for a
protocol the size of pind[i]. Copy pind[i] to pind[i-1].
This is not necessary if the protocols have the same number of times,
but will be required if protocols are allowed to have different numbers of
sampling times */
				PROTOC_destroy(&pop->pind[i-idec]);
	PROTOC_alloc(&pop->pind[i-idec],pop->pind[i].ntps,pop->pind[i].ndim);
				PROTOC_copy(&pop->pind[i],&pop->pind[i-idec]); 
				
			}
		}
	}
/* Free the remaining protocols, that have all been copied*/
	for(i=np;i<pop->np;i++) PROTOC_destroy(&pop->pind[i]);
	pop->np=np;
/*	fprintf(stderr,"\n Dans tassement :\n");
	POPROT_print(pop,0);
*/	return;
}

int ajout(PROTOC *allprot,POPROT *pop,int nprot)
{	double cri,xmul,xtest,ifin=0,xdim,rcalc;
	int iprot,j,deja,qajout=-1;

	cri=pop->det;
	xtest=(double)(pop->ndim);
	xdim=xtest;
	for(iprot=0;iprot<nprot;iprot++)
	{
		deja=0;
		for(j=0;j<pop->np;j++)
		{
			if (iprot==pop->num[j]) deja=1;
		}
		if(fabs(deja)<EPSILON)
		{
			xmul=matrace(iprot,allprot,pop);
			if(xmul>=xtest) 
			{
				xtest=xmul;
				qajout=iprot;
			}
		}
	}
	if(fabs(xtest-xdim)<EPSILON)
	{
		ifin=1;
	} 
	else
	{
		rcalc=(xtest-xdim)/(xdim*(xtest-1));
/*		printf("rcalc=%lf\n",rcalc);
*/
		PROTOC_alloc(&pop->pind[pop->np],allprot[qajout].ntps,allprot[qajout].ndim);
		PROTOC_copy(&allprot[qajout],&pop->pind[pop->np]);
		pop->num[pop->np]=qajout;
		pop->freq[pop->np]=rcalc;
		for(j=0;j<pop->np;j++) pop->freq[j]=pop->freq[j]*(1-rcalc);
		pop->np++;
	}
/*	fprintf(stderr,"\n Dans ajout :\n");
	POPROT_print(pop,0);
*/	cri=lik(pop);
	return ifin;
}

int takeout_k(POPROT *pop,PROTOC *allprot,int *nnul)
{
	int ir=-1,j;
	double xmul;
	xmul=lik(pop); /*update Mf-1 */
	for(j=0;j<pop->np;j++)
	{
		xmul=0;
		if(pop->freq[j]<=FREQMIN)
		{
			xmul=matrace(pop->num[j],allprot,pop);
			if(xmul>(double)pop->ndim)
			{
				ir=j;
				*nnul=*nnul-1;
				return ir;
			}
		}
	}	
	return ir;
}

int project_grad(POPROT *pop, PROTOC *allprot, double *gal,int ir,double
		*vmgal,int nnul)
{
/* Computes the projection of the gradient of ln(V)=ln(det(H(pop))) on the
surface sum(frequencies)=1 and freq[j]=0 for all j such as pop->freq[j]<FREQMIN
gal[j] is the jth component of the projected gradient (the direction)
gal[j]=0 if pop->freq[j]<FREQMIN and j<>ir
ir is -1 or the index removed from the active set at the previous step
vmgal is the largest abs(gal[j])
returns in, the last index j such as gal[j]<0
*/
	int j,in=-1;
	double sg,sgal;
		
	sg=lik(pop); /* Check that Mf-1 has been computed */
	sg=0;
	for(j=0;j<pop->np;j++)
	{
		gal[j]=0;
		if((pop->freq[j]>=FREQMIN) | (j==ir))
			gal[j]=matrace(pop->num[j],allprot,pop);
		sg=sg+gal[j];
	}
	sgal=0;*vmgal=0;
	for(j=0;j<pop->np;j++)
	{
		if((pop->freq[j]>=FREQMIN) | (j==ir))
		{
			gal[j]=gal[j]-sg/(double)(pop->np-nnul);
			if(fabs(gal[j])>*vmgal) *vmgal=(double)fabs(gal[j]);
			if(gal[j]<0) in=j;
			sgal=sgal+gal[j];
		}
	}
	return in;
}

int optim_lambda(POPROT *pop, double *gal, double *alkeep, double *crit,
	double dmax, double crit2, double vmgal)
{
/* dmax is too large a move along the direction
=> move along the direction by dmax/2**k until
-the move is too small : no new frequencies have been found (return 1)
-a better likelihood is found => update pop with the new frequencies
	-if the criteria has changed by more than FTOL, return 0
	-else return 1; the frequencies correspond to the last (smallest) move
*/
	int j;
	double ro=1,crit1;
	do
	{
		ro=ro*0.5;
		for(j=0;j<pop->np;j++) pop->freq[j]=alkeep[j]+ro*dmax*gal[j];
		crit1=lik(pop); /* Updates Mf-1 for new frequencies */
		if(crit1>*crit)
		{
			for(j=0;j<pop->np;j++) alkeep[j]=pop->freq[j];
			*crit=crit1;
			if(fabs((crit2-*crit)/crit2)>FTOL) {return 0;} 
			else {return 1;}
		}
	} while((dmax*ro*vmgal)>1e-4);
	for(j=0;j<pop->np;j++) pop->freq[j]=alkeep[j];
	return 1;
}

int doptimal(POPROT *pop,PROTOC *allprot)
{	
	int nnul=0, ir=-1,in=-1,nit=0;
	int j,imax,iopt;
	double crit,crit1,crit2,vmgal,dmax;
	double *gal,*alkeep;
	double sal,sal1;
	
	gal=(double *)calloc(pop->np,sizeof(double));
	alkeep=(double *)calloc(pop->np,sizeof(double));
	crit=lik(pop); /*Calcul de Mf-1 */

	for(nit=0;nit<500;nit++)
	{
/* alkeep stores the initial frequencies of the elementary protocols
*/
		for(j=0;j<pop->np;j++)
			alkeep[j]=pop->freq[j]; 
		crit2=crit;
		if(ir==(-1)) /* First time through */
		{
			nnul=0; /* Computing the number of nul frequencies */
			for(j=0;j<pop->np;j++)
			{
				if(pop->freq[j]<FREQMIN) nnul++;
			}
		}
/* Compute the projected gradient (direction)*/
		in=project_grad(pop,allprot,gal,ir,&vmgal,nnul);
		if((in<0) | (gal[in]==0)) 
		{	
			free(gal);
			free(alkeep);
			return nnul;
		}
/* Compute the optimal change along the gradient dmax (lambda)*/
		dmax=-pop->freq[in]/gal[in];
		imax=in;
		if((ir!=(-1)) & (gal[ir]<=0)) 				
		{	
			free(gal);
			free(alkeep);
			return nnul;
		}
		ir=-1;
		for(j=0;j<pop->np;j++)
		{
			if((gal[j]<0) & (fabs(pop->freq[j]/gal[j])<dmax))
			{
				dmax=-pop->freq[j]/gal[j];
				imax=j;
			}
		}
/* Updates the frequencies as beta*(t,j)=beta(t,j)+dmax*gal(j)
Compute the likelihood for the updated vector of frequencies
*/
		sal=0;
		for(j=0;j<pop->np;j++)
		{
			pop->freq[j]=pop->freq[j]+dmax*gal[j];
			sal=sal+pop->freq[j];
		}
		sal1=sal-pop->freq[imax];
		pop->freq[imax]=0;
		crit1=lik(pop);
/* If better likelihood, update alkeep (pop->freq has been updated already)
*/
		if(crit1>crit)
		{
			for(j=0;j<pop->np;j++) alkeep[j]=pop->freq[j];
			crit=crit1;
			nnul++;
		}
		else
		{
			for(j=0;j<pop->np;j++) pop->freq[j]=alkeep[j];
/* If not, try dmax*ro for ro=1/2**k until dmax*ro too small
*/		
			iopt=optim_lambda(pop,gal,alkeep,&crit,dmax,crit2,vmgal);
			if(iopt==1)
			{
				ir=takeout_k(pop,allprot,&nnul);
/* ir=-1 if for all j, matrace(q[j]) <= ndim*(ndim+1)/2, ie optimal protocol
=> in that case exit the subroutine
*/		
				if(ir==(-1)) 
				{	
					free(gal);
					free(alkeep);
					return nnul;
				}
			}
		}
	}
	free(gal);
	free(alkeep);
	fprintf(stderr,"Pb with convergence in doptimal \n");
	return -10;
}

void index_limit(POPROT *pop)
{	
	fprintf(stderr,"The number of elementary protocols has reached");
	fprintf(stderr,"the limit.\n");
	fprintf(stderr,"Please program the index limitation step...\n");
	return;
}

int main_opt(PROTOC *allprot,POPROT *pop,char *nomout,int nprot)
{	int nit=0,ifin=0,nstop,ntest;
	double cri;

	cri=lik(pop);
/*	fprintf(stderr,"Critère au point initial = %lf \n",cri);
*/	ecritprot(nomout,pop,(-1));
	for(nit=0;nit<100;nit++)
	{
/* A priori, il n'y a aucune chance d'avoir un alpha trop petit a ce stade donc
je mets l'appel a tassement a la fin
*/		nstop=0;
		do 
		{
			ifin=ajout(allprot,pop,nprot);
			nstop++;
		} while ((cri<(-1e10))& (nstop<100));
		if((nstop>99)& (cri<(-1e10))) 
		{
			fprintf(stderr,"Unable to find a new protocol to add \n");
			return -10;
		}
		if(ifin==1) break; /* Si ifin=1, on sort de la boucle sur nit */
		ntest=doptimal(pop,allprot);
		if(ntest<0)
			return -100;
		tassement(pop);
		ecritprot(nomout,pop,nit);
		if(pop->np>(pop->ndim*(pop->ndim+1)/2)) 
		{
			index_limit(pop);
			return -1;
		}
		cri=lik(pop); /* Updates matrices */
	}
	return nit;
}

/******************************************************************************
   Fonctions calculatoires
	matrace:
	lik:determinant de hh et inverse de hh
*/

double matrace(int iprot,PROTOC *prot,POPROT *pop)
{/* Computes tr(Mf(iprot)*Mf-1(poprot)) 
WARNING pop->finv must be computed beforehand (with a call to lik)*/
	int i,j,ij;
	double xcal=0;
	
	for(i=0;i<(prot->ndim);i++) 
	{
		for(j=0;j<i;j++) 
		{
/* Nb d'éléments jusqu'à la ligne i non inclue: i(i+1)/2
   Nb d'éléments jusqu'à la ligne i inclue: (i+1)(i+2)/2
   i ieme élément de la ième ligne 
*/
			ij=i*(i+1)/2+j;
			xcal=xcal+2*(prot[iprot].fisher[ij])*(pop->finv[ij]);
			}
		ij=i*(i+3)/2;
		xcal=xcal+(prot[iprot].fisher[ij])*(pop->finv[ij]);
		}
	return xcal;
}

double lik(POPROT *pop)
{/* Function to compute the inverse of the Fisher information matrix of the
population protocol pop. Returns the determinant of the matrix, also stored in
pop->det. The Fisher  information matrix of the population protocol is stored
in pop->fisher and the inverse is stored in pop->finv
*/
	int ndim=pop->ndim;
	int ncase=(int)(ndim*(ndim+1)/2);
	int i,j,jj,ifail;
	int *indx;
	double *col;
	matrix *xa;
	double cri;
	
/*	POPROT_print(pop,1);
	PROTOC_print(&prot[7]);
*/
	indx=(int *)calloc(ndim,sizeof(int));
	col=(double *)calloc(ndim,sizeof(double));
	xa=matrix_create(ndim,ndim);
	for(i=0;i<ncase;i++) 
	{
		pop->fisher[i]=0;
		for(j=0;j<pop->np;j++) 
		{
	pop->fisher[i]=pop->fisher[i]+pop->freq[j]*pop->pind[j].fisher[i];
		}
	}
	jj=0;
	for(i=0;i<ndim;i++) 
	{
		for(j=0;j<=i;j++) 
		{
			xa->m[i][j]=pop->fisher[jj];
			xa->m[j][i]=pop->fisher[jj];
			jj++;
			}
		}
	ifail=ludcmp(xa->m,ndim,indx,&cri);
	if(ifail==1) return CRI_MAX;
	for(i=0;i<ndim;i++) cri=cri*xa->m[i][i];
	pop->det=cri;
	for(i=0;i<ndim;i++) 
	{
		for(j=0;j<ndim;j++) col[j]=0.0;
		col[i]=1.0;
		lubksb(xa->m,ndim,indx,col);
		for(j=0;j<ndim;j++) 
		{
			jj=i*(i+1)/2+j;
			pop->finv[jj]=col[j]; /*=xb[j][i] using another matrix*/
			}
		}
	/* desallocation */
	matrix_destroy(xa);
	free(indx);
	return cri;
}

void lubksb(double **a, int n, int *indx, double b[]) 
{ 
	int i,ii=-1,ip,j; 
	double sum; 
 
	for (i=0;i<n;i++) { 
		ip=indx[i]; 
		sum=b[ip]; 
		b[ip]=b[i]; 
		if (ii>=0) 
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j]; 
		else if (sum) ii=i; 
		b[i]=sum; 
	} 
	for (i=n-1;i>=0;i--) { 
		sum=b[i]; 
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j]; 
		b[i]=sum/a[i][i]; 
	} 
} 

int ludcmp(double **a, int n,int *indx, double *d) 
{ 
	int i,imax,j,k; 
	double big,dum,sum,temp; 
	double *vv;
	
	vv=(double *)calloc(n,sizeof(double)); 
	*d=1.0; 
	for (i=0;i<n;i++) { 
		big=0.0; 
		for (j=0;j<n;j++) 
			if ((temp=fabs(a[i][j])) > big) big=temp; 
		if (big == 0.0) {
			fprintf(stderr,"Singular matrix in routine ludcmp \n"); 
			return 1;}
		vv[i]=1.0/big; 
	} 
	for (j=0;j<n;j++) { 
		for (i=0;i<j;i++) { 
			sum=a[i][j]; 
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j]; 
			a[i][j]=sum; 
		} 
		big=0.0; 
		for (i=j;i<n;i++) { 
			sum=a[i][j]; 
			for (k=0;k<j;k++) 
				sum -= a[i][k]*a[k][j]; 
			a[i][j]=sum; 
			if ( (dum=vv[i]*fabs(sum)) >= big) { 
				big=dum; 
				imax=i; 
			} 
		} 
		if (j != imax) { 
			for (k=0;k<n;k++) { 
				dum=a[imax][k]; 
				a[imax][k]=a[j][k]; 
				a[j][k]=dum; 
			} 
			*d = -(*d); 
			vv[imax]=vv[j]; 
		} 
		indx[j]=imax; 
		if (a[j][j] == 0.0) a[j][j]=TINY; 
		if (j != (n-1)) { 
			dum=1.0/(a[j][j]); 
			for (i=j+1;i<n;i++) a[i][j] *= dum; 
		} 
	}
	free(vv);
	return 0;
} 

/* Random number generator */
double randu(double x)
{	double xrand;
	xrand=rand();
	return xrand;
}
double gauss1(long *idum)
{
/* 	float ran1(long *idum); 
*/	static int iset=0; 
	static float gset; 
	float fac,rsq,v1,v2; 
 
	if  (iset == 0) { 
		do { 
			v1=2.0*ran1(idum)-1.0; 
			v2=2.0*ran1(idum)-1.0; 
			rsq=v1*v1+v2*v2; 
		} while (rsq >= 1.0 || rsq == 0.0); 
		fac=sqrt(-2.0*log(rsq)/rsq); 
		gset=v1*fac; 
		iset=1; 
		return v2*fac; 
	} else { 
		iset=0; 
/*		printf("Valeur de idum %ld \n",*idum);
*/		return (double)gset; 
	} 
} 
 
float ran1(long *idum) 
{ 
	int j; 
	long k; 
	static long iy=0; 
	static long iv[NTAB]; 
	float temp; 
 
	if (*idum <= 0 || !iy) { 
		if (-(*idum) < 1) *idum=1; 
		else *idum = -(*idum); 
		for (j=NTAB+7;j>=0;j--) { 
			k=(*idum)/IQ; 
			*idum=IA*(*idum-k*IQ)-IR*k; 
			if (*idum < 0) *idum += IM; 
			if (j < NTAB) iv[j] = *idum; 
		} 
		iy=iv[0]; 
	} 
	k=(*idum)/IQ; 
	*idum=IA*(*idum-k*IQ)-IR*k; 
	if (*idum < 0) *idum += IM; 
	j=iy/NDIV; 
	iy=iv[j]; 
	iv[j] = *idum; 
	if ((temp=AM*iy) > RNMX) return RNMX; 
	else return temp; 
} 

/******************************************************************************
   Main program
*/

void FedorovInit_R(char **namfich,int *ndimen,int *nbprot,int *numprot,
	double *freq,int *nbdata, double *vectps, double *fisher,int *nok, 
		int *protdep,double *freqdep)
/* Ajout d'une initialisation par l'utilisateur du protocole de population
* définis dans R et envoyé au programme C
        namfich fichiers avec les noms des 2 fichiers (nomini=matrices de Fisher
	        nomout=nom du fichier de résultats)
	ndimen 	dimensions du problème (nb de prot elementaires, dimension de la
		matrice de Fisher, cout total)
* renvoyés par le programme C
	nbprot	nb de protocoles dans le protocole de population final P
	numprot	liste des protocoles elementaires dans P
	freq	liste des fréquences des protocoles élementaires
	nbdata	liste du nb de temps dans chaque protocole élementaires
	vectps	vecteur des temps 
	fisher	matrice d'information de Fisher du protocole de population
	noik	nb of iterations if optimisation success, pb encountered if not
* défini dans R et envoyé au programme C, pour l'initialisation du protocole de
population
	protdep	premier chiffre=nb de protocoles de départ, chiffres
		suivants=numéros des protocoles élémentaires
* défini dans R, envoyé à C, modifié par C et renvoyé à R
	freqdep	fréquence des protocoles élémentaires constituant le protocole
		initial
		renvoyé : fréquences pondérées (par les coûts)
*/
{	char nomini[120],nomout[120];
	int nprot,ndim,Ntot;
	int i,j,ij;
	double sumf;
	PROTOC *allprot; /*Vector of individual protocols */
	POPROT *pop; /* Population protocol */

	*nok=0;
	sprintf(nomini,"%s\0",namfich[0]);
	sprintf(nomout,"%s\0",namfich[1]);
	nprot=ndimen[0];
	ndim=ndimen[1];
	Ntot=ndimen[2];
/*	fprintf(stderr,"Lecture des protocoles élémentaires \n");
	fprintf(stderr,"%d %d %d \n",nprot,ndim,Ntot);
*/	allprot=readfich_R(nomini,nomout,nprot,ndim,Ntot);
	if(allprot==NULL) {
		*nok=-500;
		fprintf(stderr,"Check file names.\n");		
	}
/*	PROTOC_print(&allprot[285]); 
	fprintf(stderr,"Initialisation du protocole de population \n");
	fprintf(stderr,"Taille d'un protocole individuel=%d ",sizeof(PROTOC));
	fprintf(stderr,"Taille d'un protocole de popoulation=%d ",sizeof(POPROT));
*/	if(*nok>=0) {
		pop=initprot(allprot,protdep,freqdep);
		assert(pop!=NULL);
/*	fprintf(stderr,"Initialisation à %d protocoles \n",pop->np);
*/		*nok=main_opt(allprot,pop,nomout,nprot);
		if(*nok==(-1))
		fprintf(stderr,"Number of protocols in population protocol too large\n");
		if(*nok==(-10))
		fprintf(stderr,"Could not find a new protocol to add\n");
		if(*nok==(-100))
		fprintf(stderr,"Convergence problem in the frequency optimisation step\n");	
/*	fprintf(stderr,"-------------------------------------------- \n");
	fprintf(stderr,"Results of the protocol optimisation\n");
	fprintf(stderr,"-------------------------------------------- \n");
	POPROT_print(pop,1);
	fprintf(stderr,"-------------------------------------------- \n");
*/
/* Dimensionnement des tableaux de sortie doit etre effectué avant l'appel!!!
	numprot=(int *)calloc(pop->np,sizeof(int));
*/		
		if(*nok>=0) {
			*nbprot=pop->np;
			sumf=0;
			for(i=0;i<pop->np;i++) {
				numprot[i]=pop->num[i]+1;
				nbdata[i]=pop->pind[i].ntps;
/* Calcul de la fréquence à partir de la fréquence pondérée par le coût
  calcul de la somme des fréquences pour normaliser Somme_i freq[i]=1
*/
		freq[i]=pop->freq[i]*(double)(Ntot)/(double)(nbdata[i]);
				sumf=sumf+freq[i];
				freqdep[i]=pop->freq[i];
/*				ij=ij+pop->pind[i].ntps;
*/			}
			ij=0;
			for(i=0;i<pop->np;i++) {
				freq[i]=freq[i]/sumf;
				for(j=0;j<pop->pind[i].ntps;j++) {
					vectps[ij]=pop->pind[i].tps[j];
					ij=ij+1;
				}

			}
			for(i=0;i<(pop->ndim*(pop->ndim+1)/2);i++) {
				fisher[i]=pop->fisher[i];
			}
		}
	}
	return;
}

/*----------------------------------------------------------------
*/
#undef FREE_ARG
#undef NMAX
#undef IA 
#undef IM 
#undef AM 
#undef IQ 
#undef IR 
#undef NTAB 
#undef NDIV 
#undef EPS
#undef RNMX 
#undef SMALL
#undef TINY
#undef SWAP 
#undef tol
#undef NRANSI 
#undef CRI_MAX
#undef EPSILON
#undef FREQMIN
#undef FTOL
