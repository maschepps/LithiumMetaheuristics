
/******************************************************************************
   Defining structures
*/

typedef struct PROTOC /* Structure defining an individual or elementary protocol */
{
	int	ntps; /* Number of sampling times */
	double	*tps; /* Vector of the sampling times */
	int	ndim; /* Number of random effects (=> size of matrices)*/
	double	*fisher; /* Fisher information matrix of the individual protocol */
} PROTOC;

typedef struct POPROT /* Structure defining a population protocol */
{
	int 	maxnp; /* Max number of individual protocols, needed for memory allocation*/
	int	np; /* Number of individual protocols WARNING np<=maxnp !*/
	int	ndim; /* Number of random effects (=> size of matrices)*/
	int	*num; /* In this program we keep trace of the index of individual protocols*/
	PROTOC	*pind; /* In NPML we would only need the individual protocols themselves*/
	double	*freq; /* Frequency of each individual protocols */
	double 	*fisher; /* Fisher information matrix of population protocol */
	double 	*finv; /* Inverse of the Fisher information matrix */
	double	det; /* Determinant of the matrix */
} POPROT;

typedef struct matrix
{
	int 	nrow;
	int	ncol;
	double	**m;
} matrix;


/******************************************************************************
   Construction and destruction of structures
*/

void PROTOC_alloc( PROTOC *p, int ntps, int ndim)
{	
	p->ntps = ntps;
	p->tps = (double *)calloc( p->ntps, sizeof(double));
	assert(p->tps!=NULL);	
	p->ndim=ndim;
	p->fisher = (double *)calloc( p->ndim*(p->ndim+1)/2, sizeof(double));
	assert(p->fisher!=NULL);
	return;
}

POPROT *POPROT_alloc(int np,int ndim,int maxnp)
{	POPROT *p;
	p = (POPROT *)malloc(sizeof(POPROT));
	assert(p!=NULL);
	
	p->maxnp=maxnp;
	p->np = np;
	p->ndim=ndim;
/* Keep trace of the index of the individual protocols */
	p->num = (int *)calloc( p->maxnp, sizeof(int));
	assert(p->num!=NULL);
	p->pind = (PROTOC *)calloc( p->maxnp, sizeof(PROTOC));
	assert(p->pind!=NULL);
	p->freq = (double *)calloc( p->maxnp, sizeof(double));
	assert(p->freq!=NULL);
	p->fisher = (double *)calloc( p->ndim*(p->ndim+1)/2, sizeof(double));
	assert(p->fisher!=NULL);
	p->finv = (double *)calloc( p->ndim*(p->ndim+1)/2, sizeof(double));
	assert(p->finv!=NULL);
	p->det=1.0;	
	return p;
}

void PROTOC_copy(PROTOC *p,PROTOC *p1)
{
	int i,ndat;

	for(i=0;i<p->ntps;i++) {p1->tps[i]=p->tps[i];}
	ndat=(int)(p->ndim*(p->ndim+1))/2;
	for(i=0;i<ndat;i++) {p1->fisher[i]=p->fisher[i];
/*		fprintf(stderr,"%lf ",p->fisher[i]); */
		}
	return;
}

matrix *matrix_create(int nrow, int ncol)
{	int i,j;
	matrix *mat;
	mat = (matrix *)malloc(sizeof(matrix));
	assert(mat!=NULL);
	mat->nrow=nrow;
	mat->ncol=ncol;
	mat->m = (double **)calloc( nrow, sizeof(double *));
	assert(mat->m!=NULL);	
	for( i=0 ; i<nrow ; i++ )
	{
	mat->m[i] = (double *)calloc( ncol, sizeof(double));
	assert(mat->m[i]!=NULL);	
	for( j=0 ; j<ncol ; j++ )
		mat->m[i][j]=0.0;
	}
	return mat;
	
}

void PROTOC_print(PROTOC *p)
{	int i,icas=0,j;
	fprintf(stderr,"Temps du protocole :");
	for (i=0;i<p->ntps;i++) {
			fprintf(stderr,"%lf  ",p->tps[i]);}
		fprintf(stderr,"\n ");
	fprintf(stderr,"Matrice d'information de Fisher :\n");
	for (i=0;i<p->ndim;i++) {
		for (j=0;j<=i;j++) {
			fprintf(stderr,"%lf  ", p->fisher[icas]);
			icas++;}
		fprintf(stderr,"\n ");
	}
	return;
}

void POPROT_print(POPROT *pop,int nofish)
{	int i,icas=0,j,k;

	fprintf(stderr,"Population protocol with %d individual protocols :\n",pop->np);
	for(i=0;i<pop->np;i++) 
	{
		fprintf(stderr,"protocol %d with frequency %lf\n",pop->num[i],pop->freq[i]);
		fprintf(stderr,"Sampling times :");
		for (j=0;j<pop->pind[i].ntps;j++) 
		{
			fprintf(stderr,"%lf  ",pop->pind[i].tps[j]);
		}
		fprintf(stderr,"\n ");
		if(nofish==1) {
		fprintf(stderr,"Fisher information matrix :\n");
		icas=0;
		for (j=0;j<pop->pind[i].ndim;j++) 
		{
			for (k=0;k<=j;k++) 
			{
				fprintf(stderr,"%lf  ",pop->pind[i].fisher[icas]);
				icas++;
			}
			fprintf(stderr,"\n ");
		}
		fprintf(stderr,"\n ");
		}
	}
	return;
}

void POPROT_printfisher(POPROT *pop)
{	int icas,j,k;

	fprintf(stderr,"Population protocol with %d individual protocols :\n",pop->np);
	fprintf(stderr,"Fisher information matrix :\n");
	icas=0;
	for (j=0;j<pop->ndim;j++) 
	{
		for (k=0;k<=j;k++) 
		{
		fprintf(stderr,"%lf  ",pop->fisher[icas]);
		icas++;
		}
		fprintf(stderr,"\n ");
	}
	fprintf(stderr,"Inverse matrix :\n");
	icas=0;
	for (j=0;j<pop->ndim;j++) 
	{
		for (k=0;k<=j;k++) 
		{
		fprintf(stderr,"%lf  ",pop->finv[icas]);
		icas++;
		}
		fprintf(stderr,"\n ");
	}
	return;
}

void POPROT_calculefisher(POPROT *pop)
{	int i,iprot;

	for(i=0;i<pop->ndim*(pop->ndim+1)/2;i++)
	{
		pop->fisher[i]=0.0;
		for(iprot=0;iprot<pop->np;iprot++)
		{
	pop->fisher[i]=pop->fisher[i]+pop->pind[iprot].fisher[i]*pop->freq[iprot];
		}
	}
}

void matrix_print(matrix *mat)
{	int i,j;

	for (i=0;i<mat->nrow;i++) {
		for (j=0;j<mat->ncol;j++) {
			fprintf(stderr,"%lf  ",mat->m[i][j]);}
		fprintf(stderr,"\n ");}
	return;
}

void PROTOC_destroy(PROTOC *p)
{
	free(p->fisher);
	free(p->tps);	
}

void POPROT_destroy(POPROT *p)
{
	free(p->finv);	
	free(p->fisher);
	free(p->freq);	
	free(p->pind);	
	free(p->num);	
}

void matrix_destroy(matrix *mat)
{
	int i;
	for( i=0 ; i<mat->nrow ; i++ )
	{
	free(mat->m[i]);
	}
	free(mat->m);	
}
