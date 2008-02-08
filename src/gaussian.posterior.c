/* Copyright (C) 2002  Susanne Gammelgaard Bøttcher, Claus Dethlefsen */
/* Licenced under "GPLv2 or later version" licence. */

#include <R.h>
#include <Rmath.h>

int *ivector(int nl, int nh)
{
   int *v;

   v=(int *) R_alloc((unsigned) (nh-nl+1)*sizeof(int),sizeof(int));
   if ( v == NULL ){
      error("memory allocation failure in ivector()"); return(NULL);
   }
   return v-nl;
}

void free_ivector(int *v, int nl, int nh) { free((char*) (v+nl)); }

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
   int i;
   double **m;

   m=(double **) R_alloc((unsigned) (nrh-nrl+1)*sizeof(double*),sizeof(double*));
   if ( m == NULL ){
      error("memory allocation failure 1 in dmatrix()"); return(NULL);
   }
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
	   m[i]=(double *) R_alloc((unsigned) (nch-ncl+1)*sizeof(double),sizeof(double));
      if ( m[i] == NULL ){
         error("memory allocation failure 2 in dmatrix()"); return(NULL);
      }
      m[i] -= ncl;
   }
   return m;
}

void printmat(double **mat, int nr, int nc) {
	int i,j;
	for (i=1; i<=nr; i++) {
		for (j=1; j<=nc; j++)
			Rprintf("%f\t",mat[i][j]);
		Rprintf("\n");
	}
}

void asmatrix(double *vek, double **mat, int nr, int nc) {
	int i,j;
	for (i=1; i<=nr; i++) {
		for (j=1; j<=nc; j++) {
			mat[i][j] = vek[j-1+(i-1)*nc];
		}
	}

}

double** matcopy(double **mat, int nr, int nc) {
	/* copy mat[i][j] into nat[i][j] */
	int i,j;
	double **nat;
	nat = dmatrix(1,nr,1,nc);

	for (i=1; i<=nr; i++) {
		for (j=1; j<=nc; j++) {
			 nat[i][j] = mat[i][j];
		}
	}
	return(nat);
}

double** matmult(double **a, double **b, int nra, int nca, int ncb) {
	double **c;
	int i,j,k;
	c = dmatrix(1,nra,1,ncb);
	for (i=1; i<=nra; i++)
		for (j=1; j<=ncb; j++)
			c[i][j] = 0.0;

	for (i=1; i<=nra; i++) 
		for (k=1; k<=ncb; k++)
			for (j=1; j<=nca; j++)
				c[i][k] += a[i][j]*b[j][k];
	return(c);
}

double** matsum(double **a, double **b, int nr, int nc) {
	double **c;
	int i,j;
	c = dmatrix(1,nr,1,nc);

	for (i=1; i<=nr; i++)
		for (j=1; j<=nc; j++)
			c[i][j] = a[i][j] + b[i][j];
	return(c);
}

double** matminus(double **a, double **b, int nr, int nc) {
	double **c;
	int i,j;
	c = dmatrix(1,nr,1,nc);

	for (i=1; i<=nr; i++)
		for (j=1; j<=nc; j++)
			c[i][j] = a[i][j] - b[i][j];
	return(c);
}

double** transp (double **a, int n, int m) {
	double **b;
	int i,j;
	b = dmatrix(1,m,1,n);
	for (i=1; i<=n; i++)
		for (j=1; j<=m; j++)
			b[j][i] = a[i][j];
	return(b);
}

int invers(double **a, int n, double **b, int m)
{
   int *indxc,*indxr,*ipiv;
   int i,icol=1,irow=1,j,k,l,ll;
   double big,dum,pivinv;

   if( (indxc = ivector(1,n)) == NULL){ return(-1); }
   if( (indxr = ivector(1,n)) == NULL){ return(-1); }
   if( (ipiv  = ivector(1,n)) == NULL){ return(-1); }
   for (j=1;j<=n;j++) ipiv[j]=0;
   for (i=1;i<=n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
         if (ipiv[j] != 1)
            for (k=1;k<=n;k++) {
               if (ipiv[k] == 0) {
                  if (fabs(a[j][k]) >= big) {
                     big=fabs(a[j][k]);
                     irow=j;
                     icol=k;
                  }
               } else if (ipiv[k] > 1){
                  error("Invers: Singular Matrix-1");
                  return(-1);
               }
            }
      ++(ipiv[icol]);
      if (irow != icol) {
         for (l=1;l<=n;l++){
            double temp=a[irow][l]; a[irow][l]=a[icol][l]; a[icol][l]=temp;
         }
         for (l=1;l<=m;l++){
            double temp=b[irow][l]; b[irow][l]=b[icol][l]; b[icol][l]=temp;
         }
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if (a[icol][icol] == 0.0){
         error("Invers: Singular Matrix-2");
         return(-1);
      }
      pivinv=1.0/a[icol][icol];
      a[icol][icol]=1.0;
      for (l=1;l<=n;l++) a[icol][l] *= pivinv;
      for (l=1;l<=m;l++) b[icol][l] *= pivinv;
      for (ll=1;ll<=n;ll++)
         if (ll != icol) {
            dum=a[ll][icol];
            a[ll][icol]=0.0;
            for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
            for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
         }
   }
   for (l=n;l>=1;l--) {
      if (indxr[l] != indxc[l]){
         for (k=1;k<=n;k++){
            double temp     = a[k][indxr[l]];
            a[k][indxr[l]] = a[k][indxc[l]];
            a[k][indxc[l]] = temp;
         }
      }
   }
   return(0);
}

void postc0(double *mu, double *tau, double *rho, double *phi, double
	    *loglik, double *y, int *n)
{
	int i;
	double logscale,logk,mscore;
	double oldtau,oldmu;

	for(i = 0; i < *n; i++) {
		
		logscale = log(*phi)+log1p(1.0/(*tau));
		logk = lgammafn( 0.5*(1.0+*rho) ) - lgammafn(*rho*0.5);
		logk -= 0.5*(logscale + log(M_PI));
		mscore = logk - 0.5*(*rho+1.0)*log1p( (y[i]-*mu)*(y[i]-*mu)/exp(logscale));
		*loglik += mscore;

		oldtau = *tau;
		oldmu  = *mu;

		(*tau)++;
		(*rho)++;
		*mu = (oldtau*(*mu)+y[i])/(*tau);
		*phi+= (y[i]-(*mu))*y[i] + (oldmu-(*mu))*oldtau*oldmu;
	}
}

void postc(double *mu, double *tau, double *rho, double *phi, double
	    *loglik, double *y, double *z, int *n, int *d)
{

	int i,j;
	double logscale,logk,mscore;
	double **oldtau=0, **oldmu=0, **mtau, **mmu, **tauinv=0;
	double **zero, **zi, **ziy;

	/* allocate space for matrices */
	mtau   = dmatrix(1,*d,1,*d);
	zi     = dmatrix(1,*d,1,1);
	ziy    = dmatrix(1,*d,1,1);
	mmu    = dmatrix(1,*d,1,1);
	zero   = dmatrix(1,*d,1,1);

	/* copy arguments into the matrices */
	asmatrix(mu,mmu,*d,1);
	asmatrix(tau,mtau,*d,*d);
	
	for(i = 1; i <= *n; i++) {

		tauinv = matcopy(mtau,*d,*d);
		invers(tauinv, *d, zero, 1);

		for (j=1; j<=*d; j++) {
			zi[j][1] = z[j-1+(i-1)*(*d)];
		}
		
		logscale = log(*phi) +
			log1p(
				matmult(
					transp(zi,*d,1),
					matmult(tauinv,zi,*d,*d,1),
					1,*d,1
					)[1][1]
				);
		
		logk = lgammafn( 0.5*(1.0+*rho) ) - lgammafn(*rho*0.5);
		logk -= 0.5*(logscale + log(M_PI));
		
		mscore =  logk - 0.5*(*rho+1)*
			log1p(
				(y[i-1] - matmult(
					transp(zi,*d,1),
					mmu,1,*d,1
					)[1][1]
					)
				*
				(y[i-1] - matmult(
					transp(zi,*d,1),
					mmu,1,*d,1
					)[1][1])
				/exp(logscale)
				);
		
		*loglik += mscore;
		oldtau = matcopy(mtau,*d,*d);
		oldmu  = matcopy(mmu,*d,1);
		
		mtau = matsum(mtau, 
			      matmult(zi,transp(zi,*d,1),*d,1,*d)
			      , *d, *d
			);
		tauinv = matcopy(mtau,*d,*d);
		invers(tauinv, *d, zero, 1);

		for (j=1;j<=*d;j++)
			ziy[j][1] = zi[j][1]*y[i-1];

		mmu = matmult(tauinv,
			      matsum(
				      matmult(oldtau,mmu,*d,*d,1),
				      ziy,
				      *d,1)
			      ,*d,*d,1);

		(*rho)++;
		(*phi) += (y[i-1]-
			 matmult(
				 transp(zi,*d,1),
				 mmu,1,*d,1)[1][1])*y[i-1]
			+
			matmult(
				transp(
					matminus(oldmu,mmu,*d,1),
					*d,1
					),
				matmult(
					oldtau,
					oldmu,
					*d,*d,1
					),
				1,*d,1
				)[1][1];
	} 
	
	for (i=1; i<=*d;i++)
		mu[i-1] = mmu[i][1];
	for (i=1; i<=*d; i++)
		for (j=1; j<=*d; j++)
			tau[(*d)*(j-1)+i-1] = mtau[i][j];
	
} 

