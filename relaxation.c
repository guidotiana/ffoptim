#include "ffoptim.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

int my_gsl_matrix_fprintf(FILE *stream,gsl_matrix *m, int *h);


int ShouldRelax(struct parm_s *parm)
{
	int i;

	for (i=0;i<parm->nres;i++)
		if ( parm->restype[i] == 3 ) return 1; 

	return 0;
}

void FullRelaxation(struct NOErelax_s *x, double **res, struct parm_s *parm, int nframes)
{
	int i,j,k,h1,h2,iframe,a1,a2,chk;
	double rate,k1,k2,itmp;
	gsl_matrix *m,*tmp,*v0,*em;
		
	if (parm->verbose) fprintf(stderr,"Back-calculate by full relaxation.\n");

	tmp = gsl_matrix_alloc(x->nh, x->nh);
	v0 =  gsl_matrix_alloc(x->nh, x->nh);
	em =  gsl_matrix_alloc(x->nh, x->nh);
	m =   gsl_matrix_alloc(x->nh, x->nh);

	for (i=0;i<x->nh;i++)
	{
		for (j=0;j<x->nh;j++)  gsl_matrix_set(v0,i,j,0.);

		gsl_matrix_set(v0,i,i,1.);
		for (j=0;j<parm->nv0;j++)
		{
			if ( parm->iv0[j] == x->Hlist[i] ) gsl_matrix_set(v0,i,i,parm->v0[j]);
		}
	}
	if (parm->verbose>1) fprintf(stderr," Set diagonal peaks.\n");
	if (parm->verbose>1) { for (i=0;i<x->nh;i++) fprintf(stderr," v0[%d]=%f\n",x->Hlist[i],gsl_matrix_get(v0,i,i) );  fprintf(stderr,"\n"); }

	// loop over frames, which are independent
	for (iframe=0;iframe<nframes;iframe++)
	{
        for (j=0;j<parm->nres;j++) res[j][iframe] = 0;
        
		// fill matrix m with relaxation rates
		for (i=0;i<x->nh;i++)
			for (j=i;j<x->nh;j++)
			{
				// diagonal elements
				if (i==j)
				{
					if ( parm->ntaucatoms == 0 ) k1 = x->k1;	// prefactor
					else k1 = x->kAtoms[i][j];

					rate = 0;
					for (k=0;k<x->nh;k++)
						if (k!=i) rate += k1 / pow((double)x->dist[i][k][iframe],6.);

					gsl_matrix_set(m,i,i,- parm->mixingt * rate);
				}
				// off-diagonal elements 
				else 
				{
					if ( parm->ntaucatoms == 0 ) k2 = x->k2;	// prefactor
					else k2 = x->kAtoms[i][j];

					rate = k2 / pow((double)x->dist[i][j][iframe],6.);
					gsl_matrix_set(m,i,j,- parm->mixingt * rate);
					gsl_matrix_set(m,j,i,- parm->mixingt * rate);
				}
			}
		if (parm->verbose>1) fprintf(stderr," Relaxation matrix calculated:\n");
		if (parm->verbose>1) my_gsl_matrix_fprintf(stderr, m, x->Hlist);

		// exponential
		gsl_linalg_exponential_ss(m, em, GSL_PREC_DOUBLE);
		if (parm->verbose>1) fprintf(stderr," Exponential matrix:\n");
		if (parm->verbose>1) my_gsl_matrix_fprintf(stderr, em, x->Hlist);


		gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1.0,em,v0,0.0,tmp);	// tmp = em * v0
		if (parm->verbose>1) fprintf(stderr," NOE spectrum:\n");
		if (parm->verbose>1) my_gsl_matrix_fprintf(stderr, tmp, x->Hlist);


		for (i=0;i<parm->nres;i++)
        {
			if (parm->restype[i]==3)
			{ 
				h1 = x->backH[ parm->ires1[i] ];
				h2 = x->backH[ parm->ires2[i] ];
				if (h1==-1 || h2==-1) Error("Atom in restrain list with type 3 is not an hydrogen");
				res[i][iframe] = parm->noe0 * gsl_matrix_get(tmp,h1,h2);
				if (parm->lowerMatrix) res[i][iframe] = parm->noe0 * gsl_matrix_get(tmp,h2,h1);
				if (parm->absnoe && res[i][iframe]<0) res[i][iframe] = -res[i][iframe];
			}
            if (parm->restype[i]==5)
            {
                chk=0;
                for (j=0;j<parm->nagroup;j++)
                    for (k=0;k<parm->nagroup;k++)
                    {
                        if ( !strcmp(parm->alabel1[i], parm->agrouplabel[j]) &&
                            !strcmp(parm->alabel2[i], parm->agrouplabel[k]) )
                        {
                            a1 = parm->agroup[j];
                            a2 = parm->agroup[k];

                            h1 = x->backH[ a1 ];
                            h2 = x->backH[ a2 ];
                            if (h1==-1 || h2==-1) Error("Atom in restrain list with type 3 is not an hydrogen");
                            itmp = parm->noe0 * gsl_matrix_get(tmp,h1,h2);
                            if (parm->lowerMatrix) itmp = parm->noe0 * gsl_matrix_get(tmp,h2,h1);
                            if (parm->absnoe && itmp<0) itmp = -itmp;
                            res[i][iframe] += itmp;
                            chk=1;
                        }
                    }
                if (chk==0)
                {
                    fprintf(stderr,"WARNING: label (%s %s) not defined\n",parm->alabel1[i],parm->alabel2[i]);
                    Error("Label not defined");
                }
            }
        }
	}
	if (parm->verbose) fprintf(stderr," Done.\n");

  	gsl_matrix_free(m);
  	gsl_matrix_free(tmp);
	gsl_matrix_free(v0);
}


int LookbackHydrogens(struct top_s *top, int *h, struct parm_s *parm, int *backH, int *iv0, int nv0)
{
	int i,j,ok;
	char message[30];

	// loop over hydrogens on list of diagonal NOE
	for (j=0;j<nv0;j++)
	{
		// set hydrogen list from diagonal NOE
		h[j] = iv0[j];
		
		// check they are really hydrogens
		for (i=0;i<top->natoms;i++)
			if ( h[j] == (top->top_atoms[i]).nr )
				if ( (top->top_atoms[i]).type[0] != 'H' && (top->top_atoms[i]).type[0] != 'h' )
				{
					sprintf(message,"Atom %d is of type %s which is not regarded as a hydrogen.",h[j],(top->top_atoms[i]).type);
					Error(message);
				}

		// set lookback
		backH[h[j]] = j;
	}


	// check that every off-diagonal signal has the diagonal
	for (j=0;j<parm->nres;j++)
		if (parm->restype[j]==3)
		{
			ok=0;
			for (i=0;i<nv0;i++)
			{
				if ( parm->ires1[j] == iv0[i]) ok++;
				if ( parm->ires2[j] == iv0[i]) ok++;
			}
			if ( ok != 2 ) 
			{
				sprintf(message,"Restrain %d (%d-%d) is of type 3 but does not have a diagonal intensity.",j, parm->ires1[j], parm->ires2[j]);
				Error(message);
			}
		}


	if (parm->verbose) fprintf(stderr," Checked %d hydrogens.\n",nv0);
	return nv0;
}



int SelectHydrogens(struct top_s *top, int *h, struct parm_s *parm, int *backH)
{
	int i,nh=0;

	if (parm->verbose) fprintf(stderr,"Select hydrogens for full relaxation:\n");

	for (i=0;i<top->natoms;i++) backH[i]=-1;

	// select hydrogens
	for (i=0;i<top->natoms;i++)
		if ( (top->top_atoms[i]).type[0] == 'H' || (top->top_atoms[i]).type[0] == 'h' )
		{
			// list of hydrogens
			h[nh] = (top->top_atoms[i]).nr;
			backH[ (top->top_atoms[i]).nr ] = nh;

			if (parm->verbose) fprintf(stderr," %d:\t%d\t%s\n",nh,h[nh],(top->top_atoms[i]).type);
			nh++;	
		}


	if (parm->verbose) fprintf(stderr,"Selected %d hydrogens.\n",nh);
	return nh;
}

struct NOErelax_s *AlloNOErelax(int natoms,int nh, int verbose)
{
	int i;
	struct NOErelax_s *x;
	
	x = (struct NOErelax_s *) calloc(1,sizeof(struct NOErelax_s));
	if (!x) Error("Cannot allocate NOErelax structure");

	x->Hlist = (int *) calloc(nh,sizeof(int));					
	if (!x->Hlist) Error("Cannot allocate Hlist in NOErelax structure");

	x->backH = (int *) calloc(natoms+1,sizeof(int));		// numbering of atoms starts from 1
	if (!x->backH) Error("Cannot allocate backH in NOErelax structure");
	for (i=1;i<=natoms;i++) x->backH[i] = -1;

	

	if (verbose) fprintf(stderr," Allocated NOErelax lists (#h=%d, #atoms=%d, size=%f MB)\n",nh,natoms,(float)sizeof(int) * (natoms+nh) /1024/1024);


	return x;
}

void AlloNOErelaxDist(struct NOErelax_s *x, int nframes, int verbose)
{
	int i,j;
        
	x->dist = (float ***) calloc(x->nh,sizeof(float **)); 	
	if (!x->dist) Error("Cannot allocate dist in NOErelax structure");

	for (i=0;i<x->nh;i++)
	{
		x->dist[i] = (float **) calloc(x->nh,sizeof(float *));
		if (!x->dist[i]) Error("Cannot allocate dist in NOErelax structure");

		for (j=0;j<x->nh;j++)
		{
			x->dist[i][j] = (float *) calloc(nframes,sizeof(float));
			if (!x->dist[i][j]) Error("Cannot allocate dist in NOErelax structure");
		}

	}

	if (verbose) fprintf(stderr," Allocated distances in NOErelaxDist structure (#h=%d, #frames=%d. size=%f MB)\n",x->nh,nframes,(float)sizeof(float) * x->nh * x->nh * nframes/1024/1024);
}

double Prefactor(int diag, double tauc, double omega, int verbose)
{
	double j0,jw,j2w,pp;

	j0 = tauc;
	jw = tauc / ( 1. + omega * omega * tauc * tauc );
	j2w = tauc / ( 1. + 4. * omega * omega * tauc * tauc );

	if ( diag == 1) 
	{
		pp =  KRELAX * ( 0.3 * jw + 0.1 * j0 + 0.6 * j2w );
		if (verbose) fprintf(stderr," J(0) = %e ps J(w)=%e ps J(2w)=%e ps\n",j0,jw,j2w);
		if (verbose) fprintf(stderr," Prefactor for diagonal elements     = %e ns^-1\n",pp);
		return pp / 1E9;	// input in ps, I want the result in ms
	}
	else 
	{	pp = KRELAX * ( 0.6 * j2w - 0.1 * j0 );
		if (verbose) fprintf(stderr," J(0) = %e ps J(2w)=%e ps\n",j0,j2w);
		if (verbose) fprintf(stderr," Prefactor for non-diagonal elements = %e ns^-1\n",pp);


		if ( verbose && pp > 0 )
		{
			fprintf(stderr,"+----------------------------+\n");			
			fprintf(stderr,"| Cross peaks are negative.  |\n");			
			fprintf(stderr,"+----------------------------+\n");			
		}
	
		return	pp / 1E9;	// input in ps, I want the result in ms
	}
}

void PrefactorAtom(double **x, struct parm_s *parm, int nh, int *Hlist) 
{
	int i,j,i1,i2,a1,a2,k,chk=0;

	if (parm->verbose) fprintf(stderr,"Calculate prefactors for each atom\n");

			// x is labelled with the H id (not the atom id)

	// defaults
	for (i=0;i<nh;i++)
		for (j=0;j<nh;j++) 
		{
			a1 = Hlist[i];	// atom id
			a2 = Hlist[j];	// atom id
			if (a1==a2) x[i][j] = Prefactor(1,parm->tauc,parm->omega,0);
			else x[i][j] = Prefactor(0,parm->tauc,parm->omega,0);
		}

	// fill with values read from input
	for (i=0;i<nh;i++)
		for (j=0;j<nh;j++) 
		{
			a1 = Hlist[i];		// atom id
			a2 = Hlist[j];		// atom id
			for (k=0;k<parm->ntaucatoms;k++)	// look into the table
			{
				i1 = parm->taucatoms1[k];		// atom id in tauc table
				i2 = parm->taucatoms2[k];		// atom id in tauc table
				if ( ( a1 == i1 && a2 == i2 ) || ( a1 == i2 && a2 == i1 ) ) 
				{
					if (a1==a2)  x[i][j] = Prefactor(1,parm->taucatomsv[k],parm->omega,0);
					else {  x[i][j] = Prefactor(0,parm->taucatomsv[k],parm->omega,0); x[j][i] = x[i][j]; }
					chk++; 
				}
			}
		
		}

	if (parm->verbose)  fprintf(stderr,"Found %d definitions of NOE prefactors depending on atom-based tau_c\n",chk);
	if (parm->verbose>1) 
	{
		fprintf(stderr," Prefactors with varying tauc: (default=%lf)\n",parm->tauc);
		for (i=0;i<nh;i++)
		{
			for (j=0;j<nh;j++) fprintf(stderr," %4e",x[i][j]);
			fprintf(stderr,"\n");
		}
	}
        
}

int my_gsl_matrix_fprintf(FILE *stream,gsl_matrix *m, int *h)
{
        size_t rows=m->size1;
        size_t cols=m->size2;
	int i,j;

	if (cols<15)
	{
		fprintf(stream,"    ");
		for (i=0;i<cols;i++)  fprintf(stream,"%12d ",h[i]);
		fprintf(stream,"\n");
		for (i=0;i<rows;i++)
		{
			fprintf(stream,"%4d ",h[i]);
			for (j=0;j<cols;j++) fprintf(stream,"%12g ",gsl_matrix_get(m,i,j));
			fprintf(stream,"\n");
		}
	}
	else
	{
		for (i=0;i<rows;i++)
			for (j=0;j<cols;j++) fprintf(stream,"[%4d,%4d]=%g\n",h[i],h[j],gsl_matrix_get(m,i,j));

	}
	fprintf(stream,"\n");
                return 0;
}

void MoveV0(struct top_s *x,  struct parm_s *p)
{
	int i,j,ok,n=0;

	if (p->verbose) fprintf(stderr,"Activate moves on pairs parameters based on V0\n");

	for (i=0;i<x->npairs;i++)
	{	
		ok = 0;
		for (j=0;j<p->nv0;j++) 
		{
			if ( ((x->top_pairs)+i)->i1 == p->iv0[j] ) ok++;
			if ( ((x->top_pairs)+i)->i2 == p->iv0[j] ) ok++;
		}

		if (ok == 2 && ((x->top_pairs)+i)->move == 0 )  
		{
			((x->top_pairs)+i)->move = 1;
			n++;
		}
	}

	if (p->verbose) fprintf(stderr," Activated %d pairs\n",n);
}

void MoveA0(struct top_s *x,  struct parm_s *p)
{
	// activate not the H in the V0 list, but the heavy atoms bound to them

	int i,j,ok,n=0,w,iv,k;
	struct network_s net;

	if (p->verbose) fprintf(stderr,"Activate moves on pairs parameters involving heavy atom based on V0\n");

	net = BuildNetwork(x, p->verbose);

	for (i=0;i<x->npairs;i++)					// check all pairs in top
	{	
		ok = 0;
		iv = ((x->top_pairs)+i)->i1;				// iv is the first atom of the pair
		for (k=0;k< net.nvic[iv];k++)				// loop over its neigbours
			for (j=0;j<p->nv0;j++) 				// look in the table of v0
			{
				w =  net.vic[iv][k];			// w is the neigbour of iv
				if ( w == p->iv0[j] ) ok=1;		// is it in the table of v0?
			}

		iv = ((x->top_pairs)+i)->i2;				// same with the second atom of the pair
		for (k=0;k< net.nvic[iv];k++)
			for (j=0;j<p->nv0;j++) 
			{
				w =  net.vic[iv][k];
				if ( w == p->iv0[j] && ok == 1) ok=2;
			}

		if (ok == 2 && ((x->top_pairs)+i)->move == 0 )  		// check if the two iv of the pair have neigbhours in v0
		{
			((x->top_pairs)+i)->move = 1;
			n++;
		}
	}

	if (p->verbose) fprintf(stderr," Activated %d pairs\n",n);
}

float OptimizeTauc(int nframes, struct parm_s *parm, double **restrain, double *resnew, struct NOErelax_s *NOErelax, double *eold, double *enew, double *resnew2)
{
	int iacc=0,it;
	double chi2old,chi2,k=1.,taucold=parm->tauc,taucini=parm->tauc,deltat;

	fprintf(stderr,"Optimize tauc\n");

	// calculate initial restrains
	GetRestrain(parm,restrain,nframes,NOErelax);
	ReweightRestrain(resnew,restrain,eold,enew,parm,nframes,parm->iskip,resnew2);
			
	if  (!parm->correl)
	{
		if(!parm->rescalenoe0)	
			chi2old = Chi2(parm->resexp,resnew,parm->ressigma,parm->nres,NULL);
		else
			chi2old = Chi2(parm->resexp,resnew,parm->ressigma,parm->nres,&k);
	}
	else
		chi2old =  1. - Correlation(parm->resexp,resnew,parm->nres);

	if (parm->verbose) fprintf(stderr,"initial chi2 = %lf\n",chi2old/parm->nres);

	// optimize
	for (it=0;it<parm->noptauc;it++)
	{
		if (parm->verbose) fprintf(stderr,"Optimization step tauc n. %d\n",it);
		deltat = (frand()-0.5)*taucini/10;
		parm->tauc += deltat;		
		if (parm->verbose) fprintf(stderr," attempted tauc = %lf\n",parm->tauc);

		if ( parm->ntaucatoms == 0 )		// if same tauc for all noe
		{
			NOErelax->k1 = Prefactor(1,parm->tauc,parm->omega,0);
			NOErelax->k2 = Prefactor(0,parm->tauc,parm->omega,0);	
		}
		else					// if each noe has different tauc
		{
			NOErelax->kAtoms = AlloDoubleMatrix(NOErelax->nh,NOErelax->nh);
			PrefactorAtom(NOErelax->kAtoms, parm, NOErelax->nh, NOErelax->Hlist);
		}

	                                                                                           
		// calculate restrains
		FullRelaxation(NOErelax,restrain, parm, nframes);
		ReweightRestrain(resnew,restrain,eold,enew,parm,nframes,parm->iskip,resnew2);
		
		// optimize chi2
		if  (!parm->correl)
		{
			if(!parm->rescalenoe0)	
				chi2 = Chi2(parm->resexp,resnew,parm->ressigma,parm->nres,NULL);
			else
				chi2 = Chi2(parm->resexp,resnew,parm->ressigma,parm->nres,&k);
		}
		else
			chi2 =  1. - Correlation(parm->resexp,resnew,parm->nres);
	
		if (parm->verbose) fprintf(stderr," new chi2 = %lf\n",chi2/parm->nres);

		if (chi2 < chi2old)
		{
			chi2old = chi2;
			taucold = parm->tauc;
			if (parm->verbose) fprintf(stderr,"Accepted.\n");
			iacc++;
		}
		else   
		{
			parm->tauc -= deltat;
			if (parm->verbose) fprintf(stderr,"Rejected.\n");
		}

		if  (!parm->correl)
		{
			fprintf(stderr,"\t%d\tchi2 = %lf",it,chi2old/parm->nres);
			fprintf(stderr,"\n");

		}
		else
			fprintf(stderr,"\t%d\tr = %lf\n",it,1.-chi2old);

	}

	if  (!parm->correl)
		fprintf(stderr,"\tNew chi2 = %lf\n\tAccepted changes = %d / %d\n",chi2old/parm->nres,iacc,parm->noptim);
	else
		fprintf(stderr,"\tNew r = %lf\n\tAccepted changes = %d / %d\n",1.-chi2old,iacc,parm->noptim);

	return chi2old;
exit(0);
}


