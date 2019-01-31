#include "ffoptim.h"

int main(int argc, char *argv[])
{
  struct parm_s *parm;
  struct top_s *top, *topnew, *toporig;
  struct NOErelax_s *NOErelax=NULL;
  int iter,nframes,recalculateEnergies=0;
  char topi[100],command[1000];
  double *eold,*enew,**restrain,*resnew,*resnew2;
  float chi2;
  FILE *fchi2;

  Welcome(stderr);

  if (argc!=2) Error("No paramater file. For help use: ffoptim.x --help");

  if (!strcmp(argv[1],"-h") || !strcmp(argv[1],"--help")) Help();

  // Read parameters
  parm = ReadParms(argv[1]);
  
  // Allocate
  if (parm->verbose) fprintf(stderr,"Allocate memory\n");
  eold = AlloDouble(NFRAMESMAX);
  enew = AlloDouble(NFRAMESMAX);
  if (!ShouldRelax(parm)) restrain = AlloDoubleMatrix(parm->nres,NFRAMESMAX);
  topnew = AlloTop(NATOMSMAX,NBONDSMAX,NPAIRSMAX,NANGLESMAX,
			NDIHEDRALSMAX,NIMPROPERSMAX,NEXCLUSIONSMAX,parm->verbose); 
  resnew = AlloDouble(NRESMAX);
  resnew2 = AlloDouble(NRESMAX);

  if (strcmp(parm->origtopfile,"none"))
	toporig =  ReadTop(parm->origtopfile,parm->itpdir,parm->verbose);

  CheckFiles(parm);


  // if noSOL selected, make a mdp for non-solvent atoms
  if (parm->noSOL)
  {
	if (!strcmp(parm->mdpprotegrps,"none")) Error("If noSOL active, you must define mdpprotegrps");
	if (!strcmp(parm-> mdpprotxtcgrps,"none")) Error("If noSOL active, you must define mdpprotxtcgrps");
        UpdateMdp(parm->mdpfile,"tmp.mdp",parm->mdpprotxtcgrps,parm->mdpprotegrps);
  }

  // read original top
  if (!parm->restart)
  {
	  if (parm->verbose) fprintf(stderr,"Read initial topology\n");
	  top = ReadTop(parm->topfile,parm->itpdir,parm->verbose);
	  if (parm->filltop) FillTop(top,parm->combrule,parm->nochangeH);
	  if ( strcmp( parm->nonbondedFile,"none") ) InsertNonbondedTop(top,parm->combrule,parm->verbose,parm->nonbondedFile);
	  if ( strcmp( parm->bondedFile,"none") ) InsertDihedralsTop(top,parm->verbose,parm->bondedFile);
  	  sprintf(topi,"topol.0.top");
 	  PrintTop(topi, top, parm->verbose,0);
  }
  else
  {
	  if (parm->verbose) fprintf(stderr,"Read restart topology\n");
	  sprintf(topi,"topol.%d.top",parm->restart);	
	  if (parm->verbose) fprintf(stderr,"RESTARTING from topology %s\n",topi);
	  top = ReadTop(topi,parm->itpdir,parm->verbose);
  }

  // relax NOE
  if (ShouldRelax(parm)) 
  {
	if (parm->verbose) fprintf(stderr,"Setting stuff for Full Relaxation\n");
	NOErelax = AlloNOErelax(top->natoms, parm->nv0, parm->verbose);
	NOErelax->nh = LookbackHydrogens(top,NOErelax->Hlist,parm,NOErelax->backH,parm->iv0,parm->nv0);		// nh=nv0
	AlloNOErelaxDist(NOErelax,parm->nframesmax, parm->verbose);

	if (parm->uNOE>0) AdduNOE(parm,NOErelax);
  	restrain = AlloDoubleMatrix(parm->nres,NFRAMESMAX);

	if ( parm->ntaucatoms == 0 )		// if same tauc for all noe
	{
		NOErelax->k1 = Prefactor(1,parm->tauc,parm->omega,parm->verbose);
		NOErelax->k2 = Prefactor(0,parm->tauc,parm->omega,parm->verbose);
	}
	else					// if each noe has different tauc
	{
		NOErelax->kAtoms = AlloDoubleMatrix(NOErelax->nh,NOErelax->nh);
		PrefactorAtom(NOErelax->kAtoms, parm, NOErelax->nh, NOErelax->Hlist);
	}
	if ( !parm->restart && parm->movev0) 	// move only atoms in v0 list
	{
		MoveV0(top,parm);
 		sprintf(topi,"topol.0.top");
 	  	PrintTop(topi, top, parm->verbose,0);
	}
	if ( !parm->restart && parm->movea0) 	// mov only heavy atoms linked to v0 list
	{
		MoveA0(top,parm);
 		sprintf(topi,"topol.0.top");
 	  	PrintTop(topi, top, parm->verbose,0);
	}

  }

  // make dumb index file
  sprintf(command,"printf 'splitat 0 \n q \n' | %s/make_ndx -f %s >> gro.log 2>&1",parm->grobindir,parm->groprotfile);
  system(command);
  if (parm->verbose) fprintf(stderr,"%s\n",command);
  if( access( "index.ndx", F_OK ) == -1 ) Error("Cannot create index.ndx");
  RenameIndex();

  // Initial restrains
  if (parm->verbose) fprintf(stderr,"Calculating initial restrains.\n");
  PrintTop("tmp.top", top, parm->verbose, parm->noSOL);
  PrintInitialRestrains(parm,restrain,stderr,NOErelax);

  // open chi2 file
  if (parm->restart) fchi2 = fopen(parm->chi2file,"a");
  else fchi2 = fopen(parm->chi2file,"w");
  if (!fchi2) Error("Cannot open chi2 file for writing");

  // iterations
  if (parm->verbose) fprintf(stderr,"Starting iterations.\n");
  for (iter=parm->restart; iter<parm->niter; iter++)
  {
	fprintf(stderr,"\n*** ITERATION %d ***\n",iter);
	system("rm -f topol.tpr traj.xtc ener.edr energy.xvg gro.log confout.gro dist.xvg /#* step*");

  	sprintf(topi,"topol.%d.top",iter);

	// simulation
	if (iter==0 && strcmp(parm->testtraj,"none"))
	{
		TestTraj(parm);
		recalculateEnergies = 1;
	}
	else
		SimulateGromacs(iter, parm, topi);
	CheckExplosion();

	// get energies
	if (parm->noSOL) recalculateEnergies = 1;
	nframes = GetGromacsEnergies(parm, top, eold, NOErelax, recalculateEnergies);
       
	// get restrains
	if (parm->verbose) fprintf(stderr,"Gets restrains\n");
	GetRestrain(parm,restrain,nframes,NOErelax);
         
	// optimize potential
	if (parm->verbose) fprintf(stderr,"Optimizes the potential\n");
	chi2 = OptimizePotential(top,topnew,nframes,parm,restrain,resnew,eold,enew,resnew2,toporig);

	// if needed, optimize tauc
	if (NOErelax!=NULL && parm->noptauc) 
	{
		if (parm->verbose) fprintf(stderr,"Optimizes tauc\n");
		OptimizeTauc(nframes,parm,restrain,resnew,NOErelax,eold,enew,resnew2);
	}

	// write output
	PrintGromacsOutput(iter, chi2, parm, fchi2, resnew, resnew2);

	// write topology
	sprintf(topi,"topol.%d.top",iter+1);
 	PrintTop(topi, top, parm->verbose,0);

	// save trajectory
	sprintf(command,"mv traj.xtc traj_%d.xtc",iter);
	if (parm->verbose) fprintf(stderr,"%s\n",command);	
	system(command);

  }

  fclose (fchi2);
  fprintf(stderr,"Finished regularly.\n\n");
  exit(0);	
}

/*********************************************************************************************************/


void TestTraj(struct parm_s *parm)
{
	char command[500];
	
	if (parm->verbose) fprintf(stderr,"** TESTING MODE **\n Use as trajectory %s instead of doing the simulation.\n",parm->testtraj);
	sprintf(command,"cp %s traj.xtc",parm->testtraj);	
	system(command);
	if (parm->verbose) fprintf(stderr,"%s\n",command);
	sprintf(command,"cp %s tmp.mdp",parm->mdpfile);	
	system(command);
	if (parm->verbose) fprintf(stderr,"%s\n",command);
    	sprintf(command,"cp %s end_0.gro",parm->grofile);
	system(command);

}

void PrintGromacsOutput(int iter, double chi2, struct parm_s *parm, FILE *fchi2, double *resnew, double *resnew2)
{
	int i;
	char nfres[500];
	FILE *frout;

	fprintf(fchi2,"%d\t%f\n",iter,chi2);
	sprintf(nfres,"restrains.%d.dat",iter);
	frout = fopen(nfres,"w");
	if (!frout) Error("Cannot open file for writing restrains.");
	for (i=0;i<parm->nres;i++)
		fprintf(frout,"%3d\t%3d %3d\t%10lf %10lf\t%10lf %10lf\n",
			i,parm->ires1[i],parm->ires2[i],resnew[i],sqrt(resnew2[i]-resnew[i]*resnew[i]),
			parm->resexp[i],parm->ressigma[i]);	
	fclose(frout);
}

void SimulateGromacs(int iter, struct parm_s *parm, char *topi)
{
	int i;	
  	char command[1000], dumb[500];

	// grompp (em)
  	if (parm->verbose) fprintf(stderr,"Energy minimization.\n");
	if (iter==0)
		sprintf(command,"%s/grompp -f %s -p %s -c %s  >> gro.log 2>&1",parm->grobindir,parm->emdpfile,topi,parm->grofile);
	else
		sprintf(command,"%s/grompp -f %s -p %s -c end_%d.gro  >> gro.log 2>&1",parm->grobindir,parm->emdpfile,topi,iter-1);
	if (parm->verbose) fprintf(stderr,"%s\n",command);
	system(command);
	if( access( "topol.tpr", F_OK ) == -1 ) Error("Cannot create topol.tpr");

	// mdrun (em)
	sprintf(command,"%s/mdrun -v -nt %d -nice %d >> gro.log 2>&1",parm->grobindir,parm->nthreads,parm->nice);
	if (parm->verbose) fprintf(stderr,"%s\n",command);
	system(command);
	if( access( "ener.edr", F_OK ) == -1 ) Error("Cannot create ener.edr");

	// plain MD
  	if (parm->verbose) fprintf(stderr,"Sampling.\n");
	if (!parm->replex)
	{
		// grompp (md)
		sprintf(command,"mv confout.gro start_%d.gro",iter);
		if (parm->verbose) fprintf(stderr,"%s\n",command);	
		system(command);
		system("rm -f topol.tpr traj.xtc ener.edr energy.xvg gro.log confout.gro");
		sprintf(command,"%s/grompp -f %s -p %s -c start_%d.gro >> gro.log 2>&1",parm->grobindir,parm->mdpfile,topi,iter);
		if (parm->verbose) fprintf(stderr,"%s\n",command);
		system(command);
		if( access( "topol.tpr", F_OK ) == -1 ) Error("Cannot create topol.tpr");

		// mdrun (md)
		sprintf(command,"%s/mdrun -v -nt %d -nice %d >> gro.log 2>&1",parm->grobindir,parm->nthreads,parm->nice);
		if (parm->verbose) fprintf(stderr,"%s\n",command);
		system(command);
		if( access( "traj.xtc", F_OK ) == -1 ) Error("Cannot create traj.xtc");
		if( access( "ener.edr", F_OK ) == -1 ) Error("Cannot create ener.edr");
		if( access( "confout.gro", F_OK ) == -1 ) Error("Cannot create confout.gro");
		sprintf(command,"mv confout.gro end_%d.gro",iter);
		if (parm->verbose) fprintf(stderr,"%s\n",command);	
		system(command);
	}
	// Replica Exchange
	else
	{
		// grompp (md)
		sprintf(command,"mv confout.gro start_%d.gro",iter);
		if (parm->verbose) fprintf(stderr,"%s\n",command);	
		system(command);
		system("rm -f topol.tpr traj.xtc ener.edr energy.xvg gro.log confout.gro");
		for (i=0;i<parm->ntemp;i++)
		{
			ChangeTmdp(parm->mdpfile,i,parm->temperatures);
			sprintf(command,"%s/grompp -f tmp%d.mdp -p %s -c start_%d.gro -o topol%d.tpr  >> gro.log 2>&1",
							parm->grobindir,i,topi,iter,i);
			if (parm->verbose) fprintf(stderr,"%s\n",command);
			system(command);
			sprintf(dumb,"topol%d.tpr",i);
			if( access( dumb, F_OK ) == -1 ) Error("Cannot create topol.tpr");	
		}

		// mdrun (md)
		sprintf(command,"mpirun -np %d ",parm->ntemp);
		if (strcmp(parm->hostfile,"")) sprintf(command,"%s --hostfile %s ",command,parm->hostfile);
		sprintf(command,"%s %s/mdrun_mpi -v -nt %d -multi %d -replex %d -nice %d >> gro.log 2>&1",	
			command,parm->grobindir,parm->nthreads,parm->ntemp,parm->replex,parm->nice);
		if (parm->verbose) fprintf(stderr,"%s\n",command);
		system(command);
		if( access( "traj0.xtc", F_OK ) == -1 ) Error("Cannot create traj0.xtc");
		if( access( "ener0.edr", F_OK ) == -1 ) Error("Cannot create ener0.edr");
		if( access( "confout0.gro", F_OK ) == -1 ) Error("Cannot create confout0.gro");
		sprintf(command,"mv confout0.gro end_%d.gro",iter);
		if (parm->verbose) fprintf(stderr,"%s\n",command);	
		system(command);
		sprintf(command,"mv ener0.edr ener.edr; mv traj0.xtc traj.xtc; cp topol0.tpr topol.tpr");
		if (parm->verbose) fprintf(stderr,"%s\n",command);	
		system(command);
	}


        if (parm->verbose) fprintf(stderr,"Sampling finished.\n");
}

int GetGromacsEnergies(struct parm_s *parm, struct top_s *top, double *eold, struct NOErelax_s *NOErelax, int recalculate)
{
  	char command[1000];
	int nframes;
	
	if (parm->verbose) fprintf(stderr,"Gets energies.\n");

	if (recalculate)
	{
		system("rm -f topol.tpr ener.edr gro.log");
		PrintTop("tmp.top", top, parm->verbose, parm->noSOL);
		sprintf(command,"%s/grompp -f tmp.mdp -p tmp.top -c %s -n -maxwarn 99 >> gro.log 2>&1",
			parm->grobindir,parm->groprotfile);
		if (parm->verbose) fprintf(stderr,"%s\n",command);
		system(command);
		if( access( "topol.tpr", F_OK ) == -1 ) Error("Cannot create topol.tpr for rerun");
	
		sprintf(command,"%s/mdrun -s topol.tpr -rerun traj.xtc -nt %d -nice %d >> gro.log 2>&1",parm->grobindir,parm->nthreads,parm->nice);
		if (parm->verbose) fprintf(stderr,"%s\n",command);
		system(command);
	}

	sprintf(command,"echo '%s' | %s/g_energy -f ener.edr >> gro.log 2>&1",parm->energytype,parm->grobindir);
	if (parm->verbose) fprintf(stderr,"%s\n",command);
	system(command);
	if( access( "energy.xvg", F_OK ) == -1 ) Error("Cannot create energy.xvg");
	nframes = ParseXVG("energy.xvg",eold);
	if (NOErelax && nframes>parm->nframesmax) Error("nframesmax too small in input file");
	if (parm->verbose) fprintf(stderr,"Read energies in %d frames\n",nframes);

	return nframes;
} 

void AverRestrain(double **restrain, double *average, int nframes, int nres)
{
	int i,j;

	for (i=0;i<nres;i++)
	{
		average[i] = 0;
		for (j=0;j<nframes;j++)
			average[i] += restrain[i][j]/nframes;
	}
}

double Chi2(double *x, double *xexp, double *sigma, int n, double *k)
{
	int i;
	double chi2=0,si=0,si2=0,k2;

	// if you want to always use the best k to rescale noe
	if (k!=NULL)
	{
		for (i=0;i<n;i++)
		{
			si += x[i] * xexp[i] / sigma[i] / sigma[i];
			si2 += x[i] * x[i] / sigma[i] / sigma[i];
		}
		
		k2 = si / si2;
		*k = k2;
	}
	else k2 = 1.;

	// chi2
	for (i=0;i<n;i++)
		chi2 += ( k2 * x[i] - xexp[i])*( k2 * x[i] - xexp[i]) / ( sigma[i] * sigma[i] );
        
	return chi2;
}

double Correlation(double *x, double *xexp, int n)
{
	int i;
	double r=0,xav=0,xexpav=0,d1=0,d2=0;

	for (i=0;i<n;i++)
	{
		xav += x[i]/n;
		xexpav += xexp[i]/n;
	}

	for (i=0;i<n;i++)
	{
		r += (x[i]-xav) * (xexp[i]-xexpav);
		d1 += (x[i]-xav) *  (x[i]-xav);
		d2 += (xexp[i]-xexpav) *  (xexp[i]-xexpav);
	}	
	return r/sqrt(d1)/sqrt(d2);
}


float OptimizePotential(struct top_s *top, struct top_s *topnew, int nframes, struct parm_s *parm, 
			double **restrain, double *resnew, double *eold, double *enew, double *resnew2, struct top_s *origtop)
{
	int iacc=0,it,n,i;
	double chi2old,chi2,k=1.;
	char command[1000];

	fprintf(stderr,"Optimize potential\n");


	// initial chi2
	ReweightRestrain(resnew,restrain,eold,eold,parm,nframes,parm->iskip,resnew2);

	if (!parm->correl)
	{
		if(!parm->rescalenoe0)	
			chi2old = Chi2(parm->resexp,resnew,parm->ressigma,parm->nres,NULL);
		else
		{
			chi2old = Chi2(parm->resexp,resnew,parm->ressigma,parm->nres,&k);
			fprintf(stderr,"Noe0=%lf (rescale=%lf)",k*parm->noe0,k);

		}
		fprintf(stderr,"Initial chi2 = %lf\n",chi2old/parm->nres);
	}
	else
	{
		chi2old = 1. - Correlation(parm->resexp,resnew,parm->nres);
		fprintf(stderr,"Initial r = %lf\n",1.-chi2old);
	}

	// optimize
	for (it=0;it<parm->noptim;it++)
	{
		if (parm->verbose) fprintf(stderr,"Optimization step n. %d\n",it);
		CopyTop(top,topnew);
		RandomChangeEnergy(parm,topnew,origtop);

		// regrompp
		system("rm -f topol.tpr ener.edr gro.log");
		if (parm->noSOL)
			sprintf(command,"%s/grompp -f tmp.mdp -p tmp.top -c %s -n -maxwarn 99 >> gro.log 2>&1",		
						parm->grobindir,parm->groprotfile);
		else
			sprintf(command,"%s/grompp -f %s -p tmp.top -c %s -n -maxwarn 99 >> gro.log 2>&1",
						parm->grobindir,parm->mdpfile,parm->grofile);
			
		if (parm->verbose) fprintf(stderr,"%s\n",command);
		system(command);
		if( access( "topol.tpr", F_OK ) == -1 ) Error("Cannot create topol.tpr for rerun");

		// rerun
		sprintf(command,"%s/mdrun -s topol.tpr -rerun traj.xtc -nt %d >> gro.log 2>&1",parm->grobindir,parm->nthreads);
		if (parm->verbose) fprintf(stderr,"%s\n",command);
		system(command);

		// get energies
		sprintf(command,"echo '%s' | %s/g_energy -f ener.edr >> gro.log 2>&1",parm->energytype,parm->grobindir);
		if (parm->verbose) fprintf(stderr,"%s\n",command);
		system(command);
		if( access( "energy.xvg", F_OK ) == -1 ) Error("Cannot create energy.xvg");
		n = ParseXVG("energy.xvg",enew);
		if (n != nframes) Error("n != nframes when calculating energies");
                                                                                           
		// reweight restrains
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

		if (parm->verbose)
		{
			fprintf(stderr,"Reweigthed restrains:\n\tid\t\tsimulated +/- sigma\t\texp +/- sigma\t\t\tchi\n");
			for (i=0;i<parm->nres;i++)
                if ( parm->restype[i] != 5 )
                   fprintf(stderr,"\t%3d\t%8d %8d\t%10lf +/- %10lf\t%10lf +/- %10lf\t%f\n",
								i,parm->ires1[i],parm->ires2[i],k*resnew[i],
								k*sqrt(resnew2[i]-resnew[i]*resnew[i]),
								parm->resexp[i],parm->ressigma[i],
								(float)(k*resnew[i]-parm->resexp[i])/parm->ressigma[i]);
                else
                    fprintf(stderr,"\t%3d\t%8s %8s\t%10lf +/- %10lf\t%10lf +/- %10lf\t%f\n",
                            i,parm->alabel1[i],parm->alabel2[i],k*resnew[i],
                            k*sqrt(resnew2[i]-resnew[i]*resnew[i]),
                            parm->resexp[i],parm->ressigma[i],
                            (float)(k*resnew[i]-parm->resexp[i])/parm->ressigma[i]);
                    
                
			if  (!parm->correl)
			{
				fprintf(stderr,"attempted chi2=%lf",chi2/parm->nres);
				if(parm->rescalenoe0) fprintf(stderr,"\tnoe0=%lf (rescale=%lf)",k*parm->noe0,k);
				fprintf(stderr,"\n");
			}
			else
				fprintf(stderr,"attempted r=%lf\n",1.-chi2);
		}
		if (chi2 < chi2old)
		{
			chi2old = chi2;
			CopyTop(topnew,top);
			if (parm->verbose) fprintf(stderr,"Accepted.\n");
			iacc++;
		}
		else if  (parm->verbose) fprintf(stderr,"Rejected.\n");

		if  (!parm->correl)
		{
			fprintf(stderr,"\t%d\tchi2 = %lf",it,chi2old/parm->nres);
			if(parm->rescalenoe0) fprintf(stderr,"\tnoe0=%lf (rescale=%lf)",k*parm->noe0,k);
			fprintf(stderr,"\n");

		}
		else
			fprintf(stderr,"\t%d\tr = %lf\n",it,1.-chi2old);

		system("rm -f \\#*");
	}

	if  (!parm->correl)
		fprintf(stderr,"\tNew chi2 = %lf\n\tAccepted changes = %d / %d\n",chi2old/parm->nres,iacc,parm->noptim);
	else
		fprintf(stderr,"\tNew r = %lf\n\tAccepted changes = %d / %d\n",1.-chi2old,iacc,parm->noptim);

	return chi2old;
}

void ReweightRestrain(double *resnew, double **restrain, double *eold, double *enew, struct parm_s *parm, int nframes, int iskip, double *resnew2)
{
	int i,ir;
	double boltz,z=0,emax=-99999;

	for (ir=0;ir<parm->nres;ir++) { resnew[ir] = 0; resnew2[ir] = 0; }

	// calculate emax
	for (i=iskip;i<nframes;i++)
		if ( (-enew[i]/parm->temp + eold[i]/parm->temp) > emax) 
					emax = -enew[i]/parm->temp + eold[i]/parm->temp;

	// reweight	
	for (i=iskip;i<nframes;i++)
	{
		boltz = exp(-enew[i]/parm->temp + eold[i]/parm->temp - emax);

		for (ir=0;ir<parm->nres;ir++) 	
		{
			resnew[ir] += restrain[ir][i] * boltz;
			resnew2[ir] += restrain[ir][i] * restrain[ir][i] * boltz;
                }                                                                                                           
		z += boltz;
	}

	for (i=0;i<parm->nres;i++) 
	{
		resnew[i] /= z;
		resnew2[i] /= z;
	}
}

void AverageRestrain(double *resnew, double **restrain, struct parm_s *parm, int nframes)
{
	int i,ir;
	double z=0;

	for (i=0;i<nframes;i++)
	{
		for (ir=0;ir<parm->nres;ir++) 	
		{
			resnew[ir] += restrain[ir][i];
                }                                                                                                           
		z ++;
	}

	for (i=0;i<parm->nres;i++) 
	{
		resnew[i] /= z;
	}

}

int IsHydrogen(char *x)
{
	if (x[0] == 'H') return 1;
	if (x[0] == 'h') return 1;
	return 0;
}

void RandomChangeEnergy(struct parm_s *parm, struct top_s *top, struct top_s *origtop)
{
	int change=0,ic,iw,i1,i2,i3,i4,isH,cnt=0,funct;
	double delta1,delta2,c6,c12,u,r,c6new,c12new,unew,rnew,c0,c1,c2,c3,c4,c5;

	do 
	{
		ic = irand(2);
		change = 0;
		if (ic==0 && parm->ch_pairs)		// change pairs
		{
			iw = irand(top->npairs); 
			i1 = ((top->top_pairs)+iw)->i1 ;
			i2 = ((top->top_pairs)+iw)->i2 ;
			isH=0;
			if ( IsHydrogen(((top->top_atoms)+i1 -1 )->type) || IsHydrogen(((top->top_atoms)+ i2 - 1 )->type) ) isH=1;	

			if (parm->combrule==1)
			{
				if ( !parm->nochangeH || !isH ) 	// maybe you don't want to change H
				if ( ((top->top_pairs)+iw)->move == 1 || ( parm->moveAll && ((top->top_pairs)+iw)->move > -1 )
								 || ( parm->moveAllPairs && ((top->top_pairs)+iw)->move > -1 ) )
				{
					c6 = ((top->top_pairs)+iw)->c0;
					c12 = ((top->top_pairs)+iw)->c1;
					if ( IsZero(c6) || IsZero(c12) ) { u=0; r=0.1; }			// if c6==0 or c12==0 r and u not defined (pair energy is 0)
					else
					{
						u = -0.25 * c6 * c6 / c12;
						if ( !parm->RnochangeH || !isH) r = pow(2.*c12/c6,1./6); 			// if RnochangeH defined, change only energy
					}
					do { delta1 = (frand()-0.5) * parm->sigma_pairs_u;} while (u+delta1 > 0 );		// u cannot be positive
					do { delta2 = (frand()-0.5) * parm->sigma_pairs_r;} while ( r+delta2 <= 0 );		// r cannot be negative

					unew = u + delta1;
					rnew = r + delta2;
					c6new = -2. * pow(rnew,6.) * unew;
					c12new = - pow(rnew,12.) * unew;

					if (parm->onlyc6) c12new = c12;

					if (parm->elimit0>0) 
						if ( fabs( c6 - ((origtop->top_pairs)+iw)->c0 ) > parm->elimit0 ) change = -1;
					if (parm->elimit1>0) 
						if ( fabs( c12 - ((origtop->top_pairs)+iw)->c1 ) > parm->elimit1 ) change = -1;
							
					if (change>-1)
					{
						((top->top_pairs)+iw)->c0 = c6new;
						((top->top_pairs)+iw)->c1 = c12new;
						if (parm->verbose) fprintf(stderr,"change pair energy %d (%d-%d;%s-%s) u=%e r=%e (c6=%e c12=%e) -> u=%e r=%e (c6=%e c12=%e)\n"
								,iw,i1,
								i2, ((top->top_atoms)+ i1 - 1 )->type, ((top->top_atoms)+ i2 - 1 )->type,u,r,c6,c12,unew,rnew,c6new,c12new);
					}
					change = 1;
				}
			}
			else if  (parm->combrule==2)
			{
				if ( !parm->nochangeH || !isH ) 	// maybe you don't want to change H
				if ( ((top->top_pairs)+iw)->move == 1 || ( parm->moveAll && ((top->top_pairs)+iw)->move > -1 )
								||	 ( parm->moveAllPairs && ((top->top_pairs)+iw)->move > -1 )	 )
				{
					r = ((top->top_pairs)+iw)->c0;
					u = ((top->top_pairs)+iw)->c1;
					do { delta1 = (frand()-0.5) * parm->sigma_pairs_u;} while ( u+delta1 < 0 );		// u cannot be positive
					if ( !parm->RnochangeH || !isH)								// if RnochangeH defined, change only energy
						do { delta2 = (frand()-0.5) * parm->sigma_pairs_r;} while ( r+delta2 <= 0 );	// r cannot be negative
					else delta2 = 0;
				
					if (parm->onlyc6) delta2 = r * pow(u/(u+delta1),0.08333) - r;			// c'_12=c_12 --> (r')^12*u'=r^12u --> r'=r(u/u')^1/12

	
					if (parm->elimit0>0) 
						if ( fabs( r + delta2 - ((origtop->top_pairs)+iw)->c0 ) > parm->elimit0 ) change = -1;
					if (parm->elimit1>0) 
						if ( fabs( u + delta1 - ((origtop->top_pairs)+iw)->c1 ) > parm->elimit1 ) change = -1;

					if (change>-1)
					{
						((top->top_pairs)+iw)->c0 = r + delta2;
						((top->top_pairs)+iw)->c1 = u + delta1;
						if (parm->verbose) fprintf(stderr,"change pair energy %d (%d-%d;%s-%s) u=%e r=%e -> u=%e r=%e \n",iw,i1,
							i2, ((top->top_atoms)+ i1 - 1 )->type, ((top->top_atoms)+ i2 - 1 )->type,u,r,u + delta1,r+delta2);
					}
					change = 1;
				}
					
			}
			else Error("Unknown combination rule for LJ");


		}

		if (ic==1 && parm->ch_dihedrals)	// change dihedrals
		{
			iw = irand(top->ndihedrals);
			i1 = ((top->top_dihedrals)+iw)->i1;
			i2 = ((top->top_dihedrals)+iw)->i2;
			i3 = ((top->top_dihedrals)+iw)->i3;
			i4 = ((top->top_dihedrals)+iw)->i4;
			c0 = ((top->top_dihedrals)+iw)->c0;
			c1 = ((top->top_dihedrals)+iw)->c1;
			c2 = ((top->top_dihedrals)+iw)->c2;
			c3 = ((top->top_dihedrals)+iw)->c3;
			c4 = ((top->top_dihedrals)+iw)->c4;
			c5 = ((top->top_dihedrals)+iw)->c5;
			funct = ((top->top_dihedrals)+iw)->funct;

			if ( ((top->top_dihedrals)+iw)->move == 1 || ( parm->moveAll && ((top->top_dihedrals)+iw)->move > -1 ) 	
							||	 ( parm->moveAllDih && ((top->top_dihedrals)+iw)->move > -1 )		)
			if ( funct == 3 || funct == 1 || funct == 9 )	// Rychaert-Bellemans or periodic
			{
				if ( !(parm->nodih0change) || c0 < -EPSILON || c0 > EPSILON )
				{ 
					delta1 = (frand()-0.5) * parm->sigma_dihe;
					((top->top_dihedrals)+iw)->c0 += delta1;
					change = 1;
				}
				if ( !(parm->nodih0change) || c1< -EPSILON || c1 > EPSILON )
				{ 
					delta1 = (frand()-0.5) * parm->sigma_dihe;
					((top->top_dihedrals)+iw)->c1 += delta1;
					change = 1;
				}

				if ( funct == 3 )  							//  Rychaert-Bellemans
				{
					if ( !(parm->nodih0change) || c2< -EPSILON || c2 > EPSILON )
					{ 
						delta1 = (frand()-0.5) * parm->sigma_dihe;
						((top->top_dihedrals)+iw)->c2 += delta1;
						change = 1;
					}
					if ( !(parm->nodih0change) || c3< -EPSILON || c3 > EPSILON )
					{ 
						delta1 = (frand()-0.5) * parm->sigma_dihe;
						((top->top_dihedrals)+iw)->c3 += delta1;
						change = 1;
					}
					if ( !(parm->nodih0change) || c4< -EPSILON || c4 > EPSILON )
					{ 
						delta1 = (frand()-0.5) * parm->sigma_dihe;
						((top->top_dihedrals)+iw)->c4 += delta1;
						change = 1;
					}
					if ( !(parm->nodih0change) || c5< -EPSILON || c5 > EPSILON )
					{ 
						delta1 = (frand()-0.5) * parm->sigma_dihe;
						((top->top_dihedrals)+iw)->c5 += delta1;
						change = 1;
					}
				}

				if (change==1) 
				{
					if (parm->verbose && ( funct==3 || funct==9) ) 
						fprintf(stderr,"change dihedral energy %d (%d-%d-%d-%d;%s-%s-%s-%s) %e %e %e %e %e %e -> %e %e %e %e %e %e\n",
								iw,i1,i2,i3,i4,((top->top_atoms)+i1-1)->type,((top->top_atoms)+i2-1)->type,
								((top->top_atoms)+i3-1)->type,((top->top_atoms)+i4-1)->type, c0,c1,c2,c3,c4,c5,
								((top->top_dihedrals)+iw)->c0,((top->top_dihedrals)+iw)->c1,((top->top_dihedrals)+iw)->c2, 
								((top->top_dihedrals)+iw)->c3,((top->top_dihedrals)+iw)->c4,((top->top_dihedrals)+iw)->c5);
					if (parm->verbose && funct==1)
						 fprintf(stderr,"change dihedral energy %d (%d-%d-%d-%d;%s-%s-%s-%s) %e %e -> %e %e \n",
								iw,i1,i2,i3,i4,((top->top_atoms)+i1-1)->type,((top->top_atoms)+i2-1)->type,
								((top->top_atoms)+i3-1)->type,((top->top_atoms)+i4-1)->type, c0,c1,
								((top->top_dihedrals)+iw)->c0,((top->top_dihedrals)+iw)->c1);

				}
			}
		}

		cnt++;
		if (cnt>9999) Error("Cannot make energy change");

	} while (change<1);

	// write new top to file
	PrintTop("tmp.top", top, parm->verbose, parm->noSOL);
}

void GetRestrain(struct parm_s *parm, double **restrain, int nframes, struct NOErelax_s *NOErelax)
{
	int i,i1,i2,n,j,k,chk,l;
	double av,av2,dmin;
	char command[500];

	if (parm->verbose) fprintf(stderr,"Calculate restraints.\n");

    for (i=0;i<nframes;i++)
        for (j=0;j<parm->nres;j++) restrain[j][i] = 0;
    
	// get restrains
	for (i=0;i<parm->nres;i++) 
	{
		// = distances
		if ( parm->restype[i] == 1 )
		{
			i1 = parm->ires1[i];
			i2 = parm->ires2[i];
			sprintf(command,"echo 'A%d A%d' | %s/g_dist -f traj.xtc -s topol.tpr -n index.ndx >> gro.log 2>&1",
						i1,i2,parm->grobindir);
			system(command);
			if (parm->verbose) fprintf(stderr,"%s\n",command);
			if( access( "dist.xvg", F_OK ) == -1 ) Error("Cannot create dist.xvg");
			n = ParseXVG("dist.xvg",restrain[i]);
			if (n != nframes) Error("n != nframes when calculating restrains");
		}
		// = noe0 / r^6
		else if ( parm->restype[i] == 2 )
		{
			i1 = parm->ires1[i];
			i2 = parm->ires2[i];
			sprintf(command,"echo 'A%d A%d' | %s/g_dist -f traj.xtc -s topol.tpr -n index.ndx >> gro.log 2>&1",
						i1,i2,parm->grobindir);
			system(command);
			if (parm->verbose) fprintf(stderr,"%s\n",command);
			if( access( "dist.xvg", F_OK ) == -1 ) Error("Cannot create dist.xvg");
			n = ParseXVG("dist.xvg",restrain[i]);
			av=0; av2=0; dmin=9999.;
			for (j=0;j<n;j++) 
			{
				av += restrain[i][j] / n;
				av2 += restrain[i][j]*restrain[i][j] / n;
				if (restrain[i][j]<dmin) dmin=restrain[i][j];
				restrain[i][j] = parm->noe0 / pow(restrain[i][j], 6. );
			}
			if (parm->verbose) fprintf(stderr," Average = %lf sigma = %lf min=%lf\n",av,sqrt(av2-av*av),dmin );
			if (n != nframes) Error("n != nframes when calculating restrains");
		}
        // = noe0 / r^6 with groups
        else if ( parm->restype[i] == 4 )
        {
            chk=0;
            for (l=0;l<parm->nagroup;l++)
                for (k=0;k<parm->nagroup;k++)
                {
                    if ( !strcmp(parm->alabel1[i], parm->agrouplabel[l]) &&
                        !strcmp(parm->alabel2[i], parm->agrouplabel[k]) )
                    {
                        i1 = parm->agroup[l];
                        i2 = parm->agroup[k];
                        sprintf(command,"echo 'A%d A%d' | %s/g_dist -f traj.xtc -s topol.tpr -n index.ndx >> gro.log 2>&1",
                                i1,i2,parm->grobindir);
                        system(command);
                        if (parm->verbose) fprintf(stderr,"%s\n",command);
                        if( access( "dist.xvg", F_OK ) == -1 ) Error("Cannot create dist.xvg");
                        n = ParseXVG("dist.xvg",restrain[i]);
                        
                        av=0; av2=0; dmin=9999.;
                        for (j=0;j<n;j++)
                        {
                            av += restrain[i][j] / n;
                            av2 += restrain[i][j]*restrain[i][j] / n;
                            if (restrain[i][j]<dmin) dmin=restrain[i][j];
                            restrain[i][j] += parm->noe0 / pow(restrain[i][j], 6. );
                        }
                        if (parm->verbose) fprintf(stderr," Average = %lf sigma = %lf min=%lf\n",av,sqrt(av2-av*av),dmin );
                        if (n != nframes) Error("n != nframes when calculating restrains");
                       
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
                    
	// full relaxation NOE
	if (NOErelax)
	{
		if (parm->verbose) fprintf(stderr,"Calculate restraints through full optimization.\n");
		for (i=0;i<NOErelax->nh;i++)
			for (j=i+1;j<NOErelax->nh;j++)
			{
				i1 = NOErelax->Hlist[i];
				i2 = NOErelax->Hlist[j];
				sprintf(command,"echo 'A%d A%d' | %s/g_dist -f traj.xtc -s topol.tpr -n index.ndx >> gro.log 2>&1",
														i1,i2,parm->grobindir);
				system(command);
				if (parm->verbose) fprintf(stderr,"%s\n",command);
				if( access( "dist.xvg", F_OK ) == -1 ) Error("Cannot create dist.xvg");
				n = ParseXVGf("dist.xvg",NOErelax->dist[i][j]);
				for (k=0;k<nframes;k++) NOErelax->dist[j][i][k] = NOErelax->dist[i][j][k];
				if (n != nframes) Error("n != nframes when calculating restrains");
			}

		FullRelaxation(NOErelax,restrain, parm, nframes);
	}


}

void PrintInitialRestrains(struct parm_s *parm, double **restrain, FILE *fout,struct NOErelax_s *NOErelax)
{
	int i;
	char command[500];

	if (parm->noSOL)
		sprintf(command,"%s/grompp -f tmp.mdp -p tmp.top -c %s -n >> gro.log 2>&1",parm->grobindir,parm->groprotfile);
	else
		sprintf(command,"%s/grompp -f %s -p tmp.top -c %s -n >> gro.log 2>&1",parm->grobindir,parm->mdpfile,parm->grofile);
	
	if (parm->verbose) fprintf(stderr,"%s\n",command);
	system(command);

	sprintf(command,"%s/trjconv -f %s -o traj.xtc  >> gro.log 2>&1",parm->grobindir,parm->grofile);
	system(command);
	if (parm->verbose) fprintf(stderr,"%s\n",command);

	GetRestrain(parm,restrain,1,NOErelax);

	fprintf(fout,"Initial restrains:\n");
	for (i=0;i<parm->nres;i++) 
		fprintf(fout,"%d\t%lf\t%lf %lf\n",i,restrain[i][0],parm->resexp[i],parm->ressigma[i]);	
	
}

void Welcome(FILE *fp)
{
	fprintf(fp,"\n\n***************************************\n");
	fprintf(fp,"*        ffoptim    v. %d.%d            *\n",NVER,NSUBVER);
	fprintf(fp,"***************************************\n");
	fprintf(fp,"  (c) G. Tiana, 2013\n\n");

	fprintf(fp,"launched: ");
	PrintTime(fp);	
	fprintf(fp,"pid = %d\n",getpid());
}


void PrintTime(FILE *fout)
{
	char buff[100];
   	time_t now = time (0);
   	strftime (buff, 100, "%d/%m/%Y %H:%M:%S", localtime (&now));
   	fprintf (fout,"%s\n", buff);

}

void Help(void)
{
	fprintf(stderr,"\n");
	fprintf(stderr,"Compulsory in the input file:\n");
	fprintf(stderr,"topfile <string>\tThe name of the file containing the initial topology. Comment with MOVE terms to be updated.\n");
	fprintf(stderr,"grobindir <string>\tThe directory containing gromacs binaries.\n");
	fprintf(stderr,"grofile <string>\tThe name of the file containing the initial gro conformation.\n");
	fprintf(stderr,"mdpfile <string>\tThe name of the mpd file for the MD simulation\n");
	fprintf(stderr,"emdpfile <string>\tThe name of the mpd file for the energy minimization\n");
	fprintf(stderr,"origtopfile <string>\tOriginal top file to implement elimit0 and elimit1\n");
	fprintf(stderr,"ch_pairs | ch_dihedrals\tOptimize the potential by changing pairs or dihedral energies\n");
	fprintf(stderr,"niter <int>\t\tNumber of optimization iterations\n");
	fprintf(stderr,"noptim <int>\t\tNumber of optimization steps after each MD iteration\n");
	fprintf(stderr,"temp <double>\t\tTemperature associated with the experimental data (this is not the\n");
	fprintf(stderr,"\t\t\tsimulation temperature, which is in the mdp file).\n");

	fprintf(stderr,"Optional:\n");
	fprintf(stderr,"itpdir <string>\t\tThe directory containing gromacs itp topologies\n");
	fprintf(stderr,"groprotfile <string>\tIf solvent is present, gro file of the solute alone\n");
	fprintf(stderr,"chi2file <string>\tOutput file to record chi2.\n");
	fprintf(stderr,"verbose\t\t\tPrint all commands\n");
	fprintf(stderr,"boltz <double>\t\tBoltzmann's constant (default is 0.00833 kJ/mol/K)\n");
	fprintf(stderr,"sigma_pairs_u <double>\tMaximum change for single pair energy\n");
	fprintf(stderr,"sigma_pairs_r <double>\tMaximum change for single pair range\n");
	fprintf(stderr,"sigma_dihe <double>\tMaximum change for dihedral parameter\n");
	fprintf(stderr,"noe0 <double>\t\tScaling parameter for NoEs (i.e., I = noe0 / r^6) [r in nm]\n");
	fprintf(stderr,"nthreads <int>\t\tNumber of threads for mdrun\n");
	fprintf(stderr,"correl\t\t\tOptimize correlation function instead of chi2.\n");
	fprintf(stderr,"replex <int>\t\tMake a replica exchange instead of plain MD (integer is exchange time).\n");
	fprintf(stderr,"combrule <int>\t\tUse combination rule in topol.top\n");
	fprintf(stderr,"restart <int>\t\tRestart from given iteration.\n");
	fprintf(stderr,"tauc <double>\t\tCorrelation time in picoseconds for isotropic thumbling.\n");
	fprintf(stderr,"omega <double>\t\tSpectrometer frequency in [ps]^-1.\n");
	fprintf(stderr,"mixingt <double>\tMixing time [ms].\n");
	fprintf(stderr,"nframesmax <int>\tmaximum number of frames for complete relaxation\n");
	fprintf(stderr,"nochangeH\t\tdo not update interaction if one of the involved atoms is H.\n");
	fprintf(stderr,"RnochangeH\t\tdo not update R if one of the involved atoms is H.\n");
	fprintf(stderr,"iskip <int>\t\tIn each simulation skip the first frames\n");
	fprintf(stderr,"elimit0 <double>\tStrict boundaries on variatoon of c0 in top\n");
	fprintf(stderr,"elimit1 <double>\tStrict boundaries on variatoon of c1 in top\n");
	fprintf(stderr,"energytype <string>\tWhich energy item of g_energy to use for the rescaling (default: Potential).\n");
	fprintf(stderr,"rescalenoe0\t\tcalculate chi2 rescaling intensities by the best factor.\n");
	fprintf(stderr,"onlyc6\t\t\tChange only the c6 term of the energy.\n");
	fprintf(stderr,"noSOL\t\t\tExlcude solvent from calculating energies (needs following 3 keys)\n");
	fprintf(stderr,"mdpprotxtcgrps\t\tname of the molecule group to write xtc\n");
	fprintf(stderr,"mdpprotegrps\t\tname of the molecule group to write energy\t\n");
	fprintf(stderr,"groprotfile\t\tgro file of the molecule only\n");
	fprintf(stderr,"nodih0change\t\tDo not change dihedral terms which have 0 energy\n");
	fprintf(stderr,"testtraj <string>\tRead a trajectory in 1st iteration instead of making a simulation\n");
	fprintf(stderr,"moveall\t\t\tmove all energy terms except those marked with a NOC\n");
	fprintf(stderr,"movePall\t\t\tmove all pairs energy terms except those marked with a NOC\n");
	fprintf(stderr,"moveDall\t\t\tmove all dihedrals energy terms except those marked with a NOC\n");
	fprintf(stderr,"nonbondedfile <string>\tFill the pair parameters in the top with values from itp\n");
	fprintf(stderr,"movev0\t\t\tmove all energy terms listed in the v0 array\n");
	fprintf(stderr,"movea0\t\t\tmove all energy terms of heavy atoms linked to those listed in the v0 array\n");
	fprintf(stderr,"lowermatrix\t\tif NOE relaxation enabled and NOE matrix not symmetric, use lower half\n");
	fprintf(stderr,"hostfile <string>\tFor parallel tempering, provides hostfile to mpirun\n");
	fprintf(stderr,"uNOE <double>\t\t\tAdd unobserved NOEs as zero restraints for each H\n");
	fprintf(stderr,"\n");

	exit(0);
}


