#include "ffoptim.h"
#include <dirent.h> 

struct parm_s *ReadParms(char *fname)
{
  struct parm_s *x;
  int getr=0,i1,i2,i3,nres=0,i,gett=0,getv0=0,gettau=0,getg=0;
  double f,f2,kboltz;
  char aux[1000],keyword[100],label[30],label2[30];
  FILE *fp;

  fprintf(stderr,"\nRead parameters:\n");

  fp = fopen(fname,"r");
  if (fp==NULL) Error2("Cannot open input file ",fname);

  x = calloc(1,sizeof(struct parm_s));
  if (!x) Error("Cannot allocate mc_parms");

  x->ires1 = AlloInt(NRESMAX);
  x->ires2 = AlloInt(NRESMAX);
  x->restype = AlloInt(NRESMAX);
  x->resexp = AlloDouble(NRESMAX);
  x->ressigma = AlloDouble(NRESMAX);
  x->alabel1 = AlloString(NRESMAX,30);
  x->alabel2 = AlloString(NRESMAX,30);
  x->nres = 0;

  x->iv0 = AlloInt(NRESMAX);
  x->v0 = AlloDouble(NATOMSMAX);
  x->nv0 = 0;
  x->nagroup = 0;
    
  x->ch_pairs = 0;
  x->ch_dihedrals = 0;
  kboltz = 0.00833;
  x->nthreads = 1;
  x->replex = 0;
  x->ntemp = 0;
  x->correl = 0;
  x->combrule = 1;
  x->noe0 = 1.;
  x->restart = 0;
  x->sigma_pairs_u = 0.5;
  x->sigma_pairs_r = 0.05;
  x->sigma_pairs_r = 1;
  x->nframesmax = 5000;
  x->tauc = 150; 
  x->omega = 0.00376;
  x->mixingt = 500.;
  x->absnoe = 0;
  x->nochangeH = 0;
  x->RnochangeH = 0;
  x->filltop = 0;
  x->iskip = 1;
  x->noSOL = 0;
  x->elimit0 = -1;
  x->elimit1 = -1;
  x->onlyc6 = 0;
  x->rescalenoe0 = 0;
  x->nice=10;
  x->nodih0change=0;
  x->sigma_dihe = 1.;
  x->moveAll = 0;
  x->movev0 = 0;
  x->movea0 = 0;
  strcpy(x->mdpprotxtcgrps,"none");
  strcpy(x->mdpprotegrps,"none");
  strcpy(x->groprotfile,"none");
  strcpy(x->energytype,"Potential");
  strcpy(x->chi2file,"none");
  strcpy(x->origtopfile,"none");
  strcpy(x->testtraj,"none");
  strcpy(x->nonbondedFile,"none");
  strcpy(x->bondedFile,"none");
  x->lowerMatrix = 0;
  x->ntaucatoms = 0;
  x->noptauc = 0;
  strcpy(x->hostfile,"");
  x->uNOE = -1;

  while ( fgets(aux,1000,fp) )
  {
	ReadParS(aux,"topfile",x->topfile);
	ReadParS(aux,"itpdir",x->itpdir);
	ReadParS(aux,"grobindir",x->grobindir);
	ReadParS(aux,"grofile",x->grofile);
	ReadParS(aux,"groprotfile",x->groprotfile);
	ReadParS(aux,"mdpfile",x->mdpfile);
	ReadParS(aux,"chi2file",x->chi2file);
	ReadParS(aux,"origtopfile",x->origtopfile);
	ReadParS(aux,"emdpfile",x->emdpfile);
	ReadParD(aux,"verbose",&(x->verbose));
	ReadParN(aux,"ch_pairs",&(x->ch_pairs));
	ReadParN(aux,"ch_dihedrals",&(x->ch_dihedrals));
	ReadParD(aux,"niter",&(x->niter));
	ReadParD(aux,"noptim",&(x->noptim));
	ReadParD(aux,"noptauc",&(x->noptauc));
	ReadParF(aux,"temp",&(x->temp));
	ReadParF(aux,"kboltz",&kboltz);
	ReadParF(aux,"sigma_pairs_u",&(x->sigma_pairs_u));
	ReadParF(aux,"sigma_pairs_r",&(x->sigma_pairs_r));
	ReadParF(aux,"sigma_dihe",&(x->sigma_dihe));
	ReadParF(aux,"noe0",&(x->noe0));
	ReadParD(aux,"nthreads",&(x->nthreads));
	ReadParN(aux,"correl",&(x->correl));
	ReadParD(aux,"replex",&(x->replex));
	ReadParD(aux,"combrule",&(x->combrule));
	ReadParD(aux,"restart",&(x->restart));
	ReadParD(aux,"nframesmax",&(x->nframesmax));
	ReadParF(aux,"tauc",&(x->tauc));
	ReadParF(aux,"omega",&(x->omega));
	ReadParF(aux,"mixingt",&(x->mixingt));
	ReadParF(aux,"elimit0",&(x->elimit0));
	ReadParF(aux,"elimit1",&(x->elimit1));
	ReadParN(aux,"absnoe",&(x->absnoe));
	ReadParN(aux,"nochangeH",&(x->nochangeH));
	ReadParN(aux,"RnochangeH",&(x->nochangeH));
	ReadParN(aux,"filltop",&(x->filltop));
	ReadParD(aux,"iskip",&(x->iskip));
	ReadParN(aux,"noSOL",&(x->noSOL));
	ReadParS(aux,"mdpprotegrps",x->mdpprotegrps);
	ReadParS(aux,"mdpprotxtcgrps",x->mdpprotxtcgrps);
	ReadParS(aux,"energytype",x->energytype);
	ReadParN(aux,"onlyc6",&(x->onlyc6));
	ReadParN(aux,"rescalenoe0",&(x->rescalenoe0));
	ReadParN(aux,"nodih0change",&(x->nodih0change));
	ReadParD(aux,"nice",&(x->nice));
	ReadParS(aux,"testtraj",x->testtraj);
	ReadParN(aux,"movePall",&(x->moveAllPairs));
	ReadParN(aux,"moveDall",&(x->moveAllDih));
	ReadParS(aux,"nonbondedfile",x->nonbondedFile);
	ReadParS(aux,"bondedfile",x->bondedFile);
	ReadParN(aux,"movev0",&(x->movev0));
	ReadParN(aux,"movea0",&(x->movea0));
	ReadParN(aux,"lowermatrix",&(x->lowerMatrix));
	ReadParS(aux,"hostfile",x->hostfile);
	ReadParF(aux,"uNOE",&(x->uNOE));


	if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"restrains") ) getr = 1;
	if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"endrestrains") ) getr = 0;
	if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"temperatures") ) gett = 1;
	if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"endtemperatures") ) gett = 0;
	if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"v0") ) getv0 = 1;
	if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"endv0") ) getv0 = 0;
	if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"taucatoms") ) gettau = 1;
	if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"endtaucatoms") ) gettau = 0;
    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"atomgroups") ) getg = 1;
    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"endatomgroups") ) getg = 0;
	
	if (getr)
	{
		if ( sscanf(aux,"%d %d %d %lf %lf",&i1,&i2,&i3,&f,&f2) == 5 )
		{
            if (i3 == 5 || i3==4) fprintf(stderr,"WARNING: numeric atom id will be interpreted as group name\n");
            
			x->ires1[nres]=i1;	
			x->ires2[nres]=i2;	
			x->restype[nres]=i3;	
			x->resexp[nres]=f;
			x->ressigma[nres]=f2;

			if (f2<EPSILON)
            {
                fprintf(stderr,"WARNING: problem in reading %d %d ...\n",i1,i2);
                Error("Negative or zero standard deviation for restrain");
            }
			nres++;	
			if (nres>=NRESMAX) Error("NRESMAX too small");
		}
        else if ( sscanf(aux,"%s %s %d %lf %lf",label,label2,&i3,&f,&f2) == 5 )
        {
            if (i3 != 5 && i3!=4)
            {
                fprintf(stderr,"In reading line %s %s ...\n",label,label2);
                Error("String label in restraint requires type 4 or 5");
            }
            strcpy( x->alabel1[nres], label );
            strcpy( x->alabel2[nres], label2 );
            x->restype[nres]=i3;
            x->resexp[nres]=f;
            x->ressigma[nres]=f2;
            
            if (f2<EPSILON)
            {
                fprintf(stderr,"WARNING: problem in reading %s %s ...\n",label,label2);
                Error("Negative or zero standard deviation for restrain");
            }
            nres++;
            if (nres>=NRESMAX) Error("NRESMAX too small");
        }
	}
	
	if (gett)
	{
		if ( sscanf(aux,"%lf",&f) == 1 )
		{
			x->temperatures[x->ntemp] = f;
			(x->ntemp)++;
		}	

	}
	if (getv0)
	{
		if ( sscanf(aux,"%d %lf",&i1,&f) == 2 )
		{
			x->iv0[x->nv0] = i1;
			x->v0[x->nv0] = f;
			(x->nv0)++;
		}	
	}
	if (gettau)
	{
		if ( x->ntaucatoms == 0 ) 
		{
			if ( x->nv0 == 0 ) Error("Must read v0 before taucatoms"); 
			x->taucatoms1 = AlloInt(x->nv0 * x->nv0);
			x->taucatoms2 = AlloInt(x->nv0 * x->nv0);
			x->taucatomsv = AlloDouble(x->nv0 * x->nv0);
		}

		if ( sscanf(aux,"%d %d %lf",&i1,&i2,&f) == 3 )
		{
			x->taucatoms1[ x->ntaucatoms ] = i1;
			x->taucatoms2[ x->ntaucatoms ] = i2;
			x->taucatomsv[ x->ntaucatoms ] = f;
		        x->ntaucatoms ++;
			if ( x->ntaucatoms >= x->nv0 * x->nv0 ) Error("ntaucatoms too large");
		}	
	}
    if (getg)
    {
        if ( sscanf(aux,"%s %d",label,&i1) == 2 )
        {
            strcpy( x->agrouplabel[ x->nagroup ] , label);
            x->agroup[ x->nagroup ] = i1;
            x->nagroup ++;
            if ( x->nagroup > NGROUPMAX) Error("NGROUPMAX too small");
        }
    
    }


  }


  if (nres==0) Error("Cannot read experimental restrains");
  x->nres = nres;


  x->temp *= kboltz;

  if (!strcmp(x->groprotfile,"none")) strcpy(x->groprotfile,x->grofile);

  if (x->verbose)
  {
	for (i=0;i<x->nres;i++)
        if (x->restype[i]!=5)
            fprintf(stderr,"restrain %d : %d %d	%d	%lf +/- %lf\n",i,x->ires1[i],x->ires2[i],x->restype[i],x->resexp[i],x->ressigma[i]);
        else
            fprintf(stderr,"restrain %d : %s %s    %d    %lf +/- %lf\n",i,x->alabel1[i],x->alabel2[i],x->restype[i],x->resexp[i],x->ressigma[i]);

	if (x->nv0)
	{
		for (i=0;i<x->nv0;i++)
			fprintf(stderr,"diagonal peak %d:  %d     %lf\n",i,x->iv0[i],x->v0[i]);
	}

	if (x->ntaucatoms)
		for (i=0;i<x->ntaucatoms;i++)
			fprintf(stderr,"taucatoms %d: %d-%d       %lf\n",i,x->taucatoms1[i],x->taucatoms2[i],x->taucatomsv[i]);

	if (x->replex)
	{
		fprintf(stderr,"** REPLICA EXCHANGE **\n");
		for (i=0;i<x->ntemp;i++)
			fprintf(stderr,"temperature %d : %lf\n",i,x->temperatures[i]);
	}
    if (x->nagroup>0)
    {
        fprintf(stderr,"Atom groups:\n");
        for (i=0;i<x->nagroup;i++)
            fprintf(stderr,"%s\t%d\n",x->agrouplabel[i],x->agroup[i]);
    }
  }

  fprintf(stderr,"\n");
  return x;
}

void ReadParD(char *s, char key[20], int *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %d",par);
		if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file (r=%d)\n",key,r); exit(1); }
		fprintf(stderr,"%s = %d\n",key,*par);
	}

}

void ReadParL(char *s, char key[20], long *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %ld",par);
		if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file (r=%d)\n",key,r); exit(1); }
		fprintf(stderr,"%s = %ld\n",key,*par);
	}

}

void ReadParF(char *s, char key[20], double *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %lf",par);
		if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file\n",key); exit(1); }
		fprintf(stderr,"%s = %lf\n",key,*par);
	}
}

void ReadParS(char *s, char key[20], char *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %s",par);
		if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file\n",key); exit(1); }
		fprintf(stderr,"%s = %s\n",key,par);
	}
}

void ReadParN(char *s, char key[20], int *par)
{
	int l;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		*par=1;
		fprintf(stderr,"%s\n",key);
	}
}
/********************************************************************
 Given a string, returns the content of square brackets, excluding
 spaces and returns 1,
 if there are no square brackets returns 0
 ********************************************************************/
int FindKeyword(char *string, char *keyword)
{
   int i=0,key=0,k=0;
   char c;

   do {
		   c = string[i];
		   if (c=='[') key=1;
		   if (c==']') key=0;
		   if (key==1 && c!=' ' && c!='[')
		   {
			   keyword[k] = c;
			   k++;
		   }
		   i++;

	   } while ( c != '\0' && i<500 );

	keyword[k] = '\0';

	if (!strcmp(keyword,"")) return 0;
	else return 1;
}

void CheckFiles(struct parm_s *x)
{
	int a,nstenergy=-1,nstxtcout=-2;
	char aux[500];
	FILE *fp;

	if( access( x->grofile, F_OK ) == -1 ) Error2("Cannot find file gro file ",x->grofile); 
	if( access( x->topfile, F_OK ) == -1 ) Error2("Cannot find file top file ",x->topfile); 
	if( access( x->mdpfile, F_OK ) == -1 ) Error2("Cannot find file mdp file ",x->mdpfile); 
	if( access( x->emdpfile, F_OK ) == -1 ) Error2("Cannot find file emdp file ",x->emdpfile); 

	system("rm -f topol.tpr gro.log index.ndx");
	sprintf(aux,"sed 's/^ *//;s/ *$//;s/ \\{1,\\}/ /g' %s > md.tmp",x->mdpfile);        
	system(aux);

	fp = fopen("md.tmp","r");
	if (!fp) Error2("Cannot open file ",x->mdpfile);
	while ( fgets(aux,500,fp) )
	{
		if ( sscanf(aux,"nstenergy = %d",&a) == 1 ) nstenergy = a;	
		if ( sscanf(aux,"nstxtcout = %d",&a) == 1 ) nstxtcout = a;	
	}
	if ( nstenergy != nstxtcout ) Error("nstenergy must be equal to nstxtcout in mdp file");
	fclose(fp);
	system("rm -f md.tmp");
}

int ParseXVG(char *fname, double *x)
{
	int i=0;
	double a1,a2;
	FILE *fp;
	char aux[500];

	fp = fopen(fname,"r");
	if (!fp) Error2("Cannot open ",fname);

	while ( fgets(aux,500,fp) )
  	{
		if ( aux[0] != '#' && aux[0] != '@')
			if ( sscanf(aux,"%lf %lf",&a1,&a2) == 2 )
			{
				x[i] = a2;
				i++;
				if (i>=NFRAMESMAX) Error("NFRAMESMAX too small");
			}
  	}

	if (i==0) Error2("No lines read in ",fname);
	fclose(fp);
	
	sprintf(aux,"rm -f %s",fname);
	system(aux);

	return i;
}

int ParseXVGf(char *fname, float *x)
{
	int i=0;
	float a1,a2;
	FILE *fp;
	char aux[500];

	fp = fopen(fname,"r");
	if (!fp) Error2("Cannot open ",fname);

	while ( fgets(aux,500,fp) )
  	{
		if ( aux[0] != '#' && aux[0] != '@')
			if ( sscanf(aux,"%f %f",&a1,&a2) == 2 )
			{
				x[i] = a2;
				i++;
				if (i>=NFRAMESMAX) Error("NFRAMESMAX too small");
			}
  	}

	if (i==0) Error2("No lines read in ",fname);
	fclose(fp);

	sprintf(aux,"rm -f %s",fname);
	system(aux);

	return i;
}


void RenameIndex(void)
{
  system("mv index.ndx index.tmp; sed -e 's/.*_.*_/[ A/' index.tmp > index.ndx; rm -f index.tmp");
}

void ChangeTmdp(char *mdpfile, int i, double *t)
{
	char aux[500],nout[500];
	FILE *fin,*fout;

	fin=fopen(mdpfile,"r");
	if (!fin) Error("Cannot open mdp file in ChangeTmdp");
	sprintf(nout,"tmp%d.mdp",i);
	fout=fopen(nout,"w");
	if (!fout) Error("Cannot open mdp file for output in ChangeTmdp");

	while ( fgets(aux,500,fin) )
	{
		if ( !strncmp(aux,"ref_t",5) ) fprintf(fout,"ref_t\t= %lf\n",t[i]);
		else  fprintf(fout,"%s",aux);
	}


	fclose(fin);
	fclose(fout);
}


void CheckExplosion(void)
{
  DIR *d;
  struct dirent *dir;

  d = opendir(".");
  if (d)
  {
    while ((dir = readdir(d)) != NULL)
    {
      if (!strncmp(dir->d_name,"step",4))
	{	
		Error("System has exploded");
	}
    }
    closedir(d);
  }

}

int FindWord(char *string, char *word)
{
	int i,l,m,c,k;
	
	l = strlen(word);
	m = strlen(string);

	for (i=0;i<m;i++)
	{
		c=0;
		for (k=0;k<l;k++)
			if ( i+k < m )
				if ( string[i+k] == word[k] ) c++;
		if ( c == l ) return 1;
	}
	
	return 0;
}


void UpdateMdp(char *fmdp, char *omdp, char *xtcgrps, char *energygrps)
{
	char aux[800];
	FILE *fin,*fout;

	fin = fopen(fmdp,"r");
	if (!fin) Error("Cannot open mdp file in UpdateMdp");
	fout = fopen(omdp,"w");
	if (!fout) Error("Cannot open mdp file to write in UpdateMdp");
	
	
	while ( fgets(aux,500,fin) )
	{
		if (strncmp(aux,"xtc-grps",8) && strncmp(aux,"energygrps",10))
			fprintf(fout,"%s",aux);
		else if (!strncmp(aux,"xtc-grps",8))	
			fprintf(fout,"xtc-grps\t\t= %s\n",xtcgrps);
		else if (!strncmp(aux,"energygrps",10))	
			fprintf(fout,"energygrps\t\t= %s\n",energygrps);

	}

	fclose(fin);
	fclose(fout);
}

void AdduNOE(struct parm_s *parm, struct NOErelax_s *NOErelax)
{
	int i,j,k,ok,ad=0;

	fprintf(stderr,"Add uNoe:\n");

	for (i=0;i<NOErelax->nh;i++)
		for (j=i+15;j<NOErelax->nh;j++)
		{
			for (k=0;k<parm->nres;k++)
			{
				ok = 1;
				if ( ( NOErelax->Hlist[i] == parm->ires1[k] && NOErelax->Hlist[j] == parm->ires2[k] ) ||
					( NOErelax->Hlist[i] == parm->ires2[k] && NOErelax->Hlist[j] == parm->ires1[k] )) ok = 0;
			}

			if (ok == 1) 
			{
				parm->ires1[parm->nres] = NOErelax->Hlist[i];	
				parm->ires2[parm->nres] = NOErelax->Hlist[j];	
				parm->restype[parm->nres] = 3;	
				parm->resexp[parm->nres] = 0;
				parm->ressigma[parm->nres] = parm->uNOE;
				parm->nres ++;
				if ( parm->verbose) fprintf(stderr,"%4d %4d 3 0 %lf\n", NOErelax->Hlist[i], NOErelax->Hlist[j],parm->uNOE);
				ad++;
			}
			
		}

	fprintf(stderr," %d new unobserved NOEs added\n",ad);
}

