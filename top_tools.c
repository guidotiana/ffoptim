/*
 * top_tools.c
 *
 *  Created on: Apr 21, 2009
 *      Author: guido
 */
#include "ffoptim.h"

struct top_s *AlloTop(int natoms, int nbonds, int npairs, int nangles,
		int ndihedrals, int nimpropers, int nexclusions, int verbose)
{
	struct top_s *x;
	int i;

	x = (struct top_s *) calloc(1,sizeof(struct top_s));
	x->top_atoms = (struct top_atoms_s *) calloc(natoms,sizeof(struct top_atoms_s));
	if ( !(x->top_atoms) ) Error("Error in allocating top");
	x->top_bonds = (struct top_gen_s *) calloc(nbonds,sizeof(struct top_gen_s));
	if ( !(x->top_bonds) ) Error("Error in allocating top");
	x->top_pairs = (struct top_gen_s *) calloc(npairs,sizeof(struct top_gen_s));
	if ( !(x->top_pairs) ) Error("Error in allocating top");
	x->top_angles = (struct top_gen_s *) calloc(nangles,sizeof(struct top_gen_s));
	if ( !(x->top_pairs) ) Error("Error in allocating top");
	x->top_dihedrals = (struct top_gen_s *) calloc(ndihedrals,sizeof(struct top_gen_s));
	if ( !(x->top_dihedrals) ) Error("Error in allocating top");
	x->top_impropers = (struct top_gen_s *) calloc(nimpropers,sizeof(struct top_gen_s));
	if ( !(x->top_impropers) ) Error("Error in allocating top");
	x->top_exclusions = (struct top_gen_s *) calloc(nexclusions,sizeof(struct top_gen_s));
	if ( !(x->top_impropers) ) Error("Error in allocating top");

	if (verbose) fprintf(stdout," Allocated top structure (%d kB)\n",(int)(natoms*sizeof(struct top_atoms_s)+
			  (nbonds+npairs+nangles+ndihedrals+nimpropers)*sizeof(struct top_gen_s))/1024);

	for (i=0;i<nbonds;i++)
	{
		((x->top_bonds)+i)->c0 = NONDEF;		((x->top_bonds)+i)->c1 = NONDEF;
		((x->top_bonds)+i)->c2 = NONDEF;		((x->top_bonds)+i)->c3 = NONDEF;
	}
	for (i=0;i<npairs;i++)
	{
		((x->top_pairs)+i)->c0 = NONDEF;		((x->top_pairs)+i)->c1 = NONDEF;
		((x->top_pairs)+i)->c2 = NONDEF;		((x->top_pairs)+i)->c3 = NONDEF;
		((x->top_pairs)+i)->move = 0;
	}
	for (i=0;i<nangles;i++)
	{
		((x->top_angles)+i)->c0 = NONDEF;		((x->top_angles)+i)->c1 = NONDEF;
		((x->top_angles)+i)->c2 = NONDEF;		((x->top_angles)+i)->c3 = NONDEF;
	}
	for (i=0;i<ndihedrals;i++)
	{
		((x->top_dihedrals)+i)->c0 = NONDEF;		((x->top_dihedrals)+i)->c1 = NONDEF;
		((x->top_dihedrals)+i)->c2 = NONDEF;		((x->top_dihedrals)+i)->c3 = NONDEF;
		((x->top_dihedrals)+i)->c4 = NONDEF;		((x->top_dihedrals)+i)->c5 = NONDEF;
		((x->top_dihedrals)+i)->move = 0;

	}
	for (i=0;i<nimpropers;i++)
	{
		((x->top_impropers)+i)->c0 = NONDEF;		((x->top_impropers)+i)->c1 = NONDEF;
		((x->top_impropers)+i)->c2 = NONDEF;		((x->top_impropers)+i)->c3 = NONDEF;
		((x->top_impropers)+i)->c4 = NONDEF;		((x->top_impropers)+i)->c5 = NONDEF;

	}

	return x;

}

struct top_s *ReadTop(char *fname, char *itpdir, int verbose)
{
	int i1,i2,i3,i4,i5,i6,i7,chapt,j,ihead=0,ifoot=0;
	struct top_s *x;
	char aux[1000],aux2[1000],keyword[100];
	FILE *fin;

	if (verbose) printf("Reading topology %s\n",fname);

	// Allocate structure topology
	x = AlloTop(NATOMSMAX,NBONDSMAX,NPAIRSMAX,NANGLESMAX,
			NDIHEDRALSMAX,NIMPROPERSMAX,NEXCLUSIONSMAX,verbose);

	// Read #define in the topology files and store them
	//if (verbose) printf("  Reading label definitions:\n");
	//ReadDefines(fname,itpdir,x,0);
   	// ReadDefinesSimple(itpdir,x,verbose);

	// Open file
	fin = fopen(fname,"r");
	if (fin==NULL) Error2("Cannot open topology file ",fname);

	// Read lines
	chapt = 0;
	i1 = 0; i2 = 0; i3 = 0; i4 = 0; i5 = 0; i6 = 0; i7=0;

	while ( fgets(aux,1000,fin) )
	{
	   // switch which chapter of the top file is currently reading
	   if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"defaults") )  chapt = 0;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"moleculetype") )  chapt = 0;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"atomtypes") )  chapt = 0;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"atoms") )  chapt = 1;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"bonds") )  chapt = 2;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"pairs") )  chapt = 3;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"angles") )  chapt = 4;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"dihedrals") )  chapt = 5;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"impropers") )  chapt = 6;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"system") )  chapt = 8;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"molecules") )  chapt = 8;
	   else if (  FindKeyword(aux,keyword)==1 && !strcmp(keyword,"exclusions") )  chapt = 9;
	   else if ( i1>0 && ( !strncmp(aux,"[",1) || !strncmp(aux,"#include",8) ) )  chapt = 8;   
	  // else if (  FindKeyword(aux,keyword)==1  ) chapt = 8;					// footer

	   // read header
	   if (chapt == 0)
	   {
		   if (ihead>NHEADMAX) Error("Header of top too large");
		   strcpy(x->header[ihead],aux);
		   ihead++;
	   }
	   // read atoms
	   else if (chapt == 1)
	   {
		  if ( sscanf(aux,"%d %s %d %s %s %d %f %f",
				&((x->top_atoms)+i1)->nr, ((x->top_atoms)+i1)->type, &((x->top_atoms)+i1)->resnr,
				((x->top_atoms)+i1)->residue, ((x->top_atoms)+i1)->atom, &((x->top_atoms)+i1)->cgnr,
				&((x->top_atoms)+i1)->charge, &((x->top_atoms)+i1)->mass ) == 8) i1++;
		 if (i1>=NATOMSMAX) Error("NATOMSMAX too small");
	   }

	   // read bonds
	   else if (chapt == 2)
	   	   {
	   		  if ( sscanf(aux,"%d %d %d %f %f %f %f",
	   				&((x->top_bonds)+i2)->i1, &((x->top_bonds)+i2)->i2, &((x->top_bonds)+i2)->funct,
	   				&((x->top_bonds)+i2)->c0, &((x->top_bonds)+i2)->c1, &((x->top_bonds)+i2)->c2,
	   				&((x->top_bonds)+i2)->c3) > 2) i2++;
		 	  if (i2>=NBONDSMAX) Error("NBONDSMAX too small");
	   	   }

	   // read pairs
	   else if (chapt == 3)
	   	   {
	   	   	  if ( sscanf(aux,"%d %d %d %f %f %f %f",
	   	    		&((x->top_pairs)+i3)->i1, &((x->top_pairs)+i3)->i2, &((x->top_pairs)+i3)->funct,
	   	   			&((x->top_pairs)+i3)->c0, &((x->top_pairs)+i3)->c1, &((x->top_pairs)+i3)->c2,
	   				&((x->top_pairs)+i3)->c3) > 2) 
				{
			  		if(strstr(aux,"NOC") != NULL) ((x->top_pairs)+i3)->move = -1;
			  		if(strstr(aux,"MOVE") != NULL) ((x->top_pairs)+i3)->move = 1;
					i3++;
				}
			  if (i3>=NPAIRSMAX) Error("NPAIRSMAX too small");

	   	   }

	   // read angles
  	   else if (chapt == 4)
  	   {
	   	   	  if ( sscanf(aux,"%d %d %d %d %f %f %f %f",
	  				&((x->top_angles)+i4)->i1, &((x->top_angles)+i4)->i2, &((x->top_angles)+i4)->i3, &((x->top_angles)+i4)->funct,
	   	   	   		&((x->top_angles)+i4)->c0, &((x->top_angles)+i4)->c1, &((x->top_angles)+i4)->c2,
	   	   	   		&((x->top_angles)+i4)->c3) > 3) i4++;
  			  if (i4>=NANGLESMAX) Error("NANGLESMAX too small");

	   }

	   // read dihedrals
  	   else if (chapt == 5)
  	   {
		      if (DEBUG) printf("D> %s\n",aux);

  		      // if parameters are explicitely given
  		   	  if ( sscanf(aux,"%d %d %d %d %d %f %f %f %f %f %f",
  		 	  				&((x->top_dihedrals)+i5)->i1, &((x->top_dihedrals)+i5)->i2, 
							&((x->top_dihedrals)+i5)->i3, &((x->top_dihedrals)+i5)->i4,
  		 	  				&((x->top_dihedrals)+i5)->funct, &((x->top_dihedrals)+i5)->c0, 
							&((x->top_dihedrals)+i5)->c1, &((x->top_dihedrals)+i5)->c2,
  		 	   	   	   		&((x->top_dihedrals)+i5)->c3, &((x->top_dihedrals)+i5)->c4, 
							&((x->top_dihedrals)+i5)->c5) > 5) 
  		   	  {

			   	 if(strstr(aux,"NOC") != NULL) ((x->top_dihedrals)+i5)->move = -1;
			   	 if(strstr(aux,"MOVE") != NULL) ((x->top_dihedrals)+i5)->move = 1;
  		   		 i5++;
				 if (i5>=NDIHEDRALSMAX) Error("NDIHEDRALSMAX too small");
  		   	  }  		   	  
  		   	  
			  // or if parameters are given by defined label 
  		   	  else if ( sscanf(aux,"%d %d %d %d %d %s",&((x->top_dihedrals)+i5)->i1, 
					&((x->top_dihedrals)+i5)->i2, &((x->top_dihedrals)+i5)->i3, &((x->top_dihedrals)+i5)->i4,
  	  				&((x->top_dihedrals)+i5)->funct,aux2) > 5 )
  		          {
  		    	      for (j=0;j<NLABELS;j++)
  		    		  if ( !strcmp(aux2,(x->label)[j]) )
  		    		  {
  		    			((x->top_dihedrals)+i5)->c0 = (x->labelc)[j][0];
  		    			((x->top_dihedrals)+i5)->c1 = (x->labelc)[j][1];
  		    			((x->top_dihedrals)+i5)->c2 = (x->labelc)[j][2];
  		    			((x->top_dihedrals)+i5)->c3 = (x->labelc)[j][3];
  		    			((x->top_dihedrals)+i5)->c4 = (x->labelc)[j][4];
  		    			((x->top_dihedrals)+i5)->c5 = (x->labelc)[j][5];
  		    		  }
			     if(strstr(aux,"NOC") != NULL) ((x->top_dihedrals)+i5)->move = -1;
			     if(strstr(aux,"MOVE") != NULL) ((x->top_dihedrals)+i5)->move = 1;
       		    	     i5++;
			     if (i5>=NDIHEDRALSMAX) Error("NDIHEDRALSMAX too small");
  		         }
  		   	 // or not given
  		   	 else if ( sscanf(aux,"%d %d %d %d %d",&((x->top_dihedrals)+i5)->i1, &((x->top_dihedrals)+i5)->i2, 
								&((x->top_dihedrals)+i5)->i3, &((x->top_dihedrals)+i5)->i4,
  		   	  	  				&((x->top_dihedrals)+i5)->funct) == 5 ) 
			 {	
			        if(strstr(aux,"NOC") != NULL) ((x->top_dihedrals)+i5)->move = -1;
			        if(strstr(aux,"MOVE") != NULL) ((x->top_dihedrals)+i5)->move = 1;
				i5++;
			  	if (i5>=NDIHEDRALSMAX) Error("NDIHEDRALSMAX too small");
			 }
  		      
	   }

	   // read impropers
 	   else if (chapt == 6)
  	   {
		      if (DEBUG) printf("D> %s\n",aux);

  		      // if parameters are explicitely given
  		   	  if ( sscanf(aux,"%d %d %d %d %d %f %f %f %f %f %f",
  		 	  				&((x->top_impropers)+i6)->i1, &((x->top_impropers)+i6)->i2, 
							&((x->top_impropers)+i6)->i3, &((x->top_impropers)+i6)->i4,
  		 	  				&((x->top_impropers)+i6)->funct, &((x->top_impropers)+i6)->c0, 
							&((x->top_impropers)+i6)->c1, &((x->top_impropers)+i6)->c2,
  		 	   	   	   		&((x->top_impropers)+i6)->c3, &((x->top_impropers)+i6)->c4, 
							&((x->top_impropers)+i6)->c5) > 5) 
  		   	  {
  		   		    	       i6++;
					       if (i6>=NIMPROPERSMAX) Error("NIMPROPERSSMAX too small");
  		  	  }  		   	  
  		   	  
			  // or if parameters are given by defined label or not given
  		   	  else if ( sscanf(aux,"%d %d %d %d %d %s",&((x->top_impropers)+i6)->i1, &((x->top_impropers)+i6)->i2, 
					&((x->top_impropers)+i6)->i3, &((x->top_impropers)+i6)->i4,
  	  				&((x->top_impropers)+i6)->funct,aux2) > 4 )
  		          {
  		    	   	for (j=0;j<NLABELS;j++)
  		    		  if ( !strcmp(aux2,(x->label)[j]) )
  		    		  {
  		    			((x->top_impropers)+i6)->c0 = (x->labelc)[j][0];
  		    			((x->top_impropers)+i6)->c1 = (x->labelc)[j][1];
  		    			((x->top_impropers)+i6)->c2 = (x->labelc)[j][2];
  		    			((x->top_impropers)+i6)->c3 = (x->labelc)[j][3];
  		    			((x->top_impropers)+i6)->c4 = (x->labelc)[j][4];
  		    			((x->top_impropers)+i6)->c5 = (x->labelc)[j][5];
  		    		  }
  		    
  		    	       i6++;
	 	               if (i6>=NIMPROPERSMAX) Error("NIMPROPERSSMAX too small");
  		         }
 		   	  // or not given
  		   	  else if ( sscanf(aux,"%d %d %d %d %d",&((x->top_impropers)+i5)->i1, &((x->top_impropers)+i5)->i2, 
								&((x->top_impropers)+i5)->i3, &((x->top_impropers)+i5)->i4,
  		   	  	  				&((x->top_impropers)+i5)->funct) == 5 ) i6++;
			  if (i6>=NIMPROPERSMAX) Error("NIMPROPERSSMAX too small");

	   }
	  // read exclusions
	  else if  ( chapt == 9 )
	  {
		  if (DEBUG) printf("D> %s\n",aux);

		  if ( sscanf(aux,"%d %d", &((x->top_exclusions)+i7)->i1, &((x->top_exclusions)+i7)->i2 ) == 2) i7++;
	 	  if (i7>=NEXCLUSIONSMAX) Error("NEXCLUSIONSAMX too small");
		  
	  }
	  // read footer
  	  else if ( chapt == 8 )
  	  {
  		   if (ifoot>NFOOTMAX) Error("Footer of top too large");
  		   strcpy(x->footer[ifoot],aux);
  		   ifoot++;
  	  }
	   
	}

	x->natoms = i1;
	x->nbonds = i2;
	x->npairs = i3;
	x->nangles = i4;
	x->ndihedrals = i5;
	x->nimpropers = i6;
	x->nexclusions = i7;

	fclose(fin);
    	if (verbose) fprintf(stderr," Topology read (\t%d atoms, %d bonds, %d pairs,\n\t%d angles, %d dihedrals, %d impropers, %d exclusions).\n",
    		   x->natoms,x->nbonds,x->npairs,x->nangles,x->ndihedrals,x->nimpropers,x->nexclusions);

	return x;
}

void InsertDihedralsTop(struct top_s *x, int verbose, char *bondedFile)
{
	int i,j,ok=0,nr,n=0;
	float r0,r1,r2,r3,r4,r5;
	int rf;
	char ra1[5],ra2[5],ra3[5],ra4[5],t1[5],t2[5],t3[5],t4[5],aux[500],keyword[20];
	char a0[NINCLMAX][5],a1[NINCLMAX][5],a2[NINCLMAX][5],a3[NINCLMAX][5];	
	int iff[NINCLMAX];
	float c0[NINCLMAX],c1[NINCLMAX],c2[NINCLMAX],c3[NINCLMAX],c4[NINCLMAX],c5[NINCLMAX];

	FILE *fp;

	// read dihedrals
	if (verbose) fprintf(stderr,"Read dihedrals data from included file %s\n",bondedFile);

	fp = fopen(bondedFile,"r");
	if (!fp) Error("Cannot open file with included bonded parameters");
	
 	while ( fgets(aux,500,fp) )
	{
		if ( aux[0] == '[' ) ok=0;
		strcpy(keyword,"dihedraltypes");
		if ( FindKeyword(aux,keyword) ) ok=1;
		if ( ok ==1 ) 
		{	
			nr = sscanf(aux,"%s %s %s %s  %d  %f %f %f %f %f %f",ra1,ra2,ra3,ra4,&rf,&r0,&r1,&r2,&r3,&r4,&r5);
			if ( nr > 7 )
			{
				strcpy(a0[n],ra1);	
				strcpy(a1[n],ra2);	
				strcpy(a2[n],ra3);	
				strcpy(a3[n],ra4);	
				iff[n] = rf;
				c0[n]= r0;
				c1[n]= r1;
				c2[n]= r2;
				if ( nr > 8 ) c3[n]= r3; else c3[n]=NONDEF;
				if ( nr > 9 ) c4[n]= r4; else c4[n]=NONDEF;
				if ( nr > 10 ) c5[n]= r5; else c5[n]=NONDEF;
				n++;
				if ( n > NINCLMAX ) Error("NINCLMAX too small");
			}
		}
	}

	if (verbose) fprintf(stderr," %d atomtypes read.\n",n);
	if (n == 0) Error("No atomtype read in bonded itp file");
                                                                                                                                 
	// assign dihedrals
	for (i=0;i<x->ndihedrals;i++)
		if ( ((x->top_dihedrals)+i)->c0 == NONDEF )
		{
			strcpy(t1,((x->top_atoms)+ (((x->top_dihedrals)+i)->i1)-1 )->type);
			strcpy(t2,((x->top_atoms)+ (((x->top_dihedrals)+i)->i2)-1 )->type);
			strcpy(t3,((x->top_atoms)+ (((x->top_dihedrals)+i)->i3)-1 )->type);
			strcpy(t4,((x->top_atoms)+ (((x->top_dihedrals)+i)->i4)-1 )->type);
			rf =  ((x->top_dihedrals)+i)->funct;

			// find atom types
			ok = -1;
			for (j=0;j<n;j++)
			{
				
				if ( rf == iff[j] && ( !strcmp(t1,a0[j]) ||  !strcmp(a0[j],"X") )  &&  
					( !strcmp(t2,a1[j]) ||  !strcmp(a1[j],"X") ) &&  
					( !strcmp(t3,a2[j]) ||  !strcmp(a2[j],"X") ) && 
 					(  !strcmp(t4,a3[j]) ||  !strcmp(a3[j],"X"))  ) ok = j;
				if ( rf == iff[j] && ( !strcmp(t1,a3[j]) ||  !strcmp(a3[j],"X") )  &&  
					( !strcmp(t2,a2[j]) ||  !strcmp(a2[j],"X") ) &&  
					( !strcmp(t3,a1[j]) ||  !strcmp(a1[j],"X") ) && 
 					(  !strcmp(t4,a0[j]) ||  !strcmp(a0[j],"X"))  ) ok = j;

			}
			if (ok == -1)  { sprintf(aux,"Atomtypes %s-%s-%s-%s with funct %d not found in bonded file",t1,t2,t3,t4,rf); Error(aux); }

			((x->top_dihedrals)+i)->funct = iff[ok];
			((x->top_dihedrals)+i)->c0 = c0[ok];
			((x->top_dihedrals)+i)->c1 = c1[ok];
			((x->top_dihedrals)+i)->c2 = c2[ok];
			((x->top_dihedrals)+i)->c3 = c3[ok];
			((x->top_dihedrals)+i)->c4 = c4[ok];
			((x->top_dihedrals)+i)->c5 = c5[ok];
		}
	if (verbose) fprintf(stderr," done.\n");

}



void InsertNonbondedTop(struct top_s *x, int combrule, int verbose, char *nonbondedFile)
{
	int i,j,n=0,i1,i2;
	double sigma[NINCLMAX];
	double epsilon[NINCLMAX];
	double cs,ce;
	char name[NINCLMAX][8];
	char aux[500],aux2[500];
	char t1[8],t2[8];

	FILE *fp;

	// read nonbonded
	if (verbose) fprintf(stderr,"Read nonbonded data from included file %s\n",nonbondedFile);

	fp = fopen(nonbondedFile,"r");
	if (!fp) Error("Cannot open file with included non-bonded parameters");
	
 	while ( fgets(aux,500,fp) )
		if ( sscanf(aux,"%s %*d %*f %*f %*s %lf %lf",aux2,&cs,&ce) == 3 ) 
		{
			strcpy(name[n],aux2);
			sigma[n] = cs;
			epsilon[n] = ce;
			n++;
			if ( n > NINCLMAX ) Error("NINCLMAX too small");
		}

	if (verbose) fprintf(stderr," %d atomtypes read.\n",n);
	if (n == 0) Error("No atomtype read in non-bonded itp file");

	// insert in empty pairs
	for (i=0;i<x->npairs;i++)
		if ( ((x->top_pairs)+i)->c0 == NONDEF )
		{
			strcpy(t1,((x->top_atoms)+ (((x->top_pairs)+i)->i1)-1 )->type);
			strcpy(t2,((x->top_atoms)+ (((x->top_pairs)+i)->i2)-1 )->type);

			// find corresponding atomtype
			i1 = -1;
			i2 = -1;
			for (j=0;j<n;j++)
			{
				if ( !strcmp(t1, name[j]) )  i1 = j;	
				if ( !strcmp(t2, name[j]) )  i2 = j;	
			}
			if (i1 == -1) { sprintf(aux,"Atomtype %s not found in nonbonded file",t1); Error(aux); }
			if (i2 == -1) { sprintf(aux,"Atomtype %s not found in nonbonded file",t2); Error(aux); }

			// calculate current combrule
			if (combrule == 1)
			{
				cs = sqrt(sigma[i1]*sigma[i2]);
				ce = sqrt(epsilon[i1]*epsilon[i2]);
			}
			else if (combrule == 2)
			{
				cs = (sigma[i1] + sigma[i2]) / 2;
				ce = sqrt(epsilon[i1]*epsilon[i2]);
			}

			((x->top_pairs)+i)->c0 = cs;
			((x->top_pairs)+i)->c1 = ce;
		}
	
	if (verbose) fprintf(stderr," done.\n");
}


void FillTop(struct top_s *x, int combrule, int nochangeH)
{
	int i,j,k,nold=x->npairs,n,ok,i1,i2;
	double c0,c1;

	if (combrule==1)
	{
		c0 = 0.;
		c1 = 0.;
	}
	else if (combrule==2)
	{
		c0 = 0.15;
		c1 = 0.;
	}

	for (i=0;i<x->natoms;i++)
		for (j=i+1;j<x->natoms;j++)
		{
			ok=1;
			i1 = ((x->top_atoms)+i)->nr;
			i2 = ((x->top_atoms)+j)->nr;
			if ( !nochangeH || ( (((x->top_atoms)+i)->type)[0]!='H' && (((x->top_atoms)+i)->type)[0]!='h' &&
							(((x->top_atoms)+j)->type)[0]!='H' && (((x->top_atoms)+j)->type)[0]!='h') )			     {
				for (k=0;k<nold;k++)
					if ( ( ((x->top_pairs)+k)->i1 == i1 && ((x->top_pairs)+k)->i2 == i2 ) ||
						( ((x->top_pairs)+k)->i1 == i2 && ((x->top_pairs)+k)->i2 == i1 ) )
					{
						ok=0;	
						break;
					}

				if (ok == 1)
				{
					 n = x->npairs;
					 ((x->top_pairs)+n)->i1 = i1;
					 ((x->top_pairs)+n)->i2 = i2;
					 ((x->top_pairs)+n)->funct = 1;
					 ((x->top_pairs)+n)->c0 = c0;
					 ((x->top_pairs)+n)->c1 = c1;
					 (x->npairs) ++;
				}
			}
		}

	fprintf(stderr,"Filling pairs with those not in the top files.\n");
	fprintf(stderr,"New number of pairs is %d\n",x->npairs);
}



void PrintTop(char *fname, struct top_s *x, int verbose, int noSOL)
{
   int i;
   FILE *fout;

   if (!strcmp(fname,"stdout")) fout=stdout;
   else if (!strcmp(fname,"stderr")) fout=stderr;
   else 
   {
	fout=fopen(fname,"w");
   	if (!fout) Error2("Cannot open for writing ",fname);
   }

   // print header
   for (i=0;i<NHEADMAX;i++)
	   if ( strcmp(x->header[i],"") ) fprintf(fout,"%s",x->header[i]);
   
   // print atoms
   fprintf(fout,"[ atoms ]\n");
   for (i=0;i<x->natoms;i++)
	 fprintf(fout,"\t%d\t%s\t%d\t%s\t%s\t%d\t%f\t%f \n",
			                ((x->top_atoms)+i)->nr, ((x->top_atoms)+i)->type, ((x->top_atoms)+i)->resnr,
			 				((x->top_atoms)+i)->residue, ((x->top_atoms)+i)->atom, ((x->top_atoms)+i)->cgnr,
			 				((x->top_atoms)+i)->charge, ((x->top_atoms)+i)->mass );
   fprintf(fout,"\n");

   // print bonds
   fprintf(fout,"[ bonds ]\n");
   for (i=0;i<x->nbonds;i++)
   {
	  fprintf(fout,"\t%d\t%d\t%d",((x->top_bonds)+i)->i1,((x->top_bonds)+i)->i2,((x->top_bonds)+i)->funct);
	  if (((x->top_bonds)+i)->c0 != NONDEF)
		  {
		     fprintf(fout,"\t%e", ((x->top_bonds)+i)->c0);
		     if (((x->top_bonds)+i)->c1 != NONDEF)
		     		  {
		     		     fprintf(fout,"\t%e", ((x->top_bonds)+i)->c1);
		     		     if (((x->top_bonds)+i)->c2 != NONDEF)
		     		    		     		  {
		     		    		     		     fprintf(fout,"\t%e", ((x->top_bonds)+i)->c2);
		     		    		     		     if (((x->top_bonds)+i)->c3 != NONDEF)
		     		    		     		   		     		     fprintf(fout,"\t%e", ((x->top_bonds)+i)->c3);
		     		    		     		  }
		     		  }
		  }
	  
	  if (verbose) fprintf(fout,"\t; %s-%s (%s-%s)",((x->top_atoms)+ (((x->top_bonds)+i)->i1)-1 )->atom,
			  ((x->top_atoms)+ (((x->top_bonds)+i)->i2)-1 )->atom,((x->top_atoms)+ (((x->top_bonds)+i)->i1)-1 )->type,
			  ((x->top_atoms)+ (((x->top_bonds)+i)->i2)-1 )->type );
	  
	  fprintf(fout,"\n");
   }

   // print pairs
   fprintf(fout,"[ pairs ]\n");
   for (i=0;i<x->npairs;i++)
   {
	  fprintf(fout,"\t%d\t%d\t%d",((x->top_pairs)+i)->i1,((x->top_pairs)+i)->i2,((x->top_pairs)+i)->funct);
	  if (((x->top_pairs)+i)->c0 != NONDEF)
		  {
		     fprintf(fout,"\t%e", ((x->top_pairs)+i)->c0);
		     if (((x->top_pairs)+i)->c1 != NONDEF)
		     		  {
		     		     fprintf(fout,"\t%e", ((x->top_pairs)+i)->c1);
		     		     if (((x->top_pairs)+i)->c2 != NONDEF)
		     		    		     		  {
		     		    		     		     fprintf(fout,"\t%e", ((x->top_pairs)+i)->c2);
		     		    		     		     if (((x->top_pairs)+i)->c3 != NONDEF)
		     		    		     		   		     		     fprintf(fout,"\t%e", ((x->top_pairs)+i)->c3);
		     		    		     		  }
		     		  }
		  }
	 
	  if ( ((x->top_pairs)+i)->move == -1 )  fprintf(fout,"; NOC ");
	  if ( ((x->top_pairs)+i)->move == 1 )  fprintf(fout,"; MOVE ");
 
	  if (verbose) fprintf(fout,"\t; %s-%s (%s-%s)",((x->top_atoms)+ (((x->top_pairs)+i)->i1)-1 )->atom,
			  ((x->top_atoms)+ (((x->top_pairs)+i)->i2)-1 )->atom,((x->top_atoms)+ (((x->top_pairs)+i)->i1)-1 )->type,
			  ((x->top_atoms)+ (((x->top_pairs)+i)->i2)-1 )->type );
	  
	  fprintf(fout,"\n");
   }
   
   // print angles
   fprintf(fout,"[ angles ]\n");
   for (i=0;i<x->nangles;i++)
   {
	  fprintf(fout,"\t%d\t%d\t%d\t%d",((x->top_angles)+i)->i1,((x->top_angles)+i)->i2,((x->top_angles)+i)->i3,((x->top_angles)+i)->funct);
	  if (((x->top_angles)+i)->c0 != NONDEF)
		  {
		     fprintf(fout,"\t%e", ((x->top_angles)+i)->c0);
		     if (((x->top_angles)+i)->c1 != NONDEF)
		     		  {
		     		     fprintf(fout,"\t%e", ((x->top_angles)+i)->c1);
		     		     if (((x->top_angles)+i)->c2 != NONDEF)
		     		    		     		  {
		     		    		     		     fprintf(fout,"\t%e", ((x->top_angles)+i)->c2);
		     		    		     		     if (((x->top_angles)+i)->c3 != NONDEF)
		     		    		     		   		     		     fprintf(fout,"\t%e", ((x->top_angles)+i)->c3);
		     		    		     		  }
		     		  }
		  }
	  
	  if (verbose) fprintf(fout,"\t; %s-%s-%s (%s-%s-%s)",((x->top_atoms)+ (((x->top_angles)+i)->i1)-1 )->atom,
			  ((x->top_atoms)+ (((x->top_angles)+i)->i2)-1 )->atom,((x->top_atoms)+ (((x->top_angles)+i)->i3)-1 )->atom,
			  ((x->top_atoms)+ (((x->top_angles)+i)->i1)-1 )->type,
			  ((x->top_atoms)+ (((x->top_angles)+i)->i2)-1 )->type, ((x->top_atoms)+ (((x->top_angles)+i)->i3)-1 )->type );
	  
	  fprintf(fout,"\n");
   }

   // print dihedrals
   fprintf(fout,"[ dihedrals ]\n");
   for (i=0;i<x->ndihedrals;i++)
   {
 	  fprintf(fout,"\t%d\t%d\t%d\t%d\t%d",((x->top_dihedrals)+i)->i1,((x->top_dihedrals)+i)->i2,((x->top_dihedrals)+i)->i3,
 			 ((x->top_dihedrals)+i)->i4,((x->top_dihedrals)+i)->funct);
 	  if (((x->top_dihedrals)+i)->c0 != NONDEF)
 		  {
 		    fprintf(fout,"\t%e", ((x->top_dihedrals)+i)->c0);
 		    if (((x->top_dihedrals)+i)->c1 != NONDEF)
 		    {
 		       fprintf(fout,"\t%e", ((x->top_dihedrals)+i)->c1);

		       if (((x->top_dihedrals)+i)->c3 == NONDEF)	// if multiplicity defined (print as integer)
		       {
				if (((x->top_dihedrals)+i)->c2 != NONDEF) fprintf(fout,"\t%d", (int)((x->top_dihedrals)+i)->c2);
		       }
		       else 						// if more complex
		       {
 		       		if (((x->top_dihedrals)+i)->c2 != NONDEF)
 		       		{
 		     	  		fprintf(fout,"\t%e", ((x->top_dihedrals)+i)->c2);
 		          		if (((x->top_dihedrals)+i)->c3 != NONDEF)
 		     	  		{
 		     		 		fprintf(fout,"\t%e", ((x->top_dihedrals)+i)->c3);
 		     		 		if (((x->top_dihedrals)+i)->c4 != NONDEF)
 		     		 		{
 		     		    			fprintf(fout,"\t%e", ((x->top_dihedrals)+i)->c4);
 		     		    			if (((x->top_dihedrals)+i)->c5 != NONDEF)
 		     		    			fprintf(fout,"\t%e", ((x->top_dihedrals)+i)->c5);
 		     		 		}
 		     	  		}
 		       		}
			}
 		    }
 		  }
 	  
	  if ( ((x->top_dihedrals)+i)->move == -1 )  fprintf(fout,"; NOC ");
	  if ( ((x->top_dihedrals)+i)->move == 1 )  fprintf(fout,"; MOVE ");

 	  if (verbose) fprintf(fout,"\t; %s-%s-%s-%s (%s-%s-%s-%s)",((x->top_atoms)+ (((x->top_dihedrals)+i)->i1)-1 )->atom,
 				  ((x->top_atoms)+ (((x->top_dihedrals)+i)->i2)-1 )->atom,((x->top_atoms)+ (((x->top_dihedrals)+i)->i3)-1 )->atom,
 				  ((x->top_atoms)+ (((x->top_dihedrals)+i)->i4)-1 )->atom,
 				  ((x->top_atoms)+ (((x->top_dihedrals)+i)->i1)-1 )->type, ((x->top_atoms)+ (((x->top_dihedrals)+i)->i2)-1 )->type, 
 				  ((x->top_atoms)+ (((x->top_dihedrals)+i)->i3)-1 )->type, ((x->top_atoms)+ (((x->top_dihedrals)+i)->i4)-1 )->type );
 	  
 	  fprintf(fout,"\n");
   }

   // print impropers
   if (x->nimpropers>0) fprintf(fout,"[ impropers ]\n");
   for (i=0;i<x->nimpropers;i++)
   {
 	  fprintf(fout,"\t%d\t%d\t%d\t%d\t%d",((x->top_impropers)+i)->i1,((x->top_impropers)+i)->i2,((x->top_impropers)+i)->i3,
 			 ((x->top_impropers)+i)->i4,((x->top_impropers)+i)->funct);
 	  if (((x->top_impropers)+i)->c0 != NONDEF)
 		  {
 		    fprintf(fout,"\t%e", ((x->top_impropers)+i)->c0);
 		    if (((x->top_impropers)+i)->c1 != NONDEF)
 		    {
 		       fprintf(fout,"\t%e", ((x->top_impropers)+i)->c1);
 		       if (((x->top_impropers)+i)->c2 != NONDEF)
 		       {
 		     	  fprintf(fout,"\t%e", ((x->top_impropers)+i)->c2);
 		          if (((x->top_impropers)+i)->c3 != NONDEF)
 		     	  {
 		     		 fprintf(fout,"\t%e", ((x->top_impropers)+i)->c3);
 		     		 if (((x->top_impropers)+i)->c4 != NONDEF)
 		     		 {
 		     		    fprintf(fout,"\t%e", ((x->top_impropers)+i)->c4);
 		     		    if (((x->top_impropers)+i)->c5 != NONDEF)
 		     		    	fprintf(fout,"\t%e", ((x->top_impropers)+i)->c5);
 		     		 }
 		     	  }
 		       }
 		    }
 		  }
 	  
 	  if (verbose) fprintf(fout,"\t; %s-%s-%s-%s (%s-%s-%s-%s)",((x->top_atoms)+ (((x->top_impropers)+i)->i1)-1 )->atom,
 	 				  ((x->top_atoms)+ (((x->top_impropers)+i)->i2)-1 )->atom,((x->top_atoms)+ (((x->top_impropers)+i)->i3)-1 )->atom,
 	 				 ((x->top_atoms)+ (((x->top_impropers)+i)->i4)-1 )->atom,
 	 				  ((x->top_atoms)+ (((x->top_impropers)+i)->i1)-1 )->type, ((x->top_atoms)+ (((x->top_impropers)+i)->i2)-1 )->type, 
 	 				  ((x->top_atoms)+ (((x->top_impropers)+i)->i3)-1 )->type, ((x->top_atoms)+ (((x->top_impropers)+i)->i4)-1 )->type );
 	 	  
 	  
 	  fprintf(fout,"\n");
   }
   
   // print exclusions
   if (x->nexclusions>0) fprintf(fout,"[ exclusions ]\n");
   for (i=0;i<x->nexclusions;i++)
   {
	 fprintf(fout,"\t%d\t%d\n",((x->top_exclusions)+i)->i1,((x->top_exclusions)+i)->i2);
   }
  


   // print footer
    for (i=0;i<NFOOTMAX;i++)
 	   if ( strcmp(x->footer[i],"") ) 
		if (!noSOL || !FindWord(x->footer[i],"SOL") ) fprintf(fout,"%s",x->footer[i]);

  fclose(fout);
}

int ReadDefines(char *filename, char *dirname, struct top_s *x, int n)
{
	int k,j;
	float i[6];
	FILE *fp;
	char aux[1000],filename2[100],filenamec[200],s1[100];

	// open current file
	printf("%s\n",filename);
	fp = fopen(filename,"r");	// try to open in current directory
	if (!fp) 					// if cannot, try in directory "dirname"
	{
		sprintf(filenamec,"%s/%s",dirname,filename);
		fp = fopen(filenamec,"r");
		//if (!fp) Error2("Error in opening top file for searching #define; file is",filename);
	}

	// read file content
	while ( fgets(aux,1000,fp) )
	{
		// find the #define and assign content
		if (!strncmp(aux,"#define",7))
		{
			for (j=0;j<6;j++) i[j]=NONDEF;
			k = sscanf(aux,"#define %s %f %f %f %f %f %f",s1,&i[0],&i[1],&i[2],&i[3],&i[4],&i[5]);
			if (k>1)
			{
				strcpy(x->label[n],s1);
				for (j=0;j<6;j++) (x->labelc)[n][j] = i[j];
				if (DEBUG>0) printf("\t%3d  %s\n",n,x->label[n]);
				n++;
			}
		}
		// find nested includes
		else if (!strncmp(aux,"#include",8))
		{
			sscanf(aux,"#include \"%s",filename2);
			filename2[strlen(filename2)-1]='\0';
			n = ReadDefines(filename2,dirname,x,n);
		}
	}

	fclose(fp);
	return n;
}

int ReadDefinesSimple(char *filename, struct top_s *x, int verbose)
{
	int k,j,n=0;
	float i[6];
	FILE *fp;
	char aux[1000],s1[100];

	// open current file
	if (verbose) printf("  Read #define in %s\n",filename);
	fp = fopen(filename,"r");	// try to open in current directory
	if (!fp) Error2("Error in opening top file for searching #define; file is: ",filename);


	// read file content
	while ( fgets(aux,1000,fp) )
	{
		// find the #define and assign content
		if (!strncmp(aux,"#define",7))
		{
			for (j=0;j<6;j++) i[j]=NONDEF;
			k = sscanf(aux,"#define %s %f %f %f %f %f %f",s1,&i[0],&i[1],&i[2],&i[3],&i[4],&i[5]);
			if (k>1)
			{
				strcpy(x->label[n],s1);
				for (j=0;j<6;j++) (x->labelc)[n][j] = i[j];
				if (DEBUG>0) printf("\t%3d  %s\n",n,x->label[n]);
				n++;
			}
		}

	}

	if (verbose) printf("  Read %d label definitions in top.\n",n);

	fclose(fp);
	return n;
}

void ShiftAtomsTop(struct top_s *x, int satoms, int sres, int scgnr)
{
   int i;

   // shift in [ atoms ]
   for (i=0;i<x->natoms;i++)
   {
	   ((x->top_atoms)+i)->nr += satoms;
	   ((x->top_atoms)+i)->resnr += sres;
	   ((x->top_atoms)+i)->cgnr += scgnr;
   }

   // shift in [ bonds ]
   for (i=0;i<x->nbonds;i++)
    {
 	   ((x->top_bonds)+i)->i1 += satoms;
 	   ((x->top_bonds)+i)->i2 += satoms;
    }

   // shift in [ pairs ]
   for (i=0;i<x->npairs;i++)
    {
 	   ((x->top_pairs)+i)->i1 += satoms;
 	   ((x->top_pairs)+i)->i2 += satoms;
    }

   // shift in [ angles ]
   for (i=0;i<x->nangles;i++)
    {
 	   ((x->top_angles)+i)->i1 += satoms;
 	   ((x->top_angles)+i)->i2 += satoms;
 	   ((x->top_angles)+i)->i3 += satoms;
    }

   // shift in [ dihedrals ]
    for (i=0;i<x->ndihedrals;i++)
     {
  	   ((x->top_dihedrals)+i)->i1 += satoms;
  	   ((x->top_dihedrals)+i)->i2 += satoms;
  	   ((x->top_dihedrals)+i)->i3 += satoms;
  	   ((x->top_dihedrals)+i)->i4 += satoms;
     }

    // shift in [ impropers ]
      for (i=0;i<x->nimpropers;i++)
       {
    	   ((x->top_impropers)+i)->i1 += satoms;
    	   ((x->top_impropers)+i)->i2 += satoms;
    	   ((x->top_impropers)+i)->i3 += satoms;
    	   ((x->top_impropers)+i)->i4 += satoms;
       }
     // shift in [ exclusions ]
      for (i=0;i<x->nexclusions;i++)
       {
    	   ((x->top_exclusions)+i)->i1 += satoms;
    	   ((x->top_exclusions)+i)->i2 += satoms;
       }


}

void CopyTop(struct top_s *t1, struct top_s *t2)
{
	int i,j;

	CopyTopAtoms( t1->top_atoms, t2->top_atoms, t1->natoms);
	CopyTopGen( t1->top_bonds, t2->top_bonds, t1->nbonds);
	CopyTopGen( t1->top_pairs, t2->top_pairs, t1->npairs);
	CopyTopGen( t1->top_angles, t2->top_angles, t1->nangles);
	CopyTopGen( t1->top_dihedrals, t2->top_dihedrals, t1->ndihedrals);
	CopyTopGen( t1->top_impropers, t2->top_impropers, t1->nimpropers);
	CopyTopGen( t1->top_exclusions, t2->top_exclusions, t1->nexclusions);
	
	t2->natoms = t1->natoms;
	t2->nbonds = t1->nbonds;
	t2->npairs = t1->npairs;
	t2->nangles = t1->nangles;
	t2->ndihedrals = t1->ndihedrals;
	t2->nimpropers = t1->nimpropers;
	t2->nexclusions = t1->nexclusions;

	for (i=0;i<50;i++)
	 strcpy(t2->label[i],t1->label[i]);
	for (i=0;i<NLABELS;i++)
		for (j=0;j<6;j++)
	       t2->labelc[i][j] = t1->labelc[i][j];		

	for (i=0;i<NHEADMAX;i++) strcpy(t2->header[i],t1->header[i]);
	for (i=0;i<NFOOTMAX;i++) strcpy(t2->footer[i],t1->footer[i]);

}

void CopyTopGen(struct top_gen_s *t1, struct top_gen_s *t2, int n)
{
	int i;

	for (i=0;i<n;i++)
	{
		(t2+i)->i1 = (t1+i)->i1;
		(t2+i)->i2 = (t1+i)->i2;
		(t2+i)->i3 = (t1+i)->i3;
		(t2+i)->i4 = (t1+i)->i4;
		(t2+i)->funct = (t1+i)->funct;
		(t2+i)->c0 = (t1+i)->c0;
		(t2+i)->c1 = (t1+i)->c1;
		(t2+i)->c2 = (t1+i)->c2;
		(t2+i)->c3 = (t1+i)->c3;
		(t2+i)->c4 = (t1+i)->c4;
		(t2+i)->c5 = (t1+i)->c5;
		(t2+i)->move = (t1+i)->move;
	}
	
}

void CopyTopAtoms(struct top_atoms_s *t1, struct top_atoms_s *t2, int n)
{
	int i;

	for (i=0;i<n;i++)
	{
		(t2+i)->nr = (t1+i)->nr;
		strcpy((t2+i)->type, (t1+i)->type);
		(t2+i)->resnr = (t1+i)->resnr;
		strcpy((t2+i)->residue, (t1+i)->residue);
		strcpy((t2+i)->atom, (t1+i)->atom);
		(t2+i)->cgnr = (t1+i)->cgnr;
		(t2+i)->charge = (t1+i)->charge;
		(t2+i)->mass = (t1+i)->mass;
		
	}
}

void LinkTop(struct top_s *top1, struct top_s *top2, int i1, int i2, int i3, int i4, int i5, int nt1, int verbose)
{
	int i,j;
	struct network_s net;
	
	if (verbose) printf("\nPrint linking information:   (nt1=%d)\n",nt1);
	
	net = BuildNetwork(top1,verbose);

    // bonds
	if (verbose) printf("Bonds:\n");
	((top1->top_bonds)+ top1->nbonds )->i1 = i3;
	((top1->top_bonds)+ top1->nbonds )->i2 = i4 + nt1;	// nt1 is the # of atoms of the first fragment
	((top1->top_bonds)+ top1->nbonds )->funct = 1;
	if (verbose) printf("\t%d %d-%d\n",top1->nbonds,((top1->top_bonds)+ top1->nbonds )->i1,((top1->top_bonds)+ top1->nbonds )->i2);
	top1->nbonds ++;

	//pairs
	if (verbose) printf("Pairs:\n");
	for (i=0;i<net.nvic[i2];i++)						// 1-4 of i4
		if (net.vic[i2][i] != i3) 
		{
			((top1->top_pairs)+ top1->npairs )->i1 = i4 + nt1;
			((top1->top_pairs)+ top1->npairs )->i2 = net.vic[i2][i];
			((top1->top_pairs)+ top1->npairs )->funct = 1;
			if (verbose) printf("\t%d %d-%d\n",top1->npairs,((top1->top_pairs)+ top1->npairs )->i1,((top1->top_pairs)+ top1->npairs )->i2);
			top1->npairs ++;
		}
	
	for (i=0;i<net.nvic[i5+nt1];i++)						// 1-4 of i3
		if (net.vic[i5+nt1][i] != i4+nt1) 
		{
			((top1->top_pairs)+ top1->npairs )->i1 = i3;
			((top1->top_pairs)+ top1->npairs )->i2 = net.vic[i5+nt1][i];
			((top1->top_pairs)+ top1->npairs )->funct = 1;
			if (verbose) printf("\t%d %d-%d\n",top1->npairs,((top1->top_pairs)+ top1->npairs )->i1,((top1->top_pairs)+ top1->npairs )->i2);
			top1->npairs ++;
		}
	
	for (i=0;i<net.nvic[i3];i++)						// neighbours of i3 with those of i4
		if (net.vic[i3][i] != i4+nt1) 
			for (j=0;j<net.nvic[i4+nt1];j++)						
				if (net.vic[i4+nt1][j] != i3) 
				{ 
					((top1->top_pairs)+ top1->npairs )->i1 = net.vic[i3][i];
					((top1->top_pairs)+ top1->npairs )->i2 = net.vic[i4+nt1][j];
					((top1->top_pairs)+ top1->npairs )->funct = 1;
					if (verbose) printf("\t%d %d-%d\n",top1->npairs,((top1->top_pairs)+ top1->npairs )->i1,((top1->top_pairs)+ top1->npairs )->i2);
					top1->npairs ++;
				}
	
	// angles	
	if (verbose) printf("Angles:\n");	
	for (i=0;i<net.nvic[i4+nt1];i++)						// i3-i4-vic
			if (net.vic[i4+nt1][i] != i3) 
			{
				((top1->top_angles)+ top1->nangles )->i1 = i3;
				((top1->top_angles)+ top1->nangles )->i2 = i4+nt1;
				((top1->top_angles)+ top1->nangles )->i3 = net.vic[i4+nt1][i];
				((top1->top_angles)+ top1->nangles )->funct = 1;
				if (verbose) printf("\t%d %d-%d-%d\n",top1->nangles,((top1->top_angles)+ top1->nangles )->i1,
						((top1->top_angles)+ top1->nangles )->i2, ((top1->top_angles)+ top1->nangles )->i3);
				top1->nangles ++;
			}
	
	for (i=0;i<net.nvic[i3];i++)						// i4-i3-vic
			if (net.vic[i3][i] != i4+nt1) 
			{
				((top1->top_angles)+ top1->nangles )->i1 = i4+nt1;
				((top1->top_angles)+ top1->nangles )->i2 = i3;
				((top1->top_angles)+ top1->nangles )->i3 = net.vic[i3][i];
				((top1->top_angles)+ top1->nangles )->funct = 1;
				if (verbose) printf("\t%d %d-%d-%d\n",top1->nangles,((top1->top_angles)+ top1->nangles )->i1,
						((top1->top_angles)+ top1->nangles )->i2, ((top1->top_angles)+ top1->nangles )->i3);
				top1->nangles ++;
			}
	
	// dihedrals
	if (verbose) printf("Dihedrals:\n");
	for (i=0;i<net.nvic[i3];i++)						// neighbours of i3 with those of i4
		if (net.vic[i3][i] != i4+nt1) 
			for (j=0;j<net.nvic[i4+nt1];j++)						
				if (net.vic[i4+nt1][j] != i3) 
				{
					((top1->top_dihedrals)+ top1->ndihedrals )->i1 = net.vic[i3][i];
					((top1->top_dihedrals)+ top1->ndihedrals )->i2 = i3;
					((top1->top_dihedrals)+ top1->ndihedrals )->i3 = i4+nt1;
					((top1->top_dihedrals)+ top1->ndihedrals )->i4 = net.vic[i4+nt1][j];
					((top1->top_dihedrals)+ top1->ndihedrals )->funct = 3;
					if (verbose) printf("\t%d %d-%d-%d-%d\n",top1->ndihedrals,((top1->top_dihedrals)+ top1->ndihedrals )->i1,
							((top1->top_dihedrals)+ top1->ndihedrals )->i2,((top1->top_dihedrals)+ top1->ndihedrals )->i3,
							((top1->top_dihedrals)+ top1->ndihedrals )->i4);
					top1->ndihedrals ++;
				}
	
	for (i=0;i<net.nvic[i5+nt1];i++)						// i3 with neighbours of i5
		if (net.vic[i5+nt1][i] != i4+nt1) 
				{
					((top1->top_dihedrals)+ top1->ndihedrals )->i1 = i3;
					((top1->top_dihedrals)+ top1->ndihedrals )->i2 = i4+nt1;
					((top1->top_dihedrals)+ top1->ndihedrals )->i3 = i5+nt1;
					((top1->top_dihedrals)+ top1->ndihedrals )->i4 = net.vic[i5+nt1][i];
					((top1->top_dihedrals)+ top1->ndihedrals )->funct = 3;
					if (verbose) printf("\t%d %d-%d-%d-%d\n",top1->ndihedrals,((top1->top_dihedrals)+ top1->ndihedrals )->i1,
							((top1->top_dihedrals)+ top1->ndihedrals )->i2,((top1->top_dihedrals)+ top1->ndihedrals )->i3,
							((top1->top_dihedrals)+ top1->ndihedrals )->i4);
					top1->ndihedrals ++;
				}
	
	for (i=0;i<net.nvic[i2];i++)						// i4 with neighbours of i2
		if (net.vic[i2][i] != i3) 
				{
					((top1->top_dihedrals)+ top1->ndihedrals )->i1 = i4+nt1;
					((top1->top_dihedrals)+ top1->ndihedrals )->i2 = i3;
					((top1->top_dihedrals)+ top1->ndihedrals )->i3 = i2;
					((top1->top_dihedrals)+ top1->ndihedrals )->i4 = net.vic[i2][i];
					((top1->top_dihedrals)+ top1->ndihedrals )->funct = 3;
					if (verbose) printf("\t%d %d-%d-%d-%d\n",top1->ndihedrals,((top1->top_dihedrals)+ top1->ndihedrals )->i1,
							((top1->top_dihedrals)+ top1->ndihedrals )->i2,((top1->top_dihedrals)+ top1->ndihedrals )->i3,
							((top1->top_dihedrals)+ top1->ndihedrals )->i4);
					top1->ndihedrals ++;
				}
	return;
}

struct network_s BuildNetwork(struct top_s *x, int verbose)
{
   int i,j;
   struct network_s y;
   
   for (i=1;i<NATOMSMAX+1;i++) y.nvic[i]=0;
   
   for (i=0;i<x->nbonds;i++)
   {
	   y.vic[((x->top_bonds)+i)->i1][y.nvic[((x->top_bonds)+i)->i1]] = ((x->top_bonds)+i)->i2;
	   y.nvic[((x->top_bonds)+i)->i1] ++;
	   y.vic[((x->top_bonds)+i)->i2][y.nvic[((x->top_bonds)+i)->i2]] = ((x->top_bonds)+i)->i1;
	   y.nvic[((x->top_bonds)+i)->i2] ++;
   }

   if (verbose)
   {
	   printf(" Network of bonds:\n");
	   for (i=1;i<x->natoms+1;i++)
	   {
		   printf("\t%d\t",i);
		   for (j=0;j<y.nvic[i];j++)
		     printf("%d ",y.vic[i][j]);
		   printf("\n");
	   }
   }
   
   return y;
}
