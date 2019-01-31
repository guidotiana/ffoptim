#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/types.h>

#define NVER		0
#define NSUBVER		6
#define NFRAMESMAX	15000
#define NRESMAX		3000
#define NATOMSMAX	5000
#define NBONDSMAX	10000
#define NPAIRSMAX	10000
#define NANGLESMAX	10000
#define NDIHEDRALSMAX 	10000
#define NIMPROPERSMAX 	10000
#define NEXCLUSIONSMAX 	10000
#define NLABELS 	200
#define NVICMAX  	5
#define NHEADMAX 	1000
#define NFOOTMAX 	1000
#define NINCLMAX 	200
#define DEBUG 		0
#define NGROUPMAX   100

#define NONDEF -9999.0
#define EPSILON 1E-19
#define PI 3.1415
#define KRELAX 0.56		// nm^6/ms^2 


struct parm_s
{
	char topfile[500];
	char itpdir[500];
	int verbose;
	int niter;
	char grobindir[1000];
	char grofile[500];
	char groprotfile[500];
	char mdpfile[500];
	char chi2file[500];
	char emdpfile[500];
	char origtopfile[500];
	char nonbondedFile[500];
	char bondedFile[500];
	int noptim;
	int noptauc;
	int ch_pairs;
	int ch_dihedrals;
	double sigma_pairs_u;
	double sigma_pairs_r;
	double sigma_dihe;
	double temp;
	int nthreads;
	int replex;
	double temperatures[500];
	int ntemp;	
	double noe0;
	int correl;
	int combrule;
	int restart;
	double tauc;
	double omega;
	double mixingt;
	int nochangeH;
	int RnochangeH;
	int filltop;
	double elimit0;
	double elimit1;

	int *ires1;
	int *ires2;
	int *restype;
	double *resexp;
	double *ressigma;
	int nres;
        int irun0;
	int nframesmax;
    char **alabel1;
    char **alabel2;

	int *iv0;
	double *v0;
	int nv0;
	int absnoe;
	int iskip;
	int noSOL;
	char mdpprotegrps[80];
	char mdpprotxtcgrps[80];
	char energytype[200];
	int onlyc6;
	int rescalenoe0;
	int nodih0change;
	int nice;
	int moveAll;
	int moveAllPairs;
	int moveAllDih;
	int movev0;
	int movea0;
	char testtraj[500];
	int lowerMatrix;

	int *taucatoms1;
	int *taucatoms2;
	double *taucatomsv;
	int ntaucatoms;
	char hostfile[30];
	double uNOE;
    
    char agrouplabel[30][NGROUPMAX];
    int agroup[NGROUPMAX];
    int nagroup;
};

struct top_gen_s
{
	int i1;
	int i2;
	int i3;
	int i4;
	int funct;
	float c0;
	float c1;
	float c2;
	float c3;
	float c4;
	float c5;
	int move;
};

struct top_atoms_s
{
	int nr;
	char type[30];
	int resnr;
	char residue[5];
	char atom[8];
	int cgnr;
	float charge;
	float mass;
};


struct top_s
{
	struct top_atoms_s *top_atoms;
	struct top_gen_s *top_bonds;
	struct top_gen_s *top_pairs;
	struct top_gen_s *top_angles;
	struct top_gen_s *top_dihedrals;
	struct top_gen_s *top_impropers;
	struct top_gen_s *top_exclusions;
	int natoms;
	int nbonds;
	int npairs;
	int nangles;
	int ndihedrals;
	int nimpropers;
	int nexclusions;

	char label[NLABELS][50];	// parameters defined with a #define: label...
	float labelc[NLABELS][6];		// ... and corresponding values

	char header[NHEADMAX][100];
	char footer[NFOOTMAX][100];
};

struct vector
{
        double x;
        double y;
        double z;
};

struct angles {
        double alpha;
        double beta;
        double blength;
};

struct network_s
{
	int nvic[NATOMSMAX];
	int vic[NATOMSMAX][NVICMAX];
};

struct NOErelax_s
{
	int *Hlist;
	int nh;
	float ***dist;
	double k1;
	double k2;
	int *backH;			// lookback table from atom to h table
	double **kAtoms;		// prefactors of NOE atomwise
};

// ffoptim.c
void SimulateGromacs(int iter, struct parm_s *parm, char *topi);
int GetGromacsEnergies(struct parm_s *parm, struct top_s *top, double *eold, struct NOErelax_s *NOErelax, int recalculate);
void AverRestrain(double **restrain, double *average, int nframes, int nres);
double Chi2(double *x, double *xexp, double *sigma, int n, double *k);
float OptimizePotential(struct top_s *top, struct top_s *topnew, int nframes, struct parm_s *parm, 
			double **restrain, double *resnew, double *eold, double *enew, double *resnew2, struct top_s *orig);
void RandomChangeEnergy(struct parm_s *parm, struct top_s *top, struct top_s *origtop );
void ReweightRestrain(double *resnew, double **restrain, double *eold, double *enew, struct parm_s *parms, int nframes, 
			int iskip, double *resnew2);
double Correlation(double *x, double *xexp, int n);
void GetRestrain(struct parm_s *parm, double **restrains, int nframes, struct NOErelax_s *NOErelax);
void PrintInitialRestrains(struct parm_s *parm, double **restrain, FILE *fout,struct NOErelax_s *NOErelax);
void Welcome(FILE *fp);
void Help(void);
void PrintGromacsOutput(int iter, double chi2, struct parm_s *parm, FILE *fchi2, double *resnew, double *resnew2);
void TestTraj(struct parm_s *parm);
void PrintTime(FILE *fout);
void AverageRestrain(double *resnew, double **restrain, struct parm_s *parm, int nframes);


// io.c
struct parm_s *ReadParms(char *fname);
void ReadParD(char *s, char key[20], int *par);
void ReadParL(char *s, char key[20], long *par);
void ReadParF(char *s, char key[20], double *par);
void ReadParS(char *s, char key[20], char *par);
void ReadParN(char *s, char key[20], int *par);
int FindKeyword(char *string, char *keyword);
void CheckFiles(struct parm_s *x);
int ParseXVG(char *fname, double *x);
void RenameIndex(void);
char *trimwhitespace(char *str);
void ChangeTmdp(char *mdpfile, int i, double *t);
void CheckExplosion(void);
int ParseXVGf(char *fname, float *x);
int FindWord(char *string, char *word);
void UpdateMdp(char *fmdp, char *omdp, char *xtcgrps, char *energygrps);

// misc.c
void Error(char *x);
void Error2(char *x, char *y);
double **AlloDoubleMatrix(int l, int m);
int *AlloInt(int n);
double *AlloDouble(int n);
int **AlloIntMatrix(int l, int m);
int irand(int r);
double frand();
long Randomize(int n);
int IsZero(double x);
char **AlloString(int n, int strl);


// toptools.c
struct top_s *AlloTop(int natoms, int nbonds, int npairs, int nangles,
		int ndihedrals, int nimpropers, int nexclusions, int verbose);
struct top_s *ReadTop(char *fname, char *itpdir, int verbose);
struct network_s BuildNetwork(struct top_s *x, int verbose);
void PrintTop(char *, struct top_s *x, int verbose, int noSOL);
void CopyTop(struct top_s *t1, struct top_s *t2);
void CopyTopGen(struct top_gen_s *t1, struct top_gen_s *t2, int n);
void CopyTopAtoms(struct top_atoms_s *t1, struct top_atoms_s *t2, int n);
void FillTop(struct top_s *x, int combrule, int nochangeH);
void InsertNonbondedTop(struct top_s *x, int combrule, int debug, char *nonbondedFile);
void InsertDihedralsTop(struct top_s *x, int verbose, char *bondedFile);

// relaxation.c
int ShouldRelax(struct parm_s *parm);
void FullRelaxation(struct NOErelax_s *x, double **res, struct parm_s *parm, int nframes);
struct NOErelax_s *AlloNOErelax(int natoms, int nh, int verbose);
void AlloNOErelaxDist(struct NOErelax_s *x, int nframes, int verbose);
double Prefactor(int diag, double tauc, double omega, int verbose);
int SelectHydrogens(struct top_s *top, int *h, struct parm_s *parm, int *backH);
int LookbackHydrogens(struct top_s *top, int *h, struct parm_s *parm, int *backH, int *iv0, int nv0);
void MoveV0(struct top_s *x,  struct parm_s *p);
void MoveA0(struct top_s *x,  struct parm_s *p);
void PrefactorAtom(double **x, struct parm_s *parm, int nh, int *Hlist);
float OptimizeTauc(int nframes, struct parm_s *parm, double **restrain, double *resnew,struct NOErelax_s *NOErelax, double *eold, double *enew, double *resnew2);
void AdduNOE(struct parm_s *parm, struct NOErelax_s *NOErelax);

