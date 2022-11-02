/* 
   Program Evol_div
   Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa (CSIC-UAM)
   ubastolla@cbm.csic.es
   Reads a multiple alignment and computes contact divergence and other
   structure comparison measures.
   Estimates or reads outgroups from a tree and computes clock violations

   INPUT: file with multiple alignment in FASTA format (at least 2 seqs)
   Protein names must be names of PDB files.
   The first line may be PDBDIR=<directory of PDB files>
   (default: current directory)

   OUTPUT: For each protein pair, structural scores printed are
   Contact_divergence, contact overlap, TM score (default no)

*/

#include "Contact_divergence_aux.h"
#include "D_Cont.h"
#include "protein.h"
#include "cont_list.h"
#include "allocate.h"
//#include "normalization.h"
#include "tm_score.h"
#include "read_structures.h"
#include "tree.h"
#include "CV_statistics.h"
#include "align_ss.h"
#include "PC_ali.h"
//#include "consensus_msa.h"
#include "clique_msa.h"
#include "Print_pairwise.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float S0=0.06;    // for Tajima-Nei divergence
float TM0=0.167; //  TM score of unrelated proteins

//float alpha=0.7; // Exponent for computing clock violations
int PRINT_SIM=0;  // Print matrix of similarities?
int PRINT_DIV=0;  // Print matrix of divergences?
int PRINT_CV=1;   // Print clock violations?
int PRINT_PAIR=1; // Print pairwise alignments?
int PRINT_CLIQUE=1; // Print multiple alignments based on cliques of pairwise?
int PRINT_AVE=1;  // Print average values as a function of d?
int PRINT_GAP=1; // Examine relationship between Triangle Inequality and gaps?
int ALI_SS=0;     // Modify Input alignment considering sec.str.?
int ALI_CO=1;     // Modify Input alignment targetting CO?
int ALI_TM=1;     // Modify Input alignment targetting TM?
int ALI_PC=1;     // Modify Input alignment targetting PC?
int ITMAX=3; // Rounds of optimization of target similarity scores

// Parameters for secondary structure correction
int SS_MULT=0;    // Correct multiple or pairwise alignment?
int SHIFT_MAX=4;  // Maximum allowed shift when correcting for sec.str.
float HUGE=10;    // Very large divergence
//float TM0=0.17;   // TM score of unrelated proteins
int NORM, NORMA;

#define EXT_DIFF ".diff" // Extension for comparison of alignments
#define EXT_DIV ".div"   // Extension for divergence output file
#define EXT_SIM ".sim"   // Extension for similarity output file
#define EXT_CV  ".cv"    // Extension for clock-violation output file

#define IJ_MIN_DEF 3
#define CONT_TYPE_DEF 'c'
#define CONT_THR_DEF 4.5

float fun_high=0.90;
float fun_low=0.70;
int diff_thr=6;  // Minimum number of differences for joining sequences

int IJ_MIN=IJ_MIN_DEF;        // Only contacts with |i-j|>=IJ_MIN
float CONT_THR=CONT_THR_DEF;
char CONT_TYPE=CONT_TYPE_DEF;  //a=alpha b=beta c=all atoms
char CONT_STRING[80];
char CODENAME[40]="main_Evol_div.c";
int CONT_DEF=1;

struct Prot_input{
  char name[80];
  char chain;
  char *seq;
  int len;
};

// Statistics
#define IBIN 10
#define IBIN1 11
#define ITYPE 6
int ini_sum=0;
double sum_SS_SI[2], sum_SS_TM[2], sum_SS_CO[2];
double *sum_SS_SI_TM[2], *sum_SS_SI_CO[2], *sum_SS_TM_CO[2];
double // aligned:
  sum_all_ct[ITYPE][2][IBIN1],
  sum_ali_ct[ITYPE][2][IBIN1],
  sum_aaid_ct[ITYPE][2][IBIN1],
  sum_sup_ct[ITYPE][2][IBIN1],
  sum_aaid_sup_ct[ITYPE][2][IBIN1],
  sum_noid_nosup_ct[ITYPE][2][IBIN1],
  sum_ali_cont_ct[ITYPE][2][IBIN1],
  sum_id_cont_ct[ITYPE][2][IBIN1],
  // Neighbors:
  sum_neigh_ali_ct[ITYPE][2][IBIN1],
  sum_neigh_noali_ct[ITYPE][2][IBIN1],
  sum_neigh_noali_aaid_ct[ITYPE][2][IBIN1],
  sum_sup_noneigh_ct[ITYPE][2][IBIN1],
  sum_shift_sup_noneigh_ct[ITYPE][2][IBIN1],
  // contacts
  sum_ali_cont_sup[ITYPE][3][IBIN1], sum_id_cont_sup[ITYPE][3][IBIN1],
  sum_ali_cont_aaid[ITYPE][3][IBIN1], sum_id_cont_aaid[ITYPE][3][IBIN1],
  sum_ali_cont_aaid_sup[ITYPE][3][2][IBIN1],
  sum_id_cont_aaid_sup[ITYPE][3][2][IBIN1];
double c_ave; // Mean number of contacts per residue
double CO1[ITYPE], TM1[ITYPE], SI1[ITYPE], nali1[ITYPE], PC1[ITYPE];
double CO_diff1[ITYPE][ITYPE], CO_diff2[ITYPE][ITYPE];
double SI_diff1[ITYPE][ITYPE], SI_diff2[ITYPE][ITYPE];
double ali_diff1[ITYPE][ITYPE],ali_diff2[ITYPE][ITYPE];
double TM_diff1[ITYPE][ITYPE], TM_diff2[ITYPE][ITYPE];
double PC_diff1[ITYPE][ITYPE], PC_diff2[ITYPE][ITYPE];
//double PC_load[4]={0.83,0.80,0.94,0.95}; // ali, SI, TM, CO
double PC_load[4]={0.86,0.81,0.95,0.94};
double PC0_1, PC_norm;
char *ali_name[ITYPE];
int opt_PC[ITYPE];
int Diff_opt[ITYPE];
double pairs=0;
double sum_norm_c;

void help(char *pname);
void Get_input(char *file_ali, char *file_ali_str, char *file_fun,
	       char *name, char *PDBDIR, char *PDBEXT, char *OUTG,
	       int *NORM, int *ALI_SS, int *SHIFT_MAX, int *SS_MULT,
	       int *PRINT_SIM, int *PRINT_CV, int *PRINT_DIV,
	       int *PRINT_PAIR, int *PRINT_CLIQUE,
	       int argc, char **argv);
int Get_alignment_old(struct Prot_input **Prot_input, int *Nali,
		      char *PDBDIR, char *PDBEXT, int *PRINT_SIM,
		      char *file_ali);
int Get_alignment(struct Prot_input **Prot_input, int *Nali, char *file_ali,
		  char *PDBDIR, char *PDBEXT);
int *Match_alignments(struct Prot_input *Prot1,
		      struct Prot_input *Prot2, int N);
int Find_prot(char *name, struct Prot_input *Prot, int *index, int N);
float **Read_function(char *file_fun, struct Prot_input *Prot,
		      int *index, int N);
void Sec_str_AA_propensity(struct protein *pi, struct protein *pj,
			   int *ali, int nali, float PC);
void Print_propensities(char *nameout);
double ***Prop_secstr=NULL, ***Prop_AA=NULL, *DP_bin=NULL, *DP_norm=NULL;
char *SS_name[4], AA_name[22];
int SS_code(char ss);
int AA_code(char aa);

// Auxiliary
void Write_ali(int ***Ali_pair, int i, int j, int *L_seq, int *ali_PC);
char **Assign_Seq(char ***name_seq, int *len_seq,
		  struct Prot_input *Prot_in, int N_seq,
		  int *rep_str, int *i_seq, int N_ali);
int ***Select_alis(int ***Ali_pair, int N_seq, int *len_seq,int *rep_str);
void Set_contact_type();
int Count_AA(char *seq, int N_ali);
int Seq_differences(int *id, char *seq1, char *seq2, int N_ali);
int Count_gaps(char *seq1, char *seq2, int N_ali);
float Min_dist(int i, int j, int **conformation, int *N_conf, float **div);
float Min_CV(int i, int j, int k, int **conformation, int *N_conf, float **div);
void Change_conformations(int **conformation, int N_seq, int *N_conf,
			  int *ali_str);
void Get_file_name(char *name, char *file);
float Divergence(float SI, float S0, float HUGE);
void Print_seq(int *ali, char *seq, int N_ali);
void Write_identity(FILE *file_id, int type,
		    struct protein *proti, struct protein *protj,
		    int *ali_ij, int *id_aa, int *id_sup, int *shift,
		    int *neigh_ali,int *neigh_noali,
		    int *neigh_noali_aaid,
		    int *ali_cont_ct, int *id_cont_ct,
		    int *ali_cont_sup, int *id_cont_sup,
		    int *ali_cont_aaid, int *id_cont_aaid);
void Summary_identical(FILE *file_id);
void Count_noali(int *noali, int *ali,  int *ncont, int n, float c_ave);
float Mean_freq(double *se, double *re2, double sum, double tot);
int Count_ali(int *ali, int len);
double Sum_bins(double *bin);
float Ave_se(double *se, double sum1, double sum2, int n, float nind);

// Eliminate pairs with SI<SI_thr=2*S0=0.12
float SI_thr=0.10; //0.12

void Print_ali_ss(int *ali_i, char *ss_i, int *ali_j, char *ss_j, int N_ali);

/*********************************************************************
                          MAIN routine
**********************************************************************/
int main(int argc, char **argv)
{
  // INPUT
  NORM=0; // Normalize with minimum (0) maximum (1) or geometric mean (2)?
  NORMA=0; // Normalize with nali (1) or with NORM (0)?
  char PDB_DIR[100]="./", PDB_EXT[10]="", OUTG[10]="TN";
  char file_ali[200], file_ali_str[200]="\0", file_fun[200], name_in[80]="";
  Get_input(file_ali, file_ali_str, file_fun, name_in, PDB_DIR, PDB_EXT,
	    OUTG, &NORM, &ALI_SS, &SHIFT_MAX, &SS_MULT,
	    &PRINT_SIM, &PRINT_CV, &PRINT_DIV, &PRINT_PAIR, &PRINT_CLIQUE,
	    argc, argv);

  // Alignments
  int N_ali=0, N_ali_str=0, i, j;
  struct Prot_input *Prot_in, *Prot2;
  int N_prot=Get_alignment(&Prot_in, &N_ali, file_ali, PDB_DIR, PDB_EXT);

  // Read structure alignment, if any
  int Np2=Get_alignment(&Prot2, &N_ali_str, file_ali_str, PDB_DIR, PDB_EXT);
  if(Np2 && (Np2!=N_prot)){
    printf("WARNING, structure alignment contains %d proteins instead of %d\n",
	   Np2, N_prot);
    printf("Discarding structure alignment\n"); Np2=0;
  }
  int *ali_str=NULL;
  if(Np2){
    printf("Matching multiple sequence and multiple structure alignments\n");
    ali_str=Match_alignments(Prot_in, Prot2, N_prot);
  }

  // What to print ? 
  //if((PRINT_CV==0)&&(PRINT_SIM==0))PRINT_DIV=1;


  /**************   READ PROTEIN STRUCTURES  ******************/
  // Read PDB files and compute contact matrices
  Set_contact_type();
  printf("Contact type: %c Threshold: %.2f A |i-j|>%d\n",
	 CONT_TYPE, CONT_THR,IJ_MIN);
  printf("Looking for PDB files in directory %s\n", PDB_DIR);

  // Alignments
  int **Prot_ali=Allocate_mat2_i(N_prot, N_ali);
  int **Prot_ali_str=NULL;
  if(ali_str)Prot_ali_str=Allocate_mat2_i(N_prot, N_ali_str);
  int Max_ali=N_ali; if(N_ali_str>Max_ali)Max_ali=N_ali_str;

  // Proteins prots
  struct protein prots[N_prot], *prot=prots;
  int N_pdb=0, i_seq[N_prot], L_seq[N_prot];
  for(i=0; i<N_prot; i++){
    char pdb[90];
    sprintf(pdb, "%s%s", Prot_in[i].name, PDB_EXT);
    if(Read_PDB_compress(prot, pdb, &(Prot_in[i].chain), PDB_DIR)>0){
      if(Align_seq(Prot_ali[N_pdb], N_ali,
		   Prot_in[i].seq, prot->aseq, prot->len)<0)continue;
      if(Prot_ali_str && 
	 (Align_seq(Prot_ali_str[N_pdb], N_ali_str,
		    Prot2[ali_str[i]].seq, prot->aseq, prot->len)<0)){
	ali_str[i]=-1;
      }
      int NC=Compute_contact_list(prot, CONT_TYPE, CONT_THR, IJ_MIN);
      printf("%d contacts\n", NC);
      L_seq[N_pdb]=Count_AA(Prot_in[i].seq, N_ali); //i_str[i]=N_pdb;
      Prot_in[i].len=L_seq[N_pdb];
      if(L_seq[N_pdb]!=prot->len){
	printf("WARNING, different n.residues in ali (%d) and PDB (%d)\n",
	       L_seq[N_pdb], prot->len);
	for(j=0; j<N_ali; j++)
	  if(Prot_in[i].seq[j]!='-')printf("%c",Prot_in[i].seq[j]);
	printf("\n");
	for(j=0; j<prot->len; j++){printf("%c",prot->aseq[j]);} printf("\n");
	//exit(8);
      }
      i_seq[N_pdb]=i; N_pdb++; prot++;
    }
  }
  printf("%d proteins read out of %d listed in %s\n", N_pdb, N_prot, file_ali);
  if(N_pdb<2){
    printf("ERROR, fewer than 2 proteins found\n"); exit(8);
  }

  // Average number of contacts per residue
  c_ave=0; double c_norm=0;
  for(i=0; i<N_pdb; i++){
    struct protein *proti=prots+i;
    for(j=0; j<proti->len; j++)c_ave+=proti->ncont[j];
    c_norm+=proti->len;
  }
  c_ave/=c_norm;

  /**************************************
      Secondary structure based alignment
   *******************************************/
  int **Prot_ali_ss=NULL, **Prot_ali_ij=NULL, NS=0;
  int PRINT_ALI=1; struct protein **prot_ij=NULL;
  if(SS_MULT){NS=N_pdb;}else{NS=2;}
  if(ALI_SS){
    if(Set_sec_str(prots, N_pdb)<0){ // Change - into c
      printf("WARNING, secondary structure information not found\n");
      ALI_SS=0; goto end_ss;
    }
    int s_not=Test_notali(Prot_ali,prots,N_pdb,N_ali);
    if(s_not)printf("WARNING, %d proteins had not aligned residues\n",s_not);
    //if(PRINT_ALI)Write_ss_ali(Prot_ali,prots,N_pdb,N_ali,name_in,"SeqAli");
    Prot_ali_ss=malloc(NS*sizeof(int *));
    for(i=0; i<NS; i++)Prot_ali_ss[i]=malloc(N_ali*sizeof(int));
    prot_ij=malloc(NS*sizeof(struct protein *));
    if(SS_MULT){
      for(i=0; i<NS; i++)prot_ij[i]=prots+i;
      Align_ss(Prot_ali_ss, Prot_ali, N_ali, prot_ij, N_pdb, SHIFT_MAX, 1);
      if(PRINT_ALI){
	Write_ali_ss(Prot_ali_ss, prots, N_pdb, N_ali, name_in);
	Write_ss_ali(Prot_ali_ss, prots, N_pdb, N_ali, name_in, "SecStrAli");
      }
    }else{
      Prot_ali_ij=malloc(NS*sizeof(int *));
    }
  }
  
  // Prepare output
 end_ss:
  printf("ALI_SS= %d SS_MULT= %d NS= %d SHIFT_MAX=%d\n",
	 ALI_SS,SS_MULT,NS,SHIFT_MAX);
  char ss_def[100];
  if(ALI_SS){
    sprintf(ss_def, "MSA corrected for secondary structure.");
    if(SS_MULT){strcat(ss_def, " Multiple alignment correction.");}
    else{strcat(ss_def, " Pairwise correction.");}
    char tmp[30]; sprintf(tmp, " SHIFT_MAX= %d\n",SHIFT_MAX);
    strcat(ss_def, tmp);
  }
  char norm_def[200];
  sprintf(norm_def,"# Normalization of nali: ");
  if(NORM==0){strcat(norm_def, " Minimum length\n");}
  else if(NORM==1){strcat(norm_def, " Maximum length\n");}
  else{strcat(norm_def, " Geometric mean\n");}
  if(NORMA){
    strcat(norm_def,"# Normalization of seqid, TM-score and CO: ");
    strcat(norm_def, " aligned residues\n");
  }

  char name_sim[100]; FILE *file_sim=NULL;
  if(PRINT_SIM && 0){
    Change_ext(name_sim, name_in, EXT_SIM);
    file_sim=fopen(name_sim, "w");
    fprintf(file_sim, "### Prot1 Prot2 Seq_Id Cont_Overlap TM_Score align ");
    if(ALI_SS)fprintf(file_sim, " Seq_Id_SS Cont_Ov_SS TM_SS align_SS");
    fprintf(file_sim, " Seq_Id_PC Cont_Ov_PC TM_PC align_PC\n");
    fprintf(file_sim, "### 0 0  1 1 1 1");
    if(ALI_SS)fprintf(file_sim, "  0 0 0 0");
    fprintf(file_sim, "  0 0 0 0\n");
  }

  // Conditional probabilities
  int PRINT_ID=1;
  char name_id[200]; FILE *file_id=NULL;
  int *id_aa[ITYPE], *id_sup[ITYPE], *shift[ITYPE],
    *neigh_ali[ITYPE], *neigh_noali[ITYPE], *neigh_noali_aaid[ITYPE],
    *ali_cont_ct[ITYPE], *id_cont_ct[ITYPE], 
    *ali_cont_sup[ITYPE], *id_cont_sup[ITYPE],
    *ali_cont_aaid[ITYPE], *id_cont_aaid[ITYPE];
  if(PRINT_ID){
    Change_ext(name_id, name_in, ".id");
    file_id=fopen(name_id, "w");
    for(i=0; i<ITYPE; i++){
      id_aa[i]=malloc(N_ali*sizeof(int));
      id_sup[i]=malloc(N_ali*sizeof(int));
      neigh_ali[i]=malloc(N_ali*sizeof(int));
      neigh_noali[i]=malloc(N_ali*sizeof(int));
      neigh_noali_aaid[i]=malloc(N_ali*sizeof(int));
      shift[i]=malloc(N_ali*sizeof(int));
      id_cont_ct[i]=malloc(2*sizeof(int));
      ali_cont_ct[i]=malloc(2*sizeof(int));
      id_cont_sup[i]=malloc(3*sizeof(int));
      ali_cont_sup[i]=malloc(3*sizeof(int));
      id_cont_aaid[i]=malloc(3*sizeof(int));
      ali_cont_aaid[i]=malloc(3*sizeof(int));
    }
    fprintf(file_id,"# Conservation properties of aligned residues\n");
  }

  // Allocate pairwise computations only for i>j

  ///////////////////////////////////
  printf("Computing pairwise similarities\n");
  int ali_ij[Max_ali];  // Input MSA

  float **Seq_diff_SqA=Allocate_mat2_f(N_pdb, N_pdb);

  int *ali_all[ITYPE];
  float **na_all[ITYPE];
  float **SI_all[ITYPE];
  float **TM_all[ITYPE];
  float **CO_all[ITYPE];
  float **PC_all[ITYPE];

  int *ali_tmp=malloc(Max_ali*sizeof(int));
  for(int it=0; it<ITYPE; it++){
    ali_all[it]=NULL; na_all[it]=NULL; SI_all[it]=NULL;
    TM_all[it]=NULL;  CO_all[it]=NULL; PC_all[it]=NULL;
    if(it==1){if(ALI_SS==0)continue;}
    else if(it==2){if(ALI_TM==0)continue;}
    else if(it==3){if(ALI_CO==0)continue;}
    else if(it==4){if(ali_str==NULL)continue;}
    else if(it==5){if(ALI_PC==0)continue;}
    ali_all[it]=malloc(Max_ali*sizeof(int));

    na_all[it]=Allocate_mat2_f(N_pdb, N_pdb);
    SI_all[it]=Allocate_mat2_f(N_pdb, N_pdb);
    TM_all[it]=Allocate_mat2_f(N_pdb, N_pdb);
    CO_all[it]=Allocate_mat2_f(N_pdb, N_pdb);
    PC_all[it]=Allocate_mat2_f(N_pdb, N_pdb);
  }

  int al2i=-1, al2j=-1;
  float CO[ITYPE], TM[ITYPE], SI[ITYPE], PC[ITYPE], nali[ITYPE];
  for(i=0; i<ITYPE; i++){
    opt_PC[i]=0; Diff_opt[i]=0;
    CO[i]=0; TM[i]=0; SI[i]=0; nali[i]=0; PC[i]=0;
    CO1[i]=0; TM1[i]=0; SI1[i]=0; nali1[i]=0; PC1[i]=0;
    for(j=0; j<ITYPE; j++){
      SI_diff1[i][j]=0; SI_diff2[i][j]=0;
      CO_diff1[i][j]=0; CO_diff2[i][j]=0;
      TM_diff1[i][j]=0; TM_diff2[i][j]=0;
      PC_diff1[i][j]=0; PC_diff2[i][j]=0;
      ali_diff1[i][j]=0;ali_diff2[i][j]=0;
    }
    ali_name[i]=malloc(40*sizeof(char));
  }
  strcpy(ali_name[0], "Input");
  strcpy(ali_name[1], "SS");
  strcpy(ali_name[2], "TM");
  strcpy(ali_name[3], "CO");
  strcpy(ali_name[4], "Struc");
  strcpy(ali_name[5], "PC");

  PC0_1=PC_load[0]*0.5+PC_load[1]*S0+PC_load[2]*TM0;
  PC_norm=PC_load[0]+PC_load[1]+PC_load[2]+PC_load[3];
  // PC value for random pairs

  // Store PCA alignments
  //int **PC_opt_pair=Allocate_mat2_i(N_pdb, N_pdb);
  int **Ali_pair[N_pdb]; // ali[i][j][site_i]=site_j i<j
  for(i=0; i<N_pdb; i++){
    Ali_pair[i]=Allocate_mat2_i(N_pdb, L_seq[i]);
  }

  for(i=0; i<N_pdb; i++){
    int al1i=i_seq[i]; if(ali_str)al2i=ali_str[al1i];
    struct protein *proti=prots+i, *protj=prots;
    printf("Str %s %d of %d L= %d\n", proti->name, i, N_pdb, L_seq[i]);

    if(ALI_SS && SS_MULT==0){
      prot_ij[0]=prots+i; Prot_ali_ij[0]=Prot_ali[i];
    }
    float **d2=Allocate_mat2_f(proti->len, N_ali);
    int   **nc=Allocate_mat2_i(proti->len, N_ali);

    for(j=0; j<i; j++){
 
      pairs++;

      // Normalization
      float norm_ali=proti->len, norm_c=proti->N_cont;
      if(NORM==0){
	if(protj->len<norm_ali)norm_ali=protj->len;
	if(protj->N_cont<norm_c)norm_c=protj->N_cont;
      }else if(NORM==1){
	if(protj->len>norm_ali)norm_ali=protj->len;
	if(protj->N_cont>norm_c)norm_c=protj->N_cont;
      }else{
	norm_ali=sqrt(proti->len*protj->len);
	norm_c=sqrt(proti->N_cont*protj->N_cont);
      }

      int id, al1j=i_seq[j]; if(ali_str)al2j=ali_str[al1j];
      float cali=1;

      // Input alignment
      int it=0; 
      nali[it]=Pair_ali(ali_ij, N_ali, Prot_ali[i], Prot_ali[j])/norm_ali;
      for(int s=0; s<proti->len; s++)ali_all[it][s]=ali_ij[s];
      // Similarities and divergences based on sequence alignment
      // Seq. identity
      Seq_diff_SqA[i][j]=
	Seq_differences(&id, Prot_in[al1i].seq, Prot_in[al1j].seq, N_ali);
      SI[it]=Seqid(ali_ij, id_aa[it],
		   proti->aseq, proti->len,
		   protj->aseq, protj->len);

      // TM score and statistics of optimal rotation
      float d02=0;
      TM[it]=TM_score(d2, &d02, ali_ij, norm_ali,
		      proti->xca, proti->len, protj->xca, protj->len, 0);
      Examine_neighbors(d2, ali_ij, d02, shift[it], id_sup[it],
			neigh_ali[it],neigh_noali[it],neigh_noali_aaid[it],
			proti->aseq, proti->len, protj->aseq, protj->len);

      // Contact overlap
      CO[it]=
	Contact_overlap(ali_ij, &cali, id_sup[it], c_ave,
			ali_cont_ct[it],  id_cont_ct[it],
			ali_cont_sup[it], id_cont_sup[it],
			ali_cont_aaid[it],id_cont_aaid[it],
			proti->Cont_map, proti->aseq, proti->ncont, proti->len,
			protj->Cont_map, protj->aseq, protj->ncont, protj->len);

      na_all[it][i][j]=nali[it];
      SI_all[it][i][j]=SI[it];
      TM_all[it][i][j]=TM[it];
      CO_all[it][i][j]=CO[it];
      //if(NORMA){TM_all[it][i][j]/=nali[it]; CO_all[it][i][j]/=cali;}

      //Cont_Div_SqA[i][j]=
      //Compute_Dcont(&qinf, CO[it], proti->len, protj->len, &homo, NORM);

      Write_identity(file_id, it, proti, protj, ali_ij,
		     id_aa[it], id_sup[it], shift[it],
		     neigh_ali[it], neigh_noali[it], neigh_noali_aaid[it],
		     ali_cont_ct[it], id_cont_ct[it],
		     ali_cont_sup[it], id_cont_sup[it],
		     ali_cont_aaid[it], id_cont_aaid[it]);

      // Similarities and divergences based on sec.str. corrected alignments 
      if(ALI_SS){
	it=1;
	// Alignment length
	if(SS_MULT==0){
	  prot_ij[1]=prots+j; Prot_ali_ij[1]=Prot_ali[j];
	  Align_ss(Prot_ali_ss,Prot_ali_ij,N_ali,prot_ij,NS,SHIFT_MAX,0);
	  nali[it]=Pair_ali(ali_ij, N_ali, Prot_ali_ss[0], Prot_ali_ss[1]);
	}else{
	  nali[it]=Pair_ali(ali_ij, N_ali, Prot_ali_ss[i], Prot_ali_ss[j]);
	}
	nali[it]/=norm_ali;
	for(int s=0; s<proti->len; s++)ali_all[it][s]=ali_ij[s];
	
	// Seq. identity
	int SI_high=-1, TM_high=-1, CO_high=-1;
	SI[it]=Seqid(ali_ij, id_aa[it],
		     proti->aseq, proti->len,
		     protj->aseq, protj->len);
	if(SI[it]<SI[0]){SI[it]=SI[0]; SI_high=0;}
	else if(SI[it]>SI[0]){SI_high=1; }
	if(SI_high>=0)sum_SS_SI[SI_high]++;

	// TM score
	int accept_SS=0, kt;
	TM[it]=TM_score(d2, &d02, ali_ij, norm_ali,
		      proti->xca, proti->len, protj->xca, protj->len, 0);
	if(TM[it]>=TM[0]){
	  accept_SS=1; kt=it; if(TM[it]>TM[0]){TM_high=1;}
	  Examine_neighbors(d2, ali_ij, d02, shift[it], id_sup[it],
			    neigh_ali[it],neigh_noali[it],neigh_noali_aaid[it],
			    proti->aseq, proti->len, protj->aseq, protj->len);
	}else{
	  TM[it]=TM[0]; kt=0; TM_high=0;
	}
	if(TM_high>=0)sum_SS_TM[TM_high]++;

	if(0 && (TM[it]-TM[0])<-0.03){
	  int l=proti->len; if(protj->len<l)l=protj->len;
	  printf("WARNING, Decrease of TM from %.1f to %.1f\n",
		 TM[0]*l,TM[it]*l);
	  Print_ali_ss(Prot_ali[i], proti->ss,
		       Prot_ali[j], protj->ss, N_ali);
	  printf("=========================================\n");
	  int ii=0,jj=1; if(SS_MULT){ii=i; jj=j;}
	  Print_ali_ss(Prot_ali_ss[ii], proti->ss,
		       Prot_ali_ss[jj], protj->ss, N_ali);
	}

	// Contact overlap
	CO[it]=
	  Contact_overlap(ali_ij, &cali, id_sup[it], c_ave,
			  ali_cont_ct[it],  id_cont_ct[it],
			  ali_cont_sup[it], id_cont_sup[it],
			  ali_cont_aaid[it],id_cont_aaid[it],
			  proti->Cont_map,proti->aseq,
			  proti->ncont,proti->len,
			  protj->Cont_map,protj->aseq,
			  protj->ncont,protj->len);
	if(CO[it]<CO[0]){CO_high=0;}
	else if(CO[it]>CO[0]){CO_high=1;}
	if(accept_SS==0)CO[it]=CO[0];
	if(CO_high>=0)sum_SS_CO[CO_high]++;
	if(SI_high>=0 && TM_high>=0)sum_SS_SI_TM[SI_high][TM_high]++;
	if(SI_high>=0 && CO_high>=0)sum_SS_SI_CO[SI_high][CO_high]++;
	if(TM_high>=0 && CO_high>=0)sum_SS_TM_CO[TM_high][CO_high]++;

	na_all[it][i][j]=nali[it];
	SI_all[it][i][j]=SI[it];
	TM_all[it][i][j]=TM[it];
	CO_all[it][i][j]=CO[it];
	//if(NORMA){TM_all[it][i][j]/=nali[it]; CO_all[it][i][j]/=cali;}

	Write_identity(file_id, it, proti, protj, ali_ij,
		       id_aa[kt], id_sup[kt], shift[kt],
		       neigh_ali[kt], neigh_noali[kt], neigh_noali_aaid[kt],
		       ali_cont_ct[kt], id_cont_ct[kt],
		       ali_cont_sup[kt], id_cont_sup[kt],
		       ali_cont_aaid[kt], id_cont_aaid[kt]);
      }


      // TM_score based alignment
      float d2min1[proti->len];
      if(ALI_TM){
	it=2; int *ali_TM=ali_all[it];
	for(int k=0; k<proti->len; k++)ali_tmp[k]=ali_ij[k];
	for(int iter=0; iter<ITMAX; iter++){
	  Align_TM(ali_TM, d2min1, d2, d02, ali_tmp, proti->len, protj->len);
	  TM[it]=TM_score(d2, &d02, ali_TM, norm_ali,
	  		  proti->xca, proti->len, protj->xca, protj->len, 0);
	  /*TM[it]=TM_fast(d2, d02, ali_TM, norm_ali, d2min1,
	    proti->xca, proti->len, protj->xca, protj->len);*/
	  for(int k=0; k<proti->len; k++)ali_tmp[k]=ali_TM[k];
	}
	// Statistics of TM score
	Examine_neighbors(d2, ali_TM, d02, shift[it], id_sup[it],
			  neigh_ali[it],neigh_noali[it],neigh_noali_aaid[it],
			  proti->aseq, proti->len, protj->aseq, protj->len);

	nali[it]=Count_ali(ali_TM, proti->len)/norm_ali;
	SI[it]=Seqid(ali_TM, id_aa[it],
		     proti->aseq, proti->len,
		     protj->aseq, protj->len);	

	// Contact overlap
	CO[it]=
	  Contact_overlap(ali_TM, &cali, id_sup[it], c_ave,
			  ali_cont_ct[it],  id_cont_ct[it],
			  ali_cont_sup[it], id_cont_sup[it],
			  ali_cont_aaid[it],id_cont_aaid[it],
			  proti->Cont_map,proti->aseq,proti->ncont,proti->len,
			  protj->Cont_map,protj->aseq,protj->ncont,protj->len);

	Write_identity(file_id, it, proti, protj, ali_TM,
		       id_aa[it], id_sup[it], shift[it],
		       neigh_ali[it], neigh_noali[it], neigh_noali_aaid[it],
		       ali_cont_ct[it], id_cont_ct[it],
		       ali_cont_sup[it], id_cont_sup[it],
		       ali_cont_aaid[it], id_cont_aaid[it]);
      }

     // Contact Overlap based alignment
      if(ALI_CO){
	it=3; int *ali_CO=ali_all[it], *ali_ini;
	if(ALI_TM){ali_ini=ali_all[2];}else{ali_ini=ali_ij;}
	for(int k=0; k<proti->len; k++)ali_tmp[k]=ali_ini[k];
	for(int iter=0; iter<ITMAX; iter++){
	  Align_CO(ali_CO, ali_tmp, nc,
		   proti->Cont_map, proti->len, protj->Cont_map, protj->len);
	  CO[it]=
	    Contact_overlap(ali_CO, &cali, id_sup[it], c_ave,
			    ali_cont_ct[it],  id_cont_ct[it],
			    ali_cont_sup[it], id_cont_sup[it],
			    ali_cont_aaid[it],id_cont_aaid[it],
			    proti->Cont_map,proti->aseq,
			    proti->ncont,proti->len,
			    protj->Cont_map,protj->aseq,
			    protj->ncont,protj->len);
	  for(int k=0; k<proti->len; k++)ali_tmp[k]=ali_CO[k];
	}

	nali[it]=Count_ali(ali_CO, proti->len)/norm_ali;
	SI[it]=Seqid(ali_CO, id_aa[it],
		     proti->aseq, proti->len,
		     protj->aseq, protj->len);	
	// TM score
	TM[it]=TM_score(d2, &d02, ali_CO, norm_ali,
			proti->xca, proti->len, protj->xca, protj->len, 0);
	Examine_neighbors(d2, ali_CO, d02, shift[it], id_sup[it],
			  neigh_ali[it],neigh_noali[it],neigh_noali_aaid[it],
			  proti->aseq, proti->len, protj->aseq, protj->len);

	Write_identity(file_id, it, proti, protj, ali_CO,
		       id_aa[it], id_sup[it], shift[it],
		       neigh_ali[it], neigh_noali[it], neigh_noali_aaid[it],
		       ali_cont_ct[it], id_cont_ct[it],
		       ali_cont_sup[it], id_cont_sup[it],
		       ali_cont_aaid[it], id_cont_aaid[it]);
      }


     // PC based alignment
      if(ALI_PC){
	it=5; int *ali_PC=ali_all[it], *ali_ini;

	if(ALI_CO){ali_ini=ali_all[3];}
	else if(ALI_TM){ali_ini=ali_all[2];}
	else{ali_ini=ali_ij;}
	for(int k=0; k<proti->len; k++)ali_tmp[k]=ali_ini[k];

	for(int iter=0; iter<ITMAX; iter++){
	  Align_PC(ali_PC, ali_tmp, PC_load, d02, d2, nc,
		   proti->aseq, proti->len, protj->aseq, protj->len);
	  TM[it] = TM_score(d2, &d02, ali_PC, norm_ali,
			    proti->xca, proti->len, protj->xca, protj->len, 0);
	  CO[it]=
	    Contact_overlap(ali_PC, &cali, id_sup[it], c_ave,
			    ali_cont_ct[it],  id_cont_ct[it],
			    ali_cont_sup[it], id_cont_sup[it],
			    ali_cont_aaid[it],id_cont_aaid[it],
			    proti->Cont_map,proti->aseq,
			    proti->ncont,proti->len,
			    protj->Cont_map,protj->aseq,
			    protj->ncont,protj->len);
	  for(int k=0; k<proti->len; k++)ali_tmp[k]=ali_PC[k];
	}

	nali[it]=Count_ali(ali_PC, proti->len)/norm_ali;
	SI[it]=Seqid(ali_PC, id_aa[it],
		     proti->aseq, proti->len,
		     protj->aseq, protj->len);	
	Examine_neighbors(d2, ali_PC, d02, shift[it], id_sup[it],
			  neigh_ali[it],neigh_noali[it],neigh_noali_aaid[it],
			  proti->aseq, proti->len, protj->aseq, protj->len);
      }


      // Similarities and divergences based on structure alignments
      if(ali_str && (al2i>=0) && (al2j>=0)){
	it=4;
	nali[it]=Pair_ali(ali_ij, N_ali_str, Prot_ali_str[i], Prot_ali_str[j])
	  /norm_ali;
	SI[it]=Seqid(ali_ij, id_aa[it],
		     proti->aseq, proti->len,
		     protj->aseq, protj->len);

	// TM score and statistics of optimal rotation
	TM[it]=TM_score(d2, &d02, ali_ij, norm_ali,
			proti->xca, proti->len, protj->xca, protj->len, 0);
	Examine_neighbors(d2, ali_ij, d02, shift[it], id_sup[it],
			  neigh_ali[it],neigh_noali[it],neigh_noali_aaid[it],
			  proti->aseq, proti->len, protj->aseq, protj->len);

	// Conctact overlap
	CO[it]=
	  Contact_overlap(ali_ij, &cali, id_sup[it], c_ave,
			  ali_cont_ct[it],  id_cont_ct[it],
			  ali_cont_sup[it], id_cont_sup[it],
			  ali_cont_aaid[it],id_cont_aaid[it],
			  proti->Cont_map,proti->aseq,proti->ncont,proti->len,
			  protj->Cont_map,protj->aseq,protj->ncont,protj->len);

	na_all[it][i][j]=nali[it];
	SI_all[it][i][j]=SI[it];
	TM_all[it][i][j]=TM[it];
	CO_all[it][i][j]=CO[it];
	//if(NORMA){TM_all[it][i][j]/=nali[it]; CO_all[it][i][j]/=cali;}

	Write_identity(file_id, it, proti, protj, ali_str,
		       id_aa[it], id_sup[it], shift[it],
		       neigh_ali[it], neigh_noali[it], neigh_noali_aaid[it],
		       ali_cont_ct[it], id_cont_ct[it],
		       ali_cont_sup[it], id_cont_sup[it],
		       ali_cont_aaid[it], id_cont_aaid[it]);
      }

      // Compute PC 
      int kt=5;
      for(int it=0; it<ITYPE; it++){
	if(nali[it]==0)continue;
	PC[it]=(PC_load[0]*nali[it]+PC_load[1]*SI[it]+
		PC_load[2]*TM[it]+PC_load[3]*CO[it])/PC_norm;
	PC_all[it][i][j]=PC[it];
      }

      // Optimal alignments
      for(it=5; it>=2; it--){
	int ko=it, kk=0; // Optimal alignment
	if(it==4){ // SS_ali
	  continue;
	}else if(it==5){ // PC_ali
	  for(int k=1; k<(ITYPE-1); k++)if(PC[k]>PC[kk])kk=k;
	  if(PC[5]>=PC[kk]){ko=5;}else{ko=kk;}
	}else if(it==3){ // CO_ali
	  for(int k=0; k<ITYPE; k++)if(CO[k]>CO[ko])ko=k;
	}else if(it==2){ // TM_ali
	  for(int k=0; k<ITYPE; k++)if(TM[k]>TM[ko])ko=k;
	}
	if(ko!=it){
	  PC[it]=PC[ko]; TM[it]=TM[ko]; CO[it]=CO[ko]; 
	  SI[it]=SI[ko]; nali[it]=nali[ko]; Diff_opt[it]++;
	}
	if(it==5){ // PC_ali
	  opt_PC[kk]++; kt=ko;
	  Write_ali(Ali_pair, i, j, L_seq, ali_all[ko]);
	  //PC_opt_pair[i][j]=ko;
	} 

	na_all[it][i][j]=nali[it];
	SI_all[it][i][j]=SI[it];
	PC_all[it][i][j]=PC[it];
	TM_all[it][i][j]=TM[it];
	CO_all[it][i][j]=CO[it];
	//if(NORMA){TM_all[it][i][j]/=nali[it]; CO_all[it][i][j]/=cali;}
      }

      // Propensities from PC alignment
      it=5;
      Write_identity(file_id, it, proti, protj, ali_all[kt],
		     id_aa[kt], id_sup[kt], shift[kt],
		     neigh_ali[kt], neigh_noali[kt], neigh_noali_aaid[kt],
		     ali_cont_ct[kt], id_cont_ct[kt],
		     ali_cont_sup[kt], id_cont_sup[kt],
		     ali_cont_aaid[kt], id_cont_aaid[kt]);

      Sec_str_AA_propensity(proti, protj, ali_all[kt], nali[kt], PC[kt]);

      // Compare alignments
      for(int it=0; it<ITYPE; it++){
	if(nali[it]==0)continue;
	nali1[it]+=nali[it]; SI1[it]+=SI[it];
	TM1[it]+=TM[it]; CO1[it]+=CO[it]; PC1[it]+=PC[it];
	for(int jt=it+1; jt<ITYPE; jt++){
	  if(nali[jt]==0){continue;} float d;
	  d=SI[jt]-SI[it]; SI_diff1[it][jt]+=d; SI_diff2[it][jt]+=d*d;
	  d=CO[jt]-CO[it]; CO_diff1[it][jt]+=d; CO_diff2[it][jt]+=d*d;
	  d=TM[jt]-TM[it]; TM_diff1[it][jt]+=d; TM_diff2[it][jt]+=d*d;
	  d=PC[jt]-PC[it]; PC_diff1[it][jt]+=d; PC_diff2[it][jt]+=d*d;
	  d=nali[jt]-nali[it]; ali_diff1[it][jt]+=d; ali_diff2[it][jt]+=d*d;
	}
      }

      // Print
      if(file_sim){
	fprintf(file_sim, "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f",
		proti->name, protj->name, SI[0], CO[0], TM[0], nali[0]);
	if(ALI_SS)fprintf(file_sim,"\t%.3f\t%.3f\t%.3f\t%.3f\n",
			  SI[1],CO[1],TM[1], nali[1]);
	if(ALI_TM)fprintf(file_sim,"\t%.3f\t%.3f\t%.3f\t%.3f\n",
			  SI[2],CO[2],TM[2], nali[2]);
	if(ALI_CO)fprintf(file_sim,"\t%.3f\t%.3f\t%.3f\t%.3f\n",
			  SI[3],CO[3],TM[3], nali[3]);
	if(ALI_PC)fprintf(file_sim,"\t%.3f\t%.3f\t%.3f\t%.3f\n",
			  SI[5], CO[5], TM[5], nali[5]);
      }
      protj++;
    }
    Empty_matrix_i(nc, proti->len);
    Empty_matrix_f(d2, proti->len);
  } // end pairs

  // Symmetrize
  int it;
  for(it=0; it<ITYPE; it++){
    if(ali_all[it]==NULL)continue;
    for(i=0; i<N_pdb; i++){
      for(j=0; j<i; j++){
	na_all[it][j][i]=na_all[it][i][j];
	SI_all[it][j][i]=SI_all[it][i][j];
	CO_all[it][j][i]=CO_all[it][i][j];
	TM_all[it][j][i]=TM_all[it][i][j];
	PC_all[it][j][i]=PC_all[it][i][j];
      }
    }
  }
  printf("End pairwise computations\n");

  if(file_sim){
    printf("Similarities written in file %s for all PDB pairs\n", name_sim);
    fclose(file_sim);
  }
  if(file_id){
    printf("Probabilities of identity written in %s for all pairs\n",name_id);
    Summary_identical(file_id);
    //fclose(file_id);
  }

  char name_prop[100];
  Change_ext(name_prop, name_in, ".prop");
  Print_propensities(name_prop);

  if(PRINT_ID){
    for(i=0; i<ITYPE; i++){
      free(id_aa[i]);
      free(id_sup[i]);
      free(neigh_ali[i]);
      free(neigh_noali[i]);
      free(neigh_noali_aaid[i]);
      free(shift[i]);
      free(id_cont_ct[i]);
      free(ali_cont_ct[i]);
      free(id_cont_sup[i]);
      free(ali_cont_sup[i]);
      free(id_cont_aaid[i]);
      free(ali_cont_aaid[i]);
    }
  }


  /*************************************************************************
       Minimum str. divergence over all conformations of the same sequence
  **************************************************************************/
  // Group identical sequences
  // Structural divergence: minimum among all conformations 
  printf("Grouping identical sequences by single linkage\n");
  int **conformation, *N_conf, *rep_str, *seq_clus;
  int N_seq=
    Single_linkage(&conformation, &N_conf, &rep_str, &seq_clus,
		   Seq_diff_SqA, N_pdb, diff_thr);
  printf("%d conformations grouped into %d sequences\n", N_pdb, N_seq);
  printf("Number of conformations per sequence: ");
  //for(i=0; i<N_seq; i++){printf("%d (rep=%d) ",N_conf[i], rep_str[i]);}
  printf("Average: %.1f\n", (float)N_pdb/N_seq);

  // Print MSA
  if(ALI_PC && (PRINT_PAIR||PRINT_CLIQUE)){
    char name2[85]; sprintf(name2, "%s_PC", name_in);
    int len_seq[N_seq], i;
    char **name_seq, **Seq=
      Assign_Seq(&name_seq, len_seq, Prot_in,N_seq,rep_str,i_seq,N_ali);
    int ***Ali_pair_seq=Select_alis(Ali_pair, N_seq, len_seq, rep_str);
    for(i=0; i<N_pdb; i++)Empty_matrix_i(Ali_pair[i], N_pdb);
    if(PRINT_PAIR)
      Print_pairwise(Ali_pair_seq, len_seq, N_seq, Seq, name_seq, name2);
    if(PRINT_CLIQUE)
      Clique_MSA(Ali_pair_seq, len_seq, N_seq, Seq, name_seq, name2, N_ali);
    for(i=0; i<N_seq; i++)Empty_matrix_i(Ali_pair_seq[i], N_seq);
  }


  /****************************************************************
              Function similarity 
  *****************************************************************/
  // Functional similarity (if file_fun is present)
  float **fun_sim=Read_function(file_fun, Prot_in, i_seq, N_pdb);
  float **fun_sim_Seq=NULL;
  int n_low=0, n_high=0;
  if(fun_sim){
    fun_sim_Seq=malloc(N_seq*sizeof(float *));
    for(i=0; i<N_seq; i++){
      fun_sim_Seq[i]=malloc(N_seq*sizeof(float));
      for(j=0; j<N_seq; j++)fun_sim_Seq[i][j]=-1;
    }
    for(i=0; i<N_pdb; i++){
      int i1=seq_clus[i];
      for(j=0; j<i; j++){
	if(fun_sim[i][j]<0)continue;
	int j1=seq_clus[j];
	float s=fun_sim[i][j];
	fun_sim_Seq[i1][j1]=s;
	fun_sim_Seq[j1][i1]=s;
	if(s<=fun_low)n_low++;
	if(s>=fun_high)n_high++;
      }
    }
    Empty_matrix_f(fun_sim, N_pdb);
    printf("%d pairs with function similarity <= %.3f\n",
	   n_low, fun_low);
    printf("%d pairs with function similarity >= %.3f\n",
	   n_high, fun_high);
  }


  /***************************************************************************
    Minimum structure divergence over all conformations of the same protein
  ****************************************************************************/
  printf("Computing maximum structure similarity "
	 "over all conformations of the same protein\n");
  char name_div[100];
  Change_ext(name_div, name_in, ".prot.div");
  printf("Printing grouped structure divergence in %s\n", name_div);
  FILE *file_div=fopen(name_div, "w");
  Change_ext(name_sim, name_in, ".prot.sim");
  printf("Printing grouped structure similarity in %s\n", name_sim);
  file_sim=fopen(name_sim, "w");

  char head[1000];
  sprintf(head, "# Structure similarity maximized over all different "
	  "conformations of the same protein and different alignments\n");
  if(ALI_SS){strcat(head, "# "); strcat(head, ss_def);}
  strcat(head, norm_def); 
  fprintf(file_div, "%s", head);
  fprintf(file_sim, "%s", head);

  fprintf(file_div, "###Prot1 Prot2 ");
  fprintf(file_div,
	  " 3=TN_Div_seq 4=Cont_Div_seq 5=TM_Div_seq 6=PC_Div_seq");
  int k=7;
  if(ALI_SS){
    fprintf(file_div,
	    " %d=TN_Div_SS %d=Cont_Div_SS %d=TM_Div_SS %d=PC_Div_SS",
	    k, k+1, k+2, k+3); k+=4;
  }
  if(ALI_TM){
    fprintf(file_div,
	    " %d=TN_Div_TM %d=Cont_Div_TM %d=TM_Div_TM %d=PC_Div_TM",
	    k, k+1, k+2, k+3); k+=4;
  }
  if(ALI_CO){
    fprintf(file_div,
	    " %d=TN_Div_CO %d=Cont_Div_CO %d=TM_Div_CO %d=PC_Div_CO",
	    k, k+1, k+2, k+3); k+=4;
  }
  if(ALI_PC){
    fprintf(file_div,
	    " %d=TN_Div_PC %d=Cont_Div_PC %d=TM_Div_PC %d=PC_Div_PC",
	    k, k+1, k+2, k+3); k+=4;
  }
  fprintf(file_div,"\n");
  fprintf(file_div, "### 0 0");
  fprintf(file_div, "  1 1 1 1");
  if(ALI_SS)fprintf(file_div, "  0 0 0 0");
  if(ALI_TM)fprintf(file_div, "  0 0 0 0");
  if(ALI_CO)fprintf(file_div, "  0 0 0 0");
  if(ALI_PC)fprintf(file_div, "  0 0 0 0");
  fprintf(file_div,"\n");

  fprintf(file_sim, "### Prot1 Prot2 ");
  fprintf(file_sim, " nali_seq SI_seq CO_seq TM_seq PC_seq");
  if(ALI_SS)fprintf(file_sim,"  nali_SS SI_SS CO_SS TM_SS PC_SS");
  if(ALI_TM)fprintf(file_sim,"  nali_TM SI_TM CO_TM TM_TM PC_TM");
  if(ALI_CO)fprintf(file_sim,"  nali_CO SI_CO CO_CO TM_CO PC_CO");
  if(ALI_PC)fprintf(file_sim,"  nali_PC SI_PC CO_PC TM_PC PC_PC");
  if(fun_sim_Seq)fprintf(file_sim, "  Function_similarity");
  fprintf(file_sim,"\n");
  fprintf(file_sim, "### 0 0");
  fprintf(file_sim, "  1 1 1 1 0");
  if(ALI_SS)fprintf(file_sim, "  0 0 0 0 0");
  if(ALI_TM)fprintf(file_sim, "  0 0 0 0 0");
  if(ALI_CO)fprintf(file_sim, "  0 0 0 0 0");
  if(ALI_PC)fprintf(file_sim, "  0 0 0 0 0");
  if(fun_sim_Seq)fprintf(file_sim, " 0");
  fprintf(file_sim,"\n");

  //if(ali_str)Change_conformations(conformation, N_seq, N_conf, ali_str);
  float **p_Div_Seq=NULL, **Poiss_Div_Seq=NULL; int POISS=0;
  if(POISS){
    p_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
    Poiss_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
  }
  float **Seq_Id_Seq=Allocate_mat2_f(N_seq, N_seq);
  float **TN_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
  float **TM_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
  float **PC_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
  float **Cont_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
  int o=0; // Scores based on PC_ALI, o=5
  if(ALI_PC){o=5;}else if(ali_str){o=4;}else if(ALI_CO){o=3;}
  else if(ALI_SS){o=1;}else if(ALI_TM){o=1;}

  float SI_max[ITYPE], TM_max[ITYPE], CO_max[ITYPE], PC_max[ITYPE],
    na_max[ITYPE];
  for(i=0; i<N_seq; i++){
    struct protein *pi=prots+conformation[i][0];
    for(j=0; j<i; j++){
      struct protein *pj=prots+conformation[j][0];
      // Look for minimum among conformations

      for(int it=0; it<ITYPE; it++){
	SI_max[it]=0; TM_max[it]=0; CO_max[it]=0; PC_max[it]=0; na_max[it]=0;
      }

      for(int ki=0; ki<N_conf[i]; ki++){
	int ci=conformation[i][ki];
	for(int kj=0; kj<N_conf[j]; kj++){
	  int cj=conformation[j][kj];
	  // Best over alignment it
	  for(int it=0; it<ITYPE; it++){
	    if(ali_all[it]==NULL)continue;
	    if(SI_all[it][ci][cj]>SI_max[it]){
	      SI_max[it]=SI_all[it][ci][cj]; na_max[it]=na_all[it][ci][cj];
	    }
	    if(TM_all[it][ci][cj]>TM_max[it])TM_max[it]=TM_all[it][ci][cj];
	    if(CO_all[it][ci][cj]>CO_max[it])CO_max[it]=CO_all[it][ci][cj];
	    if(PC_all[it][ci][cj]>PC_max[it])PC_max[it]=PC_all[it][ci][cj];
	  }
	}
      }

      if(file_sim){
	fprintf(file_sim, "%s\t%s", pi->name, pj->name);
	fprintf(file_sim, "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
		na_max[0], SI_max[0],CO_max[0],TM_max[0], PC_max[0]);
	if(ALI_SS){int a=1;
	  fprintf(file_sim, "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
		  na_max[a], SI_max[a],CO_max[a],TM_max[a], PC_max[a]);
	}
	if(ALI_TM){int a=2;
	  fprintf(file_sim, "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
		  na_max[a], SI_max[a],CO_max[a],TM_max[a], PC_max[a]);
	}
	if(ALI_CO){int a=3;
	  fprintf(file_sim, "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
		  na_max[a], SI_max[a],CO_max[a],TM_max[a], PC_max[a]);
	}
	if(ALI_PC){int a=5;
	  fprintf(file_sim, "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
		  na_max[a], SI_max[a],CO_max[a],TM_max[a], PC_max[a]);
	}
	if(fun_sim_Seq)fprintf(file_sim,"\t%.3f",fun_sim_Seq[i][j]);
	fprintf(file_sim,"\n");
      }

      // Print divergence
      int homo; double qinf=0, PC0=0;
      fprintf(file_div, "%s\t%s", pi->name, pj->name);
      for(int it=0; it<ITYPE; it++){
	if(ali_all[it]==NULL)continue;
	double DS=Divergence(SI_max[it], S0, HUGE),
	  DT=Divergence(TM_max[it], TM0, HUGE),
	  DC=Compute_Dcont(&qinf,CO_max[it],pi->len,pj->len,&homo,NORM);
	if(it==0){
	  PC0=(PC0_1+PC_load[3]*qinf)/PC_norm;
	  //fprintf(file_div, "\t%.3f", PC0);
	}
	double DP=Divergence(PC_max[it], PC0, HUGE);
	fprintf(file_div, "\t%.3f\t%.3f\t%.3f\t%.3f", DS, DC, DT, DP);
	if(it==o){
	  Cont_Div_Seq[i][j]=DC;
	  TN_Div_Seq[i][j]=DS;
	  TM_Div_Seq[i][j]=DT;
	  PC_Div_Seq[i][j]=DP;
	  if(POISS){
	    p_Div_Seq[i][j]=p_Div_Seq[j][i]=1-SI_max[o];
	    Poiss_Div_Seq[i][j]=Poiss_Div_Seq[j][i]=-log(SI_max[o]);
	  }
	}		
      }
      fprintf(file_div,"\n");
    }
  }
  if(file_div){
    printf("Divergences between proteins written in %s\n",name_div);
    fclose(file_div);
  }
  if(file_sim){
    printf("Similarities between proteins written in %s\n",name_sim);
    fclose(file_sim);
  }


  Empty_matrix_f(Seq_diff_SqA, N_pdb);
  for(int it=0; it<ITYPE; it++){
    if(it==o)continue;
    if(na_all[it])Empty_matrix_f(na_all[it], N_pdb);
    if(SI_all[it])Empty_matrix_f(SI_all[it], N_pdb);
    if(TM_all[it])Empty_matrix_f(TM_all[it], N_pdb);
    if(CO_all[it])Empty_matrix_f(CO_all[it], N_pdb);
    if(PC_all[it])Empty_matrix_f(PC_all[it], N_pdb);
  }

  // Structural divergences with PC_ali
  float **TN_Div_PDB=Allocate_mat2_f(N_pdb, N_pdb);
  float **TM_Div_PDB=Allocate_mat2_f(N_pdb, N_pdb);
  float **Cont_Div_PDB=Allocate_mat2_f(N_pdb, N_pdb);
  //float **PC_Div_PDB=Allocate_mat2_f(N_pdb, N_pdb);
  for(i=0; i<N_pdb; i++){
    struct protein *pi=prots+i;
    for(j=0; j<i; j++){
      struct protein *pj=prots+j;
      int homo; double qinf=0;
      Cont_Div_PDB[i][j]=
	Compute_Dcont(&qinf,CO_all[o][i][j],pi->len,pj->len,&homo,NORM);
      TN_Div_PDB[i][j]=Divergence(SI_all[o][i][j], S0, HUGE);
      TM_Div_PDB[i][j]=Divergence(TM_all[o][i][j], TM0, HUGE);
      //double PC0=(PC0_1+PC_load[3]*qinf)/PC_norm;
      //PC_Div_PDB[i][j]=Divergence(PC_all[o][i][j], PC0, HUGE);
    }
  }
  if(na_all[o])Empty_matrix_f(na_all[o], N_pdb);
  if(SI_all[o])Empty_matrix_f(SI_all[o], N_pdb);
  if(TM_all[o])Empty_matrix_f(TM_all[o], N_pdb);
  if(CO_all[o])Empty_matrix_f(CO_all[o], N_pdb);
  if(PC_all[o])Empty_matrix_f(PC_all[o], N_pdb);

  printf("%s", ss_def);
  printf("Mean alignment length:  %.3f", nali1[0]/pairs);
  if(ALI_SS)printf(" sec.str. corrected: %.3f", nali1[1]/pairs);
  if(ALI_TM)printf(" TM_ali: %.3f", nali1[2]/pairs);
  if(ALI_CO)printf(" CO_ali: %.3f", nali1[3]/pairs);
  if(nali1[4])printf(" Str_ali: %.3f\n", nali1[4]/pairs);
  if(ALI_PC)printf(" PC_ali: %.3f", nali1[5]/pairs);
  printf("\n");

  printf("Mean contact overlap: %.4f",CO1[0]/pairs);
  if(ALI_SS)printf(" sec.str.corrected: %.4f",CO1[1]/pairs);
  if(ALI_TM)printf(" TM_ali: %.4f",CO1[2]/pairs);
  if(ALI_CO)printf(" CO_ali: %.4f",CO1[3]/pairs);
  if(ALI_PC)printf(" PC_ali: %.4f",CO1[5]/pairs);
  printf("\n");

  printf("Mean TM score: %.4f",TM1[0]/pairs);
  if(ALI_SS)printf(" sec.str.corrected: %.4f",TM1[1]/pairs);
  if(ALI_TM)printf(" TM_ali: %.4f",TM1[2]/pairs);
  if(ALI_CO)printf(" CO_ali: %.4f",TM1[3]/pairs);
  if(ALI_PC)printf(" PC_ali: %.4f",TM1[5]/pairs);
  printf("\n");

  printf("Mean Principal Component: %.4f",PC1[0]/pairs);
  if(ALI_SS)printf(" sec.str.corrected: %.4f",PC1[1]/pairs);
  if(ALI_TM)printf(" TM_ali: %.4f",PC1[2]/pairs);
  if(ALI_CO)printf(" CO_ali: %.4f",PC1[3]/pairs);
  if(ALI_CO)printf(" PC_ali: %.4f",PC1[5]/pairs);
  printf("\n\n");

  /************* Exit if clock violations are not needed *****************/
  if(PRINT_CV==0){
    printf("Clock violation not required, exiting the program\n");
    return(0);
  }

  // Neighbor Joining and outgroups
  /* struct treenode *nodes=Neighbor_Joining(TajNei_Div, N_seq);
  int **N_out, ***outgroup=Build_outgroups(&N_out, nodes, N_seq); */
  float **Out_Div=TN_Div_Seq;
  if(strcmp(OUTG,"CD")==0){Out_Div=Cont_Div_Seq;}
  else if(strcmp(OUTG,"TM")==0){Out_Div=TM_Div_Seq;}
  printf("Determining outgroups through neighbor joining ");
  printf("with divergence %s\n", OUTG);
  int **N_out, ***outgroup=Outgroups_NJ(&N_out, Out_Div, N_seq);

  /*********************** Clock violations ***********************/
  int Nd=3, dist; // types of distances
  if(POISS){Nd=5;}
  char *name_dist[Nd]; float **Div[Nd], **Div_PDB[Nd];
  for(dist=0; dist<Nd; dist++)name_dist[dist]=malloc(80*sizeof(char));
  name_dist[0]="Tajima-Nei"; Div[0]=TN_Div_Seq;   Div_PDB[0]=TN_Div_PDB;   
  name_dist[1]="Cont_Div";   Div[1]=Cont_Div_Seq; Div_PDB[1]=Cont_Div_PDB;
  name_dist[2]="TM_Div";     Div[2]=TM_Div_Seq;   Div_PDB[2]=TM_Div_PDB;
  if(POISS){
    name_dist[3]="p_Div"; Div[3]=p_Div_Seq;
    name_dist[4]="Poiss_Div"; Div[4]=Poiss_Div_Seq;
  }

  // Prepare output
  char head2[5000], ALI[100];
  sprintf(head2,
	  "# %d sequences aligned with %d PDB files\n"
	  "# %d groups of sequences with < %d intragroup substitutions\n"
	  "# Average number of structures per group: %.1f\n"
	  "# Outgroups assigned with NJ using divergence %s "
	  " constructed with alignment %s\n"
	  "# Multiple sequence alignment: %s\n", 
	  N_prot, N_pdb, N_seq, diff_thr, N_pdb/(float)N_seq, file_ali,
	  OUTG, ali_name[o]);
  if(ali_str){
    strcat(head2, "Multiple structure alignment: ");
    strcat(head2, file_ali_str);
    strcpy(ALI,"StAli");
  }else{
    strcpy(ALI,"SqAli");
  }

  char name_out[100];
  Change_ext(name_out, name_in, EXT_CV);
  FILE *file_out=fopen(name_out, "w");
  fprintf(file_out, "# Output of the program %s\n", argv[0]);
  fprintf(file_out, "%s", head);
  fprintf(file_out, "# Prot1, Prot2, number_of_outgroups_C\n");
  fprintf(file_out, "# For each divergence measure d, plot:\n");
  fprintf(file_out, "# d(A,B)\n");
  fprintf(file_out, "# CV(A,B)=sum_C(d(A,C)-d(B,C))/nC*d(A,B)\n"); //^%.2f,alpha
  fprintf(file_out, "# t=|CV(A,B)|/S.E.M.(CV) (Standard Error of Mean)\n");
  fprintf(file_out, "# n_TI number of violations of Triangle Inequality\n");
  fprintf(file_out, "# number of outgroups with minority sign\n");
  fprintf(file_out, "# Note: used outgroups = num_out-num_TI\n");
  fprintf(file_out, "#Prot1 Prot2 L1-L2 num_out");
  k=5;
  for(dist=0; dist<Nd; dist++){
    fprintf(file_out, " %s:", name_dist[dist]);
    fprintf(file_out, " %d=d %d=CV %d=t_CV %d=n_TI %d=n_sign",
	    k, k+1, k+2, k+3, k+4); k+=5;
  }
  fprintf(file_out, "\n");

  // Compute CV
  float **CV_dist[Nd], **t_dist[Nd];
  int **nout_dist[Nd], **nsign_dist[Nd], **nTIV_dist[Nd];
  if(PRINT_AVE){
    for(dist=0; dist<Nd; dist++){
      CV_dist[dist]=Allocate_mat2_f(N_seq, N_seq);
      t_dist[dist]=Allocate_mat2_f(N_seq, N_seq);
      nout_dist[dist]=Allocate_mat2_i(N_seq, N_seq);
      nTIV_dist[dist]=Allocate_mat2_i(N_seq, N_seq);
      nsign_dist[dist]=Allocate_mat2_i(N_seq, N_seq);
    }
  }

  int kgap_max=40, bingap=5;
  int **N_gaps=NULL, **All_gaps=NULL, **TIV_gaps=NULL;
  if(PRINT_GAP){
    N_gaps=Allocate_mat2_i(N_seq, N_seq);
    All_gaps=Allocate_mat2_i(Nd, kgap_max+1);
    TIV_gaps=Allocate_mat2_i(Nd, kgap_max+1);
    for(i=0; i<N_seq; i++){
      int i1=i_seq[rep_str[i]];
      for(j=0; j<i; j++){
	int j1=i_seq[rep_str[j]];
	N_gaps[i][j]=Count_gaps(Prot_in[i1].seq, Prot_in[j1].seq, N_ali);
	N_gaps[j][i]=N_gaps[i][j];
      }
    }
  }

  int SI_low=0, No_out[Nd]; long Num_out=0;
  for(i=0; i<Nd; i++)No_out[i]=0;

  // Sum over pairs of sequences
  for(i=0; i<N_seq; i++){
    for(j=0; j<i; j++){
      Num_out+=N_out[i][j];
      if(Seq_Id_Seq[i][j]<SI_thr){SI_low++; continue;}
      fprintf(file_out, "%s\t%s\t%d\t%d",
	      prots[rep_str[i]].name, prots[rep_str[j]].name,
	      L_seq[rep_str[i]]-L_seq[rep_str[j]], 
	      N_out[i][j]);

      int *outg=outgroup[i][j], used[N_seq];
      for(dist=0; dist<Nd; dist++){
	float **diver=Div[dist], d=diver[i][j];
	// outgroups, triangle inequality, diff.sign
	int nout=0, nTIV=0, nplus=0, kgap=-1, TIV=0; 
	double CV1=0, CV2=0, t=0, nindep=0; // norm=pow(d,alpha);
	for(int k1=0; k1<N_out[i][j]; k1++){
	  // Eliminate outgroups that are very far away
	  int k=outg[k1]; float CV;
	  if((Seq_Id_Seq[i][k]<SI_thr)||(Seq_Id_Seq[j][k]<SI_thr))continue;
	  if(dist<3){
	    CV=Min_CV(i, j, k, conformation, N_conf, Div_PDB[dist]);
	  }else{
	    CV=diver[i][k]-diver[j][k];
	  }
	  if((CV>d)||(CV<-d)){ // check triangle inequality TI
	    used[k1]=0; TIV=1; nTIV++; 
	  }else{
	    used[k1]=1; if(TIV)TIV=0;
	    CV1+=CV; CV2+=CV*CV; nout++; if(CV>0)nplus++;
	    float s_max=0, *Sk=Seq_Id_Seq[k];
	    for(int k2=0; k2<k1; k2++){
	      if((used[k2])&&(Sk[outg[k2]]>s_max))s_max=Sk[outg[k2]];
	    }
	    nindep+=s_max;
	  }
	  if(N_gaps){
	    kgap=N_gaps[i][j]+N_gaps[i][k]+N_gaps[j][k];
	    kgap/=bingap; if(kgap>kgap_max)kgap=kgap_max;
	    All_gaps[dist][kgap]++;
	    if(TIV)TIV_gaps[dist][kgap]++;
	  }
	} // end sum over outgroups

	if(nout){CV1/=nout;}
	else{No_out[dist]++;}
	nindep=nout-nindep;
	if(nout <= 1){t=1;}
	else{
	  CV2=(CV2-nout*CV1*CV1)/(nout-1);
	  if(CV2<=0){t=1;}
	  else{t=fabs(CV1)/sqrt(CV2/nindep);}
	}
	CV1/=d;

	int nm=nout-nplus; if(nm<nplus)nplus=nm;
	fprintf(file_out, "\t%.3f\t%.3f\t%.1f\t%d\t%d",
		d, CV1, t, nTIV, nplus); //norm
	if(PRINT_AVE){
	  CV_dist[dist][i][j]=CV1;
	  t_dist[dist][i][j]=t;
	  nout_dist[dist][i][j]=nout;
	  nsign_dist[dist][i][j]=nplus;
	  nTIV_dist[dist][i][j]=nTIV;
	}
      } // end dists
      fprintf(file_out, "\n");
    }
  } // end pairs
  fclose(file_out);

  if(PRINT_AVE){
    char tmp[200];
    sprintf(tmp, 
	    "# Total number of outgroups: %ld\n# %d pairs with SI < %.2f\n",
	    Num_out, SI_low, SI_thr);
    strcat(head, tmp);
    if(Nd){
      strcat(head, "# pairs with zero outgroups: ");
      for(i=0; i<Nd; i++){
	sprintf(tmp, " %d (%s)", No_out[i], name_dist[i]);
	strcat(head, tmp);
      }
      strcat(head, "%s\n");
    }

    for(dist=0; dist<Nd; dist++){
      CV_statistics(name_in, head, OUTG, ALI,
		    name_dist[dist], Div[dist], CV_dist[dist],
		    t_dist[dist], nout_dist[dist],
		    nsign_dist[dist], nTIV_dist[dist], Seq_Id_Seq, SI_thr,
		    NULL, 0.00, 1.00, N_seq, dist, "AllFun");
      if(fun_sim_Seq){
	CV_statistics(name_in, head, OUTG, ALI,
		      name_dist[dist], Div[dist], CV_dist[dist],
		      t_dist[dist], nout_dist[dist],
		      nsign_dist[dist], nTIV_dist[dist], Seq_Id_Seq, SI_thr,
		      fun_sim_Seq, fun_high, 1.00, N_seq, dist, "SameFun");
	CV_statistics(name_in, head, OUTG, ALI,
		      name_dist[dist], Div[dist], CV_dist[dist],
		      t_dist[dist], nout_dist[dist],
		      nsign_dist[dist], nTIV_dist[dist], Seq_Id_Seq, SI_thr,
		      fun_sim_Seq, 0.00, fun_low, N_seq, dist, "DiffFun");
      }

      Empty_matrix_f(CV_dist[dist], N_seq);
      Empty_matrix_f(t_dist[dist], N_seq);
      Empty_matrix_i(nout_dist[dist], N_seq);
      Empty_matrix_i(nTIV_dist[dist], N_seq);
      Empty_matrix_i(nsign_dist[dist], N_seq);
    }
  }
  if(fun_sim_Seq)Empty_matrix_f(fun_sim_Seq, N_seq);

  /***************************************************************************/
  /* Relation between gaps and triangle inequality */
  if(PRINT_GAP){
    Change_ext(name_out, name_in, ".gaps");
    file_out=fopen(name_out, "w");
    for(int dist=0; dist<Nd; dist++){
      fprintf(file_out, "# dist=%s\n", name_dist[dist]);
      fprintf(file_out, "#ngap P(TIV) s.e. num\n");
      double norm=0;
      for(i=0; i<=kgap_max; i++)norm+=All_gaps[dist][i];
      for(i=0; i<=kgap_max; i++){
	if(All_gaps[dist][i]==0)continue;
	float p=(float)TIV_gaps[dist][i]/All_gaps[dist][i];
	fprintf(file_out, "%.1f %.3f %.3f %.3f\n", (i+0.5)*bingap,
		p, sqrt(p*(1-p)/All_gaps[dist][i]),
		All_gaps[dist][i]/norm);
      }
    }
    fclose(file_out);
    printf("Writing %s\n", name_out);
  }

  /***************************************************************************/
  return(0);
}

void help(char *pname){
  printf("Program %s\n", pname);
  printf("Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa "
	 "(CSIC-UAM), Madrid, Spain\nEmail: <ubastolla@cbm.csic.es>\n\n");
  printf("Given a multiple sequence alignment (MSA) of proteins with known structures and optionally a multiple structure alignment, it computes and prints sequence and structural similarity measures (sequence identity SI, contact overlap CO, TM-score TM, Principal component similarity PC that integrates all the above similarity measures) and the corresponding divergence measures (Tajima-Nei divergence TN=-log((SI-S0)/(1-S0)) with S0=%.2f, Contact_divergence CD=-log((CO-CO(L))/(1-CO(L)), TM_divergence=-log((TM-TM0)/(1-TM0), PC divergence=-log((PC-PC0)/(1-PC0)) for all aligned pairs.\n", S0);
  printf("For all four similarity measures SI, CO, TM, PC, the program modifies the input alignment in such a way that the similarity measures are maximized with minimal modification of the input alignment. The PC alignment, which is the most accurate one, is printed as pairwise alignments in file <>_PC.pwa if PRINT_PAIR=1 is set in the configuration file, and as multiple alignment based on the maximal cliques of the pairwise alignments in file <>_PC.msa if PRINT_CLIQUE=1 is set in the configuration file (WARNING, this computation scales as (n*L)^3, it is slow for more than 40 proteins 200 residues long).\n");
  printf("Optionally, if PRINT_CV=1 is set, the program computes and prints for all four divergence measures the violation of the molecular clock averaged over all possible outgroups identified with the Neighbor-Joining criterion, and the corresponding significance score\n");
  printf("==========================================================\n");
  printf("RUN: %s <Config file>\n", pname);
  printf("==========================================================\n");
  printf("Config file:\n");
  printf("ALI=<MSA file in FASTA format, indicating names of PDB files>\n");
  printf("Optional parameters:\n");
  printf("STR_ALI=<MSA file in FASTA format, with names of PDB files>\n");
  printf("# (multiple structure alignment, optional).\n");
  printf("FUN_SIM=<file with function similarity for pairs of PDB files, optional>\n");
  printf("OUTGROUP=<Method to assign outgroups> ");
  printf("Allowed: TN (Tajima-Nei, default) CD (Contact divergence) TM (TM-score)\n"); 
  printf("NAME= <Name of output files> (default: alignment file)\n");
  printf("PDBDIR=<directory of pdb files>  (default: current directory)\n");
  printf("PDBEXT=<extension of pdb files>  (default: none)\n");
  printf("NORM=MIN  Normalization of seq.id, TM, overlap, MIN MAX or MEAN\n");
  printf("ALI_SS=<0,1>  Correct MSA with sec.structure?\n");
  printf("SHIFT_MAX=8   Maximum shift for correcting MSA\n");
  printf("SS_MULT=0     Sec,str. correction pairwise (0) or multiple (1)?\n");
  printf("The protein name is the name of a PDB file, optionally followed\n");
  printf("by the chain index (Ex: >1opd.pdb A or >1opdA or >1opd_A)\n\n");
  printf("PRINT_SIM=<0,1>   Print similarity measures? (default: 0)\n");
  printf("PRINT_DIV=<0,1>   Print divergence measures? (default: 1)\n");
  printf("PRINT_CV=<0,1>    Print clock violations? (default: 1)\n");
  printf("PRINT_PAIR=1      Print pairwise alignments made with PC\n");
  printf("PRINT_CLIQUE=1    Print MSA based on cliques of pairwise alis");
  printf("==========================================================\n");
  printf("OUTPUT (for each pair of proteins):\n");
  printf("File with extension .sim (similarity):\n");
  printf("Sequence identity SI (sequence alignment)\n");
  printf("Contact overlap q (structural, str.ali if present)\n");
  printf("TM-score TM (structural, str.ali if present), Zhang & Skolnick Proteins 2004 57:702\n");
  printf("\n");
  printf("File with extension .div (divergence):\n");
  printf("Tajima-Nei divergence TN=-log((SI-S0)/(1-S0) S0=%.2f (Tajima F & Nei 1984, Mol Biol Evol 1:269).\n", S0);
  printf("Contact divergence  CD=-log((q-q0(L))/(1-q0(L)) (Pascual-García et al Proteins 2010 78:181-96)\n");
  printf("TM_divergence=-log(TM)\n");
  printf("\n");
  printf("File with extension .cv (clock violations):\n");
  printf("CV=sum_c (D(a,c)-D(b,c))/(n_c*D(a,b)) (Pascual-Garcia, Arenas & Bastolla, submitted).\n");
  printf("\n");
  exit(8);
}

float Mean_freq(double *se, double *relerr2, double sum, double tot){
  if(tot==0)return(0);
  float p=sum/tot;
  if(sum){
    *relerr2=((1-p)/sum); //(sum*sum)
    *se=p*sqrt(*relerr2);
  }else{
    *relerr2=1; *se=0;
  }
  return(p);
}

void Print_ali_ss(int *ali_i, char *ss_i, int *ali_j, char *ss_j, int N_ali){
  for(int i=0; i<N_ali; i++){
    if(ali_i[i]<0){if(ali_j[i]>=0)printf("-");}
    else{printf("%c", ss_i[ali_i[i]]);}
  }
  printf("\n");
  for(int i=0; i<N_ali; i++){
    if(ali_j[i]<0){if(ali_i[i]>=0)printf("-");}
    else{printf("%c", ss_j[ali_j[i]]);}
  }
  printf("\n");
}

int Get_alignment_old(struct Prot_input **Prot_input, int *Nali,
		      char *PDB_DIR, char *PDB_EXT, int *PRINT_SIM,
		      char *file_ali)
{
  // Open file
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("ERROR, alignment file %s does not exist\n", file_ali); exit(8);
  }
  // Count proteins and read path
  char string[1000]; int dir=0, n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      n++;
    }else if(n==0){
      if(strncmp(string, "PDBDIR", 6)==0){
	sscanf(string+7,"%s", PDB_DIR);
	printf("Directory for PDB files: %s\n", PDB_DIR); dir=1;
      }else if(strncmp(string, "PDBEXT", 6)==0){
	sscanf(string+7, "%s", PDB_EXT);
      }else if(strncmp(string, "PRINT_SIM", 8)==0){
	sscanf(string+9, "%d", PRINT_SIM);
      }
    }
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no sequence found in file %s\n", file_ali); exit(8);
  }
  printf("%d sequences found in %s\n", n, file_ali);

  // Allocate and read sequences
  int LMAX=10000, l=0, i;
  char chain[10], dumm[40];
  char *Seq=malloc(LMAX*sizeof(char)), *s=NULL;
  *Prot_input=malloc(n*sizeof(struct Prot_input));
  n=-1; *Nali=0;
  file_in=fopen(file_ali, "r");
  if(dir)fgets(string, sizeof(string), file_in);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(string[0]=='>'){
      n++;
      sscanf(string+1, "%s", (*Prot_input)[n].name);
      for(i=0; i<10; i++)chain[i]='\0';
      int c=sscanf(string, "%s%s\n", dumm, chain);
      if((c>1)&&(chain[0]!='\0')&&(chain[0]!='\n')){
	(*Prot_input)[n].chain=chain[0];
      }else if(string[5]=='_'){
	printf("Getting chain after _\n");
	(*Prot_input)[n].name[4]='\0';
	(*Prot_input)[n].chain=string[6];
      }else if((string[5]>=65)&&(string[5]<=90)){ // Maiuscule
	printf("Getting chain after character 4\n");
	(*Prot_input)[n].name[4]='\0';
	(*Prot_input)[n].chain=string[5];
      }else{
	(*Prot_input)[n].chain=' ';
      }
      int ic=(*Prot_input)[n].chain;
      if((*Prot_input)[n].chain=='\n' ||
	 (*Prot_input)[n].chain=='\r'||
	 (ic<48) || (ic>122) ){
	(*Prot_input)[n].chain=' ';
      }
      printf("%s %c %d\n", (*Prot_input)[n].name, (*Prot_input)[n].chain, ic);
      if((*Nali==0)&&(l)){
	*Nali=l;
	(*Prot_input)[0].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[0].seq;
	for(l=0; l<*Nali; l++){*s=Seq[l]; s++;}
      }
      if(*Nali){
	(*Prot_input)[n].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[n].seq; l=0;
      }else{
	s=Seq; l=0;
      }
    }else if(n>=0){
      char *c=string;
      while(*c!='\n'){*s=*c; l++; s++; c++;}
      if(l > LMAX){
	printf("ERROR, alignment length larger than maximum allowed %d\n", l);
	printf("Increase LMAX in code %s\n", CODENAME); exit(8);
      }
      if((*Nali)&&(l>*Nali)){
	printf("ERROR, too many column in alignment %d.",n+1);
	printf(" Expected %d, found >= %d\n", *Nali, l); exit(8); 
      }
    }
  }
  n++;
  fclose(file_in);
  printf("%d sequences with %d columns found in MSA %s\n",
	 n, *Nali, file_ali);
  return(n);
}

int Get_alignment(struct Prot_input **Prot_input, int *Nali,
		  char *file_ali, char *PDB_DIR, char *PDB_EXT)
{
  // Open file
  if(file_ali[0]=='\0')return(0);
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("ERROR, alignment file %s does not exist\n", file_ali); exit(8);
  }
  // Count proteins and read path
  char string[1000]; int n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>')n++;
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no sequence found in file %s\n", file_ali); exit(8);
  }
  printf("%d sequences found in %s\n", n, file_ali);

  // Allocate and read sequences
  int LMAX=10000, l=0, i;
  char chain[10], dumm[40];
  char Seq[LMAX], *s=NULL;
  *Prot_input=malloc(n*sizeof(struct Prot_input));
  n=-1; *Nali=0;
  file_in=fopen(file_ali, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string,"PDBDIR",6)==0){
      if(PDB_DIR[0]=='\0'){
	sscanf(string+7, "%s", PDB_DIR);
      }else{
	printf("WARNING reading line %sPDB_DIR is already set.",string);
	printf(" Keeping previous value %s\n",PDB_DIR);
      }
    }else if(string[0]=='>'){
      n++;
      sscanf(string+1, "%s", (*Prot_input)[n].name);
      for(i=0; i<10; i++)chain[i]='\0';
      int c=sscanf(string, "%s%s\n", dumm, chain);
      if((c>1)&&(chain[0]!='\0')&&(chain[0]!='\n')){
	(*Prot_input)[n].chain=chain[0];
      }else if(string[5]=='_'){
	printf("Getting chain after _\n");
	(*Prot_input)[n].name[4]='\0';
	(*Prot_input)[n].chain=string[6];
      }else if((string[5]>=65)&&(string[5]<=90)){ // Maiuscule
	printf("Getting chain after character 4\n");
	(*Prot_input)[n].name[4]='\0';
	(*Prot_input)[n].chain=string[5];
      }else{
	(*Prot_input)[n].chain=' ';
      }
      if((*Prot_input)[n].chain=='\n' ||
	 (*Prot_input)[n].chain=='\r'){
	(*Prot_input)[n].chain=' ';
      }
      printf("%s %c\n", (*Prot_input)[n].name, (*Prot_input)[n].chain);
      if((*Nali==0)&&(l)){
	*Nali=l;
	(*Prot_input)[0].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[0].seq;
	for(l=0; l<*Nali; l++){*s=Seq[l]; s++;}
      }
      if(*Nali){
	(*Prot_input)[n].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[n].seq; l=0;
      }else{
	s=Seq; l=0;
      }
    }else{
      char *c=string;
      while((*c!='\n')&&(*c!='\r')&&(*c!='\0')&&(*c!=' ')){
	if(l > LMAX){
	  printf("ERROR, alignment %s length %d > maximum allowed %d\n",
		 (*Prot_input)[n].name, l, LMAX);
	  printf("Increase LMAX in code %s\n", CODENAME); exit(8);
	}
	if((l>*Nali)&&(*Nali)){
	  // Check if all additional characters are gaps
	  int allgaps=1; char *d=c; int i=l;
	  while((*d!='\n')&&(*d!='\r')&&(*d!='\0')&&(*d!=' ')){
	    if(*d!='-'){allgaps=0;} d++; i++;
	  }
	  if(allgaps==0){
	    printf("ERROR, too many columns in sequence %s",
		   (*Prot_input)[n].name);
	    printf(" Expected %d, found %d\n", *Nali, i);
	    for(i=0; i<l; i++)printf("%c", (*Prot_input)[n].seq[i]);
	    printf("+");
	    char *d=c; int k=n-1;
	    while((*d!='\n')&&(*d!='\r')&&(*d!='\0')&&(*d!=' ')){
	      printf("%c",*d); d++;
	    }
	    printf("\n");
	    printf("Previous seq %s:\n",(*Prot_input)[k].name);
	    for(i=0; i<*Nali; i++)printf("%c", (*Prot_input)[k].seq[i]);
	    printf("\n");
	    exit(8); // end error
	  }else{
	    printf("WARNING, length %d >= Nali=%d but all gaps\n",l,*Nali);
	    break;
	  } // end test all gaps
	}// end l>Nali
	*s=*c; l++; s++; c++;
      }// end read line
    } // end type of line
  } // end read file
  n++;
  fclose(file_in);
  printf("%d sequences with %d columns found in MSA %s\n",
	 n, *Nali, file_ali);
  return(n);
}

void Write_identity(FILE *file_id, int type,
		    struct protein *proti, struct protein *protj,
		    int *ali_ij, int *id_aa, int *id_sup, int *shift,
		    int *neigh_ali,int *neigh_noali, int *neigh_noali_aaid,
		    int *ali_cont_ct, int *id_cont_ct,
		    int *ali_cont_sup, int *id_cont_sup,
		    int *ali_cont_aaid, int *id_cont_aaid)
{
  if(ini_sum==0){
    ini_sum=1;
    for(int i=0; i<ITYPE; i++){
      for(int j=0; j<=IBIN; j++){
	for(int k=0; k<2; k++){
	  sum_all_ct[i][k][j]=0;
	  sum_ali_ct[i][k][j]=0;
	  sum_aaid_ct[i][k][j]=0; 
	  sum_sup_ct[i][k][j]=0; sum_aaid_sup_ct[i][k][j]=0;
	  sum_noid_nosup_ct[i][k][j]=0;
	  sum_ali_cont_ct[i][k][j]=0; sum_id_cont_ct[i][k][j]=0;
	  sum_neigh_ali_ct[i][k][j]=0; sum_neigh_noali_ct[i][k][j]=0;
	  sum_neigh_noali_aaid_ct[i][k][j]=0;
	  sum_sup_noneigh_ct[i][k][j]=0;
	  sum_shift_sup_noneigh_ct[i][k][j]=0;
	}
	for(int k=0; k<3; k++){
	  sum_ali_cont_sup[i][k][j]=0; sum_id_cont_sup[i][k][j]=0;
	  sum_ali_cont_aaid[i][k][j]=0; sum_id_cont_aaid[i][k][j]=0;
	}
      }
    }
    for(int k=0; k<2; k++){
      sum_SS_SI[k]=0; sum_SS_TM[k]=0; sum_SS_CO[k]=0;
      sum_SS_SI_TM[k]=malloc(2*sizeof(double));
      sum_SS_SI_CO[k]=malloc(2*sizeof(double));
      sum_SS_TM_CO[k]=malloc(2*sizeof(double));
      for(int j=0; j<2; j++){
	sum_SS_SI_TM[k][j]=0; sum_SS_SI_CO[k][j]=0; sum_SS_TM_CO[k][j]=0;
      }
    }
  }

  int k, noali[2], nali[2], n_sup[2], n_aaid[2], n_aaid_sup[2], n_noid_nosup[2],
    n_neigh_ali[2], n_neigh_noali[2], n_neigh_noali_aaid[2],
    n_sup_noneigh[2], shift_sup_noneigh[2];
  for(k=0; k<2; k++){
    noali[k]=0; nali[k]=0; n_aaid[k]=0; n_sup[k]=0;
    n_aaid_sup[k]=0; n_noid_nosup[k]=0; 
    n_neigh_ali[k]=0; n_neigh_noali[k]=0; n_neigh_noali_aaid[k]=0;
    n_sup_noneigh[k]=0; shift_sup_noneigh[k]=0;
  }
  // Not aligned sites
  if(proti->len<protj->len){
    Count_noali(noali, ali_ij,  proti->ncont, proti->len, c_ave);
  }else{
    int ali2[protj->len]; Invert_ali(ali2, protj->len, ali_ij, proti->len);
    Count_noali(noali, ali2, protj->ncont, protj->len, c_ave);
  }
  for(int i=0; i<proti->len; i++){
    if(ali_ij[i]<0){
      continue; // These are not errors in alignments!
      if(neigh_noali[i]==0)continue;
      if(proti->ncont[i]<=c_ave){k=0;}// few contacts
      else{k=1;} // many contacts
      n_neigh_noali[k]++;
      if(neigh_noali_aaid[i])n_neigh_noali_aaid[k]++;
      continue;
    }
    float c=sqrt(proti->ncont[i]*protj->ncont[ali_ij[i]]);
    if(c<=c_ave){k=0;} // few contacts
    else{k=1;}  // many contacts
    nali[k]++;
    if(id_aa[i]){n_aaid[k]++; if(id_sup[i])n_aaid_sup[k]++;}
    else{if(id_sup[i]==0)n_noid_nosup[k]++;}
    if(id_sup[i]){
      n_sup[k]++;
      if(neigh_ali[i]==0){
	n_sup_noneigh[k]++; shift_sup_noneigh[k]+=shift[i];
      }
    }
    if(neigh_ali[i]){n_neigh_ali[k]++;}
    else if(neigh_noali[i]){n_neigh_noali[k]++;
      if(neigh_noali_aaid[i])n_neigh_noali_aaid[k]++;
    }
  }

  //int j=(n_aaid[0]+n_aaid[1])*IBIN/(float)(nali[0]+nali[1]);
  int j=(n_aaid[0]+n_aaid[1]), ntot=nali[0]+nali[1];
  if(ntot){j*=(IBIN/(float)ntot);}
  for(k=0; k<2; k++){
    sum_all_ct[type][k][j]+=(noali[k]+nali[k]);
    sum_ali_ct[type][k][j]+=nali[k];
    sum_sup_ct[type][k][j]+=n_sup[k];
    sum_aaid_ct[type][k][j]+=n_aaid[k];
    sum_aaid_sup_ct[type][k][j]+=n_aaid_sup[k];
    sum_noid_nosup_ct[type][k][j]+=n_noid_nosup[k];
    sum_neigh_ali_ct[type][k][j]+=n_neigh_ali[k];
    sum_neigh_noali_ct[type][k][j]+=n_neigh_noali[k];
    sum_neigh_noali_aaid_ct[type][k][j]+=n_neigh_noali_aaid[k];
    sum_sup_noneigh_ct[type][k][j]+=n_sup_noneigh[k];
    sum_shift_sup_noneigh_ct[type][k][j]+=shift_sup_noneigh[k];
    //
    sum_ali_cont_ct[type][k][j]+=ali_cont_ct[k];
    sum_id_cont_ct[type][k][j]+=id_cont_ct[k];
  }
  for(k=0; k<3; k++){
    sum_ali_cont_sup[type][k][j]+=ali_cont_sup[k];
    sum_id_cont_sup[type][k][j]+=id_cont_sup[k];
    sum_ali_cont_aaid[type][k][j]+=ali_cont_aaid[k];
    sum_id_cont_aaid[type][k][j]+=id_cont_aaid[k];
  }
}


void Get_file_name(char *name, char *file){
  char *s=file, *t=name;
  while(*s!='\0'){
    if(*s=='/'){s++; t=name;}
    else if(*s=='.'){break;}
    *t=*s; s++; t++;
  }
}

void Set_contact_type(){

  if(CONT_TYPE=='a'){strcpy(CONT_STRING, "Alpha");}
  else if(CONT_TYPE=='b'){strcpy(CONT_STRING, "Beta");}
  else if(CONT_TYPE=='c'){strcpy(CONT_STRING, "All atoms");}
  else{
    printf("WARNING, undefined contact %c\n", CONT_TYPE);
    CONT_TYPE='c'; strcpy(CONT_STRING, "All atoms");
    printf("Using default %s\n", CONT_STRING);
  }
  // Default type of contacts?
  if((CONT_TYPE!=CONT_TYPE_DEF)||(CONT_THR!=CONT_THR_DEF)||
     (IJ_MIN!=IJ_MIN_DEF))CONT_DEF=0;
}

int Count_AA(char *seq, int N_ali){
  int L=0; char *s=seq;
  for(int i=0; i<N_ali; i++){if(*s!='-')L++; s++;}
  return(L);
}

int Seq_differences(int *id, char *seq1, char *seq2, int N_ali){
  // Consider only aligned positions
  int d=0; (*id)=0; char *s1=seq1, *s2=seq2;
  for(int i=0; i<N_ali; i++){
    if((*s1!='-')&&(*s2!='-')){ // ||
      if(*s1!=*s2){d++;}else{(*id)++;}
    }
    s1++; s2++;
  }
  return(d);
}


int Count_gaps(char *seq1, char *seq2, int N_ali){
  int ngap=0, open=0; char *s1=seq1, *s2=seq2;
  for(int i=0; i<N_ali; i++){
    if((*s1!='-')&&(*s2!='-')){open=0;}
    else if(open==0){ngap++; open=1;}
    s1++; s2++;
  }
  return(ngap);
}

float Min_CV(int i, int j, int k, int **conformation, int *N_conf, float **div)
{
  float CV=10000;
  for(int k1=0; k1<N_conf[k]; k1++){
    float *dk=div[conformation[k][k1]];
    for(int i1=0; i1<N_conf[i]; i1++){
      float d_ki=dk[conformation[i][i1]];
      for(int j1=0; j1<N_conf[j]; j1++){
	float cv=d_ki-dk[conformation[j][j1]];
	if(fabs(cv)<fabs(CV)){CV=cv;}
      }
    }
  }
  return(CV);
}

void Change_conformations(int **conformation, int N_seq, int *N_conf, int *ali_str)
{
  for(int i=0; i<N_seq; i++){
    for(int k=0; k<N_conf[i]; k++){
      int l=conformation[i][k];
      conformation[i][k]=ali_str[l];
    }
  }
}

float Min_dist(int i, int j, int **conformation, int *N_conf, float **div)
{
  float d=1000; int ini=1;
  for(int i1=0; i1<N_conf[i]; i1++){
    float *d_i=div[conformation[i][i1]];
    for(int j1=0; j1<N_conf[j]; j1++){
      float d_il=d_i[conformation[j][j1]];
      if(ini){d=d_il; ini=0;}
      else if(d_il<d){d=d_il;}
    }
  }
  return(d);
}

int *Match_alignments(struct Prot_input *Prot1,
		      struct Prot_input *Prot2, int N)
{
  int *ali_str=malloc(N*sizeof(int));
  struct Prot_input *P1=Prot1;
  for(int i=0; i<N; i++){
    ali_str[i]=-1;
    for(int j=0; j<N; j++){
      if((strcmp(Prot2[j].name, P1->name)==0)&&
	 (Prot2[j].chain==P1->chain)){
	ali_str[i]=j; break;
      }
    }
    if(ali_str[i]<0){
      printf("WARNING, protein %s%c i=%d N=%d not matched\n",
	     P1->name, P1->chain, i, N);
      free(ali_str); return(NULL);
    } 
    P1++;
  }

  return(ali_str);
}

void Get_input(char *file_ali, char *file_ali_str, char *file_fun,
	       char *name_in,  char *PDB_DIR, char *PDB_EXT, char *OUTG,
	       int *NORM, int *ALI_SS, int *SHIFT_MAX, int *SS_MULT,
	       int *PRINT_SIM, int *PRINT_CV, int *PRINT_DIV,
	       int *PRINT_PAIR, int *PRINT_CLIQUE,
	       int argc, char **argv)
{
  printf("Starting %s\n", argv[0]);
  if((argc<2)||(strncmp(argv[1], "-h", 2)==0))help(argv[0]);

  FILE *file_in=fopen(argv[1], "r");
  if(file_in==NULL){
    printf("ERROR, input file %s does not exist\n\n", argv[1]);
    exit(8);
  }
  char string[100], read[80]; int ali=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "ALI=", 4)==0){
      sscanf(string+4,"%s", file_ali);
      printf("File with multiple sequence alignment: %s\n", file_ali);
      ali=1;
    }else if(strncmp(string, "STR_ALI=", 8)==0){
      sscanf(string+8,"%s", file_ali_str);
      printf("File with multiple structure alignment: %s\n", file_ali_str);
    }else if(strncmp(string, "NAME", 4)==0){
      sscanf(string+5, "%s", name_in);
    }else if(strncmp(string, "FUN_SIM", 7)==0){
      sscanf(string+8, "%s", file_fun);
    }else if(strncmp(string, "OUTGROUP", 8)==0){
      sscanf(string+9, "%s", read);
      if((strcmp(read,"TN")!=0) && 
	 (strcmp(read,"CD")!=0) && 
	 (strcmp(read,"TM")!=0)){
	printf("WARNING, outgroup method %s not implemented, using default\n",
	       read);
      }else{
	strcpy(OUTG, read);
      }
      printf("Outgroup method set to %s\n", OUTG);
    }else if(strncmp(string, "PDBDIR", 6)==0){
      sscanf(string+7,"%s", PDB_DIR);
      printf("Directory for PDB files: %s\n", PDB_DIR);
   }else if(strncmp(string, "NORM", 4)==0){
      char tmp[80]; sscanf(string+5,"%s", tmp);
      if((strcmp(tmp,"MIN")==0)||
	 (strcmp(tmp,"Min")==0)||
	 (strcmp(tmp,"min")==0)){*NORM=0;}
      else if((strcmp(tmp,"MAX")==0)||
	 (strcmp(tmp,"Max")==0)||
	 (strcmp(tmp,"max")==0)){*NORM=1;}
      else if((strcmp(tmp,"MEAN")==0)||
	 (strcmp(tmp,"Mean")==0)||
	 (strcmp(tmp,"mean")==0)){*NORM=2;}
      else{
	printf("WARNING, %s is not an allowed normalization type\n",tmp);
	printf("Keeping default MIN\n");
      }
    }else if(strncmp(string, "ALI_SS", 6)==0){
      sscanf(string+7,"%d", ALI_SS);
    }else if(strncmp(string, "SS_MULT", 7)==0){
      sscanf(string+8,"%d", SS_MULT);
    }else if(strncmp(string, "SHIFT_MAX", 9)==0){
      int tmp;
      sscanf(string+10,"%d", &tmp);
      if(tmp>=0){*SHIFT_MAX=tmp;}
      else{printf("WARNING, %d not allowed value for SHIFT_MAX\n", tmp);}
    }else if(strncmp(string, "PDBEXT", 6)==0){
      sscanf(string+7, "%s", PDB_EXT);
    }else if(strncmp(string, "PRINT_SIM", 9)==0){
      sscanf(string+10, "%d", PRINT_SIM);
    }else if(strncmp(string, "PRINT_CV", 8)==0){
      sscanf(string+9, "%d", PRINT_CV);
    }else if(strncmp(string, "PRINT_DIV", 9)==0){
      sscanf(string+10, "%d", PRINT_DIV);
    }else if(strncmp(string, "PRINT_PAIR", 10)==0){
      sscanf(string+11, "%d", PRINT_PAIR);
    }else if(strncmp(string, "PRINT_CLIQUE", 12)==0){
      sscanf(string+13, "%d", PRINT_CLIQUE);
    }else{
      printf("WARNING, unknown command %s", string);
    }
  }
  if(ali==0){
    printf("ERROR, file with alignment not provided in %s\n", argv[1]);
    exit(8);
  }
  printf("Modifying MSA with sec.str.information: ");
  if(ALI_SS){printf("YES\n");}else{printf("NO\n");}
  if(name_in[0]=='\0')Get_file_name(name_in, file_ali);
}

float **Read_function(char *file_fun, struct Prot_input *Prot,
		      int *index, int N)
{
  if(file_fun[0]=='\0')return(NULL);
  FILE *file_in=fopen(file_fun, "r");
  if(file_in==NULL)return(NULL);

  int n=0, i, j;
  float **f_sim=malloc(N*sizeof(float *)), sim;
  for(i=0; i<N; i++){
    f_sim[i]=malloc(N*sizeof(float));
    for(j=0; j<N; j++)f_sim[i][j]=-1;
  }
  char string[1000], name1[40], name2[40];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    sscanf(string, "%s%s%f", name1, name2, &sim);
    i=Find_prot(name1, Prot, index, N);
    j=Find_prot(name2, Prot, index, N);
    if((i<0)||(j<0)){
      printf("WARNING, pair %s %s not identified in function file %s\n",
	     name1, name2, file_fun);
    }else{
      f_sim[i][j]=sim; f_sim[j][i]=sim; n++;
    }
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no pair was identified in function file %s\n",
	   file_fun);
    Empty_matrix_f(f_sim, N);
    return(NULL);
  }else{
    printf("%d protein pairs identified in %s\n", n, file_fun);
    return(f_sim);
  }
}

int Find_prot(char *name, struct Prot_input *Prot, int *index, int N){
  for(int i=0; i<N; i++){
    int i1=index[i];
    if((strncmp(name, Prot[i1].name, 4)==0)&&(name[4]==Prot[i1].chain))
      return(i);
  }
  return(-1);
}

float Divergence(float SI, float S0, float HUGE)
{
  float Div;
  if(SI>S0){
    Div=-log((SI-S0)/(1.-S0)); if(Div>HUGE)Div=HUGE;
  }else{
    Div=HUGE;
  }
  return(Div);
}
void Print_seq(int *ali, char *seq, int N_ali){
  for(int k=0; k<N_ali; k++){
    if(ali[k]<0){printf("-");}
    else{printf("%c",seq[ali[k]]);}
  }
  printf("\n");
}
int Count_ali(int *ali, int len){
  int n=0, k; for(k=0; k<len; k++)if(ali[k]>=0)n++;
  return(n);
}

void Summary_identical(FILE *file_id)
{
  char *name_c[2], *name_ct[2]; int i, j, k;
  name_c[0]=malloc(9*sizeof(char)); strcpy(name_c[0], "fct");
  name_c[1]=malloc(9*sizeof(char)); strcpy(name_c[1], "mct");
  name_ct[0]=malloc(20*sizeof(char)); strcpy(name_ct[0], "few contacts");
  name_ct[1]=malloc(20*sizeof(char)); strcpy(name_ct[1], "many contacts");
  fprintf(file_id,"# ali=aligned, ide=identical aa,");
  fprintf(file_id," sup=superimposed in space, d<d0(TM)/2");
  fprintf(file_id," cons cont=conserved contact\n");
  fprintf(file_id,"# fct=residues with few structural contacts,");
  fprintf(file_id," mct= many contacts\n");
  fprintf(file_id,"# npairs= %.0f\n", pairs);
  for(i=2; i<=5; i++){
    char what[4];
    if(i==4){continue;}
    else if(i==2){strcpy(what, "TM");}
    else if(i==3){strcpy(what, "CO");}
    else if(i==5){strcpy(what, "PC");}
    fprintf(file_id,
	    "# Fraction of pairs for which %s_ali has not max %s: %.3f\n",
	    what, what, Diff_opt[i]/pairs);
  }
  fprintf(file_id,"# Fraction of pairs where the named alignment has max PC:");
  for(i=0; i<(ITYPE-1); i++)
    if(nali1[i])fprintf(file_id,"  %s %.3f", ali_name[i], opt_PC[i]/pairs);
  fprintf(file_id,"\n");

  //
  fprintf(file_id,"# Secstr_corrected_alignment:\n");
  char id[3][4]; strcpy(id[0],"SI"); strcpy(id[1],"TM"); strcpy(id[2],"CO");
  char sign[2]; sign[0]='<'; sign[1]='>';
  double *count[3];
  count[0]=sum_SS_SI;count[1]=sum_SS_TM;count[2]=sum_SS_CO;
  double ave=0, se=0, re2A=0, re2B=0, ratio=0, ser=0;
  for(i=0; i<3; i++){
    fprintf(file_id,"# ");
    for(j=0; j<2; j++){
      ave=Mean_freq(&se, &re2A, count[i][j], pairs);
      fprintf(file_id," P(%s%c%s(Input)= %.4f s.e.=%.4f",
	      id[i],sign[j],id[i],ave, se);
    }
    fprintf(file_id,"\n");
  }
  double **count2[3]; count2[0]=sum_SS_SI_TM;
  count2[1]=sum_SS_SI_CO;count2[2]=sum_SS_TM_CO; int i2=0;
  for(i=0; i<3; i++){
    for(k=i+1; k<3; k++){
      fprintf(file_id,"# ");
      for(j=0; j<2; j++)
	for(int j2=0; j2<2; j2++){
	  float PA=Mean_freq(&se, &re2A, count2[i2][j][j2], count[i][j]);
	  float PB=Mean_freq(&se, &re2B, count[k][j2], pairs);
	  double ratio=PA/PB, ser=ratio*sqrt(re2A+re2B);
	  /*fprintf(file_id," P(%s%c,%s%c)/P(%s%c)P(%s%c)",
		  id[i],sign[j],id[k],sign[j2],
		  id[i],sign[j],id[k],sign[j2]);*/
	  fprintf(file_id," Prop(%s%c,%s%c)", id[i],sign[j],id[k],sign[j2]);
	  fprintf(file_id,"= %.2f se= %.2f",ratio, ser);
	}
      fprintf(file_id,"\n"); i2++;
    }
  }

  for(i=0; i<ITYPE; i++){
    double sum, tot, sum_all[2], tot_all=0, re2;
    for(k=0; k<2; k++){
      sum_all[k]=Sum_bins(sum_all_ct[i][k]);
      tot_all+=sum_all[k];
    }
    if(tot_all==0)continue;
    float nind=sqrt(2*pairs+0.25)+0.5;
    for(j=i+1; j<ITYPE; j++){
      if(SI1[j]==0){continue;} double d1, d2;
      fprintf(file_id, "# Diff_ali_%s_vs_%s:", ali_name[j],ali_name[i]);
      d1=Ave_se(&d2, SI_diff1[i][j], SI_diff2[i][j], pairs, nind);
      fprintf(file_id, "\tSI: %.4f %.4f", d1, d2);
      d1=Ave_se(&d2, TM_diff1[i][j], TM_diff2[i][j], pairs, nind);
      fprintf(file_id, "\tTM: %.4f %.4f", d1, d2);
      d1=Ave_se(&d2, CO_diff1[i][j], CO_diff2[i][j], pairs, nind);
      fprintf(file_id, "\tCO: %.4f %.4f", d1, d2);
      d1=Ave_se(&d2, PC_diff1[i][j], PC_diff2[i][j], pairs, nind);
      fprintf(file_id, "\tPC: %.4f %.4f", d1, d2);
      d1=Ave_se(&d2, ali_diff1[i][j], ali_diff2[i][j], pairs, nind);
      fprintf(file_id, "\tali: %.4f %.4f\n", d1, d2);
    }

    fprintf(file_id,"######## %s alignment\n", ali_name[i]);
    fprintf(file_id, "### Pairwise scores:\n");
    fprintf(file_id, "# Mean ali: %.4g\n", nali1[i]/pairs);
    fprintf(file_id, "# Mean SI:  %.4f\n", SI1[i]/pairs);
    fprintf(file_id, "# Mean TM:  %.4f\n", TM1[i]/pairs);
    fprintf(file_id, "# Mean CO:  %.4f\n", CO1[i]/pairs);
    fprintf(file_id, "# Mean PC:  %.4g\n", PC1[i]/pairs);
    fprintf(file_id, "### Conservation across all residues:\n");
    double sum_ali[2];
    for(k=0; k<2; k++)sum_ali[k]=Sum_bins(sum_ali_ct[i][k]);
    ave=Mean_freq(&se, &re2, sum_ali[0]+sum_ali[1], tot_all);
    fprintf(file_id,"# Frac_ali_vs_all:\t%.4f se %.4f n %.0f\t",
	    ave, se, tot_all);
    fprintf(file_id,"P(aligned)\n");  

    double sum_ide[2];
    for(k=0; k<2; k++)sum_ide[k]=Sum_bins(sum_aaid_ct[i][k]);
    ave=Mean_freq(&se, &re2, sum_ide[0]+sum_ide[1], tot_all);
    fprintf(file_id,"# Frac_ide_vs_all:\t%.4f se %.4f n %.0f\t",
	    ave, se, tot_all);
    fprintf(file_id,"P(identical)\n");  
    
    double sum_sup[2];
    for(k=0; k<2; k++)sum_sup[k]=Sum_bins(sum_sup_ct[i][k]);
    ave=Mean_freq(&se, &re2, sum_sup[0]+sum_sup[1], tot_all);
    fprintf(file_id,"# Frac_sup_vs_all:\t%.4f se %.4f n %.0f\t",
	    ave, se, tot_all);
    fprintf(file_id,"P(superimposed)\n");  
    
    double sum_id_cont[2], sum_ali_cont[2]; sum=0; tot=0;
    for(k=0; k<2; k++){
      sum_id_cont[k]=Sum_bins(sum_id_cont_ct[i][k]); 
      sum_ali_cont[k]=Sum_bins(sum_ali_cont_ct[i][k]);
      sum+=sum_id_cont[k]; tot+=sum_ali_cont[k];
    }
    ave=Mean_freq(&se, &re2, sum, tot);
    fprintf(file_id,"# Frac_cons_cont:\t%.4f se %.4f n %.0f\t",
	    ave, se, tot);
    fprintf(file_id, "P(conserved cont) norm. by all cont\n");

    fprintf(file_id, "### Conservation many contacts vs few contacts:\n");
    // As a function of the number of contacts
    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_ali[k], sum_all[k]);
      fprintf(file_id,"# Frac_ali_vs_all_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_all[k]);
      fprintf(file_id,"P(aligned | %s)\n", name_ct[k]);  
    }
    
    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_ide[k], sum_ali[k]);
      fprintf(file_id,"# Frac_ide_vs_ali_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_ali[k]);
      fprintf(file_id,"P(identical | aligned, %s)\n", name_ct[k]);  
    }

    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_sup[k], sum_ali[k]);
      fprintf(file_id,"# Frac_sup_vs_ali_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_ali[k]);
      fprintf(file_id,"P(superimposed | aligned, %s)\n", name_ct[k]);  
    }
 
    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_id_cont[k], sum_ali_cont[k]);
      fprintf(file_id,
	      "# Frac_cons_cont_vs_ali_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_ali_cont[k]);
      fprintf(file_id,"P(conserved cont | aligned cont, %s)\n", name_ct[k]);  
    }       

    // Neighbors in space
    fprintf(file_id, "### Spatial neighbors even if not aligned:\n");
    double sum_neigh_ali[2], sum_neigh_noali[2], tot_neigh[2]; sum=0;
    for(k=0; k<2; k++){
      sum_neigh_ali[k]=Sum_bins(sum_neigh_ali_ct[i][k]);
      sum_neigh_noali[k]=Sum_bins(sum_neigh_noali_ct[i][k]);
      tot_neigh[k]=sum_neigh_ali[k]+sum_neigh_noali[k];
      ave=Mean_freq(&se, &re2, tot_neigh[k], sum_ali[k]);
      fprintf(file_id,"# Frac_neigh_vs_ali_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_ali[k]);
      fprintf(file_id,"P(neigh| %s) aligned or not\n", name_ct[k]);  
    }
    
    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_neigh_noali[k], tot_neigh[k]);
      fprintf(file_id,"# Frac_neigh_noali_vs_neigh_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, tot_neigh[k]);
      fprintf(file_id,"P(neigh, not ali| neigh, %s)\n", name_ct[k]);
      /*sum=0; for(j=0; j<=IBIN; j++)sum+=sum_neigh_noali_aaid_ct[i][k][j];
      ave=Mean_freq(&se, &re2, sum, sum_neigh_noali[k]);
      fprintf(file_id,
	    "# Frac_neigh_noali_id_ali_vs_neigh_noali:\t%.4f se %.4f n %.0f\t",
	    ave, se, sum);
	    fprintf(file_id,"P(id ali | neighb not ali)\n");  */
    }

    for(k=0; k<2; k++){
      sum=Sum_bins(sum_sup_noneigh_ct[i][k]);
      ave=Mean_freq(&se, &re2, sum, sum_sup[k]);
      fprintf(file_id,"# Frac_sup_noneigh_vs_sup_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_sup[k]);
      double sum_shift=Sum_bins(sum_shift_sup_noneigh_ct[i][k]);
      if(sum)sum_shift/=sum;
      fprintf(file_id,"shift=%.1f\t",sum_shift);
      fprintf(file_id,"P(no neigh | ali sup, %s)\n", name_ct[k]);
    }

    // Propensity superimposed-identical
    fprintf(file_id, "### Propensity superimposed - identical:\n");
    double sum_id_sup[2], sum_noid_nosup[2];
    for(k=0; k<2; k++){
      sum_id_sup[k]=Sum_bins(sum_aaid_sup_ct[i][k]);
      sum_noid_nosup[k]=Sum_bins(sum_noid_nosup_ct[i][k]);
    }

    for(k=0; k<2; k++){
      double P_sup_id=Mean_freq(&se, &re2A, sum_id_sup[k], sum_ide[k]);
      double sup_noid=sum_sup[k]-sum_id_sup[k];   
      double noid=sum_ali[k]-sum_ide[k];
      double P_sup_noid=Mean_freq(&se, &re2B, sup_noid, noid);
      ratio=P_sup_id/P_sup_noid; ser=ratio*sqrt(re2A+re2B);
      fprintf(file_id,"# P(sup|ide)/P(sup|no id)_%s:", name_c[k]);
      fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, sum_ide[k], noid);

      double P_nosup_noid=Mean_freq(&se, &re2A, sum_noid_nosup[k], noid);
      double nosup=sum_ali[k]-sum_sup[k];
      double nosup_id=nosup-sum_noid_nosup[k];   
      double P_nosup_id=Mean_freq(&se, &re2B, nosup_id, sum_ide[k]);
      ratio=P_nosup_noid/P_nosup_id; ser=ratio*sqrt(re2A+re2B);
      fprintf(file_id,"# P(no sup|no id)/P(no sup|id)_%s:", name_c[k]);
      fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, noid, sum_ide[k]);

      double P_id_sup=Mean_freq(&se, &re2A, sum_id_sup[k], sum_sup[k]);
      double id_nosup=sum_ide[k]-sum_id_sup[k];
      double P_id_nosup=Mean_freq(&se, &re2B, id_nosup, nosup);
      ratio=P_id_sup/P_id_nosup; ser=ratio*sqrt(re2A+re2B);
      fprintf(file_id,"# P(ide|sup)/P(ide|no sup)_%s:", name_c[k]);
      fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, sum_sup[k], nosup);

      double P_noid_nosup=Mean_freq(&se, &re2A, sum_noid_nosup[k], nosup);
      double noid_sup=noid-sum_noid_nosup[k];
      double P_noid_sup=Mean_freq(&se, &re2B, noid_sup, sum_sup[k]);
      ratio=P_noid_nosup/P_noid_sup; ser=ratio*sqrt(re2A+re2B);
      fprintf(file_id,"# P(no id|no sup)/P(no id|sup)_%s:", name_c[k]);
      fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, nosup, sum_sup[k]);

    } 

    // Conserved contacts
    double id_cont_sup[3], ali_cont_sup[3];
    double id_cont_aaid[3], ali_cont_aaid[3];
    double id_cont=0, ali_cont=0;
    for(k=0; k<3; k++){
      id_cont_sup[k]=0; ali_cont_sup[k]=0; 
      id_cont_aaid[k]=0; ali_cont_aaid[k]=0; 
      for(j=0; j<=IBIN; j++){
	id_cont_sup[k]+=sum_id_cont_sup[i][k][j];
	ali_cont_sup[k]+=sum_ali_cont_sup[i][k][j];
	id_cont_aaid[k]+=sum_id_cont_aaid[i][k][j];
	ali_cont_aaid[k]+=sum_ali_cont_aaid[i][k][j];
      }
      id_cont+=id_cont_sup[k];
      ali_cont+=ali_cont_sup[k];
    }
    // Conditional Prob conserved contact | superimposed
    fprintf(file_id, "### Propensity ide or sup | ide or sup in contact\n");

    double sum_a=sum_ali[0]+sum_ali[1]; //0=few cont 1=many cont
    double sum_i=sum_ide[0]+sum_ide[1];
    double sum_s=sum_sup[0]+sum_sup[1];
    double re2i,Pi=Mean_freq(&se, &re2i, sum_i, sum_a);//re2noi=re2i*(1-Pi)/Pi;
    double re2s,Ps=Mean_freq(&se, &re2s, sum_s, sum_a);//re2nos=re2s*(1-Ps)/Ps;
    
    sum_i=ali_cont_aaid[1]+2*ali_cont_aaid[2];
    double sum_a2=2*(ali_cont_aaid[0]+ali_cont_aaid[1]+ali_cont_aaid[2]);
    double re2ic=0, Pic=Mean_freq(&se, &re2ic, sum_i, sum_a2),
      re2noic=re2ic*(1-Pic)/Pic;
    fprintf(file_id, "P(id|cont)/P(id)= %.3f\n", Pic/Pi);
    sum_s=ali_cont_sup[1]+2*ali_cont_sup[2];
    sum_a2=2*(ali_cont_sup[0]+ali_cont_sup[1]+ali_cont_sup[2]);
    double re2sc=0, Psc=Mean_freq(&se, &re2sc, sum_s, sum_a2),
      re2nosc=re2sc*(1-Psc)/Psc;
    fprintf(file_id, "P(sup|cont)/P(sup)= %.3f\n", Psc/Ps);

    // Correlated substitutions
    sum=2*ali_cont_aaid[0]+ali_cont_aaid[1];
    double PA=Mean_freq(&se, &re2A, 2*ali_cont_aaid[0], sum);
    ratio=PA/(1-Pic); ser=ratio*sqrt(re2A+re2noic);
    fprintf(file_id,"# P(no id|cont no id)/P(no id):");
    fprintf(file_id,"\t%.2f se %.2f n %.0f %.0f\n",ratio,ser, sum, sum_a);

    //sum=2*ali_cont_aaid[0]+ali_cont_aaid[1];
    PA=Mean_freq(&se, &re2A, ali_cont_aaid[1], sum);
    ratio=PA/Pic; ser=ratio*sqrt(re2A+re2ic);
    fprintf(file_id,"# P(id|cont no id)/P(id):");
    fprintf(file_id,"\t%.2f se %.2f n %.0f %.0f\n", ratio,ser, sum, sum_a);

    sum=ali_cont_aaid[1]+2*ali_cont_aaid[2];
    PA=Mean_freq(&se, &re2A, 2*ali_cont_aaid[2], sum);
    ratio=PA/Pic; ser=ratio*sqrt(re2A+re2ic);
    fprintf(file_id,"# P(id|cont id)/P(id):");
    fprintf(file_id,"\t%.2f se %.2f n %.0f %.0f\n", ratio,ser, sum, sum_a);

    //sum=ali_cont_aaid[1]+2*ali_cont_aaid[2];
    PA=Mean_freq(&se, &re2A, ali_cont_aaid[1], sum);
    ratio=PA/(1-Pic); ser=ratio*sqrt(re2A+re2noic);
    fprintf(file_id,"# P(no id|cont id)/P(no id):");
    fprintf(file_id,"\t%.2f se %.2f n %.0f %.0f\n", ratio,ser, sum, sum_a);

    sum=2*ali_cont_sup[0]+ali_cont_sup[1];
    PA=Mean_freq(&se, &re2A, 2*ali_cont_sup[0], sum);
    ratio=PA/(1-Psc); ser=ratio*sqrt(re2A+re2nosc);
    fprintf(file_id,"# P(no sup|cont no sup)/P(no sup):");
    fprintf(file_id,"\t%.2f se %.2f n %.0f %.0f\n", ratio,ser, sum, sum_a);

    //sum=2*ali_cont_sup[0]+ali_cont_sup[1];
    PA=Mean_freq(&se, &re2A, ali_cont_sup[1], sum);
    ratio=PA/Psc; ser=ratio*sqrt(re2A+re2sc);
    fprintf(file_id,"# P(sup|cont no sup)/P(sup):");
    fprintf(file_id,"\t%.2f se %.2f n %.0f %.0f\n", ratio,ser, sum, sum_a);

    sum=ali_cont_sup[1]+2*ali_cont_sup[2];
    PA=Mean_freq(&se, &re2A, 2*ali_cont_sup[2], sum);
    ratio=PA/Psc; ser=ratio*sqrt(re2A+re2sc);
    fprintf(file_id,"# P(sup|cont sup)/P(sup):");
    fprintf(file_id,"\t%.2f se %.2f n %.0f %.0f\n", ratio,ser, sum, sum_a);

    //sum=ali_cont_sup[1]+2*ali_cont_sup[2];
    PA=Mean_freq(&se, &re2A, ali_cont_sup[1], sum);
    ratio=PA/(1-Psc); ser=ratio*sqrt(re2A+re2nosc);
    fprintf(file_id,"# P(no sup|cont sup)/P(no sup):");
    fprintf(file_id,"\t%.2f se %.2f n %.0f %.0f\n", ratio,ser, sum, sum_a);

    fprintf(file_id, "### Conditional probability conserved contact | sup:\n");
    double P_cons_cont_sup[3], re2ccs[3];
    for(k=0; k<3; k++){
      P_cons_cont_sup[k]=
	Mean_freq(&se, re2ccs+k, id_cont_sup[k], ali_cont_sup[k]);
      fprintf(file_id,"# Frac_cons_cont_vs_%d_sup:\t%.4f se %.4f n %.0f\t",
	      k, P_cons_cont_sup[k], se, ali_cont_sup[k]);
      fprintf(file_id, "P(conserved cont | %d sup res)\n",k);
    }

    // Conditional Prob conserved contact | identical 
    fprintf(file_id,
	    "### Conditional probability conserved contact | id aa:\n");
    double P_cons_cont_aaid[3], re2cca[3];
    for(k=0; k<3; k++){
      P_cons_cont_aaid[k]=
	Mean_freq(&se, re2cca+k, id_cont_aaid[k], ali_cont_aaid[k]);
      fprintf(file_id,"# Frac_cons_cont_vs_%d_id:\t%.4f se %.4f n %.0f\t",
	      k, P_cons_cont_aaid[k], se, ali_cont_aaid[k]);
      fprintf(file_id, "P(conserved cont | %d identical aa)\n",k);
    }


    // Propensity conserved contacts - identical
    fprintf(file_id, "### Propensity conserved contact - identical:\n");
    ratio=P_cons_cont_aaid[2]/(P_cons_cont_aaid[0]);
    ser=ratio*sqrt(re2cca[2]+re2cca[0]);
    fprintf(file_id,"# P(cons cont|2 ide)/P(cons cont|0 ide):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_aaid[2], ali_cont_aaid[0]+ali_cont_aaid[1]);
    //
    //
    ratio=(1-P_cons_cont_aaid[1])/(1-P_cons_cont_aaid[2]);
    ser=ratio*sqrt(re2cca[2]+re2cca[0]);
    fprintf(file_id,"# P(chg cont|1 ide)/P(chg cont|2 ide):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_aaid[1], ali_cont_aaid[2]);
    //
    ratio=(1-P_cons_cont_aaid[0])/(1-P_cons_cont_aaid[2]);
    ser=ratio*sqrt(re2cca[2]+re2cca[0]);
    fprintf(file_id,"# P(chg cont|0 ide)/P(chg cont|2 ide):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_aaid[0], ali_cont_aaid[2]);
    //
    double P_2aaid_id_cont=Mean_freq(&se, &re2A, id_cont_aaid[2], id_cont);
    double P_2aaid_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_aaid[2]-id_cont_aaid[2], ali_cont-id_cont);
    ratio=P_2aaid_id_cont/P_2aaid_ch_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"# P(2 ide|cons cont)/P(2 ide|chg cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, id_cont, ali_cont-id_cont);
    //
    double P_1aaid_id_cont=Mean_freq(&se, &re2A, id_cont_aaid[1], id_cont);
    double P_1aaid_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_aaid[1]-id_cont_aaid[1], ali_cont-id_cont);
    ratio=P_1aaid_ch_cont/P_1aaid_id_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"# P(1 ide|chg cont)/P(1 ide|cons cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
    ratio, ser, ali_cont-id_cont, id_cont);
    //
    double P_0aaid_id_cont=Mean_freq(&se, &re2A, id_cont_aaid[0], id_cont);
    double P_0aaid_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_aaid[0]-id_cont_aaid[0], ali_cont-id_cont);
    ratio=P_0aaid_ch_cont/P_0aaid_id_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"# P(0 ide|chg cont)/P(0 ide|cons cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont-id_cont, id_cont);

    // Propensity conserved contacts - superimposed
    fprintf(file_id, "### Propensity conserved contact - superimposed:\n");
    ratio=P_cons_cont_sup[2]/P_cons_cont_sup[0];
    ser=ratio*sqrt(re2ccs[2]+re2ccs[0]);
    fprintf(file_id,"# P(cons cont|2 sup)/P(cons cont|0 sup):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_sup[2], ali_cont_sup[0]);
    //
    ratio=(1-P_cons_cont_sup[1])/(1-P_cons_cont_sup[2]);
    ser=ratio*sqrt(re2ccs[2]+re2ccs[1]);
    fprintf(file_id,"# P(chg cont|1 sup)/P(chg cont|2 sup):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_sup[1], ali_cont_sup[2]);
    //
    ratio=(1-P_cons_cont_sup[0])/(1-P_cons_cont_sup[2]);
    ser=ratio*sqrt(re2ccs[2]+re2ccs[0]);
    fprintf(file_id,"# P(chg cont|0 sup)/P(chg cont|2 sup):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_sup[0], ali_cont_sup[2]);
    //
    double P_2sup_id_cont=Mean_freq(&se, &re2A, id_cont_sup[2], id_cont);
    double P_2sup_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_sup[2]-id_cont_sup[2], ali_cont-id_cont);
    ratio=P_2sup_id_cont/P_2sup_ch_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"# P(2 sup|cons cont)/P(2 sup|chg cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, id_cont, ali_cont-id_cont);
    //
    double P_1sup_id_cont=Mean_freq(&se, &re2A, id_cont_sup[1], id_cont);
    double P_1sup_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_sup[1]-id_cont_sup[1], ali_cont-id_cont);
    ratio=P_1sup_ch_cont/P_1sup_id_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"# P(1 sup|chg cont)/P(1 sup|cons cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont-id_cont, id_cont);
    //
    double P_0sup_id_cont=Mean_freq(&se, &re2A, id_cont_sup[0], id_cont);
    double P_0sup_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_sup[0]-id_cont_sup[0], ali_cont-id_cont);
    ratio=P_0sup_ch_cont/P_0sup_id_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"# P(0 sup|chg cont)/P(0 sup|cons cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont-id_cont, id_cont);

    fprintf(file_id, "### As a function of sequence identity:\n");
    fprintf(file_id, "#1=f_id s.e. 3=f_ali s.e. ");
    fprintf(file_id, " 5=f_sup s.e. 7=f_cons_cont s.e.");
    fprintf(file_id, " 9=f_id_mct/f_id_fct s.e.");
    fprintf(file_id, " 11=f_sup_mct/f_sup_fct s.e.");
    fprintf(file_id, " 13=f_cons_cont_mct/f_cons_cont_fct s.e.");
    fprintf(file_id, " 15=f_neigh_noali s.e.");
    fprintf(file_id, " 17=f_id_sup/f_id*f_sup s.e.");
    fprintf(file_id, " 19=f_noid_nosup/f_noid*f_nosup s.e.");
    fprintf(file_id, " 21=f_id_cons_cont/f_id*f_cons_cont s.e.");
    fprintf(file_id, " 23=f_noid_nocons_cont/f_noid*f_nocons_cont s.e.");
    fprintf(file_id, " 25=f_sup_cons_cont/f_sup*f_cons_cont s.e.");
    fprintf(file_id, " 27=f_nosup_nocons_cont/f_nosup*f_nocons_cont s.e.");
    /*for(int k=0; k<3; k++){
      fprintf(file_id, " %d=f_cont_%dsup s.e.", 23+2*k, k);
      }*/
    fprintf(file_id, " 29=n\n");

    float P_A, P_B;
    for(j=0; j<=IBIN; j++){
      double ali=sum_ali_ct[i][0][j]+sum_ali_ct[i][1][j];
      double all=sum_all_ct[i][0][j]+sum_all_ct[i][1][j];
      if(ali<100)continue;
      double aaid=sum_aaid_ct[i][0][j]+sum_aaid_ct[i][1][j];
      double re2id=0, P_aaid=Mean_freq(&se, &re2id, aaid, ali);
      fprintf(file_id, "%.3f\t%.1g", P_aaid, se);
      ave=Mean_freq(&se, &re2, ali, all);
      fprintf(file_id, "\t%.3f\t%.1g", ave, se);
      double sup=sum_sup_ct[i][0][j]+sum_sup_ct[i][1][j];
      double re2sup, P_sup=Mean_freq(&se, &re2sup, sup, ali);
      fprintf(file_id, "\t%.3f\t%.1g", P_sup, se);
      double cons_cont=sum_id_cont_ct[i][0][j]+sum_id_cont_ct[i][1][j];
      double ali_cont=sum_ali_cont_ct[i][0][j]+sum_ali_cont_ct[i][1][j];
      double re2cont=0, P_cont=Mean_freq(&se, &re2cont, cons_cont, ali_cont);
      fprintf(file_id, "\t%.3f\t%.1g", P_cont, se);

      P_A=Mean_freq(&se, &re2A, sum_aaid_ct[i][1][j], sum_ali_ct[i][1][j]);
      P_B=Mean_freq(&se, &re2B, sum_aaid_ct[i][0][j], sum_ali_ct[i][0][j]);
      ratio=P_A/P_B; se=ratio*sqrt(re2A+re2B);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);     
      P_A=Mean_freq(&se, &re2A, sum_sup_ct[i][1][j], sum_ali_ct[i][1][j]);
      P_B=Mean_freq(&se, &re2B, sum_sup_ct[i][0][j], sum_ali_ct[i][0][j]);
      ratio=P_A/P_B; se=ratio*sqrt(re2A+re2B);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      P_A=Mean_freq(&se,&re2A,sum_id_cont_ct[i][1][j],sum_ali_cont_ct[i][1][j]);
      P_B=Mean_freq(&se,&re2B,sum_id_cont_ct[i][0][j],sum_ali_cont_ct[i][0][j]);
      ratio=P_A/P_B; se=ratio*sqrt(re2A+re2B);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      P_A=Mean_freq(&se,&re2A,sum_neigh_noali_ct[i][1][j],sum_all_ct[i][1][j]);
      P_B=Mean_freq(&se,&re2B,sum_neigh_noali_ct[i][0][j],sum_all_ct[i][0][j]);
      ratio=P_A/P_B; se=ratio*sqrt(re2A+re2B);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);

      // Propensities id-sup
      double id_sup=sum_aaid_sup_ct[i][0][j]+sum_aaid_sup_ct[i][1][j];
      double P_id_sup=Mean_freq(&se, &re2A, id_sup, sup);
      ratio=P_id_sup/P_aaid; se=ratio*sqrt(re2A+re2id);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      double noid_nosup=sum_noid_nosup_ct[i][0][j]+sum_noid_nosup_ct[i][1][j];
      double P_noid_nosup=Mean_freq(&se, &re2A, noid_nosup, ali-sup);
      ratio=P_noid_nosup/(1-P_aaid); se=ratio*sqrt(re2A+re2id);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      // Propensities id-cons cont
      ratio=P_cons_cont_aaid[2]/P_cont; se=ratio*sqrt(re2cont+re2cca[2]);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      ratio=(1-P_cons_cont_aaid[0])/(1-P_cont);
      se=ratio*sqrt(re2cont+re2cca[0]);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      // Propensities sup-cons cont
      ratio=P_cons_cont_sup[2]/P_cont; se=ratio*sqrt(re2cont+re2ccs[2]);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      ratio=(1-P_cons_cont_sup[0])/(1-P_cont);
      se=ratio*sqrt(re2cont+re2ccs[0]);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);


      /*float sum_aaid_sup=sum_aaid_sup_ct[i][0][j]+sum_aaid_sup_ct[i][1][j];
      ave=Mean_freq(&se, &re2, sum_aaid_sup, sum_aaid);
      fprintf(file_id, "\t%.3f\t%.1g", ave, se);
      float sup_noaa=sum_sup_ct[i][0][j]+sum_sup_ct[i][1][j]-sum_aaid_sup;
      float noaa=sum_ali-sum_aaid;
      ave=Mean_freq(&se, &re2, sup_noaa, noaa);
      fprintf(file_id, "\t%.3f\t%.1g", ave, se);
      for(k=0; k<3; k++){
	ave=Mean_freq(&se,&re2,sum_id_cont_sup[i][k][j],
		      sum_ali_cont_sup[i][k][j]);
	fprintf(file_id, "\t%.3f\t%.1g", ave, se);
      }*/
      fprintf(file_id, "\t%.0f\n", ali);
    }
      
  }	    

  fclose(file_id);
  for(int k=0; k<2; k++){free(name_c[k]); free(name_ct[k]);}
}

void Count_noali(int *noali, int *ali,  int *ncont, int n, float c_ave){
  noali[0]=0; noali[1]=0;
  for(int i=0; i<n; i++){
    if(ali[i]<0){if(ncont[i]<=c_ave){noali[0]++;}else{noali[1]++;}}
  }
}
double Sum_bins(double *bin){
  double sum=0; for(int j=0; j<=IBIN; j++)sum+=bin[j];
  return(sum);
}

float Ave_se(double *se, double sum1, double sum2, int n, float nind){
  double ave=sum1/n;
  *se=sqrt((sum2/n-ave*ave)/nind);
  return(ave);
}


void Sec_str_AA_propensity(struct protein *pi, struct protein *pj,
			   int *ali, int nali, float PC)
{
  int ns=4, na=22, gap=3, gap2=21, nbin=5, nmax=nbin-1;
  float lbin=0.25;
  if(Prop_secstr==NULL){
    Prop_secstr=malloc(nbin*sizeof(double **));
    Prop_AA=malloc(nbin*sizeof(double **));
    DP_bin=malloc(nbin*sizeof(double));
    DP_norm=malloc(nbin*sizeof(double));
    for(int i=0; i<nbin; i++){
      Prop_secstr[i]=Allocate_mat2_d(ns, ns);
      Prop_AA[i]=Allocate_mat2_d(na, na);
      DP_bin[i]=0; DP_norm[i]=0;
    }
    for(int i=0; i<ns; i++)SS_name[i]=malloc(10*sizeof(char));
    strcpy(SS_name[0], "Coil");
    strcpy(SS_name[1], "Helix");
    strcpy(SS_name[2], "Strand");
    strcpy(SS_name[3], "Gap");
    for(int i=0; i<20; i++)AA_name[i]=AANAME1[i];
    AA_name[20]='-'; AA_name[21]='X'; 
  }

  double qinf=Compute_qinf(pi->len,pj->len,NORM);
  double PC0=(PC0_1+PC_load[3]*qinf)/PC_norm;
  double DP=Divergence(PC, PC0, HUGE);
  int bin=DP/lbin; if(bin>nmax)bin=nmax;
  DP_bin[bin]+=DP; DP_norm[bin]++;
  double **Prop_bin=Prop_secstr[bin];
  double **Prop_AA_bin=Prop_AA[bin];

  int i, j_old=0;
  for(i=0; i<pi->len; i++){
    int iss1=SS_code(pi->ss[i]);
    int aa1=AA_code(pi->aseq[i]);
    int j=ali[i];

    if(j<0){
      Prop_bin[iss1][gap]++;
      Prop_AA_bin[aa1][gap2]++;
    }else{
      if(j<j_old){
	for(int jj=j_old; jj<j; jj++){
	  int iss2=SS_code(pj->ss[jj]);
	  int aa2=SS_code(pj->aseq[jj]);
	  Prop_bin[iss2][gap]++;
	  Prop_AA_bin[aa2][gap2]++;
	}
      }
      j_old=j;
      int iss2=SS_code(pj->ss[j]);
      int aa2=SS_code(pj->aseq[j]);
      Prop_bin[iss1][iss2]++;
      Prop_AA_bin[aa2][gap2]++;
    }
  }
}

int SS_code(char ss){
  if(ss=='H'){return(1);}
  if(ss=='E'){return(2);}
  return(0);
}

int AA_code(char aa){
  for(int i=0; i<20; i++)if(aa==AANAME1[i])return(i);
  if(aa=='-')return(21);
  return(20);
}

void Print_propensities(char *nameout)
{
  int ns=4, nbin=5; // na=22, gap=3, gap2=21, nmax=nbin-1;
  //float lbin=0.25;
  FILE *file_out=fopen(nameout, "w");
  printf("Writing Sec.str. and AA propensities of aligned sites in %s\n",
	 nameout);
  char name_ss[200], tmp[20]; int a;
  sprintf(name_ss, "%s", SS_name[0]);
  for(a=1; a<ns; a++){sprintf(tmp, "\t%s", SS_name[a]); strcat(name_ss, tmp);}
  char name_AA[200];
  sprintf(name_AA, "%c", AA_name[0]);
  for(a=1; a<20; a++){sprintf(tmp, "\t%c", AA_name[a]); strcat(name_AA, tmp);}
  strcat(name_AA, "\t-\n");

  for(int i=0; i<nbin; i++){
    fprintf(file_out, "# Bin %d: <PC_div>= %.3g norm= %.3g\n",
	    i+1, DP_bin[i]/DP_norm[i], DP_norm[i]);
    if(DP_norm[i]<4)continue;
    int n=ns, a, b; double **Prop=Prop_secstr[i], F[n];
    // Symmetrize
    for(a=0; a<n; a++){ 
      for(b=0; b<a; b++){
	Prop[a][b]+=Prop[b][a];
	Prop[a][b]/=2;
	Prop[b][a]=Prop[a][b];
      }
    }
    int zero=0; double norm=0;
    for(a=0; a<n; a++){
      F[a]=0; for(b=0; b<n; b++)F[a]+=Prop[a][b];
      norm+=F[a]; if(F[a]==0)zero=1;
    }

    float Low=-log(10);
    fprintf(file_out, "# Secondary structure: Frequency\n");
    fprintf(file_out, "# %s\n", name_ss);
    fprintf(file_out, "# %.3f", F[0]/norm);
    for(a=1; a<n; a++)fprintf(file_out, "\t%.3f", F[a]/norm);
    fprintf(file_out, "\n");
    if(zero)continue;
    fprintf(file_out, "# Secondary structure: Propensities\n");
    fprintf(file_out, "#SS\t%s\n", name_ss);
    for(a=0; a<n; a++){
      fprintf(file_out, "%s", SS_name[a]);
      for(b=0; b<n; b++){
	double p=norm*Prop[a][b]/(F[a]*F[b]);
	if(p>0){p=log(p);}else{p=Low;}
	fprintf(file_out, "\t%.4g", p);
      }
      fprintf(file_out, "\n");
    }
  }
  fclose(file_out);
}

char **Assign_Seq(char ***name_seq, int *len_seq,
		  struct Prot_input *Prot_in, int N_seq,
		  int *rep_str, int *i_seq, int N_ali)
{
  *name_seq=malloc(N_seq*sizeof(char *));
  char **Seq=malloc(N_seq*sizeof(char *));
  for(int i=0; i<N_seq; i++){
    struct Prot_input *proti=Prot_in+i_seq[rep_str[i]];
    int Li=proti->len, j;
    len_seq[i]=Li;
    Seq[i]=malloc(Li* sizeof(char)); int k=0;
    for(j=0; j<N_ali; j++){
      if(proti->seq[j]!='-'){Seq[i][k]=proti->seq[j]; k++;}
    }
    if(k!=Li){
      printf("ERROR in Assign_seq, L= %d %d\n", Li, k); exit(8);
    }
    (*name_seq)[i]=malloc(80*sizeof(char));
    sprintf((*name_seq)[i], "%s%c", proti->name, proti->chain);
  }
  return(Seq);
}

void Write_ali(int ***Ali_pair, int i, int j, int *L_seq, int *ali_PC)
{
  int *ali_ij=Ali_pair[i][j], *ali_ji=Ali_pair[j][i], *ali=ali_PC, u;
  int error=0;
  for(u=0; u<L_seq[j]; u++){ali_ji[u]=-1;}
  for(u=0; u<L_seq[i]; u++){
    ali_ij[u]=*ali;
    if(*ali>=L_seq[j]){
      ali_ij[u]=-1;
      if(error)continue;
      printf("ERROR in alignment %d %d, res %d > %d\n",
	     i, j, *ali+1, L_seq[j]);
      for(int v=0; v<L_seq[i]; v++)printf("%d ",ali_PC[v]);
      printf("\n"); error=1; //exit(8);
    }else if(*ali>=0){
      ali_ji[*ali]=u;
      ali_ij[u]=*ali;
    }
    ali++;
  }
}

int ***Select_alis(int ***Ali_pair, int N_seq, int *len_seq, int *rep_str)
{
  int ***Ali_pair_seq=malloc(N_seq*sizeof(int **));
  for(int i=0; i<N_seq; i++){
    Ali_pair_seq[i]=malloc(N_seq*sizeof(int *));
    int Li=len_seq[i], i_str=rep_str[i];
    for(int j=0; j<N_seq; j++){
      Ali_pair_seq[i][j]=malloc(Li*sizeof(int));
      if(j==i)continue;
      int *ali_sel=Ali_pair_seq[i][j];
      int *ali=Ali_pair[i_str][rep_str[j]];
      for(int k=0; k<Li; k++)ali_sel[k]=ali[k];
    }
  }
  return(Ali_pair_seq);
}
