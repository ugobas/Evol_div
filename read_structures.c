//#include "mammoth.h"
#include "protein.h"
#include "cont_list.h"
#include "Contact_divergence_aux.h"
#include "read_structures.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#define T_SEC "C"
//#define S_SEC "C"
//#define G_SEC "H"
#define LEN_SS_THR 4  /* Minimal length of a secondary structure element */
#define N_ELE_MAX 1000

/********************************/

int Read_pdb_names(struct protein **prots, char *PDB_PATH, char *file_name);
int Get_pdbid(char *pdbid, char *file_name, char *chain);
int Remove_char(char *string, char c);

int Read_processed_proteins(struct protein *prot, int N_prot,
			    char *file_cm, char *file_seq, int IJ_MIN);
char Convert_ss(char ss, char ss_old);
int Set_ele(char ss, int *len, int *n_sec_el, short *ss_num, int i,
	    struct sec_el *sec_el);
int Read_PDB_list(struct protein *prots, int N_pdb, char *PDB_PATH,
		  int CM_FILE, char CONT_TYPE, float CONT_THR, int IJ_MIN);
int Print_proteins(struct protein *prots, int N_pdb,
		   char *file_cm, char *file_seq,
		   char CONT_TYPE, float CONT_THR, int IJ_MIN);

#define NCHAR 200
char FILE_CM[NCHAR], FILE_SEQ[NCHAR];

/***********************************************************
             READING STRUCTURES
************************************************************/
int Read_PDB_compress(struct protein *prot,
		      char *Prot_name, char *chain, char *PDB_PATH)
{
  char FILE_TMP[100]="pdb.tmp"; int L;
  char filename[200];
  sprintf(filename,"%s%s", PDB_PATH, Prot_name);
  if(Get_compression(prot->name_file)==0){
    L=Read_pdb(prot, filename, chain);
  }else{
    char command[400];
    sprintf(command,"gunzip -c %s > %s\n", filename, FILE_TMP);
    system(command);
    L=Read_pdb(prot, FILE_TMP, chain);
    sprintf(command,"rm -f %s\n", FILE_TMP);
    system(command);
  }
  //printf("%d residues\n", L);
  return(L);
}

int Read_structures(struct protein **prots, char *PDB_PATH, char *FILE_LIST,
		    char *EXT_CM, char *EXT_SEQ, int CONT_DEF,
		    char CONT_TYPE, float CONT_THR, int IJ_MIN)
{
  int CM_FILE=0;
  int N_pdb=Read_pdb_names(prots, PDB_PATH, FILE_LIST);
  if(N_pdb<=0){
    printf("WARNING, no protein read in %s\n", FILE_LIST); return(0);
  }

  // If default contact definition, look for processed contact matrices
  // in FILE_CM and FILE_SEQ. If not found, CM_FILE=0
  if(CONT_DEF){
    Change_ext(FILE_CM, FILE_LIST, EXT_CM);
    Change_ext(FILE_SEQ,FILE_LIST, EXT_SEQ);
    CM_FILE=Read_processed_proteins(*prots, N_pdb, FILE_CM, FILE_SEQ,IJ_MIN);
  }

  // Read list of PDB files
  Read_PDB_list(*prots, N_pdb, PDB_PATH, CM_FILE, CONT_TYPE, CONT_THR, IJ_MIN);
  if((CM_FILE==0)&&(CONT_DEF)){
    Print_proteins(*prots, N_pdb, FILE_CM, FILE_SEQ,CONT_TYPE,CONT_THR,IJ_MIN);
  }
  return(N_pdb);
}

int Read_pdb_names(struct protein **prots, char *PDB_PATH, char *file_name)
{
  FILE *file_in=fopen(file_name, "r");
  int n=0, i; char string[5000], chain[9];

  PDB_PATH[0]='\0'; *prots=NULL;
  if(file_in==NULL){
    printf("WARNING, file %s does not exist", file_name); return(0);
  }

  // Count
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string, "DIR=", 4)==0 ||
       strncmp(string, "MAMMOTH",7)==0)continue;
//     if((string[0]!='#')&&(string[0]!='/')&&(string[0]!='\n'))n++;
    if((string[0]!='#')&&(string[0]!='\n'))n++;
  }
  fclose(file_in);
  if(n==0)return(0);

  // Allocate and read
  *prots=malloc(n*sizeof(struct protein));
  struct protein *prot=*prots;
  file_in=fopen(file_name, "r"); n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if((string[0]=='#')||
    (string[0]=='\n')  ||
    strncmp(string, "MAMMOTH",7)==0)continue;
    if(strncmp(string, "DIR=", 4)==0){
      for(i=0; i<4; i++)string[i]=' ';
      strcpy(PDB_PATH, string); continue;
    }
    chain[0]='\0';
    sscanf(string, "%s %s", prot->name_file, chain);
    Get_pdbid(prot->name, prot->name_file, chain);
    prot->chain=chain[0];
    printf("file=%s PDBID= %s chain=%c\n",
	   prot->name_file, prot->name, prot->chain);
    n++; prot++;
  }
  fclose(file_in);
  Remove_char(PDB_PATH, '\n');
  printf("DIR: %s\n", PDB_PATH);
  return(n);
}

int Remove_char(char *string, char c){
  char *ptr=string;
  while(*ptr!='\0'){if(*ptr==c)*ptr='\0'; ptr++;}
  return(0);
}


int Read_PDB_list(struct protein *prots, int N_pdb, char *PDB_PATH,
		  int CM_FILE, char CONT_TYPE, float CONT_THR, int IJ_MIN)
{
  int i; char FILE_TMP[100];
  strcpy(FILE_TMP, "pdb.tmp");
  for(i=0; i<N_pdb; i++){
    char filename[200], command[400], chain=prots[i].chain;
    sprintf(filename,"%s%s", PDB_PATH, prots[i].name_file);
    if(Get_compression(prots[i].name_file)==0){
      Read_pdb(prots+i, filename, &chain);
    }else{
      sprintf(command,"gunzip -c %s > %s\n", filename, FILE_TMP);
      system(command);
      Read_pdb(prots+i, FILE_TMP, &chain);
    }
    if(CM_FILE==0)
      Compute_contact_list(prots+i, CONT_TYPE, CONT_THR, IJ_MIN);
    if(i&&(((int)(0.001*i)*1000)==i))printf("%4d PDB files read\n", i);
  }
  if(Get_compression(prots[N_pdb-1].name_file)){
    char command[300];
    sprintf(command,"rm -f %s\n", FILE_TMP); system(command);
  }
  printf("%4d PDB files read\n", i);
  return(0);
}

int Print_proteins(struct protein *prots, int N_prot,
		   char *file_cm, char *file_seq,
		   char CONT_TYPE, float CONT_THR, int IJ_MIN)
{
  FILE *file_out=fopen(file_seq, "w"); int i, j;
  char CONT_STRING[40];

  // Print sequence
  for(i=0; i<N_prot; i++){
    if((prots[i].len==0)||(prots[i].Cont_map==NULL))continue;
    fprintf(file_out, "# %s %d\n", prots[i].name, prots[i].len);
    char *seq=prots[i].aseq, *ss=prots[i].ss;
    char **pdbres=prots[i].pdbres;
    int *ncont=prots[i].ncont;
    for(j=0; j<prots[i].len; j++){
      fprintf(file_out, "%c %5s %c %2d\n", seq[j],
	      pdbres[j], ss[j], ncont[j]);
    }
  }
  fclose(file_out);

  // Print contact list
  file_out=fopen(file_cm, "w");
  if(CONT_TYPE=='a'){strcpy(CONT_STRING, "Alpha");}
  else if(CONT_TYPE=='b'){strcpy(CONT_STRING, "Beta");}
  else if(CONT_TYPE=='c'){strcpy(CONT_STRING, "All atoms");}
  fprintf(file_out,"# %s contacts. |i-j| >= %d. Threshold distance: %.3f A\n",
	  CONT_STRING, IJ_MIN, CONT_THR);
  for(i=0; i<N_prot; i++){
    if((prots[i].len==0)||(prots[i].Cont_map==NULL))continue;
    fprintf(file_out, "# %d %d  %s %c\n",
	    prots[i].len, prots[i].N_cont, prots[i].name, prots[i].exp_meth);
    short **Cont_list=prots[i].Cont_map;
    for(j=0; j<prots[i].len; j++){
      short *res=Cont_list[j];
      while(*res>=0){
	fprintf(file_out, "%d %d\n", j, *res); res++;
      }
    }
  }
  fclose(file_out);
  return(0);
}

int Get_pdbid(char *pdbid, char *file_name, char *chain){ //, char chain
  char *ptr=file_name; int i=0;
  while((*ptr!='\0')&&(*ptr!='\n')){
    if(*ptr=='.')break;
    if(*ptr=='/'){i=0; ptr++;}
    pdbid[i]=*ptr; i++; ptr++;
  }
  if(chain[0]!='\0'){pdbid[i]=chain[0];}
  else{pdbid[i]='\0';}
  pdbid[i+1]='\0';
  return(0);
}


int Read_processed_proteins(struct protein *prot, int N_prot,
			    char *file_cm, char *file_seq, int IJ_MIN)
{
  FILE *file_in;
  int ip, i, L, Nc, nc, nc_max=0, Np=0, i1, i2;
  short **Cont_map; int *n_cont;
  char string[1000], dum[40], dum1[40], dum2[40], name2[100], exp_meth[4];
  char ss, ss_old;
  struct sec_el sec_el[N_ELE_MAX], *sec_el_2;
  char **names;
  
  // Initialize
  for(ip=0; ip<N_prot; ip++){
    prot[ip].len=0;
    prot[ip].N_cont=0;
    prot[ip].ss=NULL;
    prot[ip].ss_num=NULL;
    prot[ip].aseq=NULL;
    prot[ip].ncont=NULL;
    prot[ip].Cont_map=NULL;
  }

  // Test if files exist
  file_in=fopen(file_cm, "r");
  if(file_in==NULL)return(0);
  fclose(file_in);
  file_in=fopen(file_seq, "r");
  if(file_in==NULL)return(0);

  names=malloc(N_prot*sizeof(char *));
  for(ip=0; ip<N_prot; ip++){
    names[ip]=malloc(200*sizeof(char));
    strcpy(names[ip], prot[ip].name);
  }

  /* Read sequences and secondary structures */
  printf("Reading sequences and sec str in %s\n", file_seq);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    sscanf(string, "%s%s%d%s", dum, name2, &L, exp_meth);
    ip=Find_name(name2, names, N_prot, Np);
    if(ip>=0){
      struct protein *P=prot+ip; Np++;
      P->len=L; 
      P->aseq=malloc(L*sizeof(char));
      P->ss=malloc(L*sizeof(char));
      P->exp_meth=exp_meth[0];
      ss_old=' ';
      for(i=0; i<L; i++){
	fgets(string, sizeof(string), file_in);
	if(string[0]=='*')string[0]='X';
	P->aseq[i]=string[0];
	P->ss[i]= Convert_ss(string[8], ss_old);
	ss_old=string[8]; string[8]='A';
	sscanf(string, "%s%s%s%d", dum, dum1, dum2, &nc);
	if(nc>nc_max)nc_max=nc;
      }
    }else{
      printf("Protein %s not in list\n", name2);
      for(i=0; i<L; i++)fgets(string, sizeof(string), file_in);
    }
  }
  fclose(file_in);
  printf("%d proteins found\n", Np);

  /* Secondary structure elements */
  for(ip=0; ip<N_prot; ip++){
    int n_sec_el=0, len=0; char sec_old=' ';
    if(prot[ip].len==0)continue;
    prot[ip].ss_num=malloc(prot[ip].len*sizeof(short));
    for(i=0; i<prot[ip].len; i++){
      ss=prot[ip].ss[i];
      if(ss!=sec_old){
	if(sec_old!=' ')	// End element
	  Set_ele(sec_old, &len, &n_sec_el, prot[ip].ss_num, i, sec_el);
	sec_old=ss;
      }
      if(ss!=' '){len++; prot[ip].ss_num[i]=n_sec_el;}
      else{prot[ip].ss_num[i]=-1;}
    }
    // End protein, write sec str elements
    if(len)Set_ele(sec_old, &len, &n_sec_el, prot[ip].ss_num, i, sec_el);
    sec_el_2=malloc(n_sec_el*sizeof(struct sec_el));
    prot[ip].sec_el=sec_el_2; prot[ip].n_sec_el=n_sec_el;
    for(i=0; i<n_sec_el; i++)sec_el_2[i]=sec_el[i];
  }

  // Test if all proteins have been found
  // If not, all proteins are read from PDB files
  if(Np!=N_prot){
    printf("WARNING, only %d proteins instead of %d found in %s\n",
	   Np, N_prot, file_seq);
    printf("Missing proteins:\n");
    for(ip=0; ip<N_prot; ip++){
      if(prot[ip].len==0)printf("%s\n", names[ip]);
    }
    printf("Reading coordinates and contact matrices from PDB files\n");
    for(ip=0; ip<N_prot; ip++)free(names[ip]);
    free(names);
    return(0);
  }

  /* Read contact matrices */
  printf("Maximum num. of contacts per residue: %d\n", nc_max);
  nc_max++; //if(nc_max>20)exit(8);
  printf("Reading contact matrices in %s\n", file_cm);
  file_in=fopen(file_cm, "r"); Np=0;
  fgets(string, sizeof(string), file_in);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    sscanf(string, "%s%d%d%s", dum, &L, &Nc, name2);
    ip=Find_name(name2, names, N_prot, Np);
    if(ip>=0){
      prot[ip].N_cont=Nc; Np++;
      //Cont_map=Allocate_mat2_i(L, nc_max+1);
      Cont_map=malloc(L*sizeof(short *));
      prot[ip].Cont_map=Cont_map;
      n_cont=malloc(L*sizeof(int));
      prot[ip].ncont=n_cont;
      for(i=0; i<L; i++){
	n_cont[i]=0;
	prot[ip].Cont_map[i]=malloc((nc_max)*sizeof(short));
      }
      for(i=0; i<Nc; i++){
	fgets(string, sizeof(string), file_in);
	sscanf(string, "%d %d", &i1, &i2);
	if((i2-i1)<IJ_MIN)continue;
	Cont_map[i1][n_cont[i1]]=i2; n_cont[i1]++;
      }
      for(i=0; i<L; i++)Cont_map[i][n_cont[i]]=-1;
    }else{
      for(i=0; i<Nc; i++)fgets(string, sizeof(string), file_in);
    }
  }
  fclose(file_in);


  // Test if all proteins have been found
  if(Np!=N_prot){
    printf("WARNING, only %d proteins instead of %d found in %s\n",
	   Np, N_prot, file_cm);
    printf("Missing proteins:\n");
    for(ip=0; ip<N_prot; ip++){
      if(prot[ip].N_cont==0)printf("%s\n", names[ip]);
    }
    printf("Reading coordinates and contact matrices from PDB files\n");
    for(ip=0; ip<N_prot; ip++)free(names[ip]);
    free(names);
    return(0);
  }

  printf("End reading proteins\n");
  for(ip=0; ip<N_prot; ip++)free(names[ip]);
  free(names);
  return(1);
}

char Convert_ss(char ss, char ss_old){

  if((ss=='T')||(ss=='G')){
    return(' ');
  }else if((ss!=' ')&&(ss!='H')&&(ss!='E')){
    return(ss_old);
    /*if(ss=='S'){
      return(S_SEC[0]);
      }else if(ss=='T'){
      return(T_SEC[0]);
      }else if(ss=='G'){
      return(G_SEC[0]);
      }else if((ss!='H')&&(ss!='E')){
      return(' ');*/
  }else{
    return(ss);
  }
}

int Set_ele(char ss, int *len, int *n_sec_el, short *ss_num, int i,
	    struct sec_el *sec_el)
{
  int j; struct sec_el *sec;
  if(*len<LEN_SS_THR){
    for(j=i-(*len); j<i; j++)ss_num[j]=-1;
    return(0);
  }
  sec=sec_el+(*n_sec_el);
  sec->len=*len; sec->type=ss; sec->ini=i-(*len);
  (*len)=0; (*n_sec_el)++;
  return(0);
}
