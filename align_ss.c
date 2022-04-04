#include "protein.h"
#include "align_ss.h"
#include "Contact_divergence_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>

int *warn_notali=NULL;

static int SS_match(int **Ali, char **ss_prot, int m, int i,
		    int ini_res, int end_res, int move_res);

int Set_sec_str(struct protein *prots, int N_pdb)
{
  // Determine number of secondary structures
  float ss_ave=0;
  for(int i=0; i<N_pdb; i++){
    if(prots[i].ss==NULL){
      printf("ERROR, secondary structure of %s (i=%d) not present\n",
	     prots[i].name, i); return(-1);
    }
    char *ss=prots[i].ss; int nss=0; char type='-';
    for(int j=0; j<prots[i].len; j++){
      if(ss[j]=='-'){ss[j]='c';}
      else if(ss[j]!=type){nss++; type=ss[j];}
    }
    ss_ave+=nss;
  }
  printf("Mean number of secondary structure elements: %.1f\n",ss_ave/N_pdb);
  return(0);
}


void Write_ali_ss(int **Prot_ali_ss, struct protein *prots,
		  int N_pdb, int N_ali, char *name_in)
{
  char f_ali_ss[100]; sprintf(f_ali_ss,"%s_ss_ali.fasta",name_in);
  printf("Writing %s\n", f_ali_ss);
  FILE *f_ali=fopen(f_ali_ss,"w");
  for(int i=0; i<N_pdb; i++){
    fprintf(f_ali,">%s\n",prots[i].name);
    for(int j=0; j<N_ali; j++){
      if(Prot_ali_ss[i][j]<0){
	fprintf(f_ali,"-");
      }else{
	fprintf(f_ali,"%c",prots[i].aseq[Prot_ali_ss[i][j]]);
      }
    }
    fprintf(f_ali,"\n");
  }
  fclose(f_ali);
}

void Write_ss_ali(int **ali, struct protein *prots, int N_pdb, int N_ali,
		  char *name, char *what)
{
  char f_ali_ss[100]; sprintf(f_ali_ss,"%s_%s_ss.fasta",name,what);
  printf("Writing %s\n", f_ali_ss);
  FILE *f_ali=fopen(f_ali_ss,"w");
  for(int i=0; i<N_pdb; i++){
    fprintf(f_ali,">%s\n",prots[i].name);
    for(int j=0; j<N_ali; j++){
      if(ali[i][j]<0){
	fprintf(f_ali,"-");
      }else{
	char c=prots[i].ss[ali[i][j]]; if(c=='-')c='c';
	fprintf(f_ali,"%c",c);
      }
    }
    fprintf(f_ali,"\n");
  }
  fclose(f_ali);

}

int Align_ss(int **Ali,
	     int **Ali_store, int N_ali, struct protein **prots, int m,
	     int SHIFT_MAX, int verbose)
{
  // Store Alignments
  int all_gap[N_ali], i, j;
  for(j=0; j<N_ali; j++){
    all_gap[j]=1;
    for(i=0; i<m; i++){
      if(Ali_store[i][j]>=0){all_gap[j]=0; break;}
    }
  }
  int nali=0, *column[m], L[m];
  char *ss_prot[m];
  for(i=0; i<m; i++){
    int *ali=Ali[i], *col=malloc(N_ali*sizeof(int));
    ss_prot[i]=(prots[i])->ss; L[i]=(prots[i])->len;
    column[i]=col; for(j=0; j<L[i]; j++){col[j]=-1;}
    nali=0; int *al=Ali_store[i];
    for(j=0; j<N_ali; j++){
      if(all_gap[j])continue;
      ali[nali]=al[j]; if(ali[nali]>=0)col[ali[nali]]=nali; nali++;
    }
    for(j=nali; j<N_ali; j++){ali[j]=-1;}
  }

  
  /* For every column in MSA, test whether there is a gap in a seq and a residue in the core of a secondary str. on another one. In this case, determine the initial and final column of the secondary structure. Compute the score of moving the gaps in these columns that do not coincide with the tails of the secondary stucture either to the initial or to the final column. Perform the change with minimum score. Move to the final column of the secondary structure. */

  int Test_gap[m]; for(i=0; i<m; i++)Test_gap[i]=1;
  for(int icol=0; icol<nali; icol++){
    // Is there a gap?
    int gap=0;
    for(i=0; i<m; i++){
      if(Ali[i][icol]<0){if(Test_gap[i])gap=1;}
      else{Test_gap[i]=1;}
    }
    if(gap==0)continue;
    // Yes, there is a gap. Is it inside a secondary structure?
    int ini_ss=-1, end_ss=-1;
    int ini_ssi[m], end_ssi[m];
    for(i=0; i<m; i++){
      ini_ssi[i]=-1; end_ssi[i]=-1;
      int iss=Ali[i][icol]; j=icol;
      if(iss<0){
	while(j<nali){if(Ali[i][j]>=0){iss=Ali[i][j]; break;} j++;}
      }
      char *ss=ss_prot[i]+iss; 
      if((*ss!='H')&&(*ss!='E'))continue; // It is a loop
      if((j>icol)&&(*ss!=*(ss-1)))continue; // leave if gap and no ss.
      int *col=column[i], iss2=iss; char *ss2=ss;
      ini_ssi[i]=col[iss];
      while((iss2>=0)&&(*ss2==*ss)){
	if((col[iss2]<ini_ssi[i])&&(col[iss2]>=0))ini_ssi[i]=col[iss2];
	iss2--; ss2--;
      }
      iss2=iss; ss2=ss;
      end_ssi[i]=col[iss];
      while((iss2<L[i])&&(*ss2==*ss)){
	if(col[iss2]>end_ssi[i])end_ssi[i]=col[iss2];
	iss2++; ss2++;
      }
      if((ini_ss<0)||(ini_ssi[i]<ini_ss))ini_ss=ini_ssi[i];
      if((end_ss<0)||(end_ssi[i]>end_ss))end_ss=end_ssi[i];
    }
    if((ini_ss<0)||(end_ss<0)||(ini_ss==icol)||(end_ss==icol)){
      // Accept gap at extreme of secondary structure element
      for(i=0; i<m; i++){if(Ali[i][icol]<0)Test_gap[i]=0;}
      continue;
    }
    int last_col=end_ss;
    for(i=0; i<m; i++){
      if(Test_gap[i]==0)continue;
      // Gap inside secondary structure, move toward ini_ss or end_ss
      // Identify all gaps and initial and final residues
      // New: decide extremes of secondary structure excluding target seq
      if(0){ini_ss=-1; end_ss=-1;
	for(int j=0; j<m; j++){
	  if(j==i)continue;
	  if((ini_ss<0)||(ini_ssi[j]<ini_ss))ini_ss=ini_ssi[j];
	  if((end_ss<0)||(end_ssi[j]>end_ss))end_ss=end_ssi[j];
	}
	if(ini_ss<0)continue;
      }

      int *col=column[i], *ali=Ali[i], j, jj;
      int ngap=0, open=0, inigap[nali], endgap[nali];
      int ini_res=-1, end_res=-1, next_ss=-1;
      char *ss_ini=NULL, *ss=ss_prot[i];
      for(j=icol; j<=end_ss; j++){
	if(ali[j]<0){
	  if(open==0){inigap[ngap]=j; open=1;}
	  endgap[ngap]=j;
	}else{
	  if(open){open=0; if(inigap[ngap]!=ini_ss)ngap++;}
	  if(ini_res<0)ini_res=j;
	  char *ssj=ss+ali[j];
	  if(ss_ini==NULL){ss_ini=ssj;}
	  else if((*ss_ini!='H')&&(*ss_ini!='E')){ss_ini=ssj;}
	  //else if((*ssj!=*ss_ini)&&((*ssj=='H')||(*ssj=='E')))
	  //{next_ss=j;break;}
	  end_res=j;
	}
      }
      if(open){open=0; if(endgap[ngap]!=end_ss)ngap++;}
      if((next_ss>=0)&&(next_ss<last_col))last_col=next_ss;
      // Backward
      if(Ali[i][icol]<0){
	j=icol-1; while(Ali[i][j]<0){inigap[0]=j; j--;}
      }
      for(j=icol-1; j>=ini_ss; j--){
	if(ali[j]>=0){
	  char *ssj=ss+ali[j];
	  if(ss_ini==NULL){ss_ini=ssj;}
	  else if((*ss_ini!='H')&&(*ss_ini!='E')){ss_ini=ssj;}
	  //else if((*ssj!=*ss_ini)&&((*ssj=='H')||(*ssj=='E'))){break;}
	  ini_res=j;
	}
      }
      if((ngap==0)||(ini_res<0)){Test_gap[i]=0; continue;}
 
      // Sequence i has gaps inside secondary structure
      // Determine shifts to the left and to the right
      int lgap[ngap];
      for(j=0; j<ngap; j++)lgap[j]=endgap[j]-inigap[j]+1;

      int store_l=ini_res, store_r=end_res; char dir[10];
 
      // recursively move the gap with lowest shift
      int todo[ngap]; for(j=0; j<ngap; j++)todo[j]=1;
      for(jj=0; jj<ngap; jj++){

	int inimin=1000, igap1=-1, endmax=-1, igap2=-1;
	for(j=0; j<ngap; j++){
	  if(todo[j]==0)continue;
	  if(inigap[j]<inimin){inimin=inigap[j]; igap1=j;}
	  if(endgap[j]>endmax){endmax=endgap[j]; igap2=j;}
	}
	int igap=-1, res1=-1, res2=-1;
	int shift_min=1000, ss_match_max=-10, store_res, gap1;
	if(igap1>=0){
	  for(int store=store_l; store>=0; store--){
	    if(ali[store]<0)break;
	    int store_tmp=store+lgap[igap1]; // Where residues are moved
	    int ss_match=
	      SS_match(Ali,ss_prot,m,i,store,inigap[igap1]-1,store_tmp);
	    if((igap)&&(ss_match<ss_match_max))break;
	    if((igap<0)||(ss_match>ss_match_max)){
	      igap=igap1; strcpy(dir,"left");
	      ss_match_max=ss_match;
	      store_res=store_tmp;
	      gap1=store; // Where gap is moved
	      res1=ali[store];
	      res2=ali[inigap[igap]-1];
	      shift_min=inigap[igap]-store;
	    }
	  }
	}
	if(igap2>=0){
	  int store_tmp=inigap[igap2], match_max=-1;
	  for(int store=store_r; store<nali; store++){
	    if(ali[store]<0)break;
	    int ss_match=
	      SS_match(Ali,ss_prot,m,i,endgap[igap2]+1,store,store_tmp);
	    //if((igap<0)||(shift<shift_min)){
	    if((match_max>=0)&&(ss_match<=match_max))break;
	    match_max=ss_match;
	    if((ss_match >= ss_match_max)||(igap<0)){
	      igap=igap2; strcpy(dir,"right");
	      ss_match_max=ss_match;    
	      shift_min=store-endgap[igap2];
	      gap1=store-lgap[igap]+1;
	      store_res=store_tmp;
	      res1=ali[endgap[igap]+1];
	      res2=ali[store];
	    }
	  }
	}
	if(igap<0 || shift_min<=0)continue;
	if(shift_min >SHIFT_MAX){
	  if(verbose)
	    printf("WARNING, icol %d gap %d of %d shift=%d > %d, keeping ali\n",
		   icol,jj+1, ngap, shift_min, SHIFT_MAX);
	  //printf("gap: %d-%d store_l= %d store_r=%d\n",
	  //	 inigap[igap], endgap[igap], store_l, store_r);
	  break;
	}
	todo[igap]=0;

	if(res1<0 || res2<0 || (res2<res1)){
	  printf("WARNING gap %d to %s, res1=%d res2=%d exiting\n",
		 jj+1,dir,res1,res2);
	  for(j=0; j<ngap; j++){
	    printf("gap %d: %d-%d\n",j,inigap[j],endgap[j]);
	  }
	  continue;
	}

	// The move is accepted
	if(strcmp(dir,"left")==0){
	  // Shift to the left
	  store_l=store_res;
	}else{
	  store_r=gap1-1;
	}
	// Move residues and gap
	for(int k=res1; k<=res2; k++){
	  col[k]=store_res; ali[store_res]=k; store_res++;
	}
	for(j=gap1; j<(gap1+lgap[igap]); j++)ali[j]=-1;
	if(verbose){
	  printf("icol %d moving gap %d l=%d from col %d to %d shift= %d\n",
		 icol,jj+1,lgap[igap],inigap[igap],gap1,shift_min);
	}	

	// Test and print

	/*if(shift_r[igap]<shift_l[igap]){
	  // Test (right)
	  k=store_r-lgap[igap]; int ak=ali[k];
	  if((ak>=0)&&(ak<(nali-1))&&(col[ak]>=col[ak+1])){
	    printf("ERROR moving to %s (end), col[%d]=%d col[%d]=%d\n",
		   dir, ak, col[ak], ak+1, col[ak+1]); exit(8);
	  }
	  k=inigap[igap]; ak=ali[k];
	  if((ak>0)&&(col[ak]<=col[ak-1])){
	    printf("ERROR moving to %s (ini), col[%d]=%d col[%d]=%d\n",
		   dir, ak, col[ak], ak-1, col[ak-1]); exit(8);
	  }
	}else{
	  // Test left
	  k=store_l+lgap[igap]; int ak=ali[k];
	  if((ak>0)&&(col[ak]<=col[ak-1])){
	    printf("ERROR moving to %s (ini), col[%d]=%d col[%d]=%d\n",
		   dir, ak, col[ak], ak-1, col[ak-1]); exit(8);
	  }
	  k=endgap[igap]; ak=ali[k];
	  if((ak>=0)&&(ak<(nali-1))&&(col[ak]>=col[ak+1])){
	    printf("ERROR moving to %s (end), col[%d]=%d col[%d]=%d\n",
		   dir, ak, col[ak], ak+1, col[ak+1]); exit(8);
	  }
	  }*/
	int ires=ali[ini_res];
	if((ires>0) && (col[ires]<=col[ires-1])){
	  printf("ERROR moving gaps to %s  col[%d]=%d col[%d]=%d\n",
		 dir, ires, col[ires], ires-1, col[ires-1]);
	  printf("sequence= %d gap %d of %d l=%d shift= %d\n",
		 i, jj, ngap, lgap[igap], shift_min);
	  printf("col %d columns of sec.str: %d-%d sec.str. residues: %d-%d\n",
		 icol, ini_ss, end_ss, ires, ali[end_res]);
	  printf("Columns of residues: ");
	  for(j=ires; j<=ali[end_res]; j++)printf(" %d", col[j]);
	  printf("\n");
	  if(ini_res)printf("Previous column of residues: %d\n",col[ires-1]);
	  exit(8);
	}
      } // end gaps
    } // end sequences
    icol=last_col; // All columns until last_col have been changed
  }

 for(i=0; i<m; i++){
    int *col=column[i];
    for(j=0; j<L[i]; j++){
      if((col[j]>=nali)){
	printf("ERROR, sequence %d res %d column %d not in [0,%d]\n",
	       i, j, col[j], nali-1); exit(8);
      }else if(j && (col[j]<=col[j-1])&&(col[j-1]>=0)&&(col[j]>=0)){
	printf("ERROR in seq %d, column[%d]=%d<=column[%d]=%d\n",
	       i, j, col[j], j-1, col[j-1]); exit(8);
      }
    }
  }


 // Free memory
 for(i=0; i<m; i++){free(column[i]);}
 //exit(8);

 return(nali);
}

int SS_match(int **Ali, char **ss_prot, int m, int i,
	     int ini_res, int end_res, int new_res)
{
  // New match of fragment of prot i in columns ini_res-end_res minus
  // Previous match 
  int ss_match=0; int new=new_res;
  int *Ali_i=Ali[i]; char *ss_i=ss_prot[i];
  for(int j=ini_res; j<=end_res; j++){
    int res=Ali_i[j]; if(res<0)continue;
    if(ss_i[res]=='c')continue;
    for(int k=0; k<m; k++){
      if(k==i)continue;
      if((Ali[k][j]>=0)&&(ss_prot[k][Ali[k][j]]==ss_i[res]))ss_match--;
      if((Ali[k][new]>=0)&&(ss_prot[k][Ali[k][new]]==ss_i[res]))ss_match++;
    }
    new++;
  }
  return(ss_match);
}

int Test_notali(int **Ali, struct protein *prot, int n, int N_ali)
{
   // Test if some structural position is not aligned
  int n_not=0, s_not=0;
  int col[N_ali], i, j;
  for(i=0; i<n; i++){
    int L=prot[i].len;
    for(j=0; j<N_ali; j++)col[j]=-1;
    for(j=0; j<N_ali; j++)if(Ali[i][j]>=0)col[Ali[i][j]]=j;
    int not_i=0; int notali[L];
    for(j=0; j<L; j++)if(col[j]<0){notali[not_i]=j; not_i++;}
    if(not_i){
      printf("WARNING, protein %d L=%d %d columns not aligned: ",i,L,not_i);
      for(j=0; j<not_i; j++)printf(" %d",notali[j]);
      printf("\n"); n_not+=not_i; s_not++;
    }
  }
  if(n_not){
    printf("WARNING, %d residues not aligned in %d proteins\n",n_not, s_not);
    if(0){printf("Exiting\n"); exit(8);}
  }
  return(s_not);
}
