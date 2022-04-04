/* TM score = 1/L sum_i 1/[1+(di/d0)^2]
   d0=TM_coeff*(L-L_offset)^(1/3)-TM_offset
*/
#define TM_coeff 1.24
#define TM_offset 1.8
#define L_offset 15
//#include "normalization.h"
#include "McLachlan.h"
#include "tm_score.h"
#include "Contact_divergence_aux.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

//

int Set_weights(int *w1, int *w2, float *w_ali, int *ali, int n1, int n2);
int Set_coord_ali(float *xca1, float **xca1_store, 
		  float *xca2, float **xca2_store, int *ali, int n1);

int Set_coord_cm(float *xca1, float **xca1_store, int *w, int n1);

float dist2(float *x1, float *x2);
void Rotate_vector(float *x, double rot[3][3]);
void Closest_neighbor(int *neighb1, float *dmin1, float *score, float d0,
		      float **d2, int n1, int n2, int order);
void All_distances(float **d2, float *x1, int n1, float *x2, int n2);

float TM_score(float **d2, float *d02, int *ali, float ltar,
	       float **xca1_store, int n1, float **xca2_store, int n2,
	       int verbose)
{
  if(ltar <= L_offset)return(-1);
  *d02=TM_coeff*pow(ltar-L_offset, 1./3)-TM_offset;
  if(*d02 <= 0)return(-1);
  if(verbose)printf("TM score, d0= %.3f n=%3.0f\n", *d02, ltar);
  (*d02)*=(*d02);

  // Copy coordinates of aligned residues
  float xca1[3*n1], xca2[3*n2];
  int n=Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);
  if(n<4)return(0);

  double TM_max=0;
  double rot[3][3];
  float w_TM[n1], w_tmp[n1], w[n1], dd[n1];
  int i;
  // Optimizing rotations
  for(int L_ini=n; L_ini>=(n/4); L_ini*=0.5){
    if(L_ini<3)break;
    for(int ini=0; ini<n; ini+=L_ini){
      int end=ini+L_ini; if(end>n)break;
      // initial superimposed fragment
      for(i=0; i<n; i++){
	if((i>=ini)&&(i<end)){w[i]=1;}else{w[i]=0;}
      }
      int nan=0;
      for(i=0; i<(3*n); i++)if(isnan(xca1[i])||isnan(xca2[i])){nan=1; break;}
      if(nan)Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);

      float TM_tmp=0, dd_opt[n1]; int align;
      for(int round=0; round<n; round++){
	float TM_tmp2=0;
	align=0; for(i=0; i<n; i++)if(w[i])align++; if(align<4)break;
	for(int round2=0; round2<n; round2++){
	  if(align<4)break;
	  rmsd_mclachlan_f(rot, xca1, xca2, w, n);
	  double TM=0, d2_max=0; int k=0, imax=-1;
	  for(i=0; i<n; i++){
	    dd[i]=dist2(xca1+k, xca2+k); k+=3;
	    TM+=1./(1+dd[i]/(*d02));
	    if(w[i]&&(dd[i]>d2_max)){d2_max=dd[i]; imax=i;}
	  }
	  if(TM>TM_tmp2){TM_tmp2=TM; for(i=0; i<n; i++)dd_opt[i]=dd[i];}
	  else {break;}
	  // Remove most divergent pair
	  if(imax>=0){w[imax]=0; align--;}
	}
	if(TM_tmp2>TM_tmp){
	  TM_tmp=TM_tmp2; for(i=0; i<n; i++)w_tmp[i]=w[i];
	}
	else{break;}
	// Superimpose only if di<d0
	for(i=0; i<n; i++)if(dd_opt[i]<(*d02)){w[i]=1;}else{w[i]=0;}
      } // end rounds
      if(TM_tmp>TM_max){
	TM_max=TM_tmp; for(i=0; i<n; i++)w_TM[i]=w_tmp[i];
      }
      if(verbose){
	printf("remove d<d0 L=%3d ini=%3d align=%d ltar=%.0f TM= %.1f\n",
	       L_ini, ini, align, ltar, TM_tmp);
      }
    }
  }
  if(verbose){
    printf("TM= %.1f d0=%.1f ltar=%.0f\n", TM_max, sqrt(*d02), ltar);
  }
  if(d2){
    // Superimpose all residues, including not aligned and compute distances
    n=Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);
    rmsd_mclachlan_f(rot, xca1, xca2, w_TM, n);
    int w1[n1], w2[n2];
    Set_weights(w1, w2, w_TM, ali, n1, n2);
    Set_coord_cm(xca1, xca1_store, w1, n1);
    Set_coord_cm(xca2, xca2_store, w2, n2);
    for(int k=0; k<(3*n2); k+=3)Rotate_vector(xca2+k, rot);
    All_distances(d2, xca1, n1, xca2, n2);
  }
  return(TM_max/ltar);
}

float TM_fast(float **d2_out, float d02, int *ali, int ltar, float *d2min1,
	      float **xca1_store, int n1, float **xca2_store, int n2)
{
  float TM_max=0, w[n1], w_TM[n1];
  float xca1[3*n1], xca2[3*n2];
  int n=Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);
  int i, k=0, align=0;
  for(i=0; i<n1; i++){
    if(ali[i]<0)continue;
    if(d2min1[i]<d02){w[k]=1; align++;}else{w[k]=0;} k++;
  }
  double rot[3][3];
  for(int round=0; round<n; round++){
    if(align<4)break;
    rmsd_mclachlan_f(rot, xca1, xca2, w, n);
    double TM=0, d2_max=0; int k=0, imax=-1;
    for(i=0; i<n; i++){
      float d2=dist2(xca1+k, xca2+k); k+=3;
      TM+=1./(1+d2/d02);
      if(w[i]&&(d2>d2_max)){d2_max=d2; imax=i;}
    }
    if(TM>TM_max){TM_max=TM; for(i=0; i<n; i++)w_TM[i]=w[i];}
    else {break;}
    // Remove most divergent pair
    if(imax>=0){w[imax]=0; align--;}
  }
  if(d2_out){
    // Superimpose all residues, included not aligned and compute distances
    Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);
    rmsd_mclachlan_f(rot, xca1, xca2, w_TM, n);
    int w1[n1], w2[n2];
    Set_weights(w1, w2, w_TM, ali, n1, n2);
    Set_coord_cm(xca1, xca1_store, w1, n1);
    Set_coord_cm(xca2, xca2_store, w2, n2);
    for(int k=0; k<(3*n2); k+=3)Rotate_vector(xca2+k, rot);
    All_distances(d2_out, xca1, n1, xca2, n2);
  }
  return(TM_max/ltar);
}

void Examine_neighbors(float **d2, int *ali, float d02,
		       int *shift, int *id_3D,
		       int *neigh_ali, int *neigh_noali,
		       int *neigh_noali_aaid,
		       char *seq1, int n1, char *seq2, int n2)
{
  // fraction of d0 to identify superimposed residues
  float d_thr=0.5, thr=d_thr*d02; int i;

  // Identify aligned and well superimposed residues below threshold
  if(id_3D){
    for(i=0; i<n1; i++){
      if(ali[i]<0){id_3D[i]=0; continue;}
      else{
	if(d2[i][ali[i]]<thr){id_3D[i]=1;}else{id_3D[i]=0;}
      }
    }
  }

  // Find the closest atom of str.2 for each atom of str.1 and vice-versa
  int neighb1[n1]; float d2min1[n1], score1[n1];
  Closest_neighbor(neighb1, d2min1, score1, d02, d2, n1, n2, 0);
  int neighb2[n2]; float d2min2[n2], score2[n2];
  Closest_neighbor(neighb2, d2min2, score2, d02, d2, n2, n1, 1);

  // Find neighbors that are not aligned
  if(neigh_noali && neigh_ali){
    int ali2[n2]; Invert_ali(ali2, n2, ali, n1);
    for(i=0; i<n1; i++){
      int j=neighb1[i];
      if((j>=0)&&(neighb2[j]==i)&&(d2min1[i]<thr)){
	if(ali[i]==j){neigh_ali[i]=1;}else{neigh_noali[i]=1;}
      }else{
	neigh_ali[i]=0;
	neigh_noali[i]=0;
      }

      if(id_3D[i] && neigh_ali[i]==0 && j>=0){
	int s1=abs(i-neighb2[j]), s2=abs(j-ali[i]);
	if(s1>s2){shift[i]=s1;}else{shift[i]=s2;}
      }else{
	shift[i]=0;
      }

      if(neigh_noali[i]){
	// Residue i is closest to not aligned residue
	// is ali identical?
	if(seq1 && ((ali[i] >=0 && seq1[i]==seq2[ali[i]])||
		    (j>=0 && ali2[j]>=0 && seq2[j]==seq1[ali2[j]]))){
	  neigh_noali_aaid[i]=1;	  
	}else{
	  neigh_noali_aaid[i]=0;
	}
      }
    }
  }

}

void Align_TM(int *ali_new, float *d2min1, float **d2, float d02, int *ali,
	      int n1, int n2)
{
  // Change alignment, optimize the TM score of new alignment
  float tthr=1./(1+0.75);

  // Find the closest atom of str.2 for each atom of str.1 and vice-versa
  int neighb1[n1]; float score1[n1]; //d2min1[n1], 
  Closest_neighbor(neighb1, d2min1, score1, d02, d2, n1, n2, 0);
  int neighb2[n2]; float d2min2[n2], score2[n2];
  Closest_neighbor(neighb2, d2min2, score2, d02, d2, n2, n1, 1);
  Align_neighbors(ali_new, ali, tthr, neighb1,score1,n1,neighb2,score2,n2);
}


int Set_weights(int *w1, int *w2, float *w_ali, int *ali, int n1, int n2)
{
  int i, k=0;
  for(i=0; i<n1; i++)w1[i]=0;
  for(i=0; i<n2; i++)w2[i]=0;
  for(i=0; i<n1; i++){
    if(ali[i]<0)continue;
    if(w_ali[k]){w1[i]=1; w2[ali[i]]=1;} k++;
  }
  return(0);
}

void All_distances(float **d2, float *x1, int n1, float *x2, int n2)
{
  for(int i=0; i<n1; i++){
    float *x1i=x1+3*i; int k=0;
    for(int j=0; j<n2; j++){
      d2[i][j]=dist2(x1i, x2+k); k+=3;
    }
  }
}

void Closest_neighbor(int *neighb1, float *dmin1, float *score, float d0,
		      float **d2, int n1, int n2, int order)
{
  for(int i=0; i<n1; i++){
    float d_min=1000, dd; int j_min=-1;
    for(int j=0; j<n2; j++){
      if(order==0){dd=d2[i][j];}else{dd=d2[j][i];}
      if(dd<d_min){d_min=dd; j_min=j;}
    }
    neighb1[i]=j_min; dmin1[i]=d_min; score[i]=1./(1+d_min/d0);
  }
}


void Rotate_vector(float *x, double rot[3][3]){
  float v[3];
  for(int i=0; i<3; i++){
    double w=0; for(int j=0; j<3; j++)w+=rot[i][j]*x[j]; v[i]=w;
  }
  for(int i=0; i<3; i++)x[i]=v[i];
}

int Set_coord_cm(float *x, float **x_store, int *w, int n){
  int i, j, k=0, m=0;
  double xc[3]; for(i=0; i<3; i++)xc[i]=0;
  for(i=0; i<n; i++){
    if(w[i])m++;
    for(j=0; j<3; j++){x[k]=x_store[i][j]; if(w[i])xc[j]+=x[k]; k++;}
  }
  for(i=0; i<3; i++)xc[i]/=m;
  k=0;
  for(i=0; i<n; i++){
    for(j=0; j<3; j++){x[k]-=xc[j]; k++;}
  }
  return(0);
}

float dist2(float *x1, float *x2){
  float dx=x1[0]-x2[0], dy=x1[1]-x2[1], dz=x1[2]-x2[2];
  return(dx*dx+dy*dy+dz*dz);
}

int Set_coord_ali(float *xca1, float **xca1_store,
		  float *xca2, float **xca2_store, int *ali, int n1)
{
  int k=0; 
  for(int i1=0; i1<n1; i1++){
    int i2=ali[i1]; if(i2<0)continue;
    for(int a=0; a<3; a++){
      xca1[k]=xca1_store[i1][a];
      xca2[k]=xca2_store[i2][a]; k++;
    }
  }
  return(k/3);
}
