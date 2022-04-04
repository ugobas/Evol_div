float TM_score(float **d2, float *d02, int *ali, float ltar,
	       float **xca1_store, int n1, float **xca2_store, int n2,
	       int verbose);
void Align_TM(int *ali_new, float *d2min1, float **d2, float d02, int *ali,
	      int n1, int n2);
float TM_fast(float **d2_out, float d02, int *ali, int ltar, float *d2min1,
	      float **xca1_store, int n1, float **xca2_store, int n2);
void Examine_neighbors(float **d2, int *ali, float d02,
		       int *shift, int *id_3D,
		       int *neigh_ali, int *neigh_noali, int *neigh_noali_aaid,
		       char *seq1, int n1, char *seq2, int n2);
