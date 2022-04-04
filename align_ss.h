int Align_ss(int **Ali,
	     int **Ali_store, int N_ali, struct protein **prots, int m,
	     int SHIFT_MAX, int verbose);
void Write_ali_ss(int **Prot_ali_ss, struct protein *prots,
		  int N_pdb, int N_ali, char *name_in);
void Write_ss_ali(int **ali, struct protein *prots, int N_pdb, int N_ali,
		  char *name, char *what);
int Set_sec_str(struct protein *prots, int N_pdb);

int Test_notali(int **Ali, struct protein *prot, int n, int N_ali);
