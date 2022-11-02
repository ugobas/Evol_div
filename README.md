# Evol_div
Generates structurally modified multiple sequence alignments (MSA) and pairwise structural scores from an input MSA of proteins with known structure.  
Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa <ubastolla@cbm.csic.es>  
Associated paper: Piette O, Abia D, Bastolla U (2022) PC_sim: An integrated measure of protein sequence and structure similarity for improved alignments and evolutionary inference. submitted.

The code includes the needlemanwunsch aligner developed by Dr. Andrew C. R. Martin in the Profit suite of programs, (c) SciTech Software 1993-2007

## OVERVIEW:
1) Given a multiple sequence alignment (MSA) of proteins with known structure, for all aligned pairs it computes and prints (if required) sequence and structure similarity measures: aligned fraction, sequence identity, TM score, contact overlap and hybrid sequence-structure similarity PC_sim based on the main Principal Component of the four above similarity score.
2) The program computes five modified pairwise alignments (PA) that target different similarity scores: secondary structure (SS_ali), TM score (TM_ali), contact overlap (CO_ali) and PC_sim (PC_ali). Except for SS_ali, the modified alignments are based on the following algorithm: (1) The program identifies the closest residue of each residue of protein A in protein B under the targeted score; (2) It identifies neighbors (double best matches) (3) It makes frames of neighbors that are aligned in the input alignment; (4) It aligns neighbors that are consistent with the frames and are not aligned in the input alignment. In this way, the program does not have to score gaps.  
3) The program obtains a new MSA based on the PAs that target PC_sim in this way. (1) It transforms the PAs into the graph that connects the aligned residues. If all PAs are consistent with an MSA, each column of the MSA corresponds to a maximal clique of the graph, i.e. a maximal set of fully interconnected residues. (2) It determines the maximal cliques of the graph iteratively, exploiting the list of the neighbours of each residue in the native structure, through a non exhaustive strategy that reaches a good compromise between completeness and computational efficiency. (3) It assembles the cliques that are reciprocally consistent, i.e. they do not violate sequential order, starting from the largest one. (4) It assigns the residues that are not assigned to any maximal clique to the clique most connected to them if this  assignment is consistent with all other pre-existing cliques. If this is not possible, unassigned residues seed a new column. (5) It reconstructs the MSA from the set of all ordered columns (maximal cliques) and prints it for subsequent use in the .msa file.
4) For each pair and each modified alignment, the program computes the evolutionary divergence associated to each similarity score (see below). The PC_div divergence is the one that correlates best with the other divergence types.

If multiple conformations per protein are given (the program identifies them through sequence clustering), the structural similarity (divergence) between two proteins is computed as the maximum (minimum) across all the examined conformations. The output files have extension .prot.sim and .prot.dvi, respectively.
The similarity and divergence scores are computed both for the input alignment and for modified alignments that target different similarity scores (secondary structure, TM score, contact overlap and PC_sim). The best results are obtained targeting the hybrid sequence and structure similarity measure PC_sim.

Similarity measures:
1) Aligned fraction ali.  
2) Sequence identity SI.  
3) Contact overlap CO.  
4) TM-score TM (Zhang & Skolnick Proteins 2004 57:702).  
5) PC_sim, based on the main Principal Component of the four above similarity scores.  

Divergence measures:
1) Tajima-Nei divergence TN=-log((SI-S0)/(1-S0)) with S0=0.06 (Tajima F & Nei 1984, Mol Biol Evol 1:269).  
2) Contact_divergence CD=-log((q-q0(L))/(1-q0(L)) (Pascual-Garcia et al Proteins 2010 78:181-96).  
3) TM_divergence=-log((TM-TM0)/(1-TM0)).  
4) PC_divergence=-log((PC-PC0)/(1-PC0)).  

## USAGE:
### COMPILE:  
unzip Evol_div.zip  
make -f Evol_div.makefile  
cp Evol_div ~/bin/ (or whatever path directory you like)  

### RUN:  
Evol_div <Configuration file>

### EXAMPLE: 
Evol_div Input_Evol_div_Aldolase.in

The package provides the input MSA Aldolase_C1_MAFFT.fas  
You have to make a folder in which you store the PDB files named in Aldolase_C1_MAFFT.fas (characters 2-5 of ">1a5aA00" are the PDB code, 6 is chain identifier, 7-8 are not needed), name the folder in the "PDBDIR=<>" record (default is current directory) and name the extension of the files in the  "PDBDIR=<>" record (defult is .pdb)

### Configuration file: 

ALI= MSA in FASTA format, indicating names of PDB files.  
NORM=MIN   Normalize seq.id, TM, Cont.overlap based on MIN MAX or MEAN length  
PRINT_SIM=<0,1>   Print similarity measures? (default: 0)  
PRINT_CV=<0,1>    Print molecular clock violations? (default: 1)  
PRINT_PAIR=<0,1>  Print pairwise alignments that target PC_sim (default 0)  
PRINT_CLIQUE=<0,1> Print MSA based on cliques of PAs that target PC_sim  
ALI_SS=<0,1>    Correct multiple sequence alignment with secondary structure?  
SS_MULT=<0,1>   SS correction pairwise (0) or multiple (1)?  
SHIFT_MAX=9     Maximum shift when correcting MSA with sec.str.  

Parameters that can be passed either from the Input file or from the MSA file:  
PDBDIR=<directory of pdb files>  (default: current directory)  
PDBEXT=<extension of pdb files>  (default: .pdb)  

The protein name is the name of a PDB file, optionally followed by the chain index (Ex: >1opd.pdb A or >1opdA or >1opd_A);  

### OUTPUT
File with extension .sim .prot.sim (similarity)
File with extension .div .prot.div (divergence)
For each pair of conformations (optionally) and pair of proteins (clusters of conformations within 5 mutations in the sequence), the 5 similarity scores and 4 divergence scores are printed for 5 alignments: input, SS_ali (modified to target sec.str.), TM_ali (target TM score), CO_ali (target Contact overlap), PC_ali (target PC_sim).

.msa  
Multiple sequence alignment obtained with the clique method from the modified pairwise alignments that target PC_sim (PC_ali)  

.id   
Statistics of the five studied pairwise alignments:  
######## Input alignment  
######## SS alignment  
######## TM alignment  
######## CO alignment  
######## PC alignment  
