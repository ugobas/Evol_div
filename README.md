# Evol_div
Structural scores from multiple sequence alignments
Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa
Email: <ubastolla@cbm.csic.es>

It includes the needlemanwunsch aligner developed by Dr. Andrew C. R. Martin
in the Profit suite of programs, (c) SciTech Software 1993-2007

Given a multiple sequence alignment (MSA) of proteins with known structure, for all aligned pairs it computes and prints (if required) sequence and structure similarity measures (aligned fraction, sequence identity, TM score, contact overlap and hybrid sequence-structure similarity PC_sim) and estimates of the evolutionary divergence based on each of them. 
The structural similarity (divergence) between two proteins is the maximum (minimum) across all the examined conformations, if more than one is found. The output files have extension .prot.sim and .prot.dvi, respectively.
The similarity and divergence scores are computed both for the starting alignment and for modified alignments that target different similarity scores (TM score, contact overlap and PC_sim).

Similarity measures:
(1) Aligned fraction ali,
(2) Sequence identity SI,
(3) Contact overlap CO,
(4) TM-score TM (Zhang & Skolnick Proteins 2004 57:702)
(5) PC_sim, based on the main Principal Component of the four above similarity scores

Divergence measures:
(1) Tajima-Nei divergence TN=-log((SI-S0)/(1-S0)) with S0=0.06 (Tajima F & Nei 1984, Mol Biol Evol 1:269),
(2) Contact_divergence CD=-log((q-q0(L))/(1-q0(L)) (Pascual-Garcia et al Proteins 2010 78:181-96),
(3) TM_divergence=-log((TM-TM0)/(1-TM0)).
(4) PC_divergence=-log((PC-PC0)/(1-PC0)).


COMPILE:
>unzip Evol_div.zip
>make -f Evol_div.makefile
>cp Evol_div ~/bin/ (or whatever path directory you like)

RUN:
>Evol_div <alignment file>

EXAMPLE: Evol_div Input_Cont_Div_50044.aln
(you have to modify the names of the pdb files and directory in
Input_Cont_Div_50044.aln)

INPUT: 

ALI= MSA in FASTA format, indicating names of PDB files
NAME <name of output files> (default: MSA name)
PRINT_SIM=<0,1>   Print similarity measures? (default: 0)
PRINT_CV=<0,1>    Print clock violations? (default: 1)

Parameters that can be passed either from the Input file or from the MSA file:
PDBDIR=<directory of pdb files>  (default: current directory)
PDBEXT=<extension of pdb files>  (default: none)
The protein name is the name of a PDB file, optionally followed by the chain index (Ex: >1opd.pdb A or >1opdA or >1opd_A)\n\n");

OUTPUT (for each pair of proteins):

File with extension .sim (similarity, optional):
Sequence identity SI
Contact overlap CO
TM-score TM (structural), Zhang & Skolnick Proteins 2004 57:702
PC_sim

File with extension .div (divergence):
"Tajima-Nei divergence TN=-log((SI-S0)/(1-S0) S0=0.06 (Tajima F & Nei 1984, Mol Biol Evol 1:269).
Contact divergence  CD=-log((q-q0(L))/(1-q0(L)) (Pascual-Garcia et al Proteins 2010 78:181-96)
TM_divergence=-log((TM-TM0)/(1-TM0)) TM0=0.167 (Zhang & Skolnick Proteins 2004)
PC_divergenceProgram Evol_div
Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa
Email: <ubastolla@cbm.csic.es>

It includes the needlemanwunsch aligner developed by Dr. Andrew C. R. Martin
in the Profit suite of programs, (c) SciTech Software 1993-2007

Given a multiple sequence alignment (MSA) of proteins with known structure, for all aligned pairs it computes and prints (if required) sequence and structure similarity measures (aligned fraction, sequence identity, TM score, contact overlap and hybrid sequence-structure similarity PC_sim) and estimates of the evolutionary divergence based on each of them. 
The structural similarity (divergence) between two proteins is the maximum (minimum) across all the examined conformations, if more than one is found. The output files have extension .prot.sim and .prot.dvi, respectively.
The similarity and divergence scores are computed both for the starting alignment and for modified alignments that target different similarity scores (TM score, contact overlap and PC_sim).


Similarity measures:
(1) Aligned fraction ali,
(2) Sequence identity SI,
(3) Contact overlap CO,
(4) TM-score TM (Zhang & Skolnick Proteins 2004 57:702)
(5) PC_sim, based on the main Principal Component of the four above similarity scores

Divergence measures:
(1) Tajima-Nei divergence TN=-log((SI-S0)/(1-S0)) with S0=0.06 (Tajima F & Nei 1984, Mol Biol Evol 1:269),
(2) Contact_divergence CD=-log((q-q0(L))/(1-q0(L)) (Pascual-Garcia et al Proteins 2010 78:181-96),
(3) TM_divergence=-log((TM-TM0)/(1-TM0)).
(4) PC_divergence=-log((PC-PC0)/(1-PC0)).


COMPILE:
>unzip Evol_div.zip
>make -f Evol_div.makefile
>cp Evol_div ~/bin/ (or whatever path directory you like)

RUN:
>Evol_div <alignment file>

EXAMPLE: Evol_div Input_Cont_Div_50044.aln
(you have to modify the names of the pdb files and directory in
Input_Cont_Div_50044.aln)

INPUT: 

ALI= MSA in FASTA format, indicating names of PDB files
NAME <name of output files> (default: MSA name)
PRINT_SIM=<0,1>   Print similarity measures? (default: 0)
PRINT_CV=<0,1>    Print clock violations? (default: 1)

Parameters that can be passed either from the Input file or from the MSA file:
PDBDIR=<directory of pdb files>  (default: current directory)
PDBEXT=<extension of pdb files>  (default: none)
The protein name is the name of a PDB file, optionally followed by the chain index (Ex: >1opd.pdb A or >1opdA or >1opd_A)\n\n");

OUTPUT (for each pair of proteins):

File with extension .sim (similarity, optional):
aligned fraction
Sequence identity SI
Contact overlap q
TM-score TM (structural), Zhang & Skolnick Proteins 2004 57:702
PC_sim Principal component of the four above similarity scores

File with extension .div (divergence):
"Tajima-Nei divergence TN=-log((SI-S0)/(1-S0) S0=0.06 (Tajima F & Nei 1984, Mol Biol Evol 1:269).
Contact divergence  CD=-log((q-q0(L))/(1-q0(L)) (Pascual-Garcia et al Proteins 2010 78:181-96)
TM_divergence=-log((TM-TM0)/(1-TM0)) TM0=0.167 (Zhang & Skolnick Proteins 2004)
PC_divergence
They are printed both for the orignal MSA and for structurally modified MSA that target secondary structure (SS_ali), TM score (TM_ali), contact overlap (CO_ali) or PC_sim (PC_ali).
  
They are printed both for the orignal MSA and for the modified alignments.

