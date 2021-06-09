# MISPhyl (Mutual Information Supermatrix Phylogenomic Tree)
A phylogenomics pipeline which reconstructs tree based on the single core orthologous genes (SCG) / proteins (SCP) from each sample with the utilization of mutual information to select the phylogenetically optimal genes.

## CITATION
1. Alexey M Kozlov, Diego Darriba, Tomáš Flouri, Benoit Morel, Alexandros Stamatakis, RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference, Bioinformatics, Volume 35, Issue 21, 1 November 2019, Pages 4453–4455
2. Darriba, D., Posada, D., Kozlov, A. M., Stamatakis, A., Morel, B., & Flouri, T. (2020). ModelTest-NG: a new and scalable tool for the selection of DNA and protein evolutionary models. Molecular Biology and Evolution, 37(1), 291-294. doi.org/10.1093/molbev/msz189
3. Flouri T., Izquierdo-Carrasco F., Darriba D., Aberer AJ, Nguyen LT, Minh BQ, von Haeseler A., Stamatakis A. (2014) The Phylogenetic Likelihood Library. Systematic Biology, 64(2): 356-362. doi:10.1093/sysbio/syu084
4. Kazutaka Katoh, Kazuharu Misawa, Kei‐ichi Kuma, Takashi Miyata, MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform, Nucleic Acids Research, Volume 30, Issue 14, 15 July 2002, Pages 3059–3066
5. Lechner et al (2011). Proteinortho: Detection of (Co-)Orthologs in Large-Scale Analysis. BMC Bioinformatics 2011 Apr 28;12(1):124.
6. Tan, J.L., Khang, T.F., Ngeow, Y.F. et al. A phylogenomic approach to bacterial subspecies classification: proof of concept in Mycobacterium abscessus. BMC Genomics 14, 879 (2013). https://doi.org/10.1186/1471-2164-14-879

# Table of Contents
* [Software and Script](#software-and-script)
* [Prerequisites](#prerequisites)
* [Getting Started](#getting-started)
* [Options](#options)
* [Main Output Folders](#main-output-folders)

## Software and Script
All the dependencies needed for the script are included within the tarball file. 
1) ProteinOrtho : 6.0.24 (Perl: 5 version 32 , Python: 3.8.3 , BLAST: 2.9.0-2 , Diamond: 2.0.4)
2) Pal2Nal: v14.1
3) Mafft : v7.453
4) Mutual Information script: 3.6.3
5) ModelTest-NG: 0.1.7
6) RAXML-NG : 1.0.1


## Prerequisites
Linux 64-bit system is required.
To run this script, you need:
1. Perl v5.08 or higher (test by typing "perl -v" in terminal)
2. Python v3.0 or higher (test by typing "python -V" in terminal)
3. Biopython module
```
	$ sudo apt install python3-pip
	$ pip3 install biopython
```

4. R language (test by typing "Rscript --version" in terminal)
5. R seqinr and parallel library
```
	$sudo apt-get update -y
	$sudo apt-get install -y r-cran-seqinr
```
If u do not have the sudo right, please contact your system administrator.




# Getting Started
## How to run this program?
1. Users are required to change each sequence header/description in all input files corresponding to the species.
<br />Format: <br>
`>[speciesName]_[accessionID]....`
<br><br />Users could utilize the `renameInput.py` script provided in dependencies/ to rename their files in working child directory.
Ensure there is NO '_', underscore in your species name and accession ID.
```
	$python3 ./dependencies/renameInput.py
```
2.Ensure there is no invalid character in your input files. 
<br>If you require help, `removeInvalidCharacter.py` is provided in dependencies/ folder.

3.For codon based alignment, ENSURE:<br>
<br>  a) Same IDs are used in both protein and nucleotide input files
<br>  b) amino acid files in main directory inputfolder/ whereas nucleotide files in directory ntfolder/ .



### Run all steps from 1 to 4 
1.Users are required to create a folder inputfolder/.
```
	$mkdir inputfolder
```
2.Put the files in folder inputfolder/. 
3.In your current directory, run the pipeline script. Eg:
```
	python3 pipeline.py [option]
```
  a) Run in default mode which accepts input files as amino acid sequence and utilize Diamond as blast program
```
	$python3 pipeline.py -e faa -f aa
```
  b) Run with nucleotide input files, blastn program and mutual information mode ON (ensure blastn is present in your system)
```	
	$python3 pipeline.py -e fna -f nt -d blastn -u
```
  c) Run codon-based alignment with mutal information mode ON. Amino acid files with .faa file extension in inputfolder/ and nucleotide files in ntfolder/.
```
	$python3 pipeline.py -o -u -f aa -e faa
```

Note: Codon Alignment (AA and corresponding NT files MUST have same filename, file extension need not to be)

### Step 1 : ProteinOrtho & (optional) Codon Alignment
1.Users are required to create a folder inputfolder/.
```
	$mkdir inputfolder
```
2.Put the files in folder inputfolder/. 
Example:<br>
a)Run step 1 with quiet mode ON and prefix for proteinortho as "project1"
```
	$python3 pipeline.py -s 1 -q -p project1 -e fa -f aa
```
3.Slight difference if codon based alignment is ENABLED, ENSURE:
<br>  i) Same IDs are used in both protein and nucleotide input files
|    Example Condition   |          Protein           |  		   	Nucleotide				     | 
|:------------:|:---------------------------|:-------------------------------|
|      Same ID      |>H.sapiens_ACE1180<br>ACDACDACD<br>>H.sapiens_ACD12739<br>ACDDCACDDC       |	>H.sapiens_ACE1180<br>GCUUGUGAUGCUUGUGAUGCUUGUGAU<br>>H.sapiens_ACD12739<br>GCUUGUGAUGAUUGUGCUUGUGAUGAUUGU |
|       Same Tag    |>H.sapiens_ACE80_1<br>ACDACDACD<br>>H.sapiens_ACD12739_2<br>ACDDCACDDC |>H.sapiens_ACE1180_1<br>GCUUGUGAUGCUUGUGAUGCUUGUGAU<br>><br>H.sapiens_ACDS2_2<br>GCUUGUGAUGAUUGUGCUUGUGAUGAUUGU	|


<br>  ii) amino acid files in main directory inputfolder/ whereas nucleotide files in directory ntfolder/ .
<br>  iii) RUN step 1 with codon alignment ENABLED will automatically finish up until the step 2, MSA.
<br>Example:
<br>a) Run all steps with codon alignment, mutual information mode ON.
```
	$python3 pipeline.py -e fasta -f aa -u -o -c 4
```
<br>a) Run step 1 to step 2 with codon alignment, mutual information mode ON
```
	$python3 pipeline.py -e fasta -f aa -u -o -c 4 -s 1
```


### Step 2 : Multiple Sequence Alignment
1. Make a directory orthologFamily/.
```
	$mkdir orthologFamily/
```
2. Move your .fasta files into folder orthologFamily/.
3.Example:
<br>a) Run step 2 with mafft program
```
	$python3 pipeline.py -s 2
```


### Step 3 : (Optional) Mutual Information & Concatenation 
1. Make a directory msa/.
```
	$mkdir msa
```
2. Move aligned files into folder msa/.
3.Example:
a) Mutual information mode ON with 10 median ranked genes and a aligned output file named aligned.fa 
```
	$python3 pipeline.py -s 3 -u -n aligned.fa -t 10
```
b) Step 3 with codon alignment and mutual information mode ON
```
	$python3 pipeline.py -s 3 -u -o
```


### Step 4 : Model Selection & Tree Construction 
1. Put your MSA file in current directory.
2. Ensure you have your partition.txt file in your current directory. 
If you run all the steps from 1 to 4, you need not to worry for this. Partition.txt is produced in step 3.
3. Example:<br>
a) Run step 4 with input file MSA.fa, nucleotide, output files prefix "tree1", 2 cpus, boostrapping of 250 and a partition file named "partition.txt".
	```
	$python3 pipeline.py -s 4 -r tree1 -c 2 -b 250 -n MSA.fa -f nt -o partition.txt
	```


## Options
|    Argument   |             Type             | Default 									|                            Description                                           |
|:-------------:|:----------------------------:|:------------------------------------------:|--------------------------------------------------------------------------------- |
|-h             | N/A                          |N/A 					 					| show help message																   |
|-s             | int                          |0 <br>(all)	 			 					| select step to be run <br>0:all (from step 1 to 3)<br>1:proteinortho<br>2:msa(muscle/mafft)<br>3:raxml-ng |
|-q             | N/A	                       |OFF	 	 				 					| do not show progress message   						|
|-c             | int	                       |-1 <br>(all available)	 					| number of cpu / threads to be utilized   				|
|-e             | string                       |faa 	 				 					| input file extension {fasta,faa,fna,fa}   			|
|-f             | string                       |N/A 	 				 					| type of input sequences {protein:aa / nucleotide:nt}	|
|-u             | N/A                          |OFF		 				 					| mutual information : select phylogenetically optimal genes for phylogenetic interference	    |
|-t             | int	                       |myproject	 	 		 					| prefix for proteinOrtho resulting file names   		|
|-p             | string                       |faa 	 				 					| input file extension {fasta,faa,fna,fa}   			|
|-d             | string                       |protein [diamond]<br> nucleotide [blastn] 	| blast program available for proteinOrtho 				|
|-i             | string                       |./dependencies/	 							| binpath for proteinOrtho blast program selection	    |
|-n             | string	                   |MSA.fa	 	 								| multiple sequence alignment output filename in FASTA format  |	
|-a             | int                          |0		 									| number of maximum iterations in mafft	    			|
|-o             | string	                   |partition.txt	 							| partition filename   									|
|-m             | string                       |modeltest 	 								| modelTest-NG output file prefix   					|
|-r             | string                       |T1 	 										| prefix for raxml-ng output files 						|
|-v             | string                       |bic		 								    | model selection for tree construction {bic,aic,aicc}	|
|-b             | int	                       |500	 	 									| number of bootstrap replicates for raxml-ng  			|	
	       		     		
     		


## Main Output Folders

### STEP 1 FOLDERS:
1.orthologFamily: core orthologous proteins/genes
2.nt_orthologFamily: corresponding core orthologous nucleotides (codon alignment)
3.codonAlignment: codon aligned nucleotides

### STEP 2 FOLDER:
4.msa: multiple sequence alignment

### STEP 3 FILES:
5.MSA.fa
<br>6.MI_genes.csv: Mutual Information file
<br>7.partition.txt

### STEP 4 FOLDER:
8.treeConstruction
