# MISPhyl (Mutual Information Supermatrix Phylogenomic) Pipeline
Phylogenomics analyses are commonly applied to answer various research questions pertaining to relationships of species and events on Earth. Although phylogenomic tree reconstructions have been widely used in research, it is still a great challenge for many researchers to deal with its multi-step technical procedure and data handling, especially on genomic data. Herein we present MISPhyl, a user-friendly pipeline which utilises supermatrix-based procedure to yield phylogenomic tree. While a supermatrix phylogenomic tree aims to amplify phylogenetic signals, there are chances to include phylogenetic noises into the tree reconstruction. To address the issue, this automated pipeline has also implemented a Mutual Information (MI) approach to perform systematic selection of genes with optimal phylogenetic signals for phylogenomic inference. The MI approach has been previously discussed for its ability to generate a reliable phylogenomic tree and identify species-specific markers in Mycobacterium abscessus Complex (Tan et al., 2013).


# Table of Contents
* [Software and Script](#software-and-script)
* [Prerequisites](#prerequisites)
* [Getting Started](#getting-started)
* [Options](#options)
* [Main Output Folders](#main-output-folders)
* [Highlight](#highlight)

## Software and Script
All the dependencies needed for the script are included within the tarball file. 
1) ProteinOrtho : 6.0.24 (Perl: 5 version 32 , Python: 3.8.3 , BLAST: 2.9.0-2 , Diamond: 2.0.4)
2) Pal2Nal: v14.1
3) Mafft : v7.490
4) Mutual Information script: R: 3.6.3
5) ModelTest-NG-static: 0.1.7 _(Please note that modeltest-ng-static binary file relies on compatible hardware)_
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
## How to run this program? [CRITICAL]
1. Users are required to change each sequence header/description in all input files corresponding to the species.
<br />Format: <br>
`>[speciesName]_[accessionID]....`
<br><br />Users could utilize the `renameInput.py` script provided in dependencies/ to rename their files in working child directory.
Ensure there is NO '_', underscore in your species name and accession ID.
```
	$python3 ./dependencies/renameInput.py
```
2. Ensure there is no invalid character in your input files. 
<br>If you require help, `removeInvalidCharacter.py` is provided in dependencies/ folder.

3. For codon based alignment, ENSURE:
	1. Same IDs are used in both protein and nucleotide input files
	2. Amino acid files in main directory inputfolder/ whereas nucleotide files in directory ntfolder/ .
	3. Value for option -f is set to "aa".



### Run all steps from 1 to 4 
1. Users are required to create a folder inputfolder/.
```
	$mkdir inputfolder
```
2. Put the files in folder inputfolder/. 
3. In your current directory, run the pipeline script. Eg:
```
	python3 MISPhyl.py [option]
```
  a) Run in default mode which accepts input files as amino acid sequence and utilize Diamond as blast program
```
	$python3 MISPhyl.py -e faa -f aa -u
```
  b) Run with nucleotide input files, blastn program and mutual information mode ON (ensure blastn is present in your system)
```	
	$python3 MISPhyl.py -e fna -f nt -d blastn
```
  c) Run codon-based alignment with mutal information mode ON. Amino acid files with .faa file extension in inputfolder/ and nucleotide files in ntfolder/.
```
	$python3 MISPhyl.py -o -f aa -e faa
```

Note: Codon Alignment (AA and corresponding NT files MUST have same filename, file extension need not to be)

### Step 1 : ProteinOrtho & (optional) Codon Alignment
1. Users are required to create a folder inputfolder/.
```
	$mkdir inputfolder
```
2. Put the files in folder inputfolder/. 
Example:<br>
a) Run step 1 with quiet mode ON and prefix for proteinortho as "project1"
```
	$python3 MISPhyl.py -s 1 -q -p project1 -e fa -f aa -u
```
3.Slight difference if codon based alignment is ENABLED, ENSURE:
<br>  i) Same IDs are used in both protein and nucleotide input files <br>
|    Example Condition   |          Protein           |  		   	Nucleotide				     | 
|:------------:|:---------------------------|:-------------------------------|
|      Same ID      |>**H.sapiens**_ACE1180<br>ACDACDACD<br>>**H.sapiens**_ACD12739<br>ACDDCACDDC       |	>**H.sapiens**_ACE1180<br>GCUUGUGAUGCUUGUGAUGCUUGUGAU<br>>**H.sapiens**_ACD12739<br>GCUUGUGAUGAUUGUGCUUGUGAUGAUUGU |
|       Same Tag    |>H.sapiens_ACE80_**1**<br>ACDACDACD<br>>H.sapiens_ACD12739_**2**<br>ACDDCACDDC |>H.sapiens_ACE1180_**1**<br>GCUUGUGAUGCUUGUGAUGCUUGUGAU<br>><br>H.sapiens_ACDS2_**2**<br>GCUUGUGAUGAUUGUGCUUGUGAUGAUUGU	|



<br>  ii) amino acid files in main directory inputfolder/ whereas nucleotide files in directory ntfolder/ .
<br>  iii) RUN step 1 with codon alignment ENABLED will automatically finish up until the step 2, MSA.
<br>Example:
<br>a) Run all steps with codon alignment, mutual information mode ON.
```
	$python3 MISPhyl.py -e fasta -f aa -o -c 4
```
<br>b) Run step 1 to step 2 with codon alignment, mutual information mode ON
```
	$python3 MISPhyl.py -e fasta -f aa -o -c 4 -s 1
```


### Step 2 : Multiple Sequence Alignment
1. Make a directory orthologFamily/.
```
	$mkdir orthologFamily/
```
2. Move your .fasta files into folder orthologFamily/.
3.Example:
	1. Run step 2 with mafft program
```
	$python3 MISPhyl.py -s 2
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
	$python3 MISPhyl.py -s 3 -u -n aligned.fa -t 10
```
b) Step 3 with codon alignment and mutual information mode ON
```
	$python3 MISPhyl.py -s 3 -o
```


### Step 4 : Model Selection & Tree Construction 
1. Put your MSA file in current directory.
2. Ensure you have your partition.txt file in your current directory. 
If you run all the steps from 1 to 4, you need not to worry for this. Partition.txt is produced in step 3.
3. Example:<br>
a) Run step 4 with input file MSA.fa, nucleotide, output files prefix "tree1", 2 cpus, boostrapping of 250 and a partition file named "partition.txt".
	```
	$python3 MISPhyl.py -s 4 -r tree1 -c 2 -b 250 -n MSA.fa -f nt -o partition.txt
	```


## Options
|    Argument   |             Type             | Default 									|                            Description                                           |
|:-------------:|:----------------------------:|:------------------------------------------:|--------------------------------------------------------------------------------- |
|-h             | N/A                          |N/A 					 					| show help message																   |
|-s             | int                          |0 <br>(all)	 			 					| select step to be run <br>0:all (from step 1 to 3)<br>1:proteinortho<br>2:msa(muscle/mafft)<br>3:raxml-ng |
|-cpus             | int	                       |-1 <br>(all available)	 					| number of cpu / threads to be utilized   				|
|-f             | string                       |faa 	 				 					| input file extension {fasta,faa,fna,fa}   			|
|-i             | string                       |N/A 	 				 					| type of input sequences {protein:aa / nucleotide:nt}	|
|-mi             | N/A                          |ON		 				 					| mutual information : select optimal phylogenetical signal genes for phylogenomic interference	    |
|-p             | int	                       |myproject	 	 		 					| prefix for proteinOrtho resulting file names   		|
|-f             | string                       |faa 	 				 					| input file extension {fasta,faa,fna,fa}   			|
|-algo             | string                       |protein [diamond]<br> nucleotide [blastn] 	| blast program available for proteinOrtho 				|
|-path             | string                       |./dependencies/	 							| binpath for proteinOrtho blast program selection	    |
|-msa             | string	                   |MSA.fa	 	 								| multiple sequence alignment output filename in FASTA format  |	
|-maxiter             | int                          |0		 									| number of maximum iterations in mafft	    			|
|-partition             | string	                   |partition.txt	 							| partition filename   									|
|-n             | string                       |modeltest 	 								| modelTest-NG output file prefix   					|
|-x             | string                       |T1 	 										| prefix for raxml-ng output files 						|
|-model             | string                       |bic		 								    | model selection for tree construction {bic,aic,aicc}	|
|-b             | int	                       |500	 	 									| number of bootstrap replicates for raxml-ng  			|	
	       		     		
     		


## Main Output Folders

### Step 1 Folders:
1. orthologFamily/: core orthologous proteins/genes
2. nt_orthologFamily/: corresponding core orthologous nucleotides (codon alignment)
3. codonAlignment/: codon aligned nucleotides

### Step 2 Folder:
4. msa/: multiple sequence alignment files

### Step 3 Files:
5. MSA.fa: concatenated MSA file
6. MI_genes.csv: Mutual Information file
7. partition.txt

### Step 4 Folder:
8. treeConstruction/: constructed tree files

## Highlight
If unfortunately, you encounter this error when reaching tree construction step:
```
ERROR: modeltest-ng-static binary file seems to be not compatible with your hardware. But no worries.
There are two recommended ways to solve this issue: 
1. Download the source files from modeltest-ng github \'https://github.com/ddarriba/modeltest/wiki/Download-and-Install\' and run the partition file using modeltest-ng instead of modeltest-ng-static AND comment the try and except block in MISPhyl.py script (line 495 to 499). Rerun step 4, tree construction.
2.Comment the try and except block in MISPhyl.py script (line 495 to 499) AND stick to one substitution model for all the genes (kindly make use of the argument \'-model\' to provide the wanted subtitution model). Rerun tree construction, step 4.
```
Please follow the suggested ways to resolve it. Have fun!

## References
1. Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). https://doi.org/10.1038/s41592-021-01101-x
2. Darriba, Di. et.al. (2020). ModelTest-NG: A New and Scalable Tool for the Selection of DNA and Protein Evolutionary Models. Molecular Biology and Evolution, 37, 291-294. https://doi.org/10.1093/molbev/msz189
3. Katoh, K. et.al. (2002). MAFFT: A novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research, 30, 3059-3066. https://doi.org/10.1093/nar/gkf436
4. Kozlov, A. M. et.al. (2019). RAxML-NG: A fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35, 4453-4455. https://doi.org/10.1093/bioinformatics/btz305
5. Lechner, M. et.al. (2011). Proteinortho: Detection of (Co-)orthologs in large-scale analysis. BMC Bioinformatics, 12, 124. https://doi.org/10.1186/1471-2105-12-124
6. Suyama, M. et.al. (2006). PAL2NAL: Robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Research, 34, W609-W612. https://doi.org/10.1093/nar/gkl315
7. Tan, J. L. et.al. (2013). A phylogenomic approach to bacterial subspecies classification: Proof of concept in Mycobacterium abscessus. BMC Genomics, 14, 879. https://doi.org/10.1186/1471-2164-14-879
