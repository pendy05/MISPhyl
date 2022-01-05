#!/usr/bin/env python
'''
Tools & Script:  1)Proteinortho
            2)pal2nal
        	3)MAFFT
		4)Mutual Information Script        
		5)ModelTest-NG
		6)RAxML-NG
        '''
from Bio import SeqIO
#from colorama import Fore, Back, Style 
import argparse
import os
import csv
import re
import glob
import sys
import subprocess
import math
import textwrap

parser = argparse.ArgumentParser(prog='PROG',
				formatter_class=argparse.RawDescriptionHelpFormatter,
				description=textwrap.dedent('''\
                           =====================================
				Pipeline for Phylogenomics
			   =====================================
                        This is a phylogenomics pipeline that offers
 
1)Single Core Copy Orthologous Genes/Proteins Extraction and Codon Alignment
2)Multiple Sequence Alignment 		
3)Mutual Information Selection (Nucleotide only) and Concatenation
4)Model Selection and Tree Construction

For more information and details, read on readme.txt.



=========================================
                  CITATION
=========================================
1. Alexey M Kozlov, Diego Darriba, Tomáš Flouri, Benoit Morel, Alexandros Stamatakis, RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference, Bioinformatics, Volume 35, Issue 21, 1 November 2019, Pages 4453–4455
2. Darriba, D., Posada, D., Kozlov, A. M., Stamatakis, A., Morel, B., & Flouri, T. (2020). ModelTest-NG: a new and scalable tool for the selection of DNA and protein evolutionary models. Molecular Biology and Evolution, 37(1), 291-294. doi.org/10.1093/molbev/msz189
3. Flouri T., Izquierdo-Carrasco F., Darriba D., Aberer AJ, Nguyen LT, Minh BQ, von Haeseler A., Stamatakis A. (2014) The Phylogenetic Likelihood Library. Systematic Biology, 64(2): 356-362. doi:10.1093/sysbio/syu084
4. Kazutaka Katoh, Kazuharu Misawa, Kei‐ichi Kuma, Takashi Miyata, MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform, Nucleic Acids Research, Volume 30, Issue 14, 15 July 2002, Pages 3059–3066
5. Lechner et al (2011). Proteinortho: Dmetection of (Co-)Orthologs in Large-Scale Analysis. BMC Bioinformatics 2011 Apr 28;12(1):124.
6. Tan, J.L., Khang, T.F., Ngeow, Y.F. et al. A phylogenomic approach to bacterial subspecies classification: proof of concept in Mycobacterium abscessus. BMC Genomics 14, 879 (2013). https://doi.org/10.1186/1471-2164-14-879

'''))

parser.add_argument("-e", type=str, action="store", dest="fileextension",
			choices=["faa","fa","fasta","fna"], default="fa", help="input file extension [default: fa]")
parser.add_argument("-f", type=str, action="store", dest="inputtype",
			choices=["nt","aa"], help="type of input sequences: nucleotide[nt] or protein[aa] ")
parser.add_argument("-s", type=int, action="store",dest="step", default=0,
			choices=[0,1,2,3,4], help="step to be done [default:0]\
			0 = all ; 1 = proteinortho (protein extraction); 2 = mafft (multiple sequence alignment) ; 3 = mutual information & concatenation ,4 = raxml-ng (tree construction)")
parser.add_argument("-q", help="Do not show progress message", action="store_true", dest="quiet")
parser.add_argument("-w", help="Have minimum of four input files; required for tree construction (DEFAULT:ON)", action="store_true", dest="limit")
parser.add_argument("-c", type=int, action="store", dest="cpus",
			default=-1,help="number of cpu threads to be utilized [default: all available]")


#Proteinortho
parser.add_argument("-p", type=str, action="store", dest="name",
                    default="myproject", help = "prefix for proteinortho resulting file names")
parser.add_argument("-d", type=str, action="store", dest="blastprogram",
			default="diamond", help="blast program for proteinortho [default:diamond(protein) blastn(nucleotide)] \
			{autoblast|blastp|blastn|blastp_legacy|blastn_legacy|tblastx_legacy|diamond|usearch|ublast|lastp|lastn|rapsearch|topaz|blatp|blatn|mmseqsp|mmseqsn}")
parser.add_argument("-i",type= str, action="store", dest="binpath",
			default="",help="binpath for the blast program selected")

#codon alignment
parser.add_argument("-o", action="store_true", dest="codon",
			help="codon based alignment, translate protein alignments to nucleotide alignments (default:OFF)")

#mafft
parser.add_argument("-n", type=str, action="store", dest="MSAout",
                    default="MSA.fa", help = "MSA output filename in FASTA format with .fa extension")
parser.add_argument("-a", type=int, action="store", dest="mafftmaxiteration",
                    default=0, help = "number of mafft max iterations")

#mutualInformation
parser.add_argument("-u", help="select phylogenetically optimal genes for phylogenetic interference based on mutual information\
                         ,ONLY available for nucleotide (default:ON)", action="store_false", dest="mutualinfo")
parser.add_argument("-t", type=int, action="store",dest="medianrange",
			default=50, help="number of median-ranked range genes in MI_genes.csv to be concatenated (default:50)")


#concatenation
parser.add_argument("-k", type=str, action="store", dest="partition",
                    default="partition.txt", help = "partition filename [default: partition.txt]")
#modeltest-ng
parser.add_argument("-m", type=str, action="store", dest="modeltestPrefix",
                    default="modeltest", help = "prefix for modelTest-NG output file")

#raxml-ng
parser.add_argument("-r", type=str, action="store", dest="prefix",
                    default="T1" , help = "prefix for raxml-ng output file")
parser.add_argument("-v", type=str, action="store", dest="subsModel",
			choices=["aic","aicc","bic"], default="bic", help="type of substitution model to be used in tree construction ")
#parser.add_argument("-o", type=str, action="store", dest="model",
#                    default="LG+G8+F", help="raxml-ng evolutionary model {default- nucleotide:GTR+G ; protein: LG+G8+F}")
parser.add_argument("-b",type=int, action="store", dest="bootstrap", 
			default=500, help="number of bootstrap replicates for raxml-ng [default:500]")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

cmd1="" #for proteinOrtho
cmd2="" #for MSA; mafft
cmd3="" #for tree construction

if args.quiet: #quiet mode ON
	cmd1 += "-silent "
	cmd2 += "--quiet "
	cmd3 += "--log WARNING "

if args.cpus > os.cpu_count():
	sys.exit("Maximum CPU: %s"%os.cpu_count())
elif args.cpus == -1 : #all available cpus
	args.cpus = os.cpu_count()

cmd1 += ("-cpus=%s " %args.cpus)
cmd2 += ("--thread %s " %args.cpus)
cmd3 += ("--threads %s " %args.cpus)

'''
proteinortho
=>cluster the proteins into families
=>extract core SCP from protein families
'''

def proteinextraction():
	print("\nSTEP 1: Core Orthologous Protein Extraction\n")
	if args.codon: #codon based alignment mode ON
		if not os.path.exists("./ntfolder/"):
			sys.exit("ERROR:\"ntfolder\" child directory is NOT FOUND!\n")

	if not os.path.exists("./inputfolder/"):
		print("ERROR:\"inputfolder\" child directory is not found!\n")
		sys.exit()
	
	if not args.inputtype:
		sys.exit("ERROR:Please select the sequence type, nt or aa.\n")
	
	MAX= len(glob.glob1("./inputfolder/","*.%s" %(args.fileextension))) #count number of . sample file in current directory
	if args.step == 0:
		if MAX < 4 and args.limit != True: #must have >= 4 samples
			sys.exit("\nERROR: Less than 4 files with .%s extension being detected. Minimum of 4 input files are required for the tree construction. \nIf you would like to disable this feature for step 1 in all-in-one step, please include \"-w\" in your command line." %(args.fileextension))
	
	if args.blastprogram == "diamond":#default blast program
		if args.inputtype == "nt":
			os.system("perl ./dependencies/proteinortho-master/proteinortho6.pl %s -p=blastn -project=%s inputfolder/*.%s " %(cmd1, args.name, args.fileextension))
		else:
			os.system("perl ./dependencies/proteinortho-master/proteinortho6.pl %s -binpath=\"./dependencies/\" -p=%s -project=%s inputfolder/*.%s " %(cmd1, args.blastprogram, args.name, args.fileextension))
	else:
		if args.binpath == "":#search from directories in executable path
			os.system("perl ./dependencies/proteinortho-master/proteinortho6.pl %s -p=%s -project=%s inputfolder/*.%s " %(cmd1, args.blastprogram, args.name, args.fileextension))
		else:#bin path provided by user
			os.system("perl ./dependencies/proteinortho-master/proteinortho6.pl %s -binpath=\"%s\" -p=%s -project=%s inputfolder/*.%s " %(cmd1, args.binpath, args.blastprogram, args.name, args.fileextension))
	
	if os.path.isfile("%s.proteinortho.tsv" %(args.name)): #check if this file present in PWD
		filesize = os.path.getsize("%s.proteinortho.tsv" %(args.name))

		if filesize == 0:#is the file empty?
			sys.exit("\nERROR: File \"%s.proteinortho.tsv\" is EMPTY!\n" %(args.name))
		
		else: #not empty then proceed
            	#get coreSCG in .tsv file
			SCGoutput = open("coreSCP.tsv","w")
			with open("%s.proteinortho.tsv" %(args.name)) as tsv:
				fileContent = csv.reader(tsv, delimiter = '\t')
				tsvOutput = csv.writer(SCGoutput, delimiter = '\t')
				for row in fileContent:
					if row[0] == "# Species":
						tsvOutput.writerow(row)
					elif row[0] == row[1] == str(MAX):
						tsvOutput.writerow(row)
			SCGoutput.close()

			os.system("./dependencies/proteinortho-master/src/proteinortho_grab_proteins.pl -tofiles coreSCP.tsv inputfolder/*.%s" %(args.fileextension))
			for filename in os.listdir("."):
				if filename.startswith("coreSCP.tsv.OrthoGroup") and filename.endswith(".fasta"):
					os.rename(filename,filename.replace(".tsv.","_"))

			if not os.path.exists('./orthologFamily'):
				os.makedirs("./orthologFamily")
			os.system("mv coreSCP*  %s* orthologFamily/" %args.name)
			print("\nSingle Core Orthologous Genes/Protein Extraction done.")
			print("Output files could be found in directory orthologFamily/ .\n")
		print("\nCore Orthologous Protein / Gene Extraction Done\n")
		if args.codon:
			codonAlignment()
			
	else:
		sys.exit("\nERROR:File \"%s.proteinortho.tsv\" is NOT EXISTED in current directory \"%s\"\n" %(args.name, os.getcwd()))


'''
MSA
=>align the core SCP sequence of protein families
=>optional: gene selection for nucleotide sequence using mutual information(default: OFF)
=>concatenate aligned MSAs to a file
'''

#sort file based on filename from OrthoGroup1 up to OrthoGroup#
numbers = re.compile(r'(\d+)')
def numericalSort(value): 
	parts = numbers.split(value)
	parts[1::2] = map(int,parts[1::2])
	return parts

def msa():
	print("\nSTEP2: Multiple Sequence Alignment\n")
	if os.path.exists("./orthologFamily"): #if the imput directory is created
		for fl in os.listdir("./orthologFamily/"): 
			if fl.endswith(".fasta"):
				os.system("./dependencies/mafft-linux64/mafft.bat %s --maxiterate %s --auto ./orthologFamily/%s > aligned%s" % (cmd2, args.mafftmaxiteration, fl, fl))
    
		if not os.path.exists('./msa'):
			os.makedirs('./msa')
		os.system("mv aligned* msa/")
		print("\nMultiple Sequence Alignment done.")
		print("Output files could be found in directory msa/ .\n")
		print("\nMultiple Sequence Alignment Done\n")
	else:
		sys.exit("ERROR: Input directory \"orthologFamily\" is NOT FOUND in current directory.\n For MSA, please ensure this directory, \"orthologFamily\" is created with input files located inside.")
	

def codonAlignment():
	print("\nSTEP 1.5: Codon Alignment\n")
	if os.path.isfile("./orthologFamily/coreSCP.tsv"):
		filesize = os.path.getsize("./orthologFamily/coreSCP.tsv")
	else:
		sys.exit("ERROR: coreSCP.tsv is NOT FOUND!")
	
	if filesize == 0:#is the file empty?
		print("ERROR: File \" orthologFamily/coreSCP.tsv \" is EMPTY!")
		sys.exit()

	else:
		speciesNum = 0
		with open ("./orthologFamily/coreSCP.tsv") as inOrtho:
			content = csv.reader(inOrtho,delimiter="\t")
			count = 0
			col = 0
			
			flag = False
			for index,record in enumerate(content):
				if record[0] == "# Species":
					speciesNum = len(record)-3 #3 fixed col are attributes
					contentList = [[] for _ in range(speciesNum)]
					for i in range(len(contentList)):
                                                contentList[i].append(record[i+3].split(".")[0])
					flag = True
				
				if flag and index == 0:
					pass
				elif flag and index >=1:
					for i in range(len(contentList)):
                                                 contentList[i].append(record[i+3].split("_")[-1])
				else:
					sys.exit("ERROR: Check the coreSCP.tsv file!!")# shouldnt have this error

		orthoNum = len(contentList[0])-1
		ls = [i for i in os.listdir("./ntfolder/")] # need change to i.split(".")[0]
		ls.sort()
		print("Expected Orthologs Number to be extracted from corresponding files: ",orthoNum)

		absentNtOrthologs=[]#write to the absentNtOrthologs file only
		absentNtID=[] #remove the AA orthoGroup and its corresponding Nt file from the folder
		for pattern, target in zip(contentList, ls):
			count=0
			if target.find(pattern[0]) != -1:# found the corresponding nt file
				group=0
				with open("./ntfolder/%s" % target, "r") as targetfile:
					for i in range(1,len(pattern)):
						cnt = 0
						targetfile.seek(0)
						flag=True
						
						for record in SeqIO.parse(targetfile,"fasta"):
							if record.id.split("_")[-1] == pattern[i]: #aa ID/tag == nt ID/tag
								count+=1
								with open("nt_OrthoGroup%s.fasta" %str(i-1), "a") as outfile:
									outfile.write(">"+ record.description +"\n"+ str(record.seq) +"\n")
									group +=1
									if group >= orthoNum:
										group = 0
								flag=False
								break
						if flag:
							absentNtOrthologs.append("%s,%s"%(target,pattern[i]))
							absentNtID.append(i-1)
								
				print("Number of sequences found in nt file %s : %d"%(pattern[0],count))
						
			else:
				sys.exit("\n ERROR: Corresponding nt file %s is FOUND for aa file %s!" %(target,pattern))

		print("Absent Ortholog (file,accessionID): ",absentNtOrthologs)
		print("Ortholog Group with absent sequence: Group ", absentNtID)
		with open("absentNtOrthologs.csv","w") as writefile:
			writefile.write("Corresponding Nucleotide File, Not Found ID\n")
			for i in absentNtOrthologs:# orthologs id not found in corresponding nt file
				writefile.write("%s\n"%str(i))		
					
			count += 1
		if not os.path.exists("./nt_orthologFamily"):
			os.makedirs("./nt_orthologFamily")

		os.system("mv nt_OrthoGroup* nt_orthologFamily/")

		count = 0
		#msa
		msa()

		for i, j in zip(sorted(glob.glob("./msa/alignedcoreSCP_OrthoGroup*.fasta"),key=numericalSort),sorted(glob.glob("./nt_orthologFamily/nt_OrthoGroup*.fasta"),key=numericalSort)):		
			flag=True			
			for m in absentNtID:			
				if int(re.findall(r'\d+',i)[0]) == m and int(re.findall(r'\d+',j)[0]) == m: #get digit from the filename eg: 0 from coreSCP_OrthoGroup0.fasta
					count+=1
					flag = False
					break
			if flag:
				print("ortholog: ",i, "\tnt:", j)
				os.system("./dependencies/pal2nal.v14.1/pal2nal.pl %s %s > codon%d.fasta -output fasta"%(i,j,count))
				count +=1

		if not os.path.exists("./codonAlignment"):
			os.makedirs("./codonAlignment")

		os.system("mv codon* codonAlignment/")
	print("\nCodon Alignment Done\n")

def concatenation():
	print("\nSTEP 3: Concatenation\n")
	table={}
	pos, endpos = 0,0
	partition=""
	filelocation = "./msa/"
	if not args.inputtype:
		sys.exit("ERROR: Please select the sequence type, nt or aa.\n")
		
	if args.codon:
		filelocation="./codonAlignment/"

	for i in os.listdir(filelocation):# check if there is any empty file, if yes, remove it
		if os.stat("%s%s"%(filelocation,i)).st_size == 0:
			os.system("rm %s%s"%(filelocation,i))

	if len(os.listdir(filelocation)) == 0:
		sys.exit("ERROR: No file in directory %"%(filelocation))
		
	if (args.mutualinfo and args.inputtype == "nt" and filelocation == "./msa/") or (args.mutualinfo and args.codon and filelocation == "./codonAlignment/"): #mutual information mode ON
        #run Rscript MI
		command = "Rscript"
		path2script = "./dependencies/MutualInfo.R"
		cmd = [command, path2script, filelocation, str(args.cpus)]
		try:
			x = subprocess.check_output(cmd,universal_newlines=True)
		except:
			print("ERROR encountered for mutual information. Ensure there is no invalid character in your samples and the samples are located in folder msa/.")
			sys.exit(1)
		print("\nMI_genes.csv can be found in current directory.\n")

        #gene selection
		with open('MI_genes.csv', 'r', newline='') as infile:
			csv_input = csv.DictReader(infile)
			filtered = [row for row in csv_input if float(row["gappiness"]) < 50]	#filter data with > 50 gappiness
			data = sorted(filtered, key=lambda row: row['MI_genes'])
			MImedian = float(data[-1]["MI_genes"]) / 2 		#get median of MI
			if MImedian - float(data[0]["MI_genes"]) > 0:
                		min_dff = MImedian - float(data[0]["MI_genes"])
			elif float(data[0]["MI_genes"]) - MImedian > 0:
				min_dff = float(data[0]["MI_genes"]) - MImedian
			itr = 0
			#select the row closest to median
			filenames = []
			for i,j in enumerate(data):
				dff = MImedian - float(j["MI_genes"])
				dff1 = float(j["MI_genes"]) - MImedian
				filenames.append(j[""]+".fasta")
				if dff < min_dff and dff > 0:
					min_dff, itr = dff, i
				elif dff1 < min_dff and dff1 > 0:
					min_dff, itr = dff1, i
				else:
					pass
            

			half = args.medianrange/2
			if str(half).split(".")[-1] == "0":
				half=int(half)

			#slice median ranked genes
			if itr < half :
				extra = half -itr
				start, end = int(0), int(itr+half+extra) 
			elif itr + math.ceil(half) > len(data)-1:
				start, end = int(len(data)-args.medianrange), int(len(data))
			elif type(half) == int:    
				start, end = int(itr-half), int(itr+half)
			elif type(half) == float:
				start, end = int(itr-math.ceil(half))+1,int(itr+math.ceil(half))

			if os.path.isfile("%s" % (args.partition)):
				print("%s is found in this directory and will be overwritten!!" % (args.partition))
			with open(args.partition,"w") as partfile:#create partition file
		    		#loop through files for concatenation
				for filename in sorted(filenames[start:end]):
					with open(filename,'r',newline=''):
						for record in SeqIO.parse(filename,"fasta"):
							endpos = len(str(record.seq))
							seq = str(record.seq)
							newheader = record.description.split("_",1)[0]
							table[newheader] = table.get(newheader,'') + seq
					partitionName = filename.split("/")[-1].split(".")[0]						
					partfile.write("nt, "+partitionName+"="+str(pos+1)+"-"+str(endpos+pos)+"\n")#write to partition file
					pos = endpos+pos

	else:
		if os.path.exists(filelocation): #if the input directory is created
			if os.path.isfile("%s" % (args.partition)):
				print("%s is found in this directory and will be overwritten!!" % (args.partition))
			with open(args.partition,"w") as partfile:#create partition file
				for fl in sorted(os.listdir(filelocation),key=numericalSort):
					if fl.endswith(".fasta") or fl.endswith(".fa"):
						with open("%s%s"%(filelocation,fl),"r+") as infile:
							for record in SeqIO.parse(infile,"fasta"):
								endpos = len(str(record.seq))
								seq = str(record.seq)
								newheader = record.description.split("_",1)[0]
								table[newheader] = table.get(newheader,'') + seq
						partitionName = fl.split(".")[0]						
						partfile.write(args.inputtype+", "+partitionName+"="+str(pos+1)+"-"+str(endpos+pos)+"\n")#write to partition file
						pos = endpos+pos
		else:
			print("ERROR: Input directory \"%s\" is NOT FOUND in current directory.\n"%(filelocation))
			sys.exit()

	if os.path.isfile("%s" % (args.MSAout)):
		print("%s is found in this directory and will be overwritten!!" % (args.MSAout))

	with open(args.MSAout,'w') as outfile:
		for header, seq in table.items():
			outfile.write(">"+header+'\n'+seq+'\n')

	if args.mutualinfo and args.inputtype == "nt":
		print("\nMutual Information and Concatenation done.")
		print("Output files %s can be found in current directory.\n"%args.MSAout)
	else:
		print("\nConcatenation done.")
		print("Output files %s can be found in current directory.\n"%args.MSAout)

'''
Tree Reconstruction
=>observe the evolutionary relationship between the samples
'''
def treeconstruction():
	print("\nSTEP 4: Tree Construction\n")
	if not args.inputtype:
		sys.exit("ERROR: Please select the sequence type, nt or aa.\n")

	if os.path.exists("%s"%(args.MSAout)):
		if os.path.isfile("%s" % (args.partition)):#if partition.txt is provided
			if args.codon:
				os.system("./dependencies/modeltest-ng-static -i %s -d nt -p %s -T raxml -o %s -q %s"%(args.MSAout, args.cpus, args.modeltestPrefix, args.partition))
			else:
				os.system("./dependencies/modeltest-ng-static -i %s -d %s -p %s -T raxml -o %s -q %s"%(args.MSAout,args.inputtype, args.cpus, args.modeltestPrefix, args.partition))
			if os.path.isfile("%s.out"%(args.modeltestPrefix)):	
				try:	
					os.system("./dependencies/raxml-ng %s --all --msa %s --model %s.part.%s --prefix %s --bs-trees %s" % (cmd3, args.MSAout, args.modeltestPrefix,args.subsModel, args.prefix, args.bootstrap))
				except:
					sys.exit("ERROR: Refer log file of RAxML-NG for more details.")
			elif args.inputtype == "nt":
				os.system("./dependencies/raxml-ng %s --all --msa %s --model GTR+G --prefix %s --bs-trees %s" % (cmd3, args.MSAout, args.prefix, args.bootstrap))
			elif args.inputtype == "aa":
				os.system("./dependencies/raxml-ng %s --all --msa %s --model LG+G8+F --prefix %s --bs-trees %s" % (cmd3, args.MSAout, args.prefix, args.bootstrap))	

		else:
			#os.system("./dependencies/modeltest-ng-static -i %s -d %s -p %s -T raxml -o %s -q %s"%(args.MSAout,args.inputtype, args.cpus, args.modeltestPrefix, args.partition))
			modelOutfile="%s.out" % (args.modeltestPrefix)
			if os.path.isfile(modelOutfile):	
				with open(modelOutfile,"r") as infile:
					text=infile.read()
				substitutions = re.findall(r"  > raxml-ng .*",text)

				if args.subsModel=="bic":
					os.system("./dependencies/raxml-ng %s --all --msa %s --model %s --prefix %s --bs-trees %s" % (cmd3, args.MSAout, substitutions[0].split()[-1],args.prefix, args.bootstrap))
				elif args.subsModel=="aic":
					os.system("./dependencies/raxml-ng %s --all --msa %s --model %s --prefix %s --bs-trees %s" % (cmd3, args.MSAout, substitutions[1].split()[-1],args.prefix, args.bootstrap))
				elif args.subsModel=="aicc":
					os.system("./dependencies/raxml-ng %s --all --msa %s --model %s --prefix %s --bs-trees %s" % (cmd3, args.MSAout, substitutions[2].split()[-1],args.prefix, args.bootstrap))
				
			elif args.inputtype == "nt":
				os.system("./dependencies/raxml-ng %s --all --msa %s --model GTR+G --prefix %s --bs-trees %s" % (cmd3, args.MSAout, args.prefix, args.bootstrap))
			elif args.inputtype == "aa":
				os.system("./dependencies/raxml-ng %s --all --msa %s --model LG+G8+F --prefix %s --bs-trees %s" % (cmd3, args.MSAout, args.prefix, args.bootstrap))
	
	else:
		sys.exit("ERROR: %s"%(args.MSAout), " is NOT FOUND in current directory ", os.getcwd())
    
	if not os.path.exists('./treeConstruction'):
		os.makedirs('./treeConstruction')

	os.system("mv %s* treeConstruction/" %(args.prefix))
	if os.path.exists('./%s.raxml.reduced.phy'%args.MSAout):
		os.system("mv ./%s.raxml.reduced.phy treeConstruction/" % args.MSAout)
	if os.path.exists('./%s.raxml.log' % args.MSAout):
		os.system("mv ./%s.raxml.log treeConstruction/" % args.MSAout)
	print("\nTree construction done.")
	print("Output files can be found in directory treeConstruction/.\n")

#program starts here
#step selection
if __name__ == "__main__":
	#reset the permission of mafft compiled files
	os.system("tar xzvf mafft-7.490-linux.tgz")
	os.system("chmod 777 ./dependencies/mafft-linux64/mafft.bat ./dependencies/mafft-linux64/mafftdir/bin/mafft ./dependencies/raxml-ng ./dependencies/MutualInfo.R ./dependencies/modeltest-ng-static")
	
	if args.step == 0 and args.codon: #all steps are chosen with codon alignment
		proteinextraction()
		concatenation()
		treeconstruction()

	elif args.step == 0:
		proteinextraction()
		msa()
		concatenation()
		treeconstruction()

	elif args.step == 1: #step 1: single core orthologous nucleotide / proteine extraction
		proteinextraction()

	elif args.step == 2: #step 2: multiple sequence alignment
		msa()

	elif args.step == 3: #step 3: mutual information (default OFF) & concatenation
		concatenation() 

	elif args.step == 4: #step 4: tree construction
		treeconstruction()

	else:
		sys.exit("ERROR: select the step to be done or leave the option as 0 (default) to run all. \n[0:all, 1:proteinortho, 2:msa (mafft), 3:mutual information(optional) & concatenation, 4:tree cosntruction (raxml-ng) ]")


