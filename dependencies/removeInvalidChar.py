#!usr/bin/python

import multiprocessing
import logging
import os
import argparse
from Bio import SeqIO
import time
import sys
logging.basicConfig(filename='./log.log',level=logging.DEBUG)


parser = argparse.ArgumentParser()
parser.add_argument("-t", type=str, action="store", dest="inputtype", required=True,
			choices=["nt","aa"], help="input sequence type (nucleotide,nt or amino acid,aa)")
parser.add_argument("-d", type=str, action="store", dest="directory",
			default="inputfolder", help="working child directory")
parser.add_argument("-o", type=str, action ="store", dest="outfile",
			default="removedSeq.fa", help="file with removed sequence")
parser.add_argument("-c", type=int, action="store",default=-1, dest="cpus",help="number of CPUs [default: all available CPU]")
args = parser.parse_args()

def removeInvalidCharCPUs(filename):
	try:
		with open(args.outfile,"w") as removedSeq:	
			table={}
			num = 0
			with open("./%s/%s"%(args.directory, filename), 'r+', newline='') as infile:
				content=""
				removedSeq.write("\n"+filename+ "\n")
				for record in SeqIO.parse(infile,"fasta"):
					seq = str(record.seq).upper()
					header = record.description
					flag = True
					for char in seq:
						if char not in validchar:
							removedSeq.write(">"+header+"\n"+seq+"\n")
							flag = False
							num+=1
							break
					if flag:           
						content = content+ ">"+ header +"\n"+ seq +"\n"
		      
				infile.seek(0)
				infile.truncate()
				infile.write(content)
				infile.close()

			print("File ", filename," has ", num ," sequences being removed. ")
		print("\nYou may find these removed sequences in ", args.outfile)
	except:
		print(filename+ " Not successful")
		logging.warning(filename + " Not complete")

def removeInvalidChar():
	try:
		with open(args.outfile,"w") as removedSeq:
			for filename in os.listdir("./%s"% (args.directory)):	
				table={}
				num = 0
				
				with open("./%s/%s"%(args.directory, filename), 'r+', newline='') as infile:
					
					content=""
					removedSeq.write("\n"+filename+ "\n")
					for record in SeqIO.parse(infile,"fasta"):
						seq = str(record.seq).upper()
						header = record.description
						flag = True
						for char in seq:
							if char not in validchar:
								removedSeq.write(">"+header+"\n"+seq+"\n")
								flag = False
								num+=1
								break
						if flag:           
							content = content+ ">"+ header +"\n"+ seq +"\n"
			      
					infile.seek(0)
					infile.truncate()
					infile.write(content)
					infile.close()

				print("File ", filename," has ", num ," sequences being removed. ")
		print("\nYou may find these removed sequences in ", args.outfile)
	except:
		print("ERROR: Check your files.")

character=["ATCG","ACDEFGHIKLMNPQRSTVWY"]
validchar =""
if args.inputtype == "nt":
	validchar = character[0]
elif args.inputtype == "aa":
	validchar = character[1]
	
if __name__ == "__main__":

	cpu_available = os.cpu_count()
	if args.cpus > cpu_available:
		sys.exit("\nERROR: Maximum number of CPU available in your device is %s. \nERROR: You picked %s cpus!!" % (cpu_available,args.cpus))
	
	if os.path.exists("./%s"% (args.directory)):#if this child directory is present and has content
		if os.path.isfile("%s"%args.outfile):
			print("%s is found in this directory and will be overwritten!!" % (args.outfile))

		if args.cpus == -1: #all available cpu
			file_list = os.listdir("./%s"% (args.directory))	
			pool = multiprocessing.Pool(cpu_available)
			pool.map(removeInvalidCharCPUs, file_list)
		elif args.cpus <= 1:# single cpu
			removeInvalidChar()	
		elif args.cpus <= cpu_available:#determined by user
			file_list = os.listdir("./%s"% (args.directory))	
			pool = multiprocessing.Pool(args.cpus)
			pool.map(removeInvalidCharCPUs, file_list)

	else:
		print("ERROR: child directory ", args.directory, " is NOT FOUND / is EMPTY in current directory ", os.getcwd())
		sys.exit()
