#!/usr/bin/python

from Bio import SeqIO
import os
import argparse
import sys
import time
import multiprocessing
import logging
logging.basicConfig(filename='./log.log',level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument('-d', type=str, action="store", dest="directory", 
			default="inputfolder", help="child directory name which contains the files where headers needed to be renamed (default: inputfolder)")
parser.add_argument('-s',type=str, action="store", dest="headerfile",
			default="",help="double-column space-separated file with multifasta filename and corresponding species name to be added, read the README file for more details")
parser.add_argument("-c", type=int, action="store",default=-1, dest="cpus",help="number of CPUs [default: all available CPU]")
parser.add_argument("-o",action="store_true", dest="inorder",help="extra tag at header for aa and nt corresponding files with different accession ID [default:OFF]")
args = parser.parse_args()

def singleCPU(filelist, headerlist):
	for filename, header in zip(filelist, headerlist):
		with open(filename,"r+") as fasta_file:
			count=0
			content=""
			print("Header [", header, "] added to file: ", fasta_file)

			for record in SeqIO.parse(fasta_file,"fasta"):
				seq = str(record.seq)
				if args.inorder:
					heading = (">%s" %header) +"_"+record.description.split(" ",1)[0]+"_"+str(count)+" "+ record.description.split(" ",1)[1]
				else:
					heading = (">%s" %header) +"_"+record.description
				content = content + heading+"\n"+seq+"\n"
				count+= 1
			fasta_file.seek(0)
			fasta_file.truncate()
			fasta_file.write(content)
			fasta_file.close()

def multiCPUs(filename, header):
	with open(filename,"r+") as fasta_file:
		count=0
		content=""
		print("Header [", header, "] added to file: ", fasta_file)

		for record in SeqIO.parse(fasta_file,"fasta"):
			seq = str(record.seq)
			if args.inorder:
				heading = (">%s" %header) +"_"+record.description.split(" ",1)[0]+"_"+str(count)+" "+ record.description.split(" ",1)[1]
			else:
				heading = (">%s" %header) +"_"+record.description
			content = content + heading+"\n"+seq+"\n"
			count+= 1
		fasta_file.seek(0)
		fasta_file.truncate()
		fasta_file.write(content)
		fasta_file.close()

if __name__ == "__main__":
	cpu_available= os.cpu_count()

	if args.cpus > cpu_available: #more than available
		sys.exit("\nERROR: Maximum number of CPU available in your device is %s. \nERROR: You picked %s cpus!!" % (cpu_available,args.cpus))
	if os.path.exists("./%s"% (args.directory)):#if this child directory is present and has content
		folderpath = r"./%s" %(args.directory) # make sure to put the 'r' in front
		file_list  = [os.path.join(folderpath, name) for name in os.listdir(folderpath)]
	else:
		sys.exit("ERROR: child directory %s is NOT FOUND / is EMPTY in current directory %s"%(args.directory, os.getcwd()))

	header_list=[]
	if args.headerfile == "":#no file with headername being provided
		name=""
		for i in file_list:
			print("Current file: ", i.split("/")[-1])
			name = input("Species name to be inserted in the header of this file: ")
			header_list.append(name)
	elif os.path.isfile("%s"%args.headerfile): #a header file is provided
		with open(args.headerfile,"r") as fl: #read from the headerfile
			filename = []
			for line in fl:
				element = line.split()#first element: input filename, second element: species name
				filename.append(element[0])
				header_list.append(element[1])

	if args.cpus == -1: # all available cpu	
		pool = multiprocessing.Pool(cpu_available)
		pool.starmap(multiCPUs, zip(file_list,header_list))
	elif args.cpus <= 1: # single cpu
		singleCPU(file_list, header_list)
	elif args.cpus <= cpu_available:#determined by user	
		pool = multiprocessing.Pool(args.cpus)
		pool.starmap(multiCPUs, zip(file_list,header_list))