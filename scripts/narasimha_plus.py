#! /usr/bin/env python3
##############################################################################################################################
#																															 #
#      NARASIMHA Primer Designer																					 		 #
#																															 #
#																															 #
# Author: Nidhi Koundinya																									 #
# Date Created: April 30, 2021																								 #
##############################################################################################################################
# Written for Python 3.6
#
# Usage:
# $ conda env create -f ./scripts/environment.yaml
# $ conda activate narasimha-primers
# $ python3 run.py -i <path to input BED or fasta file> -g <path to reference genome>
#
# Arguments:
# 
#		-h, --help                  show this help message and exit
#		-i, --input                 Enter the path to the BED or FASTA file
#		-g, --genome                Enter the path to the reference genome
#
#
# Inputs:
# 1. A 6 column BED file containing coordinates of regions to be targeted
# 2. A fasta file of reference genome
#
# Outputs:
# 1. A directory "primer_out" at current location with directory "final" inside.
# 2. A csv file "all_primers.csv" inside the final directory with a list of primers for both rounds.
# 3. "round1primers.csv" and "round2primers.csv" inside the final directory with a list of the template,
#	  forward and reverse primers for round1 and round2 respectively.
#
#
############################################################################################################################


import os
import re
import csv
import sys
import math
import argparse
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

primer3_path="primer3_core"
bedtools_path="bedtools"

new_path = os.path.join(os.getcwd() + "/primer_out_plus")
temp_path = os.path.join(new_path + "/temp")
final_path = os.path.join(new_path + "/final")
seqs_path = os.path.join(temp_path + "/sequences.fasta")


#Deletes any folder with primers named "primers_out"
if os.path.isdir(new_path):
	rm_cmd = ["rm","-rf",new_path]
	subprocess.run(rm_cmd)

os.mkdir(new_path)
os.mkdir(temp_path)
os.mkdir(final_path)
fragment_size_max=150
spumm="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
number_of_regions=3

#Deletes temp file from "primer_out"
def delete_temp():
	if os.path.isdir(temp_path):
		rm_cmd = ["rm","-rf",temp_path]
		subprocess.run(rm_cmd)


#parsing command line arguments
def command_Parse():
	parser = argparse.ArgumentParser(prog='python3 narasimha.py', usage='%(prog)s -i <path to input BED or fasta file> -g <path to reference genome> -l <length of amplicon>',
    		description='Generates primers for NARASIMHA Assay: https://www.nature.com/articles/s41408-020-0313-6', epilog="***** NARASIMHA - PRIMER DESIGNER *****")
	def check_input(file):
		base, ext = os.path.splitext(file)
		if not os.path.exists(file):
			raise argparse.ArgumentTypeError('{} not found. Check path.'.format(file))
		elif ext.lower() not in ('.bed', '.fasta'):
			raise argparse.ArgumentTypeError('\nInput must be a BED file or a FASTA file')
		return file

	def check_ref(file):
		base, ext = os.path.splitext(file)
		if not os.path.exists(file):
			raise argparse.ArgumentTypeError('{} not found. Check path.'.format(file))
		elif ext.lower() not in ('.fasta'):
			raise argparse.ArgumentTypeError('\nInput must be a FASTA file.')
		return file
	
	parser.add_argument("-i", "--input", help='Enter the path to the BED or FASTA file', type=check_input , required=True)
	parser.add_argument("-g", "--genome", help='Enter the path to the reference genome', type=check_ref , required=True)
	parser.add_argument("-l", "--amp_length", help='Enter the length of amplicon required', type=int, required=True)
	return parser


###########################################################################
#Creating BEDfile for analysis after primers for both rounds are generated#
###########################################################################
def creating_bedfile(bed_input,input_bedfile):
	if bed_input==True:
		record_id=[]
		record_sequence=[]
		for record in SeqIO.parse(f'{temp_path}/ispcr_sequences1.fasta',"fasta"):
				record_id.append(record.id)
				record_sequence.append(str(record.seq))
		df1=pd.DataFrame({'gene':record_id,'fragment':record_sequence})
		df2 = pd.read_csv(input_bedfile,delimiter="\t",header=None)
		gene_list=df2[3]
		genes=[genes.upper() for genes in gene_list]
		g_index=[]
		for i,n in enumerate(record_id):	
			#Remove '_fragmentnumber'
			match=re.search("_\d+$",n).group()
			name = n.replace(match, "")
			g_index.append(genes.index(name))
		df3 = df2.iloc[g_index,:]
		
		df3 = df3.reset_index(drop=True)
		df_final=pd.concat([df1,df3],axis=1,ignore_index=True)
		spumm_length=len(spumm)
		
		exons=df_final[5]
		exon_list=[exo.upper() for exo in exons]
		for i,exons in enumerate(exon_list):
			if i<(len(exon_list)-1):
				if exon_list[i]==exon_list[i+1]:
					strand=df_final.loc[i,7]
					if strand=='+':
						frag_length=len(df_final.loc[i,1])-spumm_length
						df_final.loc[i,4]=df_final.loc[i,3]+frag_length #end=start+length_of_fragment
						df_final.loc[i+1,3]=df_final.loc[i,3]+frag_length #start of next fragment=end of previous
					elif strand=='-':
						frag_length=len(df_final.loc[i,1])-spumm_length
						df_final.loc[i,3]=df_final.loc[i,4]-frag_length
						df_final.loc[i+1,4]=df_final.loc[i,4]-frag_length
		
		df_primers = pd.read_csv(f'{temp_path}/ispcr_query2.txt',delimiter="\t",header=None)

		names = df_final.iloc[:,0].values
		df_new=pd.DataFrame(index=names)
		exons=df_final[0].tolist()
		names2=df_primers[0].tolist()

		fragments=df_final[1].tolist()
		chrom=df_final[2].tolist()
		start=df_final[3].tolist()
		end=df_final[4].tolist()
		genes=df_final[5].tolist()
		strand=df_final[7].tolist()
		primers=df_primers[2].tolist()

		df_new['Fragment']=pd.Series(fragments,index=exons)
		df_new['chr']=pd.Series(chrom,index=exons)
		df_new['start']=pd.Series(start,index=exons)
		df_new['end']=pd.Series(end,index=exons)
		df_new['genes']=pd.Series(genes,index=exons)
		df_new['strand']=pd.Series(strand,index=exons)
		df_new['primers']=pd.Series(primers,index=names2)

		df_new.fillna(0, inplace=True)
		df_new = df_new.reset_index(drop=True)
		df_new['genes']=names
		
		for i in range(len(df_new)):
			primer=df_new.iloc[i,6]
			if primer==0:
				next
			else:
				primer=primer.split("gtctcgtgggctcggagatgtgtataagagacag")[1]	#illumina_tag
				primerSeq=Seq(primer)
				reverse_comp_primer=str(primerSeq.reverse_complement())
				frag=df_new.iloc[i,0]
				m=re.search(reverse_comp_primer,frag)
				primer_end=m.end()
				strand=df_new.iloc[i,5]
				if strand=='+':
					extra_len=len(frag[m.end():])
					df_new.iloc[i,3]=df_new.iloc[i,3]-extra_len
				if strand=='-':
					extra_len=len(frag[m.end():])
					df_new.iloc[i,2]=df_new.iloc[i,2]+extra_len
		
		df_done=df_new.drop(['Fragment','strand','primers'],axis=1)
		df_done.to_csv(f'{final_path}/bedfile.bed',sep='\t',index=False,header=False)

	else:
		print("Cannot generate bedfile for analysis if input was a fasta file!")

###################
##Running Primer3##
###################
def runPrimer3(pcrRound):
	input_file = temp_path + "/inputs" + pcrRound + ".txt"
	
	#removing extra newline at the end of input file
	truncInput_cmd = ["truncate", "-s", "-1", input_file]
	subprocess.run(truncInput_cmd)
	
	#Running primer3
	output_file = "--output=" + temp_path + "/primerOutput" + pcrRound
	primer3_cmd = [primer3_path, input_file , output_file]
	subprocess.run(primer3_cmd)

#################
##Running isPcr##
#################
def runIsPcr(pcrRound):
	ispcr_sequences =temp_path+"/ispcr_sequences"+pcrRound+".fasta"
	ispcr_query=temp_path+"/ispcr_query"+pcrRound+".txt"
	pcrSequences=temp_path+"/pcrSequences_"+pcrRound+".fasta"
	ispcr_cmd = ["isPcr", ispcr_sequences, ispcr_query, pcrSequences]
	subprocess.run(ispcr_cmd)
	
	pattern="s/:.*//"
	cmd=["sed", "-i", pattern, pcrSequences]
	subprocess.run(cmd)
	pcrRound='2'
	return pcrRound

#########################
##Creating Final Output##
#########################
def compiling_primers():
	df2 = pd.read_csv(f'{final_path}/round1primers.csv',delimiter=",")
	names = df2.iloc[:,0].values
	df2=pd.DataFrame(index=names)
	
	seq_id_1=[]
	r1_left_primer=[]
	r1_right_primer=[]
	seq_id_2=[]
	r2_left_primer=[]
	r2_right_primer=[]
	
	with open(f'{temp_path}/ispcr_query1.txt','r') as primers:
		reader=csv.reader(primers,delimiter="\t")
		for word in reader:
			seq_id_1.append(word[0].upper())
			r1_left_primer.append(word[1])
			r1_right_primer.append(word[2])
	with open(f'{temp_path}/ispcr_query2.txt','r') as primers:
		reader = csv.reader(primers, delimiter='\t')
		for word in reader:
			seq_id_2.append(word[0].upper())
			r2_left_primer.append(word[1])
			r2_right_primer.append(word[2])
			
	df2['round1_Forward']=pd.Series(r1_left_primer,index=seq_id_1)
	df2['round1_Reverse']=pd.Series(r1_right_primer,index=seq_id_1)
	df2['round2_Forward']=pd.Series(r2_left_primer,index=seq_id_2)
	df2['round2_Reverse']=pd.Series(r2_right_primer,index=seq_id_2)
	
	
	df2.to_csv(f'{final_path}/all_primers.csv')

##########################
##Parsing primer3 output##
##########################
def parseP3_Output(pcrRound):
	def generate_p3log(line,counter):
		with open(f'{final_path}/primer3_log.txt','a') as log_file:
			if re.search('(?<=PRIMER_LEFT_EXPLAIN=).*$',line):
				log_file.write(f'{seq_ids[counter]}\n')
				counter+=1
				left=re.findall('(?<=PRIMER_LEFT_EXPLAIN=).*$',line)[0]
				log_file.write(f'Left Primers: {left}\n')
			if re.search('(?<=PRIMER_RIGHT_EXPLAIN=).*$',line):
				right=re.findall('(?<=PRIMER_RIGHT_EXPLAIN=).*$',line)[0]
				log_file.write(f'Right Primers: {right}\n')
			if re.search('(?<=PRIMER_PAIR_EXPLAIN=).*$',line):
				pair=re.findall('(?<=PRIMER_PAIR_EXPLAIN=).*$',line)[0]
				log_file.write(f'Primer Pairs: {pair}\n\n')
			if re.search(f'(?<=PRIMER_ERROR=).*$',line):
				log_file.write(f'{seq_ids[counter]}\n')
				err=re.findall('(?<=PRIMER_ERROR=).*$',line)[0]
				log_file.write(f'Error: {err}\n\n')
				counter+=1
	def parseIdAndTemplate():
		with open(f'{temp_path}/primerOutput{pcrRound}','r') as p_out:
			for line in p_out:
				if re.search('SEQUENCE_ID=',line):
					m=re.findall('(?<=SEQUENCE_ID=).*$',line)[0]
					seq_ids.append(m)
				if re.search('SEQUENCE_TEMPLATE=',line):
					m=re.findall('(?<=SEQUENCE_TEMPLATE=).*$',line)
					templates.append(m[0])
	def parsePrimers():
		with open(f'{temp_path}/primerOutput{pcrRound}','r') as p_out:
			for line in p_out:
				if re.search('SEQUENCE_ID=',line):
					next
				generate_p3log(line,counter)
				if re.search('PRIMER_PAIR_NUM_RETURNED=0',line):
					left_primers.append("NO_PRIMER_RETURNED")
					right_primers.append("NO_PRIMER_RETURNED")
					error.append("NO PRIMER FOUND")
					continue
				else:
					if re.search('(?<=PRIMER_LEFT_0_SEQUENCE=).*$',line):
						seq=re.findall('(?<=PRIMER_LEFT_0_SEQUENCE=).*$',line)[0]
						if pcrRound=='1':
							left_primers.append(seq)
						elif pcrRound=='2':
							left_primers.append(illumina_tag.lower()+seq)							
						error.append('')
						
					if re.search('(?<=PRIMER_RIGHT_0_SEQUENCE=).*$',line):
						seq=re.findall('(?<=PRIMER_RIGHT_0_SEQUENCE=).*$',line)[0]
						if pcrRound=='1':
							right_primers.append(seq)
						elif pcrRound=='2':
							#right_primers.append(illumina_tag.lower()+seq)
							right_primers.append(seq)
				if re.search('PRIMER_ERROR=',line):
					left_primers.append("NO_PRIMER_RETURNED")
					right_primers.append("NO_PRIMER_RETURNED")
					error.append(re.findall('(?<=PRIMER_ERROR=).*$',line)[0])
					continue
	def create_ispcrQuery():
		with open(f'{temp_path}/ispcr_query{pcrRound}.txt','w') as text_file:
			for i in range(0,len(seq_ids)):
				if (right_primers[i]=="NO_PRIMER_RETURNED"):
					continue
				text_file.write(f'{seq_ids[i]}\t{left_primers[i]}\t{right_primers[i]}\n')
	def createOutputFile():
		primers_out = {}
		primers_out['Sequence_ID'] = seq_ids
		primers_out['Template'] = templates
		primers_out['Left_Primers'] = left_primers
		primers_out['Right_Primers'] = right_primers
		primers_out['Errors'] = error
		data = list(zip(primers_out['Sequence_ID'],primers_out['Template'],primers_out['Left_Primers'],primers_out['Right_Primers'],primers_out['Errors']))
		df = pd.DataFrame(data=data)
		df.to_csv(f'{final_path}/round{pcrRound}primers.csv', index=False, header=fieldnames)
		print(f'round{pcrRound}primers.csv created in FINAL directory')

	seq_ids=[]
	templates=[]
	left_primers=[]
	right_primers=[]
	err_seq=[]
	err_temps=[]
	error=[]
	illumina_tag="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
	fieldnames=["Sequence_ID","Template","Left_Primers","Right_Primers","Errors"]
	parseIdAndTemplate()
	with open(f'{final_path}/primer3_log.txt','a') as log_file:
		log_file.write(f'ROUND {pcrRound}\n\n')
		counter=0
	parsePrimers()
	createOutputFile()
	create_ispcrQuery()


class DesignPrimer:
	def __init__(self, pcrRound, args):
		self.input=args.input
		self.genome=args.genome
		self.pcrRound=pcrRound
		self.bed_input=False #Default
		self.amp_length=args.amp_length
	
	def getSequences(self):
		base, ext = os.path.splitext(self.input)
		self.ext=ext
		if ext=='.bed':
			self.bed_input=True
			bed_cmd = [bedtools_path, "getfasta", "-fi", self.genome, "-bed", self.input, "-fo", seqs_path, "-name", "-s"]
			subprocess.run(bed_cmd)
			for record in SeqIO.parse(seqs_path,"fasta"):
				record_id = record.id.split('::')[0]
				cmd="sed -i 's/>"+record.id+"/>"+record_id+"/' " + seqs_path
				os.system(cmd)
		elif ext=='.fasta':
			fasta_cmd=["cp", self.input, seqs_path]
			subprocess.run(fasta_cmd)
		return self.bed_input
	
	def fragment_sequences(self,record_sequence):
		fragment_size=201
		fragments=[]
		i=1
		pos=0
		intron_length = 100     # Assumption that intronic region at the start of each fragment is 100 bases
		seq_length = len(record_sequence)
		#while (fragment_size > fragment_size_max):
		#	fragment_size = math.ceil(seq_length/i)
		#	i+=(number_of_regions - 1)
	#		if (fragment_size < fragment_size_max):
	#			break

		if seq_length <= self.amp_length:
			fragment_size = seq_length
		else:
			fragment_size = self.amp_length

		fragments.append(record_sequence[pos:intron_length])
		pos = intron_length
		while pos < seq_length:
			lower = pos
			upper = pos + fragment_size
			if upper >= seq_length:
				upper = seq_length

			fragments.append(record_sequence[lower:upper])
			pos += fragment_size

		#fragments.append(record_sequence[upper:seq_length])
		return fragments
	
	#########################################
	##Calculating target and product length##
	#########################################
	def make_p3Settings(self):
		if self.pcrRound=='1':
			for record in SeqIO.parse(seqs_path,"fasta"):
				record_id = record.id
				record_sequence = str(record.seq)
				fragments = DesignPrimer.fragment_sequences(self,record_sequence)
		
				for index,seq in enumerate(fragments):
					#print (record_id,seq)
					header = '{}_{}'.format(record_id,index+1)
					primer_values = ["PRIMER_TASK=generic","PRIMER_PICK_LEFT_PRIMER=1","PRIMER_PICK_INTERNAL_OLIGO=0","PRIMER_PICK_RIGHT_PRIMER=1",
										"PRIMER_NUM_RETURN=1","PRIMER_OPT_SIZE=20","PRIMER_MIN_SIZE=16","PRIMER_MAX_SIZE=27","PRIMER_EXPLAIN_FLAG=1"]
					if (len(seq)<27):
						print(f'{header} is too short.\n')

					targ=30
					length=len(seq)
					prod_min_size=27		#min should be >= max primer size
					seq_start = 27
					prod_max_size=210
					a=int(.5*(length))
					#a = prod_max_size - 10
					length_plus_spumm = length + len(spumm) - 1
					seq_length = length - seq_start
					#print (header, length, a)
					#for longer fragments
					if (len(seq)>61):
						target=f'SEQUENCE_TARGET={seq_start},{seq_length}'
						product_range=f'PRIMER_PRODUCT_SIZE_RANGE={a}-{prod_max_size} {prod_min_size}-{prod_max_size}'
						#product_range=f'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,,{length_plus_spumm},'
						primer_values.append(product_range)
						primer_values.append(target)
						primer_values.append("=")
						
					#for short fragments
					else:
						target=f'SEQUENCE_TARGET={seq_start},{seq_length}'
						product_range=f'PRIMER_PRODUCT_SIZE_RANGE={prod_min_size}-{prod_max_size}'
						#product_range=f'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,,{length_plus_spumm},'
						primer_values.append(product_range)
						primer_values.append(target)
						primer_values.append("=")
					self.make_p3Input(header,seq,primer_values)
		elif self.pcrRound=='2':
			for record in SeqIO.parse(f'{temp_path}/pcrSequences_1.fasta',"fasta"):
				name=str(record.id)
				sequence=str(record.seq.upper())
				with open(f'{temp_path}/inputs2.txt','a') as file:
					file.write(f'SEQUENCE_ID={name}\nSEQUENCE_TEMPLATE={sequence.upper()}\n')
					primer_values = ["PRIMER_TASK=generic","PRIMER_PICK_LEFT_PRIMER=1","PRIMER_PICK_INTERNAL_OLIGO=0","PRIMER_PICK_RIGHT_PRIMER=1",
											"PRIMER_NUM_RETURN=1","PRIMER_OPT_SIZE=20","PRIMER_MIN_SIZE=16","PRIMER_MAX_SIZE=27","PRIMER_EXPLAIN_FLAG=1"]
					
					spumm_inv = "TTCTCGCAGCACATCCCTTT"
					m=re.search(spumm_inv, sequence)
					targ_start=m.end()
					actual_len=len(sequence[m.end():])
					#print (actual_len)

					seq_start = 27
					length = len(sequence)
					seq_length = length - seq_start - len(spumm_inv)

					if (length > 45):
						target=f'SEQUENCE_TARGET={seq_start},{seq_length}'
						#Excluding last 5 bases to avoid generating same right primer in round 2
						exc = seq_start + 5
						excluded=f'SEQUENCE_EXCLUDED_REGION=0,5'
					else:
						target=f'SEQUENCE_TARGET={seq_start},{seq_length}'
						#Excluding last 2 bases to avoid generating same right primer in round 2
						exc = seq_start + 2
						excluded=f'SEQUENCE_EXCLUDED_REGION=0,2'
					primer_values.append(excluded)
					primer_values.append(target)
					product_range=f'PRIMER_PRODUCT_SIZE_RANGE=27-200'
					primer_values.append(product_range)
					
					primer_values.append("=")
					for i,items in enumerate(primer_values):
						file.write(items)
						if i==len(primer_values)-1:
							break
						file.write('\n')
					file.write('\n')
				with open(f'{temp_path}/ispcr_sequences2.fasta','a') as txt_file2:
					txt_file2.write(f'>{name}\n{sequence}\n')
				
	def make_p3Input(self,header,seq,primer_values):
		with open(f'{temp_path}/inputs{self.pcrRound}.txt','a') as txt_file:
			txt_file.write('SEQUENCE_ID={}\nSEQUENCE_TEMPLATE={}{}\n'.format(header.upper(),seq.upper(),spumm.lower()[::-1]))
			for i,items in enumerate(primer_values):
				txt_file.write(items)
				if i==len(primer_values)-1:
					break
				txt_file.write('\n')
			txt_file.write('\n')
		with open(f'{temp_path}/ispcr_sequences1.fasta','a') as txt_file2:
			txt_file2.write(f'>{header.upper()}\n{seq.upper()}{spumm.lower()[::-1]}\n')

def main(pcrRound='1'):
	if pcrRound=='1':
		#parse cmd args
		parser = command_Parse()
		args = parser.parse_args()
		
		print("Starting....")
		input_bedfile=args.input
		dp_1 = DesignPrimer(pcrRound, args)
		bed_input=dp_1.getSequences()
		print("Creating Input File for Primer3")
		dp_1.make_p3Settings()
		print("Input File Created")

		runPrimer3(pcrRound)
		parseP3_Output(pcrRound)
		pcrRound=runIsPcr(pcrRound)
		
	if pcrRound=='2':
		args.input="pcrSequences.fasta"
		dp_2 = DesignPrimer(pcrRound, args)
		dp_2.make_p3Settings()
		runPrimer3(pcrRound)
		parseP3_Output(pcrRound)
		compiling_primers()
		creating_bedfile(bed_input,input_bedfile)
		delete_temp()
		print("Done")

if __name__ == "__main__":
	main()
