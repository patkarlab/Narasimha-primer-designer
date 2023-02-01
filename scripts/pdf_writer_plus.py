#! /usr/bin/env python3

import csv
from fpdf import FPDF
import re
import sys

pdf = FPDF()
pdf.add_page()
pdf.set_font("Arial", size = 10)
intron_length = 100	# Assuming that the last 100 bases are from intron
sequences_fasta = sys.argv[1]
round1_primer_csv = sys.argv[2]
round2_primer_csv = sys.argv[3]
round1_primers = "round1_primers"
round2_primers = "round2_primers"
spumm="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
#dna_rna = 'RNA'

def reverse_compliment_dna (input_seq):
	compliment = {'A': 'T', 'T':'A', 'G':'C', 'C':'G'}
	rev_compliment_string = ''
	for bases in input_seq[::-1]:
		rev_compliment_string = rev_compliment_string + compliment[bases]

	return rev_compliment_string

fragment_names = dict()
gene_exons_list = []
seq_length_dict = dict()
fragment_dict = dict()
round1_primer_dict = dict()
region_length_dict = dict()
amp_seq = dict()
dna_rna_dict = dict()

with open(sequences_fasta, "r") as file:
	#csv_handle = csv.reader(file, delimiter = '')
	#header = next (csv_handle)	# Removing the header
	for lines in file:
		lines = re.sub ('\s+','',lines)
		if '>' in lines:
			#target_name = lines.split('>')
			lines = re.sub (">","", lines)
			#gene_exon = 
			gene_name = lines.split('_')[0].upper()
			exon_name = lines.split('_')[1].upper()
			#fragment_number = target_name[2]
			gene_exon = gene_name + '_'  + exon_name
		else:
			fragment_dict[gene_exon] = lines.upper() + spumm[::-1]
			seq_length_dict[gene_exon] = len(lines)

		#if gene_exon in gene_exons_list:
	#		fragment_dict[gene_exon] = lines[1].upper()
	#	else:
	#		gene_exons_list.append(gene_exon)
	#		print (gene_exon)
	
with open(round1_primer_csv, "r") as file:
	csv_handle = csv.reader(file)
	header = next (csv_handle)  # Removing the header
	for lines in csv_handle:
		target_name = lines[0]
		gene_name = target_name.split('_')[0]
		exon_name = target_name.split('_')[1]
		gene_exon = gene_name + '_'  + exon_name
		#target_seq = lines[1].upper()
		target_seq = fragment_dict[gene_exon]
		left_primer = lines[2].upper()
		right_primer = lines[3].upper()
		#fragment_names[target_name] = target_seq
		start_index = 0
		end_index = 0
		start_index_right_primer = len(target_seq) - 1
		end_index_right_primer = len(target_seq) - 1
		new_target_seq = ''
		sep_list = []
		rev_intron = ''
		intron = ''
		region_length = ''		
		ampli_seq = ''
		reg_start_index = 0
		reg_end_index = 0
		target_seq_del_spumm = ''
		amp_region = ''
		dna_rna = 'RNA'

		ampli_seq = re.sub (spumm[::-1],'',lines[1], flags = re.IGNORECASE)
		target_seq_del_spumm = re.sub (spumm[::-1], '', target_seq, flags = re.IGNORECASE)
		for reg_match in re.finditer(ampli_seq, target_seq_del_spumm):
			reg_start_index = reg_match.start()
			reg_end_index = reg_match.end()			

		#print (target_name,reg_start_index, reg_end_index)
		if '_' not in left_primer:
			for match in re.finditer(left_primer, target_seq):
				start_index = match.start()
				end_index = match.end()

		sep_list.append(start_index)
		sep_list.append(end_index)

		if '_' not in right_primer:
			reverse_comp_primer = reverse_compliment_dna(right_primer)
			for right_match in re.finditer(reverse_comp_primer, target_seq):
				start_index_right_primer = right_match.start()
				end_index_right_primer = right_match.end()
		else:
			reverse_comp_primer = ''
		
		sep_list.append(start_index_right_primer)
		sep_list.append(end_index_right_primer)

		for base_indices in range(len(target_seq)):
			sep = ''
			if base_indices in sep_list:
				sep = '*'
			new_target_seq = new_target_seq + sep + target_seq[base_indices]
		
		#intron = new_target_seq[-intron_length:-1].lower()
		#for all_bases in new_target_seq[::-1]:
		#	if all_bases.isalpha():
				#print (new_target_seq.index( new_target_seq [-1]))
		#		if new_target_seq.index(all_bases) < intron_length:
		#			print (all_bases, new_target_seq.index(all_bases))
		#			rev_intron = rev_intron + all_bases
		#		else: 
		#			break

		#intron = rev_intron[::-1]
		#spumm_exon = new_target_seq[0:len(target_seq) - 1 - intron_length]
		count = 0
		for all_bases in new_target_seq:
			temp_base = all_bases
			if all_bases.isalpha():
				if count < intron_length:
					temp_base =  all_bases.lower()
				count = count + 1

			intron = intron + temp_base

		if '_' not in left_primer:
			count2 = 0
			for base_pairs in new_target_seq:
				temp_base = base_pairs
				if count2 >= (reg_start_index) and count2 <= (reg_end_index):
					amp_region = amp_region + temp_base
				if base_pairs.isalpha():
					count2 = count2 + 1

			if reg_end_index <= intron_length:
				dna_rna = 'DNA'
		#spumm_exon = re.sub (intron,"", new_target_seq)
		#new_target_seq = spumm_exon + intron.lower()

		new_target_seq = intron
		string = target_name + '\t\t' + round1_primers + '\t\t' + str(reg_start_index) + '-' + str(reg_end_index) + '\t\t' + dna_rna  + '\t\t' + new_target_seq
		round1_primer_dict[target_name] = string
		seq_length_dict[target_name] = str(reg_start_index) + '-' + str(reg_end_index)
		amp_seq[target_name] = amp_region
		dna_rna_dict[target_name] = dna_rna
		#pdf.write(5, string)
		#pdf.write(5,'\n')
		#pdf.ln ()

with open(round2_primer_csv, "r") as file2:
	csv_handle = csv.reader(file2)
	csv_header = next (csv_handle)
	seq = ''
	for lines in csv_handle:
		target_name2 = lines[0]
		if seq != lines[1]:
			gene_name2 = target_name2.split('_')[0]
			exon_name2 = target_name2.split('_')[1]
			gene_exon2 = gene_name2 + '_'  + exon_name2
			target_seq2 = fragment_dict[gene_exon2]
			left_primer2 = lines[2].upper()
			right_primer2 = lines[3]
			start2_index = 0
			end2_index = 0
			start2_index_left_primer_sub = 0
			end2_index_left_primer_sub = 0
			start2_index_right_primer = len(target_seq2) - 1
			end2_index_right_primer = len(target_seq2) - 1
			new_target_seq = ''
			sep_list = []
			intron = ''
			#subset_seq = lines[1]
			subset_seq = re.sub ('\*', '', amp_seq[target_name2])
			start2_index_right_primer_sub = len(subset_seq) - 1
			end2_index_right_primer_sub = len(subset_seq) - 1
			subset_target_seq = ''

			# The right primer of round 2 has illumina motif attached written in small case
			right_primer_subset = ''
			for alphabets in right_primer2:
				if alphabets.isupper():
					right_primer_subset = right_primer_subset + alphabets
	
			#print (target_name2 , right_primer_subset, start2_index, end2_index, start2_index_right_primer, end2_index_right_primer)
			if '_' not in left_primer2:
				for match in re.finditer(left_primer2, target_seq2):
					start2_index = match.start()
					end2_index = match.end()

				for match2 in re.finditer(left_primer2, subset_seq):
					start2_index_left_primer_sub = match2.start()
					end2_index_left_primer_sub = match2.end()
	
			#print (target_name2 , right_primer_subset, start2_index, end2_index)
			sep_list.append(start2_index)
			sep_list.append(end2_index)

			if '_' not in right_primer2:
				reverse_comp_primer = reverse_compliment_dna(right_primer_subset)
				for right_match in re.finditer(reverse_comp_primer, target_seq2):
					start2_index_right_primer = right_match.start()
					end2_index_right_primer = right_match.end()

				#for right_match2 in re.finditer(reverse_comp_primer, subset_seq):
				#	start2_index_right_primer_sub = right_match2.start()
				#	end2_index_right_primer_sub = right_match2.end()		
			else:
				reverse_comp_primer = ''

			sep_list.append(start2_index_right_primer)
			sep_list.append(end2_index_right_primer)
			#print (target_name2 , right_primer_subset, start2_index, end2_index, start2_index_right_primer, end2_index_right_primer)
			for base_indices in range(len(target_seq2)):
				sep = ''
				if base_indices in sep_list:
					sep = '*'
				new_target_seq = new_target_seq + sep + target_seq2[base_indices]

			#intron = new_target_seq[-intron_length:-1].lower()
			#spumm_exon = new_target_seq[0:len(target_seq2) - 1 - intron_length]
			count = 0
			for all_bases in new_target_seq:
				temp_base = all_bases
				if all_bases.isalpha():
					if count < intron_length:
						temp_base =  all_bases.lower()
					count = count + 1
				intron = intron + temp_base

			new_target_seq = intron
			for indices in range(len(subset_seq)):
				sep = ''
				if indices == start2_index_left_primer_sub or indices == end2_index_left_primer_sub:
					sep = '*'
				subset_target_seq = subset_target_seq + sep + subset_seq[indices]	
			
			string = target_name2 + '\t\t' + round2_primers + '\t\t' + seq_length_dict[target_name2] + '\t\t' + dna_rna_dict[target_name2] + '\t\t' + new_target_seq
			pdf.write(5, round1_primer_dict[target_name2] + '\n')
			pdf.ln ()
			pdf.write(5, amp_seq[target_name2] + '\n')
			pdf.ln ()
			pdf.write(5, string + '\n')
			pdf.ln ()
			pdf.write(5, subset_target_seq + '\n')
			pdf.ln ()	
			
			seq = lines[1]
		else:
			pass

pdf.output ("test_plus.pdf")
