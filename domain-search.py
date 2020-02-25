import csv
import re
import numpy as np
import pandas as pd

with open("mgen-3-140-s002.csv") as f:
    filter_df = pd.read_csv(f)

gb_path = "K56-annotated_ftsztag2_tn.gb"
csv_path = "1mhdtm_id.csv"
min_reads = 1	# Minimum number of reads to consider one half non-essential
min_weight = 0
nucleotides = ['a','c','t','g'] # Which characters represent nucleotides
min_ratio = 40000

'''
	parse_annotation

	Input: A single line of text from the .gb file which contains a gene annotation.
	Output: A dictionary containing:
		-The start position of the gene
		-The end position of the gene
		-A boolean, "complement", which is true if the gene is on the opposite strand

	Description: Checks if the line contains the word 'complement', which indicates that the gene
	is found in the opposite strand. Then, parses the line to find the start and end of the gene.
	If the line given was not actually an annotation, then split() will not return the correct
	number of "parts", and we return None to indicate that there was no gene.
'''
def parse_annotation(line):
	if 'complement' in line:
		complement = True
		bp_range = line[11:-1]	# This range turns complement(start..end) into start..end
	else:
		complement = False
		bp_range = line

	bp_range = bp_range.split('..')
	if len(bp_range) == 2:
		start = int(re.sub(r'[<>]', '', bp_range[0]))	# We remove the '>' or '<' which indicate partial features
		end = int(re.sub(r'[<>]', '', bp_range[1]))
		return {
			'start': start,
			'end': end,
			'complement': complement
		}
	else:
		return None


'''
	parse_file

	Input: A string, the path to the file we wish to parse.
	Output: The entire genome as a string, with one character representing each nucleotide,
		The annotations as a list of dictionaries, with one dictionary per gene.

	Description: Opens the file to read. Once we see "ORIGIN" we know we are in the genome and
	can directly append each character read to our genome string. Otherwise, we check if the
	line we are reading contains the word 'gene', which generally indicates an annotation in the
	.gb file format. If so, we pass the line to parse_annotation to parse the information about
	the gene annotation.
'''
def parse_file(filename):
	genome = ''
	annotations = []
	with open(filename) as f:
		in_genome = False
		for line in f:
			line = line.strip()
			if line == "ORIGIN":
				in_genome = True
			elif in_genome:
				genome = genome + ''.join([c for c in line if c in nucleotides])
			else:
				line = line.split()
				if line[0] == 'gene':
					annotation = parse_annotation(line[1])
					if annotation is not None:
						annotations.append(annotation)
								
	return genome, annotations

'''
	load_csv
	
	Input: The filename of the csv file to load.
	Output: A list of dictionaries.

	Description: Opens the target csv file, and stores the position and count of reads in
	a dictionary with corresponding keys. Output is a list of each such dictionary.
'''
def load_csv(filename):
	rows = []
	with open(filename) as f:
		reader = csv.DictReader(f)
		for row in reader:
                    rows.append({'position': int(row['position']), 'count': int(row['count']), 'strand': row['strand']})
	return rows
			

'''
	create_empty_read_vector

	Input: The genome nucleotide sequence.
	Output: A 1D zero vector with dimension equal to the length of the genome.

	Description: Create an empty vector which can be used to record the hit count at each
	position in the genome.
'''
def create_empty_read_vector(genome):
	vec_size = len(genome)
	pos_vec = np.zeros((vec_size), dtype = np.int64)
	neg_vec = np.zeros((vec_size), dtype = np.int64)
	return pos_vec, neg_vec


'''
	fill_read_vector

	Input: A 1D zero vector with length equal to the genome, a list of dictionaries with positional counts
	Output: A 1D vector with read counts at each position in the genome.

	Description: Adds all of the "counts" in the count list to the corresponding position in the input vector.
'''
def fill_read_vector(pos_vec, neg_vec, counts):
	for c in counts:
		if c['strand'] == '1':
			pos_vec[c['position']] += c['count']
		elif c['strand'] == '-1':
			neg_vec[c['position']] += c['count']
		else:
			print(f"Unknown strand! <{c['strand']}>")

	return pos_vec, neg_vec


'''
	get_skewed_genes

	Input: A list of annotation dictionaries with gene start and ends, a vector of read counts
	Output: A list of genes with reads in only one side of the gene.

	Description: For each annotation, the start/middle/end are calculated of 80% of the gene's interior.
	Then, the read counts of the "left" and "right" side are summed. Any gene containing reads in exactly
	one half of the gene is added to the list, as long as the read count in this half exceeds some threshold.

	The threshold has a default set at the top of the file, but may be set to a value based on the mean/std dev.
'''
def get_skewed_genes(annotations, pos_vec, neg_vec):
	candidates = []
	for a in annotations:
		length = a['end'] - a['start']
		offset = length * 0.1
		start = int(a['start'] + offset)
		end = int(a['end'] - offset)

                # search 3 midpoints
		half1 = int((end + start) * 0.5) - int((end-start) * 0.1)
		half2 = int((end + start) * 0.5) + int((end-start) * 0.1)
		half3 = int((end + start) * 0.5)
		splits = [half3]

		for half in splits:
			# check strand 
			if (a['complement']):
				left_reads = np.sum(neg_vec[start:half])
				right_reads = np.sum(neg_vec[half:end])
			else:
				left_reads = np.sum(pos_vec[start:half])
				right_reads = np.sum(pos_vec[half:end])

			min_reads_r = (end - half) * 0.14
			min_reads_l = (half - start) * 0.14

			if a['start'] == 7133330 or a['start'] == 7354153:
				print("START: ", a['start'])
				print("LEFT: ", left_reads)
				print("MINL: ", min_reads_l)
				print("RIGHT: ", right_reads)
				print("MINR: ", min_reads_r)
			if left_reads != right_reads and (left_reads == 0 or (right_reads / left_reads) >= min_ratio) and half != half2:
				# final check: make sure we have a significant number of reads on the non empty end
				if right_reads > min_reads_r:
					candidates.append({'start': a['start'], 'end': a['end'],'left_reads': left_reads, 'right_reads': right_reads, 'strand': a['complement']})
					break
			if left_reads != right_reads and (right_reads == 0 or (left_reads / right_reads) >= min_ratio) and half != half1:
				if left_reads > min_reads_l:
					candidates.append({'start': a['start'], 'end': a['end'],'left_reads': left_reads, 'right_reads': right_reads, 'strand': a['complement']})
					break

	return candidates


def calc_mean_reads(annotations, vec):
	read_counts = []
	for a in annotations:
		length = a['end'] - a['start']
		offset = length * 0.1
		start = int(a['start'] + offset)
		end = int(a['end'] - offset)
		reads = np.sum(vec[start:end])
		weighted_reads = reads / length
		reads = weighted_reads
		read_counts.append(reads)

	read_counts.sort()
	read_counts = np.array(read_counts)
	mean = np.mean(read_counts)
	median = np.median(read_counts)
	std = np.std(read_counts)
	print("Maximum reads: ", np.amax(read_counts))
	print("Minimum reads: ", np.amin(read_counts))
	print("Mean reads: ", mean)
	print("Median reads: ", median)
	print("Standard deviation: ", std)
	print("80%: ", np.mean(read_counts) * 0.8)

	# plot reads
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = plt.axes()
	x = range(len(read_counts))
	y = read_counts
	ax.plot(x, y)
	plt.show()
	return (mean, median, std)


def parse_locus_line(line):
	return line.strip()[12:-1]


def parse_sequence_line(line):
	return line.strip()[9:-1]


def add_tags(gene):
	target = str(gene['start']) + '..' + str(gene['end'])
	with open(gb_path) as f:
		for line in f:
			if ('CDS' in line or 'RNA' in line) and '..' in line:
				start = line.split()[1].split('..')[0]
				end = line.split()[1].split('..')[1]
				if 'complement' in start:
					start = start[11:]
					end = end[:-1]
					comp = True
				else:
					comp = False

				start = re.sub('<', '', start)
				start = re.sub('>', '', start)
				end = re.sub('<', '', end)
				end = re.sub('>', '', end)
				start = int(start)
				end = int(end)
				if gene['start'] == start and gene['end'] == end:
					# if the strands do not match
					if comp != gene['strand']:
						print("MISMATCHED STRANDS: ", gene['start'], " : ", gene['strand'])
					# tRNA have a few extra lines before the locus tag
					if 'tRNA' in line:
						f.readline()
						f.readline()
						f.readline()
						f.readline()
						gene['ref_seq'] = 'None'
					# parse the locus tag from the next line
					locus_line = f.readline()
					# skip inference line
					f.readline()
					# get the sequence line
					sequence_line = f.readline()
					gene['locus_tag'] = parse_locus_line(locus_line)
					if not 'ref_seq' in gene:
						gene['ref_seq'] = parse_sequence_line(sequence_line)
					if '=' in gene['ref_seq']:
						gene['ref_seq'] = 'None'
					return(gene)
		print("unable to find tag for gene")
		print(target)
		gene['locus_tag'] = 'Unknown'
		gene['ref_seq'] = 'Unknown'
				
		
# Actual program starts here

print("Loading genome and annotations...")
genome, annotations = parse_file(gb_path)
print("Loading insertion site counts...")
counts = load_csv(csv_path)
print("Create the empty vector of positional reads...")
pos_read_vec, neg_read_vec = create_empty_read_vector(genome)
print("Fill vector with read counts...")
pos_read_vec, neg_read_vec = fill_read_vector(pos_read_vec, neg_read_vec, counts)

print("Look for genes with essential domains...")
candidates = get_skewed_genes(annotations, pos_read_vec, neg_read_vec)
print("Adding locus tag information to all candidates...") # This should just be done in the initial parse but I got lazy
for c in candidates:
	c = add_tags(c)
print("Filtering candidates based on locus tag and refseq...")
candidates = [c for c in candidates if 'WQ' in c['locus_tag'] or 'WP' in c['ref_seq']]
print("Found ", len(candidates), " candidates before filtering.")
tags = [tag for tag in filter_df['Locus tag']] + [tag for tag in filter_df['Locus tag2']]
candidates = [c for c in candidates if c['locus_tag'] in tags]
print("Found ", len(candidates), " candidates after filtering.")
print("Saving results to csv file...")
with open("results.csv", 'w') as f:
	writer = csv.DictWriter(f, ['start', 'end', 'strand', 'ref_seq', 'locus_tag', 'left_reads', 'right_reads'])
	writer.writeheader()
	for data in candidates:
		writer.writerow(data)



