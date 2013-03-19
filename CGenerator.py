import random

def randomGenome(filename,length):
	generated_genome = ''
	generated_filename = 'reference_genomes/'
	generated_filename += filename
	generated_file = open(generated_filename,'w')
	for num in xrange(length):
		random_nucleotide = random.choice(['a','t','c','g'])
		generated_genome += random_nucleotide
	generated_file.write(generated_genome)

def mutate(gene):
	if gene == 'a':
		return random.choice(['t','c','g'])
	elif gene == 't':
		return random.choice(['a','c','g'])
	elif gene == 'c':
		return random.choice(['a','t','g'])
	elif gene == 'g':
		return random.choice(['a','c','t'])

def makeHaplotypes(filename):
	reference_genome_path = 'reference_genomes/'
	reference_genome_path += filename
	haplotype0_path = 'run_data/'
	haplotype0_path += filename
	haplotype1_path = haplotype0_path
	haplotype0_path += '.case0'
	haplotype1_path += '.case1'
	mutation_region_size = 10
	mutation_region_number = 3
	with open(reference_genome_path,'r') as reference_content:
		reference_genome = reference_content.read()
	haplotype0 = list(reference_genome)
	haplotype1 = list(reference_genome)
	reference_genome_length = len(reference_genome)
	for genome_index in range(0,reference_genome_length,mutation_region_size):
		for range_index in random.sample(range(10),mutation_region_number):
			haplotype0[genome_index+range_index] = mutate(haplotype0[genome_index+range_index])
			haplotype1[genome_index+range_index] = mutate(haplotype1[genome_index+range_index])
	haplotype0 = "".join(haplotype0)
	haplotype1 = "".join(haplotype1)
	haplotype0_file = open(haplotype0_path,'w')
	haplotype1_file = open(haplotype1_path,'w')
	haplotype0_file.write(haplotype0)
	haplotype1_file.write(haplotype1)

def makeReads(filename, read_length):
	haplotype0_path = 'run_data/'
	haplotype0_path += filename
	haplotype1_path = haplotype0_path
	haplotype0_path += '.case0'
	haplotype1_path += '.case1'
	read_path = 'run_data/'
	read_path += filename
	read_path += '.reads'
	with open(haplotype0_path,'r') as haplotype0_content:
		haplotype0 = haplotype0_content.read()
	with open(haplotype1_path,'r') as haplotype1_content:
		haplotype1 = haplotype1_content.read()
	haplotype_length = len(haplotype0)
	reads=[]
	for genome_index in range(0,haplotype_length,read_length):
		reads.append(haplotype0[genome_index:genome_index+read_length])
		reads.append(haplotype1[genome_index:genome_index+read_length])
	for genome_index in range(read_length/2,haplotype_length-read_length,read_length):
		reads.append(haplotype0[genome_index:genome_index+read_length])
		reads.append(haplotype1[genome_index:genome_index+read_length])
	random.shuffle(reads)
	read_file = open(read_path,'w')
	for read in reads:
		read_file.write(read)
		read_file.write('\n')
