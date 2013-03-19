import time
import json
from collections import defaultdict
from itertools import izip, islice
def unique_sorted(values):
    "Return a sorted list of the given values, without duplicates."
    values = sorted(values)
    if not values:
        return []
    consecutive_pairs = izip(values, islice(values, 1, len(values)))
    result = [a for (a, b) in consecutive_pairs if a != b]
    result.append(values[-1])
    return result

def getReads(filename):
	read_path = 'run_data/'
	read_path += filename
	read_path += '.reads'
	with open(read_path,'r') as read_content:
		reads = read_content.read()
	return reads.splitlines()

def getReference(filename):
	reference_genome_path = 'reference_genomes/'
	reference_genome_path += filename
	with open(reference_genome_path,'r') as reference_content:
		reference_genome = reference_content.read()
	return reference_genome

def compareStrings(reference,read):
	if reference == read:
		return 0
	else:
		length = len(read)
		reference = list(reference)
		read = list(read)
		mismatches = 0
		for index in range(length):
			if read[index] != reference[index]:
				mismatches += 1
		return mismatches

def alignReadOld(filename,read):
	reference_genome_path = 'reference_genomes/'
	reference_genome_path += filename
	with open(reference_genome_path,'r') as reference_content:
		reference_genome = reference_content.read()
	reference_length = len(reference_genome)
	read_length = len(read)
	best_match_posistion = -1
	best_match_errors = read_length + 1
	for index in range(reference_length-read_length):
		num_mismatches = compareStrings(reference_genome[index:index+read_length],read)
		if num_mismatches < best_match_errors:
			best_match_index = index
			best_match_errors = num_mismatches
	return best_match_index

def prettySNPS(str1,str2):
	str1 = list(str1)
	str2 = list(str2)
	str3 = ''
	for index in range(len(str1)):
		if str1[index] == str2[index]:
			str3 += ' '
		else:
			str3 += '|'
	str1 = "".join(str1)
	str2 = "".join(str2)
	print str1
	print str3
	print str2

def readInRange(start,end,read):
	result = []
	for snp in read:
		if snp[0]>=start and snp[0]<=end:
			result.append(snp)
	return result

def readsOverlap(last_read, next_read):
	read_length = len(last_read)
	last_read_end = last_read[-1][0]
	next_read_start = next_read[1][0]
	if last_read_end < next_read_start:
		return False
	else:
		last_part = readInRange(next_read_start,last_read_end,last_read)
		next_part = readInRange(next_read_start,last_read_end,next_read)
		if last_part == next_part:
			return True
		else:
			return False

def readToSNPS(reference,read,offset):
	length = len(read)
	reference = list(reference)
	read = list(read)
	snps = []
	for index in range(length):
		if read[index] != reference[index]:
			snps.append( (index+offset,read[index]) )
	return snps

def SNPSToHaplotypes(read_in):
	reads = read_in[:]
	haplotype0 = reads[0]
	haplotype1 = reads[1]
	for read in reads:
		if readsOverlap(haplotype0,read):
			haplotype0 = sorted(set(haplotype0)|set(read))
		if readsOverlap(haplotype1,read):
			haplotype1 = sorted(set(haplotype1)|set(read))
	return (haplotype0,haplotype1)

def readsToSNPS(filename):
	reference_genome_path = 'reference_genomes/'
	reference_genome_path += filename
	with open(reference_genome_path,'r') as reference_content:
		reference_genome = reference_content.read()
	reference_length = len(reference_genome)
	reads = getReads(filename)
	read_length = len(reads[0])	
	snps = []
	for read in reads:
		for reference_index in range(reference_length-read_length+1):
			possible_snps = readToSNPS(reference_genome[reference_index:reference_index+read_length],read,reference_index)
			if len(possible_snps) <= (read_length/10)*3 + 1:
				snps.append(possible_snps)
				break
	snps = unique_sorted(snps)
	return snps

def makeCacheOld(filename,cacheLen):
	reference_genome_path = 'reference_genomes/'
	reference_genome_path += filename
	with open(reference_genome_path,'r') as reference_content:
		reference_genome = reference_content.read()
	reference_length = len(reference_genome)
	cache = defaultdict(list)
	for index in range(reference_length-cacheLen+1):
		substring = reference_genome[index:index+cacheLen]
		cache[substring].append(index)
	return cache

def chainsWithinOne(chain):
	results = []
	for i in range(len(chain)):
		this_str = list(chain[:])
		this_str[i] = 'a'
		results.append("".join(this_str))
		this_str[i] = 't'
		results.append("".join(this_str))
		this_str[i] = 'c'
		results.append("".join(this_str))
		this_str[i] = 'g'
		results.append("".join(this_str))
	return results

def chainsWithinK(chain,k):
	if k == 1:
		return chainsWithinOne(chain)
	else:
		results = chainsWithinOne(chain)
		final_results = []
		for result in results:
			final_results += chainsWithinK(result,k-1)
		return final_results


def withinDistance(cache,key,distance):
	poss_strings = chainsWithinK(key,distance)
	locations = []
	for string in poss_strings:
		locations += cache[string]
	locations = sorted(set(locations))
	return locations

def makeCache(reference_genome,cacheLen):
	reference_length = len(reference_genome)
	cache = defaultdict(list)
	for index in range(reference_length-cacheLen+1):
		substring = reference_genome[index:index+cacheLen]
		cache[substring].append(index)
	return cache


def rTH(filename,cache_len):
	start_time = time.time()
	reference_genome_path = 'reference_genomes/'
	reference_genome_path += filename
	output_path = 'run_data/'
	output_path += filename
	output_hap0 = output_hap1 = output_path
	output_hap0 += '.gen0'
	output_hap1 += '.gen1'
	with open(reference_genome_path,'r') as reference_content:
		reference_genome = reference_content.read()
	reference_length = len(reference_genome)
	reads = getReads(filename)
	read_length = len(reads[0])	
	cache = makeCache(reference_genome,cache_len)
	snps = []
	for read in reads:
		substring = read[:cache_len]
		indexes = withinDistance(cache,substring,3)
		for match_index in indexes:
			possible_snps = readToSNPS(reference_genome[match_index:match_index+read_length],read,match_index)
			if len(possible_snps) <= (read_length/10)*3+1:
				snps.append(possible_snps)
				break
	snps = unique_sorted(snps)
	haplotype0_snps = snps[0]
	haplotype1_snps = snps[1]
	snps.remove(snps[0])
	snps.remove(snps[0])
	for read in snps:
		if readsOverlap(haplotype0_snps,read):
			haplotype0_snps = sorted(set(haplotype0_snps)|set(read))
		if readsOverlap(haplotype1_snps,read):
			haplotype1_snps = sorted(set(haplotype1_snps)|set(read))
	reference_genome = list(reference_genome)
	hap0_file = open(output_hap0,'w')
	hap0_next = haplotype0_snps.pop(0)
	hap1_file = open(output_hap1,'w')
	hap1_next = haplotype1_snps.pop(0)
	for index in range(reference_length):
		if index == hap0_next[0]:
			hap0_file.write(hap0_next[1])
			if len(haplotype0_snps) != 0:
				hap0_next = haplotype0_snps.pop(0)
		else:
			hap0_file.write(reference_genome[index])
		if index == hap1_next[0]:
			hap1_file.write(hap1_next[1])
			if len(haplotype1_snps) != 0:
				hap1_next = haplotype1_snps.pop(0)
		else:
			hap1_file.write(reference_genome[index])
	elapsed_time = time.time() - start_time
	return elapsed_time


def readsToHaps(filename):
	start_time = time.time()
	reference_genome_path = 'reference_genomes/'
	reference_genome_path += filename
	output_path = 'run_data/'
	output_path += filename
	output_hap0 = output_hap1 = output_path
	output_hap0 += '.gen0'
	output_hap1 += '.gen1'
	with open(reference_genome_path,'r') as reference_content:
		reference_genome = reference_content.read()
	reference_length = len(reference_genome)
	reads = getReads(filename)
	read_length = len(reads[0])	
	snps = []
	for read in reads:
		for reference_index in range(reference_length-read_length+1):
			possible_snps = readToSNPS(reference_genome[reference_index:reference_index+read_length],read,reference_index)
			if len(possible_snps) <= (read_length/10)*3 + 1:
				snps.append(possible_snps)
				break
	snps = unique_sorted(snps)
	haplotype0_snps = snps[0]
	haplotype1_snps = snps[1]
	snps.remove(snps[0])
	snps.remove(snps[0])
	for read in snps:
		if readsOverlap(haplotype0_snps,read):
			haplotype0_snps = sorted(set(haplotype0_snps)|set(read))
		if readsOverlap(haplotype1_snps,read):
			haplotype1_snps = sorted(set(haplotype1_snps)|set(read))
	reference_genome = list(reference_genome)
	hap0_file = open(output_hap0,'w')
	next_snp = haplotype0_snps.pop(0)
	for index in range(reference_length):
		if index == next_snp[0]:
			hap0_file.write(next_snp[1])
			if len(haplotype0_snps) != 0:
				next_snp = haplotype0_snps.pop(0)
		else:
			hap0_file.write(reference_genome[index])
	hap1_file = open(output_hap1,'w')
	next_snp = haplotype1_snps.pop(0)
	for index in range(reference_length):
		if index == next_snp[0]:
			hap1_file.write(next_snp[1])
			if len(haplotype1_snps) != 0:
				next_snp = haplotype1_snps.pop(0)
		else:
			hap1_file.write(reference_genome[index])
	elapsed_time = time.time() - start_time
	return elapsed_time

def rTHOptOld(filename):
	start_time = time.time()
	reference_genome_path = 'reference_genomes/'
	reference_genome_path += filename
	output_path = 'run_data/'
	output_path += filename
	output_hap0 = output_hap1 = output_path
	output_hap0 += '.gen0'
	output_hap1 += '.gen1'
	with open(reference_genome_path,'r') as reference_content:
		reference_genome = reference_content.read()
	reference_length = len(reference_genome)
	reads = getReads(filename)
	read_length = len(reads[0])	
	sequences = []
	hap0 =[]
	hap1 =[]
	for read in reads:
		for reference_index in range(reference_length-read_length+1):
			possible_snps = readToSNPS(reference_genome[reference_index:reference_index+read_length],read,reference_index)
			if len(possible_snps) <= (read_length/10)*3 + 1:
				for seq in sequences:
					if readsOverlap(seq,possible_snps):
						seq = sorted(set(seq)|set(possible_snps))
						break
					elif readsOverlap(possible_snps,seq):
						seq = sorted(set(possible_snps)|set(seq))
						break
				else:
					sequences.append(possible_snps)
				break
	sequences = unique_sorted(sequences)
	hap0 = sequences[0]
	hap1 = sequences[1]
	sequences.remove(sequences[0])
	sequences.remove(sequences[0])	
	print hap0
	print hap1
	print sequences
	for seq in sequences:
		if readsOverlap(hap0,seq):
			hap0 = sorted(set(hap0)|set(seq))
		elif readsOverlap(hap1,seq):
			hap1 = sorted(set(hap1)|set(seq))
	hap0_file = open(output_hap0,'w')
	hap1_file = open(output_hap1,'w')
	hap0_next = hap0.pop(0)
	hap1_next = hap1.pop(0)
	for index in range(reference_length):
		if index == hap0_next[0]:
			hap0_file.write(hap0_next[1])
			if len(hap0) != 0:
				hap0_next = hap0.pop(0)
		else:
			hap0_file.write(reference_genome[index])
		if index == hap1_next[0]:
			hap1_file.write(hap1_next[1])
			if len(hap1) != 0:
				hap1_next = hap1.pop(0)
		else:
			hap1_file.write(reference_genome[index])
	elapsed_time = time.time() - start_time
	return elapsed_time




