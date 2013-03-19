Gemini Haplotype Assembler
==========================

Takes a bunch of reads, aligns them, and builds them back into two
haplotypes.  

Assumptions:
* There are no sequencing errors!
* There are only two haplotypes!

Files and Notes:
================

!!!IMPORTANT!!!
---------------
The algorithm uses the included directory structure for it's usuage of file-based storage.  Removing
elements from the folder directory structre (i.e. removing run_data or reference_genomes) will result
in the failure of the algorithm.

	HAssembler.py 	|  Defines the actual baseline assembler and solution assembler functions, as well
					|  some auxillary methods neccasay for their operation.
					|
	CGenerator.py   |  Generates random genomes, as well as loads existing genomes from text files in
					|  raw format.
					|
	Test.py 		|  Defines useful structure for testing and comparing the operation of the baseline
					|  and solutions methods

Useage:
=======

General Run:
------------
readsToHaps defines the baseline algorithm; rTH defines the optomized version.  Testing methods are
provided for convienance and clarity of use.  For example, to run a comparision test over a random
genome of length 10k with read length 100:

	import Test
	Test.trial(10000,100)

This gives the following example output:

	3.32129096985 True | 0.682209014893 True
	^				^         	^			^
	baseline 		base 		opt 		opt
	runtime 		correct 	runtime		correct


Generator Methods:
------------------
* randomGenome(filename,length)
* mutate(gene)
* makeHaplotypes(filename)
* makeReads(filename, read_length)

Assembler Methods:
------------------
* unique_sorted(values)
* getReads(filename)
* getReference(filename)
* compareStrings(reference,read)
* alignReadOld(filename,read)
* prettySNPS(str1,str2)
* readInRange(start,end,read)
* readsOverlap(last_read, next_read)
* readToSNPS(reference,read,offset)
* SNPSToHaplotypes(read_in)
* readsToSNPS(filename)
* makeCacheOld(filename,cacheLen)
* chainsWithinOne(chain)
* chainsWithinK(chain,k)
* withinDistance(cache,key,distance)
* makeCache(reference_genome,cacheLen)
* rTH(filename,cache_len)
* rTHOptOld(filename)
* readsToHaps(filename)

Test Methods:
-------------
* runBase(g_length,r_length)
* runOpt(g_length,r_length)
* test()
* trial(g_length,r_length)
