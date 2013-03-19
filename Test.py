import CGenerator
import HAssembler

def runBase(g_length,r_length):
	CGenerator.randomGenome('auto',g_length)
	CGenerator.makeHaplotypes('auto')
	CGenerator.makeReads('auto',r_length)
	baseline = HAssembler.readsToHaps('auto')
	return baseline

def runOpt(g_length,r_length):
	CGenerator.randomGenome('auto',g_length)
	CGenerator.makeHaplotypes('auto')
	CGenerator.makeReads('auto',r_length)
	optomized = HAssembler.rTH('auto',10)
	return optomized

def test():
		case0_path = case1_path = gen0_path = gen1_path = 'run_data/'
		case0_path += 'auto.case0'
		case1_path += 'auto.case1'
		gen0_path += 'auto.gen0'
		gen1_path += 'auto.gen1'
		with open(case0_path,'r') as case0_content:
			case0 = case0_content.read()
		with open(case1_path,'r') as case1_content:
			case1 = case1_content.read()
		with open(gen0_path,'r') as gen0_content:
			gen0 = gen0_content.read()
		with open(gen1_path,'r') as gen1_content:
			gen1 = gen1_content.read()
		if case0 != gen0:
			if case0 != gen1:
				passed = False
			else:
				if case1 != gen0:
					passed = False
				else:
					passed = True
		else:
			if case1 != gen1:
				passed = False
			else: 
				passed = True
		return passed


def trial(g_length,r_length):
	opt = runOpt(g_length,r_length)
	optresults = test()
	base = runBase(g_length,r_length)
	baseresults = test()
	opt = str(opt)
	optresults = str(optresults)
	base = str(base)
	baseresults = str(baseresults)
	print opt + " " + optresults + " | " + base + " " + baseresults


