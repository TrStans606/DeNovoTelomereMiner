import subprocess
import os
import argparse
import glob
import urllib.parse
import pandas as pd
import re
import copy
import shutil
import sys

# Command line argument set up
parser = argparse.ArgumentParser(
    prog='DenovoTelomereMiner',
    argument_default='try --simple',
    description=(
        'An all-in-one pipeline and analysis tool for de-novo Telomeres'
    ),
    epilog=(
        'The -s and -i options are mutually exclusive. If simple mode is not used options -g is required.'
    ),
)

# Command line option for separated files 
parser.add_argument(
    '-s',
    action='store',
    nargs=2,
    help=(
        'User defined name for two separated fastq files. This argument expects two file names separated by a space.'
    ),
    default='n/a'
)

# Command line option for interleaved files
parser.add_argument(
    '-i',
    action='store',
    help=(
        'User defined name for an interleaved fastq file.'
    ),
    default='n/a'
)

#command line argument for a genome file
parser.add_argument(
	'-g',
	action='store',
	help=(
		'The name of the assembled genome fasta file used for allignment'
	),
	required='n/a'
)

#command line argument for a directory
parser.add_argument(
	'-d',
	action='store',
	help=(
		'The name of the directory where all output files will be placed. The name defaults to the name of the fasta file provided. Directory names must be unique'
	),
	default='n/a'
)

#specifies the tel reapeat you want to search for
parser.add_argument(
	'-t',
	action='store',
	help=(
		'The foreward tel repeat you want to search for. Default is CCCTAA.'
	),
	default='CCCTAA'
)

#command lind argument for choosing filter files
parser.add_argument(
	'-f',
	action='store', nargs='*', 
	help=(
		'List all of the seqeunces you want removed from reads. Must be in fasta format.'    
	),
	default='n/a'
)

#command line argument for determining if a cluster is a canidate de-novo 
parser.add_argument(
	'--cut', 
	action='store', 
	type=int, 
	help= ('The number of reads in a cluster before it is labeled a canidate de-novo. The default is five.'
	),
	default=5
)

#command line argument for activating simple mode
parser.add_argument(
	'--simple', 
	action='store_true',
	help = ('Activates Simple mode, a step by step input mode'
	),
	default=False
)

#command lind argument for choosing additional blast files
parser.add_argument(
	'--add', 
	action='store', 
	nargs='*', 
	help=('List all addition seqeunces you want compared to the reads. All must be in fasta format'
	),
	default='n/a'
)
parser.add_argument(
	'--config', 
	action='store_true',
	help = ('alters the config file in program'
	),
	default=False
)

# Parse the arguments
args = parser.parse_args()

#config file maker

def config():
	print('Manual Time setup')
	readDir = input('Raw Reads Directory: ')
	genomeDir = input('Assembled genome Directory: ')
	filterDir = input('Filters directory: ')
	addDir = input('Additional Blast directory: ')
	with open('config.ini', 'w') as write:
		write.write(
		f'readDirectory={readDir}\ngenomeDirectory= \
			{genomeDir}\nfilterDirectory={filterDir}\naddDirectory= \
				{addDir}')

#simple input collection mode
def inputCollect():
	
	choiceMade = False
	cutOff = 5
	while not choiceMade:
		choice = input(
			'Please enter if you are using sperated (s) genome files or an \
				interleaved (i) file ')
		if choice == 'i' or 's':
			choiceMade = True
		else:
			print('that is not a valid input. Please enter i for interleaved\
				  or s for seperated')
	if choice == 's':
		
		genomeR1 = input(
			'Specify the name of the R1 file you want: ')
		while not filePathTest(os.path.join(readDir + genomeR1)):
			genomeR1 = input(
				'Specify the name of the R1 file you want: ')
		
		genomeR2= input(
				'Specify the name of the R2 file you want: ')
		while not filePathTest(os.path.join(readDir + genomeR2)):
			genomeR2= input(
					'Specify the name of the R2 file you want: ')
	
	elif choice == 'i':
		genomeI = input(
			'Specify the name of the interleaved fastq file ')
		while not filePathTest(os.path.join(readDir + genomeI)):
			genomeI = input(
				'Specify the name of the interleaved fastq file ')
	
	blastGenome = input(
		'Specify genome file to be blasted against the single reads')
	while not filePathTest(os.path.join(genomeDir + blastGenome)):
		blastGenome = input(
			'Specify genome file to be blasted against the single reads')
	
	directory = input(
		'Specify the output directory name. Press enter for default: ')
	if directory == '':
		directory = blastGenome
	if choice == 's':
		return [cutOff, choice, genomeR1, genomeR2, blastGenome, directory]
	elif choice == 'i':
		return [cutOff, choice, genomeI, blastGenome, directory]

#tests for the validity of a file path
def filePathTest(filename):
	if os.access(filename, os.F_OK):
		return True
	else:
		print(
			"That file doesn't exist. Please enter a valid filename")
		return False

def fileListBuilder(file):
	with open(file, 'r') as read:
		lines = read.readlines()
		lineCnt = len(lines)
	return lines, lineCnt

def teloPortPipeline():
	dirs = \
		['teloPortOut', 
		'muscleOut',
		'consOut',
		'telomereReads/telomereClusters', 
		'telomereReads/deNovoTelomeres',
		'blastOut/blastGenome',
		'blastOut/blastFilter',
		'blastOut/blastAdd',
		'blastOut/blastDict'
		]
	for path in dirs:
		os.makedirs(os.path.join('Outputs', 
			directory, 
			path))
	if choice == 's':
		telomereFinderSeperated()
	elif choice =='i':
		telomereFinderInter()
	idTagger()
	#    os.remove(
	#        f'Outputs/{directory}/teloPortOut/pairReads.fastq'
	#        )
	junctionFinder()
	sequenceQuality()
	mmseqs2Call()
	#wcdestCall()
	#wcdInterrogate()

#runs telomereFinder with the -s argument
#this is used for seperated fastq files
def telomereFinderSeperated():
	command = ['Programs/bin/telomereFinder', 
		'-s',
		os.path.join(readDir, genomeR1),
		os.path.join(readDir, genomeR1),
		'-o',
		'Outputs/' + directory + '/teloPortOut/',
		'-t',
		telRepeat]
	print(' '.join(command))
	subprocess.run(' '.join(command),
		shell=True,
		check=True)

#runs telomereFinder with the -i argument
#this is used for interleaved fastq files
def telomereFinderInter():
	command = ['Programs/TeloPort/build/apps/telomereFinder',
		'-i',
		os.path.join(readDir, genomeI),
		'-o',
		f'Outputs/{directory}/teloPortOut/',
		'-t',
		telRepeat]
	print(' '.join(command))
	subprocess.run(' '.join(command),
		shell=True,
		check=True)

#tags every read with a unique id number for tracking
def idTagger():
	pairLines, pairLineCnt = fileListBuilder(
		os.path.join('Outputs',
			directory,
			'teloPortOut','pairReads.fastq')
		)
	newFile = os.path.join('Outputs',
		directory,
		'teloPortOut','subTelReads.fastq')
	with open(newFile, 'a') as write:
		x = 0
		for line in pairLines:
			if x % 4 == 0:
				newLine = f"@{directory}_{x}\t{line.lstrip('@')}"
				write.write(newLine)
			else:
				write.write(line)
			x += 1

#calls junctionfinder
def junctionFinder():
	command = ['Programs/bin/junctionFinder',
					'-i',
					os.path.join('Outputs',
								 directory,
								 'teloPortOut','subTelReads.fastq'),
					'-o',
					os.path.join('Outputs',
								 directory,
								 'teloPortOut'),
					f'-t {telRepeat}',
					'--revc 0',
					'--splitJunc 1',
					'--splitDir 0']
	print(' '.join(command))
	subprocess.run(' '.join(command),
				   shell=True,
				   check=True)

#calls sequenceQuality
def sequenceQuality():
	command = ['Programs/bin/sequenceQuality',
					'-i',
					os.path.join('Outputs',
								 directory,
								 'teloPortOut','telAdjSeq.fastq'),
					'-o',
					os.path.join('Outputs',
								 directory,
					'teloPortOut',f'{directory}hiQualityTelAdjSeq.fasta'),
					'--ofmt',
					'fasta',
					'-c',
					'60',
					'-l',
					'40']
	print(' '.join(command))
	subprocess.run(' '.join(command),
				   shell=True,
				   check=True)

#calls mmseqs2
def mmseqs2Call():
	command = ['mmseqs easy-cluster',
			   os.path.join('Outputs',
								 directory,
					'teloPortOut',f'{directory}hiQualityTelAdjSeq.fasta'),
					os.path.join('Outputs',
								 directory,
					'teloPortOut','clusters'),
					"tmp/"]
	subprocess.run(' '.join(command),
					shell=True,
					check=True)

#processes the clustering into reads for further use
def mmseqs2_processing(cutOff_mm):
	print("Entered mmseqs2_processing")
	clusters = pd.read_csv(f"Outputs/{directory}/teloPortOut/clusters_cluster.tsv",sep="\t",header=None)
	labels = {}
	i=0
	for value in clusters.loc[:,0]:
		if value not in labels:
			labels[value] = f'cluster{i}'
			i+=1
	i=0
	for value in clusters.iloc[:,0]:
		clusters.loc[i,0] = labels[value]
		i+=1
	with open(f"Outputs/{directory}/teloPortOut/{directory}hiQualityTelAdjSeq.fasta",'r') as read:
		seqs = read.readlines()
	i=0
	for value in clusters.iloc[:,0]:
		with open(f"Outputs/{directory}/telomereReads/{value}.fasta",'a') as write:
			for j in range(len(seqs)):
				if re.search(clusters.loc[i,1],seqs[j]):
					write.write(f"{seqs[j]}{seqs[j+1]}")
			i+=1
	for value in pd.unique(clusters.iloc[:,0]):
		with open(f"Outputs/{directory}/telomereReads/{value}.fasta",'r') as read:
			lines = read.readlines()
		if len(lines) /2 <= cutOff_mm:
			shutil.move(f"Outputs/{directory}/telomereReads/{value}.fasta",f"Outputs/{directory}/telomereReads/deNovoTelomeres")
	#dry run before checking if the user wants to reset the cluster cutoff
	cluster_num = 0
	clusters_rename = glob.glob(f"Outputs/{directory}/telomereReads/*.fasta")
	for file in clusters_rename:
		cluster_num += 1
	#checks if the user wants to reset the cluster cutoff
	reset_check = False
	while not reset_check:
		print(f"You have generated {cluster_num} clusters.")
		reset_gate = input("Would you like to redo the cluster cutoff? (y/n) ")
		match reset_gate:
			case "y":
				recluster(cutOff_mm)
				reset_check = True
			case "n":
				reset_check = True
			case _:
				print("That is not a valid input. Please enter y for yes or n for no.")
				reset_check = False
	#the proper move function
	cluster_num = 0
	clusters_rename = glob.glob(f"Outputs/{directory}/telomereReads/*.fasta")
	for file in clusters_rename:
		shutil.move(file,f"Outputs/{directory}/telomereReads/telomereClusters/cluster{cluster_num}.fasta")
		cluster_num += 1
	singles = glob.glob(f"Outputs/{directory}/telomereReads/deNovoTelomeres/cluster*.fasta")
	with open(f"Outputs/{directory}/telomereReads/deNovoTelomeres/{directory}deNovos.fasta",'a') as write:
		for file in singles:
			with open(file,'r') as read:
				for line in read:
					write.write(f"{line.rstrip()}\n")
			os.remove(file)
	return cutOff_mm

def recluster(cutOff):
	new_cutOff = int(input(f"Please enter the new cluster cutoff (the previous value was {cutOff}): "))
	files = glob.glob(f"Outputs/{directory}/telomereReads/deNovoTelomeres/*.fasta")
	for file in files:
		os.remove(file)
	files = glob.glob(f"Outputs/{directory}/telomereReads/*.fasta")
	for file in files:
		os.remove(file)
	mmseqs2_processing(new_cutOff)

#runs muscle and EMBOSS cons
def autoMuscle():
	print("muscle entered")
	clusters = os.path.join('Outputs',
										directory,
										'telomereReads',
										'telomereClusters')
	cntCluster = len(glob.glob(f'{clusters}/cluster*.fasta'))
	for i in range(0,cntCluster):
		command = ['muscle',
						'-align',
						f'{clusters}/cluster{i}.fasta',
						'-output',
						os.path.join('Outputs',
									 directory,
									 'muscleOut',
									 f'{directory}cluster{i}.msa')]
		subprocess.run(' '.join(command),
					   shell=True,
					   check=True)

	allignments = os.path.join('Outputs',
							   directory, 
							   'muscleOut')
	cons = os.path.join('Outputs',
							   directory, 
							   'consOut')
	cntAllignment = len(glob.glob(f'{allignments}/{directory}cluster*.msa'))
	for i in range(0, cntAllignment):
		command = ['cons',
						'-sequence',
						f'{allignments}/{directory}cluster{i}.msa',
						'-outseq',
						f'{cons}/cons{i}.fasta',
						'-name',
						f'cluster{i}']
		subprocess.run(' '.join(command),
					   shell=True,
					   check=True)
	return cntCluster

#turns all consenus seqeunces into one file and blasts
def consBlast():
	cons = os.path.join('Outputs',
							   directory, 
							   'consOut')
	blast = os.path.join('Outputs',
							   directory, 
							   'blastOut')
	singles = os.path.join('Outputs',
							   directory, 
							   'telomereReads',
							   'deNovoTelomeres')
	cntCons = len(glob.glob(f'{cons}/cons*.fasta'))
	for i in range(0,cntCons):
		with open(f'{cons}/cons{i}.fasta', 'r') as read:
			consLines = read.readlines()
		with open(f'{blast}/{directory}catCons.fasta', 'a') as write:
			for line in consLines:
				write.write(line)
	query = f'{singles}/{directory}deNovos.fasta'
	subject = f'{blast}/{directory}catCons.fasta'
	output =  f'{blast}/{directory}blastCons.txt'
	dust = ''
	blastRun(query,
		 subject,
		  output,
		  dust)

#runs blast            
def blastRun(query, subject, output, dust):
	outfmt ='-outfmt "6 qseqid sseqid pident length mismatch \
		gapopen qstart qend sstart send evalue qlen"'
	command = ['blastn',
						'-query',
						query,
						'-subject',
						subject,
						outfmt,
						dust,
						'>>',
						output]
	subprocess.run(' '.join(command),
					   shell=True,
					   check=True)

def falsePosFilter():
	clusters = os.path.join('Outputs',
							   directory,
							   'telomereReads',
							   'telomereClusters')
	deNovoPath = os.path.join('Outputs',
							   directory, 
							   'telomereReads',
							   'deNovoTelomeres')
	trueDeNovos = 0
	falsePos = dictMaker(os.path.join('Outputs',
									  directory,
									  'blastOut',
									  f'{directory}blastCons.txt'), 
						 3)
	singles, cnt = fileListBuilder(
		os.path.join('Outputs',
					 directory,
					 'telomereReads',
					 'deNovoTelomeres',
					 f'{directory}deNovos.fasta'))
	for i in range(0,cnt):
		if singles[i][0] == '>':
			deNovo = singles[i].split('\t')[0].lstrip('>')
			if deNovo in falsePos:
				cluster = falsePos[deNovo][0]
				with open(f'{clusters}/{cluster}.fasta', 'a') as write:
					write.write(singles[i])
					write.write(singles[i+1])
			else:
				with open(f'{deNovoPath}/{directory}trueDeNovos.fasta', 'a') \
					as write:
						write.write(singles[i])
						write.write(singles[i+1])
						trueDeNovos += 1
	return trueDeNovos
											
def dictMaker(file, secondVal):               
	with open(file, 'r') as read:
		dictionary = {}
		for line in read:
			item = line.split('\t')
			if item[0] not in dictionary:
				dictionary[item[0]] = (item[1], item[secondVal])
			elif secondVal == 3:
				if item[3] > dictionary[item[0]][1]:
					dictionary[item[0]]=(item[1], item[3])
	return dictionary
					
def histogramBuilder():
	clusters = os.path.join('Outputs',
				 directory,
				 'telomereReads',
				 'telomereClusters')
	clustersCnt = len(glob.glob(f'{clusters}/cluster*.fasta'))
	with open(os.path.join('Outputs',
						   directory,
						   'clusterHistogram.txt'), 'a') as \
							  write:
								write.write('Distribution of Clusters\n')
	for x in range(0,clustersCnt):
		readCnt = 0
		with open(f'{clusters}/cluster{x}.fasta', 'r') as read:
			for line in read:
				if line[0] == '>':
					readCnt += 1
		with open(os.path.join('Outputs',
							   directory,
							   'clusterHistogram.txt'), 'a') as \
								  write:
									line = f'cluster{x}\t{readCnt}\t'
									bar = "#" * readCnt
									line = f'{line} {bar}'
									write.write(line +'\n')
				
def deNovoFilter():
	if filtersCnt > 0:
		for file in args.f:
			blastRun(os.path.join('Outputs',
							   directory,
							   'telomereReads',
							   'deNovoTelomeres',
							   f'{directory}trueDeNovos.fasta'),
				  f'{filterDir}{file}',
				  os.path.join('Outputs',
							   directory,
							   'blastOut',
							   'blastFilter',
							   f'{directory}{file.split[0]}.txt'),
				  '')
			with open( os.path.join('Outputs',
						  directory,
						  'blastOut',
						  'blastFilter',
						  f'{directory}filters.txt'), 'a') as write:
				lines, linesCnt = \
					fileListBuilder(os.path.join('Outputs',
								  directory,
								  'blastOut',
								  'blastFilter',
								  f'{directory}{file.split[0]}.txt'))
				write.write(lines)
		filterSingles = dictMaker(os.path.join('Outputs',
					  directory,
					  'blastOut',
					  'blastFilter',
					  f'{directory}filters.txt', 7))
		filteredDeNovos = 0
		lines, lineCnt = fileListBuilder(
				os.path.join('Outputs',
				  directory,
				  'telomereReads',
				  'deNovoTelomeres',
				  f'{directory}trueDeNovos.fasta'))
		
		for i in range(0,lineCnt):
			if lines[i][0] == '>':
				ID = lines[i].split('\t')[0].lstrip('>')
				if ID in filterSingles:
					with open(os.path.join('Outputs',
					  directory,
					  'telomereReads',
					  'deNovoTelomeres',
					  f'{directory}filteredDeNovos.fasta'), 'a') as write:
						
						write.write(lines[i])
				else:
					with open(os.path.join('Outputs',
					  directory,
					  'telomereReads',
					  'deNovoTelomeres',
					  f'{directory}blastDeNovos.fasta'), 'a') as write:
						filteredDeNovos += 1
						write.write(lines[i])
		return filteredDeNovos
	else:
		lines, lineCnt = fileListBuilder(
				os.path.join('Outputs',
				  directory,
				  'telomereReads',
				  'deNovoTelomeres',
				  f'{directory}trueDeNovos.fasta'))
		with open(os.path.join('Outputs',
		  directory,
		  'telomereReads',
		  'deNovoTelomeres',
		  f'{directory}blastDeNovos.fasta'), 'a') as write:
			filteredDeNovos = lineCnt / 2
			for line in lines:
				write.write(line)
		return filteredDeNovos

def telContigDictMaker():
	with open(os.path.join('Outputs',
						   directory,
						   'blastOut',
						   'blastDict',
						   f'{directory}telRepeats'), 'a') as write:
		write.write(f'>{blastGenome} telomere repeats\n')
		write.write(telRepeat * 8)
	blastRun(f'{genomeDir}/{blastGenome}', 
		  os.path.join('Outputs',
					   directory,
					   'blastOut',
					   'blastDict',
					   f'{directory}telRepeats'),
		  os.path.join('Outputs',
					   directory,
					   'blastOut',
					   'blastDict',
					   f'{directory}blastDictOut6.txt'),
		  '-dust no')
	telContigs = {}
	posistionStart = 0
	posistionEnd = 0
	with open(os.path.join('Outputs',
				 directory,
				 'blastOut',
				 'blastDict',
				 f'{directory}blastDictOut6.txt'), 'r') as read:
		for line in read:
			chromosome = line.split('\t')[0]
			qEnd = line.split('\t')[6]
			qStart = line.split('\t')[5]
			contigLen = line.split('\t')[10]
			if f'{chromosome}*s' not in telContigs:
				posistionStart = 0
			if f'{chromosome}*e' not in telContigs:
				posistionEnd = 0
			#this determines if a telomere match is at the beginning or ending
			#of a contig by seeing if it is in the first or last third
			if int(qEnd) < (int(contigLen)/3):
				if int(posistionStart) < int(qEnd):
					telContigs[f'{chromosome}*s'] = \
						('start', qEnd)
					posistionStart = qEnd
			if int(qEnd) > (int(contigLen) / 3):
				if int(posistionEnd) < int(qEnd):
					telContigs[f'{chromosome}*e'] = \
						('end', qStart)
					posistionEnd = qStart
	return telContigs
		
def blastInterrogate():

	readLines, linesCnt = fileListBuilder(os.path.join( \
						'Outputs',
						directory,
						'blastOut',
						'blastGenome',
						f'{directory}blastGenomeOut6.txt'))
	for x in range(0, linesCnt):
		telContigCheck = False
		blastAn = ''
	   # queryName = readLines[x].split('\t')[0]
		subjectName = readLines[x].split('\t')[1]
		telDisS = 0
		
		if (subjectName + '*s') in telContigs:
			blastAn += "\tTelYs"
			telContigCheck = True
		else:
			blastAn += "\tTelNs"
		
		if telContigCheck:
			posistion = telContigs[readLines[x].split("\t")[1]+ '*s'][1]
			telDisS = abs(int(readLines[x].split('\t')[7]) - int(posistion))
			
		telContigCheck = False
		telDisE = 0
		
		if f'{subjectName}*e' in telContigs:
			blastAn += "\tTelYe"
			telContigCheck = True
		
		if telContigCheck:
			posistion = telContigs[readLines[x].split('\t')[1] + '*e'][1]
			telDisE =  abs(int(readLines[x].split('\t')[8]) - int(posistion))
		
		telContigCheck = False
		
		if min(telDisS, telDisE) > 0:
			blastAn += f'\t{min(telDisS, telDisE)}\t'
			
		elif min(telDisS, telDisE) == 0:
			blastAn += f'\t{max(telDisS, telDisE)}\t'
		
		revS = False

		revQ = False
		
		if int(readLines[x].split('\t')[5]) \
			> int(readLines[x].split('\t')[6]):
				revS = True
				
		if int(readLines[x].split('\t')[7]) \
			> int(readLines[x].split('\t')[8]):
				revQ = True
				
		if not revS and revQ:
			rev = "revT"
		else:
			rev = "revF"
			
		blastAn += f'{rev}\n'
		
		blastAn = readLines[x].rstrip('\n') + blastAn
		
		with open(os.path.join('Outputs',
							   directory,
							   f'{directory}blastGenomeOut6Annotated.txt'),\
				  'a') as append:
			append.write(blastAn)
			
def gffBuilder():
	lines, lineCnt = fileListBuilder(
		os.path.join('Outputs',
					 directory,
					 f'{directory}blastGenomeOut6Annotated.txt'))
	with open(os.path.join('Outputs',
					directory,
					f'{directory}denovos.gff'),\
			  'a') as append:
		append.write('##gff-version 3.1.26\n')
		cnt = '0'
		for line in lines:
			if line.split('\t')[14].rstrip() == 'revF':
				strand = '+'
			elif line.split('\t')[14].rstrip() == 'revT':
				strand = '-'
			name = line.split("\t")[1]
			length = \
				abs(int(line.split('\t')[8]) - int(line.split('\t')[7]))
			
			
			gffLine = '\t'.join(
							[urllib.parse.quote(name),
							'teloPortWrapper',
							'telomere',
							line.split('\t')[7],
							line.split('\t')[8],
							line.split('\t')[9],
							strand,
							'.'])
			gffLine += f'\tID=Telomere{cnt};Name={name};Length={str(length)}'
			telDist = line.split('\t')[12]
			gffLine += f';Distance from Telomere={telDist}\n'
			
			append.write(gffLine)
			
def resultsBuilder():
	resultsLine = f'{directory}\n'
	
	rawReadCnt = 0
	
	if choice == 's':
		with open(f'{readDir}/{genomeR1}', 'r') as read:
			for line in read:
				if line[0] == '@':
					rawReadCnt += 1
		with open(f'{readDir}/{genomeR2}', 'r') as read:
			for line in read:
				if line[0] == '@':
					rawReadCnt += 1
	elif choice == 'i':
		with open(f'{readDir}/{genomeI}', 'r') as read:
			for line in read:
				if line[0] == '@':
					rawReadCnt += 1
					
	resultsLine += f'Number of raw reads processed: {str(rawReadCnt)}\n'
	
	teloCnt = 0 
	
	with open(os.path.join('Outputs',
			  directory,
			  'teloPortOut',
			  'telReads.fastq'), 'r') as read:
		for line in read:
			if line[0] == '@':
				teloCnt += 1
				
	resultsLine += f'Number of Telomere reads processed: {teloCnt}\n'
	
	deNovoCnt = 0
	
	with open(os.path.join('Outputs',
			  directory,
			  'telomereReads',
			  'deNovoTelomeres',
			  f'{directory}deNovos.fasta'), 'r') as read:
		for line in read:
			if line[0] == '>':
				deNovoCnt += 1
				
	resultsLine += f'Number of de-novo telomeres found by TeloPort: \
		{deNovoCnt}\n'
				
	resultsLine += f'Number of true de-novo telomeres: {trueDeNovos}\n'
	
	resultsLine += f'Number of clusters: {trueClusterCnt}\n'
	
	filteredOut = trueDeNovos - filteredDeNovos
	
	resultsLine += f'Number of true de-novos filtered out: {filteredOut}\n'
	
	resultsLine += \
		f'Number of true de-novo compared to the genome: {filteredDeNovos}\n'
	
	genomesCnt = 0
	
	with open(os.path.join('Outputs',
				  directory,
				  'blastOut',
				  'blastGenome',
				  f'{directory}blastGenomeOut6.txt'), 'r') as read:
		genomesCnt = len(read.readlines())
	
	resultsLine += \
		f'Number of de novo telomere matches to the genome: {genomesCnt}\n'
	
	with open(os.path.join('Outputs',
						   directory,
						   'clusterHistogram.txt'), 'r') as read:
		histogram = read.readlines()
			
	with open(os.path.join('Outputs',
						   directory,
						   f'{directory}results.txt'), 'a') as append:
		for line in histogram:
			append.write(line)

		append.write(resultsLine)
	with open(f'Outputs/{directory}/{directory}results.txt', 'a') as write:
		with open(f'Outputs/{directory}/{directory}seeds.tsv', 'r') as read:
			for line in read:
				write.write(line)
		
def path_checker(path):
	if not os.path.exists(path):
		print(f"It seems the {path} defined in config.ini doesnt exist please update the file to include the proper path.")
		print("You can either edit config.ini in a text editor or using the --config flag when running again")
		sys.exit()
#finds seed
def seed_finder():
	#creates a dictionary of IDs and their own sequences
	fullSeqs = {}
	with open(f'Outputs/{directory}/{directory}blastGenomeOut6Annotated.txt', 'r') as read:
		for line in read:
			if line.split('\t')[0] not in fullSeqs:
				fullSeqs[line.split('\t')[0]] = ""
	#seeds is intially a copy of fullSeqs but will get the seeds instead of the sequence
	seeds = copy.deepcopy(fullSeqs)
	#writes the sequences to the dictionary
	seqCheck = False
	with open(f'Outputs/{directory}/teloPortOut/subTelReads.fastq','r') as read:
		for line in read:
			if seqCheck:
				fullSeqs[read_name] = line.rstrip()
				seqCheck = False
			if line.split('\t')[0].lstrip("@") in fullSeqs:
				read_name = line.split('\t')[0].lstrip("@")
				seqCheck = True
	#print(fullSeqs)
	#this finds the DeNovo reads to go along with the fullseqs
	seqCheck = False
	with open(f'Outputs/{directory}/telomereReads/deNovoTelomeres/{directory}blastDeNovos.fasta','r') as read:
		for line in read:
			if seqCheck:
				overlap = overlap_finder(fullSeqs[read_name].upper(),line.rstrip().upper(),read_name)
				seeds[read_name] = overlap
				seqCheck = False
			if line.split('\t')[0].lstrip(">") in fullSeqs:
				read_name = line.split('\t')[0].lstrip(">")
				seqCheck = True
	for key,value in seeds.items():
		with open(f'Outputs/{directory}/{directory}seeds.tsv','a') as write:
			write.write(f"{key}\t {value}\n")

def overlap_finder(full_read, deNovo_read,id):
	#grabs the reverse version of the telomere repeat
	rev_tel = ""
	for bp in telRepeat:
		match bp.upper():
			case "A":
				rev_tel += "T"
			case "T":
				rev_tel += "A"
			case "C":
				rev_tel += "G"
			case "G":
				rev_tel += "C"
	rev_tel = rev_tel[::-1]
	#start is given value of a in case it is never initialized a will be used as an error check
	start = "a"
	global genome_read
	genome_read = ""
	#true false gates are used to grab the genome after the headers are grabbed
	genome_gate = False
	pos = 1
	with open(f'Outputs/{directory}/blastOut/blastGenome/{directory}blastGenomeOut6.txt','r') as read:
		for line in read:
			#grabs the chromosome we are reading from and start and end values of the read
			if line.split('\t')[0] == id:
				chromosome = line.split('\t')[1]
				if line.split('\t')[8] > line.split('\t')[7]:
					start_genome = int(line.split('\t')[7])
					end_genome = int(line.split('\t')[8])
				else:
					start_genome = int(line.split('\t')[8])
					end_genome = int(line.split('\t')[7])
				chr = ""
				#reads in the chromosome as a string
				with open(f'{genomeDir}/{blastGenome}','r') as read:
					for line in read:
						if genome_gate:
							if line[0] == ">":
								genome_gate = False
							if genome_gate:
								chr += line.rstrip()
						if line[0] == ">":
							if line.lstrip(">").rstrip() == chromosome:
								genome_gate = True
				#reads the chromosome to find the start and end posistion
				for bp in chr:
					if pos >= start_genome - 20 and pos <= end_genome+(len(full_read)-len(deNovo_read))+20:
						genome_read += bp.upper()
					pos+=1
				break
	for i in range(0,len(full_read)-len(deNovo_read)+1):
		if full_read[i:i+len(deNovo_read)] == deNovo_read:
			start = i
			break
	seed = ""
	if start == "a":
		return "No Overlap. Possible error"
	if start != 0:
		beggining = True
		for i in range(20,0,-1):
			if beggining:
				if genome_read[i].upper() in telRepeat:
					seed = genome_read[i].lower() + seed
				else:
					seed += "|"
					seed = genome_read[i].upper() + seed
				past_bp = genome_read[i].upper()
			else:
				past_bp += genome_read[i].upper()
				if past_bp in telRepeat:
					seed = genome_read[i].lower() + seed
				else:
					seed = genome_read[i].upper() + seed
			beggining = False
		for i in range(20,40):
			seed += full_read[i].upper()
		seed = 'foreward:' + seed
		# #seed += "t."
		# beggining = True
		# for i in range(start,start-10,-1):
		# 	if beggining:
		# 		if genome_read[i].upper() in telRepeat:
		# 			seed = genome_read[i].lower() + seed
		# 		else:
		# 			seed += "|"
		# 			seed = genome_read[i].upper() + seed
		# 		past_bp = genome_read[i].upper()
		# 	else:
		# 		past_bp += genome_read[i].upper()
		# 		if past_bp in telRepeat:
		# 			seed = genome_read[i].lower() + seed
		# 		else:
		# 			seed = genome_read[i].upper() + seed
		# 	beggining = False		
		# for i in range(start,start+10):
		# 	seed += full_read[i].upper()
		# seed = 'foreward:' + seed
	else:
		# beggining = True
		# for i in range(len(deNovo_read)-10,len(deNovo_read)):
		# 	if beggining:
		# 		if genome_read[i].upper() in rev_tel:
		# 			seed = genome_read[i].lower() + seed
		# 		else:
		# 			seed += "|"
		# 			seed = genome_read[i].upper() + seed
		# 		past_bp = genome_read[i].upper()
		# 	else:
		# 		past_bp += genome_read[i].upper()
		# 		if past_bp in rev_tel:
		# 			seed = genome_read[i].lower() + seed
		# 		else:
		# 			seed = genome_read[i].upper() + seed
		# 	beggining = False
		# #seed += "|"
		# #seed += "t."
		# for i in range(len(deNovo_read),len(deNovo_read)+10):
		# 	seed += full_read[i].upper()
		# seed = 'reverse:' + seed
		beggining = True
		length = end_genome-start_genome
		for i in range(length,length+20):
			if beggining:
				if genome_read[i].upper() in rev_tel:
					seed += genome_read[i].lower() 
				else:
					seed += "|"
					seed += genome_read[i].upper() 
				past_bp = genome_read[i].upper()
			else:
				past_bp += genome_read[i].upper()
				if past_bp in rev_tel:
					seed += genome_read[i].lower() 
				else:
					seed += genome_read[i].upper() 
			beggining = False
		for i in range(length-20,length):
			seed = full_read[i].upper() + seed
		seed = 'reverse:' + seed
	return seed
#assigns the path varibles from the config file
if args.config:
	config()
with open('config.ini','r') as read:
	configLines = read.readlines()
	readDir = configLines[0].split('=')[1].split('\n')[0]
	genomeDir = configLines[1].split('=')[1].split('\n')[0]
	filterDir = configLines[2].split('=')[1].split('\n')[0] 
	addDir = configLines[3].split('=')[1]
#makes sure the folders in config.ini exist
path_checker(readDir.strip())
path_checker(genomeDir.strip())
path_checker(filterDir.strip())
path_checker(addDir.strip())

if args.simple:
	inputResults = inputCollect()
	telRepeat = args.t
	if len(inputResults) == 6:
		cutOff = inputResults[0]
		choice = inputResults[1]
		genomeR1 = inputResults[2]
		genomeR2 = inputResults[3]
		blastGenome = inputResults[4]
		directory = inputResults[5]
	else:
		cutOff = inputResults[0]
		choice = inputResults[1]
		genomeI = inputResults[2]
		blastGenome = inputResults[3]
		directory = inputResults[4]
elif args.s == 'n/a':
	choice = 'i'
	genomeI = args.i
elif args.i == 'n/a':
	choice = 's'
	genomeR1 = args.s[0]
	genomeR2 = args.s[1]
if not args.simple:
	blastGenome = args.g
	cutOff = args.cut
	telRepeat = args.t
if args.d == 'n/a' and not args.simple:
	directory = blastGenome.split('.')[0]
else:
	directory = args.d
if args.add != 'n/a':
	addBlastCnt = len(args.add)
else:
	addBlastCnt = 0
if args.f != 'n/a':
	filtersCnt = len(args.f)
else:
	filtersCnt = 0

teloPortPipeline()
cutOff = mmseqs2_processing(cutOff)
trueClusterCnt = autoMuscle()
consBlast()
trueDeNovos = falsePosFilter()
histogramBuilder()
filteredDeNovos = deNovoFilter()
telContigs = telContigDictMaker()
blastRun(os.path.join('Outputs',
				   directory,
				   'telomereReads',
				   'deNovoTelomeres',
				   f'{directory}blastDeNovos.fasta'),
	  f'{genomeDir}/{blastGenome}',
	  os.path.join('Outputs',
				   directory,
				   'blastOut',
				   'blastGenome',
				   f'{directory}blastGenomeOut6.txt'),
	  '')
if addBlastCnt > 0:
	for file in args.add:
		blastRun(f'{addDir}file',
			  f'{genomeDir}/{blastGenome}',
			  os.path.join('Outputs',
						   directory,
						   'blastOut',
						   'blastAdd',
						   f'{directory}{file.split(".")[0]}blast.txt'),
			  '')

blastInterrogate()
gffBuilder()
seed_finder()
resultsBuilder()             
				   
#print(args.s)
