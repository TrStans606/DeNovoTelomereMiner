
import os
import time
import argparse

#command line arguments

parser = argparse.ArgumentParser(prog='TeloPortWrapper',description='An all-in-one pipeline and anaylsis tool for de-novo Telomeres',epilog='Using no arguments will cause the program to run in simple mode. The -s and -i options are mutally exclusive. If simple mode is not used options -b is required.')

parser.add_argument('-s', action='store', nargs=2, help='Allows for seperated fastq file names to be inputted on the command line. This argument exepcts two file names seperated by a space.', default=['x','x'])

parser.add_argument('-i', action='store', help='Allows for the interleaved fastq name to be inputted on the command line. This argument ecpects one file name.', default="x")

parser.add_argument('-b', action='store', help='Allows for the assembled genome file name to be inputted on the command line.', default="x")

parser.add_argument('-d', action='store', help='The name of the directory where all output files will be placed. The directory name will default to the name of the assembled genome file. Directory names must be unique', default="x")

parser.add_argument('--cut', action='store', type=int, help='The number of reads in a cluster before it becomes labeled as a set of single reads. Default is 5.', default=5)

parser.add_argument('--simple', action='store_true', help='Activates Simple mode, a step by step input mode using default values', default=False)

args = parser.parse_args()

#determines if an interleaved of seperated input is needed

def inputCollect():
	choiceMade = False
	global choice, genomeR1, genomeR2, genomeI, directory, blastD, cutoff
	cutoff = 5
	while not choiceMade:
		choice = input('Please enter if you are using sperated (s) genome files or an interleaved (i) file ')
		if choice == 'i' or choice == 's':
			choiceMade = True
		else:
			print('that is not a valid input. Please enter i for interleaved or s for seperated')
	if choice == 's':

		genomeR1 = input('Specify the name of the R1 file you want: ')
		filePathTest('Genomes/' + genomeR1)
		genomeR2= input('Specify the name of the R2 file you want: ')
		filePathTest('Genomes/' + genomeR2)
	elif choice == 'i':
		genomeI = input('Specify the name of the interleaved fastq file ')
		filePathTest('Genomes/' + genomeI)
	blastD = input('Specify genome file to be blasted against the single reads ')
	directory = input('Specify the output directory name. Press enter for default: ')
	if directory == '':
		directory = blastD
	filePathTest('Genomes/db/' + blastD)

#test to make sure the file path is valid

def filePathTest(fileName):

	try:
		test= open(fileName, 'r')
	except:
		print('We cannot find one of the files you specified please reinput')
		inputCollect()
	else:
		test.close()
#write the input file to a list varible 

def fileListBuilder(input):

	read = open(input, 'r')
	lineCnt = 0
	lines = []
	for line in read:
		lines.append(line)
		lineCnt = lineCnt + 1
	read.close()
	return lines, lineCnt
#TeloPort Pipeline
def teloPortPipeline():

	os.system('mkdir -p Programs/TeloPort/Outputs/' + directory + '/{PipelineOut,muscle_out,interrogate_out/single_reads,blast_out/{blastGenomes,blastMorepeats,blastMoters,blastRdna,blastTel},cons_out}')
	if choice == 's':
		os.system('Programs/TeloPort/build/apps/telomereFinder -s ' + 'Genomes/' + genomeR1 + ' ' + 'Genomes/' + genomeR2 + ' -o ' + 'Programs/TeloPort/Outputs/' + directory + '/PipelineOut/')
	elif choice == 'i':
		os.system('Programs/TeloPort/build/apps/telomereFinder -i ' + genomeI + ' -o ' + 'Programs/TeloPort/Outputs/' + directory + '/PipelineOut/')
	os.system('mv Programs/TeloPort/Outputs/' + directory + '/PipelineOut/pairReads.fastq Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'subtelReads.fastq')
	lines, lineCnt = fileListBuilder('Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'subtelReads.fastq')
	write = open('Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'subtelReads2.fastq', 'a')
	for x in range(lineCnt):
		if lines[x][0] == '@' and x % 4 == 0:
			newLine = '@' + directory + str(x) + '\t' + lines[x].split('@')[1]
			write.write(newLine)
		else:
			write.write(lines[x])
	write.close()
	os.system('rm Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'subtelReads.fastq')
	os.system('mv Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'subtelReads2.fastq Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'subtelReads.fastq')
	os.system('Programs/TeloPort/build/apps/junctionFinder -i Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'subtelReads.fastq -o Programs/TeloPort/Outputs/' + directory + '/PipelineOut/ --revc false --splitJunc true --splitDir false')
	os.system('mv Programs/TeloPort/Outputs/' + directory + '/PipelineOut/telAdjSeq.fastq Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'mergedTelAdjSeqs.fastq')
	os.system('rm Programs/TeloPort/Outputs/' + directory + '/PipelineOut/telSeq.fastq')
	os.system('Programs/TeloPort/build/apps/sequenceQuality -i Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'mergedTelAdjSeqs.fastq -o Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'HiQualTelAdjacent.fasta --ofmt fasta -c 0 -l 30')
	os.system('wcd Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'HiQualTelAdjacent.fasta -l 40 -T 5 -H 0 --show_clusters --histogram > Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'clusters.wcd')
	os.system('Programs/TeloPort/build/apps/wcdInterrogate -w Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'clusters.wcd -i Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'HiQualTelAdjacent.fasta -s Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'mergedTelAdjSeqs.fastq -f fastq -o Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'Clusters.out -r Programs/TeloPort/Outputs/' + directory + '/interrogate_out/ --indices --sort --size 1')
#This will turn all single clusters into single reads
def singleMaker():
	last = lastClusterCount()
	first = firstSingleCount()
	y = len(os.listdir('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/')) - 1
	for cnt in range(y):
		read = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/cluster' + str(cnt) + '.fasta', 'r')
		for line in read:
			if line[0] == '>':
				newLine = line.split('\t')[0] + '\tcluster' + str(cnt) + '\t' + line.split('\t')[1]
				write = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + 'cluster' + str(cnt) + '.fasta', 'a')
				write.write(newLine)
				write.close()
			else:
				write = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + 'cluster' + str(cnt) + '.fasta', 'a')
				write.write(line)
				write.close()
		read.close()
		os.system('rm Programs/TeloPort/Outputs/' + directory + '/interrogate_out/cluster' + str(cnt) + '.fasta')
	while first <= last:
		os.system('cat Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + 'cluster' + str(first) + '.fasta >> Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'single_reads.fasta')
		os.system('rm Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + 'cluster' + str(first) + '.fasta')
		first = first + 1
#this counts the last cluster in the directory
def lastClusterCount():
	intOut = open('Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'Clusters.out', 'r')
	cnt = 0
	for line in intOut:
		if line[0] == 'c':
			cnt = cnt + 1
			last = cnt - 1
	intOut.close()
	return last
#this will count the first single read
def firstSingleCount():
	intOut = open('Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'Clusters.out', 'r')

	global clusterDict
	clusterDict = {}
	for line in intOut:
		if line[0:2] == "cl":
			clusterDict[line.split(';')[0]] = line.split('=')[1]

	intOut = open('Programs/TeloPort/Outputs/' + directory + '/PipelineOut/' + directory + 'Clusters.out', 'r')
	cnt = 0

	for line in intOut:
		if line[0:2] == "cl":
			if int(line.split('=')[1]) <= int(cutoff):
				first = cnt
				break
			cnt = cnt + 1
	return first
#Code to automatically run muscle
def autoMuscle():
	start = 0
	end = firstSingleCount() - 1
	while start <= end:
		os.system('Programs/bin/muscle -in Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + 'cluster' + str(start) + '.fasta -out Programs/TeloPort/Outputs/' + directory + '/muscle_out/cons' + str(start) + '.msa')
		start = start + 1
	start = 0
	while start <= end:
		os.system('cons -sequence Programs/TeloPort/Outputs/' + directory + '/muscle_out/cons' + str(start) + '.msa -outseq Programs/TeloPort/Outputs/' + directory + '/cons_out/clusterCons' + str(start) + '.fasta -name cluster' + str(start))
		start = start + 1
#This will concatonate all consensus sequences into a single file to blast
def catBlast():
	fileCnt = os.listdir('Programs/TeloPort/Outputs/' + directory + '/cons_out/')
	num = len(fileCnt)

	start = 0

	command = ''

	while start < num:
		command = command + ('Programs/TeloPort/Outputs/' + directory + '/cons_out/clusterCons' + str(start) + '.fasta ')
		start = start + 1
	os.system('cat ' + command + '> Programs/TeloPort/Outputs/' + directory + '/blast_out/' + directory + 'catCons.fasta')
	os.system('rm -r Programs/TeloPort/Outputs/' + directory + '/cons_out/')
#Calls NCBI blastn
#def blast(query, subject, output):
#	os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query ' + query +  ' -subject ' + subject +  ' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qlen" >> ' + output)
#Filters false posistives
def falsePosFilter():
	global falsePos
	singles, cnt = fileListBuilder('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'single_reads.fasta')
	true_singles = 0
	for x in range(cnt):
		if singles[x][0] == '>':
			try:
				falsePos[singles[x].split('>')[1].split('\t')[0]]
			except:
				write = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSingles.fasta', 'a')
				write.write(singles[x])
				write.write(singles[x + 1])
				true_singles = true_singles + 1
				write.close()
			else:
				write = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + falsePos[singles[x].split('>')[1].split('\t')[0]][0] + '.fasta', 'a')
				write.write(singles[x])
				write.write(singles[x + 1])
				write.close()
#non telconig dictionary maker
def dictMaker(input, secondVal):
	read = open(input, 'r')
	dict = {}
	for line in read:
		try:
			dict[line.split('\t')[0]]
		except:
			dict[line.split('\t')[0]] = (line.split('\t')[1], int(line.split('\t')[secondVal]))
		else:
			if secondVal == 3:
				if line.split('\t')[3] > falsePos[line.split('\t')[0]][1]:
					dict[line.split('\t')[0]] = (line.split('\t')[1], line.split('\t')[3])
	read.close()
	return dict
#this will generate an updated collection of cluster sizes
def clusterHistogramBuilder():
	fileCnt = os.listdir('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/')
	num = len(fileCnt) - 1
	start = 0
	append = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + 'clusterHistograph.txt','a')
	while start < num:
		clusterCnt = 0
		read = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + 'cluster' + str(start) + '.fasta', 'r')
		for line in read:
			if line[0] == '>':
				clusterCnt = clusterCnt + 1
		command = 'cluster' + str(start) + '\t' + str(clusterCnt) + '\t'
		for x in range(0, clusterCnt):
			command = command + '#'
		command = command + '\n'
		append.write(command)
		read.close()
		start = start + 1
	append.close()
#This will filter out reads which match to Rdna, Moters, and Morepeats and blast the unmatched reads against the genomes
#def moFilter(dictName, writeName, attempt):
def moFilter():
	global moterSingles, RdnaSingles, MorepeatsSingles, ran, true_Singles
	true_Singles = 0
	lines, lineCnt = fileListBuilder('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSingles.fasta')
	for x in range(lineCnt):
		if lines[x][0]  == '>':
			try:
				moterSingles[lines[x].split('\t')[0].split('>')[1]]
			except:

				#filter out rdna
				try:
					RdnaSingles[lines[x].split('\t')[0].split('>')[1]]
				except:

					#filter out morepeats
					try:
						MorepeatsSingles[lines[x].split('\t')[0].split('>')[1]]
					except:
						write = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSinglesGenome.fasta', 'a')
						write.write(lines[x])
						write.write(lines[x + 1])
						write.close()
						true_Singles = true_Singles + 1
					else:
						write = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSinglesMorepeats.fasta', 'a')
						write.write(lines[x])
						write.write(lines[x + 1])
						write.close()

				else:
					write = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSinglesRdna.fasta', 'a')
					write.write(lines[x])
					write.write(lines[x + 1])
					write.close()

			else:
				write = open('Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSinglesMoter.fasta', 'a')
				write.write(lines[x])
				write.write(lines[x + 1])
				write.close()




#			try:
#				dictName[lines[x].split('\t')[0].split('>')[1]]
#			except:
#				if attempt:
#					try:
#						ran[x]
#					except:
#						write = open(writeName, 'a')
#						write.write(lines[x])
#						write.write(lines[x + 1])
#						write.close()
#						true_Singles = true_Singles + 1
#			else:
#				try:
#					ran[x]
#				except:
#					write = open(writeName, 'a')
#					write.write(lines[x])
#					write.write(lines[x + 1])
#					write.close()
#					ran[x] = True
#This will create a director of every telcontig
def telContigDictMaker():
	os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Genomes/db/' + blastD + ' -subject Genomes/db/telRepeats.fasta -dust no -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qlen" >> Programs/TeloPort/Outputs/' + directory + '/blast_out/blastTel/' + directory + 'telContigsOut6.txt')
	read = open('Programs/TeloPort/Outputs/' + directory + '/blast_out/blastTel/' + directory + 'telContigsOut6.txt', 'r')
	global telContigs
	telContigs = {}
	posistionStart = 0
	posistionEnd = 0
	for line in read:

		try:
			telContigs[line.split('\t')[0] + '*s']
		except:
			posistionStart = 0

		try:
			telContigs[line.split('\t')[0] + "*e"]
		except:
			posistionEnd = 0


		if int(line.split('\t')[7]) < (int(line.split('\t')[11]) / 3):
			if posistionStart < int(line.split('\t')[7]):
				telContigs[line.split('\t')[0] + "*s"] = ("start", int(line.split('\t')[7]))
				posistionStart = int(line.split('\t')[7])

		if int(line.split('\t')[7]) > (int(line.split('\t')[11]) / 3):
			if posistionEnd < int(line.split('\t')[7]):
				telContigs[line.split('\t')[0] + "*e"] = ("end", int(line.split('\t')[6]))
				posistionEnd = int(line.split('\t')[6])




	read.close()

	write = open('Programs/TeloPort/Outputs/' + directory + '/blast_out/blastTel/' + directory + 'TelDictionary.txt', 'a')

	write.write(str(telContigs))

	write.close()
#blast interrogation
def blastInterrogation():
	readLines, lineNum = fileListBuilder('Programs/TeloPort/Outputs/' + directory+ '/blast_out/blastGenomes/' + directory + 'blastGenomeOut6.txt')
	command = ""
	for cnt in range(lineNum):
		telContigCheck = False
		command = ""

		queryName = readLines[cnt].split('\t')[0]
		print(queryName)

		subjectName = readLines[cnt].split('\t')[1]
		print(subjectName)
		telDisS = 0

		try:
			telContigs[subjectName + "*s"]
		except:
			command = "\tTelNs"
		else:
			command = "\tTelYs"
			telContigCheck = True


		if telContigCheck:
			posistion = telContigs[readLines[cnt].split('\t')[1] + '*s'][1]
			telDisS =  abs(int(readLines[cnt].split('\t')[8]) - posistion)
		telContigCheck = False
		telDisE = 0

		try:
			telContigs[subjectName + "*e"]
		except:
			command = command + "\tTelNe"
		else:
			command = command + "\tTelYe"
			telContigCheck = True


		if telContigCheck:
			posistion = telContigs[readLines[cnt].split('\t')[1] + '*e'][1]
			telDisE =  abs(int(readLines[cnt].split('\t')[8]) - posistion)

		telContigCheck = False
		if abs(telDisS) < abs(telDisE) and telDisS != 0:
			command = command + "\t" + str(abs(telDisS)) + "\t"
		elif abs(telDisS) > abs(telDisE) and telDisE == 0:
			command = command + "\t" + str(abs(telDisS)) + "\t"
		elif abs(telDisE) < abs(telDisS) and telDisE != 0:
			command = command + "\t" + str(abs(telDisE)) + "\t"
		elif abs(telDisE) > abs(telDisS) and telDisS == 0:
			command = command + "\t" + str(abs(telDisE)) + "\t"
		else:
			command = command + "\t" + str(abs(telDisE)) + "\t"

		try:
			moterContigs[subjectName]
		except:
			command = command + "MoterN\t"
		else:
			command = command + "MoterY\t"

		revS = False

		revQ = False
		if int(readLines[cnt].split('\t')[6]) > int(readLines[cnt].split('\t')[7]):
			revS = True

		if int(readLines[cnt].split('\t')[8]) > int(readLines[cnt].split('\t')[9]):
			revQ = True
		rev = "revF"

		if not revS and revQ:
			rev = "revT"
		command = command + rev +'\n'

		command = readLines[cnt].split('\n')[0] + command

		append = open('Programs/TeloPort/Outputs/' + directory+ '/blast_out/blastGenomes/' + directory + 'blastGenomeIntOut6.txt', 'a')

		append.write(command)

		append.close()

		cnt = cnt + 1
#Builds a Gff file
def gffBuilder():
	lines, linecnt = fileListBuilder('Programs/TeloPort/Outputs/' + directory+ '/blast_out/blastGenomes/' + directory + 'blastGenomeIntOut6.txt')
	write = open('Programs/TeloPort/Outputs/' + directory+ '/' + directory + '.gff', 'a')
	command = '##gff-version 3.1.26\n'
	write.write(command)
	for x in range(linecnt):
		if int(lines[x].split('\t')[8]) < int(lines[x].split('\t')[9]):
			strand = '+'
		elif int(lines[x].split('\t')[8]) > int(lines[x].split('\t')[9]):
			strand = '-'
		length = abs(int(lines[x].split('\t')[9]) - int(lines[x].split('\t')[8]))
		command = lines[x].split('\t')[1] + '\t' + 'teloPortWrapper' + '\t' + 'telomere' + '\t' + lines[x].split('\t')[8] + '\t' + lines[x].split('\t')[9] + '\t' + lines[x].split('\t')[10] + '\t' + strand + '\t' + '.' + '\t' + 'ID=Telomere' + str(x) + ';Name=' + lines[x].split('\t')[0] + ';Length=' + str(length) + ';Distance from Telomere=' + lines[x].split('\t')[14] + '\n' 

		write.write(command)
#generates results txt file
def resultsBuilder():
	command = directory + '\n'

	rawReadCnt = 0

	if choice == 's':
		read = open('Genomes/' +genomeR1, 'r')
		for line in read:
			if line[0] == '@':
				rawReadCnt = rawReadCnt + 1
		read.close()
		read = read = open('Genomes/' + genomeR2, 'r')
		for line in read:
			if line[0] == '@':
				rawReadCnt = rawReadCnt + 1
		read.close()
	elif choice == 'i':
		read = open('Genomes/' + genomeI, 'r')
		for line in read:
			if line[0] == "@":
				rawReadCnt = rawReadCnt + 1
		read.close()
	command = command + ('Number of raw reads processed: ' + str(rawReadCnt) + '\n')

	teloCnt = 0

	read = open('Programs/TeloPort/Outputs/' + directory + '/PipelineOut/telReads.fastq', 'r')

	for line in read:
		if line[0] == '@':
			teloCnt = teloCnt + 1
	read.close()

	os.system('rm Programs/TeloPort/Outputs/' + directory + '/PipelineOut/telReads.fastq')

	command = command + ('Number of Telomere reads processed: ' + str(teloCnt) + '\n')

	command = command + ('Number of de novo telomeres: ' + str(true_Singles) + '\n')

	clusterCnt = len(os.listdir('Programs/TeloPort/Outputs/' + directory + '/interrogate_out')) - 2

	command = command + ('Number of clusters: ' + str(clusterCnt) + '\n')

	genomesCnt = 0

	read = open('Programs/TeloPort/Outputs/' + directory + '/blast_out/blastGenomes/' + directory + 'blastGenomeIntOut6.txt', 'r')

	for line in read:
		genomesCnt = genomesCnt + 1
	read.close()

	command = command + ('Number of de novo telomeres that match up with the genome: ' + str(genomesCnt) + '\n')

	moterCnt = 0

	read = open('Programs/TeloPort/Outputs/' + directory + '/blast_out/blastMoters/' + directory + 'blastMoterOut6.txt', 'r')

	for line in read:
		moterCnt = moterCnt + 1
	read.close()

	command = command + ('Number of de novo telomeres that match up with the MOTER sequences: ' + str(moterCnt) + '\n')

	comand = command + 'Time speant running ' + (str(time.time() - startTime)) + '\n'

	command = command + ('Distribution of clusters \n')

	append = open('Programs/TeloPort/Outputs/' + directory + '/' + directory + 'results.txt', 'a')

	append.write(command)

	append.close()

	os.system('cat Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + 'clusterHistograph.txt >> Programs/TeloPort/Outputs/' + directory + '/' + directory + 'results.txt')

	os.system('rm Programs/TeloPort/Outputs/' + directory + '/interrogate_out/' + directory + 'clusterHistograph.txt')

	os.system('cat Programs/TeloPort/Outputs/' + directory+ '/blast_out/blastGenomes/' + directory + 'blastGenomeIntOut6.txt >> Programs/TeloPort/Outputs/' + directory + '/' + directory + 'results.txt')




if args.simple:
	inputCollect()
elif args.b == "x" and args.i == "x":
	inputCollect()
elif args.s != 'x':
	genomeR1 = args.s[0]
	genomeR2 = args.s[1]
	if args.d == 'x':
		directory = args.b
	else:
		directory = args.d
	blastD = args.b
	cutoff = args.cut
	choice = 's'
elif args.i != 'x':
	genomeI = args.i
	if args.d == 'x':
		directory = args.b
	else:
		directory = args.d
	blastD = args.b
	cutoff = args.cut
	choice = 'i'
startTime = time.time()
teloPortPipeline()
singleMaker()
autoMuscle()
catBlast()
os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'single_reads.fasta -subject Programs/TeloPort/Outputs/' + directory + '/blast_out/' + directory + 'catCons.fasta -outfmt 6 > Programs/TeloPort/Outputs/' + directory + '/blast_out/' + directory + 'blastCons.txt')
falsePos = dictMaker('Programs/TeloPort/Outputs/' + directory + '/blast_out/' + directory + 'blastCons.txt', 3)
falsePosFilter()
clusterHistogramBuilder()
os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSingles.fasta -subject Genomes/db/MoTeRs.gb -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue slen" >> Programs/TeloPort/Outputs/' + directory + '/blast_out/blastMoters/' + directory + 'blastMoterOut6.txt')
os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSingles.fasta -subject Genomes/db/Mo_rDNA.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue slen" >> Programs/TeloPort/Outputs/' + directory + '/blast_out/blastRdna/' + directory + 'blastRdnaOut6.txt')
os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSingles.fasta -subject Genomes/db/MoRepeats.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue slen" >> Programs/TeloPort/Outputs/' + directory + '/blast_out/blastMorepeats/' + directory + 'blastMorepeatsOut6.txt')
moterSingles = dictMaker('Programs/TeloPort/Outputs/' + directory + '/blast_out/blastMoters/' + directory + 'blastMoterOut6.txt', 7)
RdnaSingles = dictMaker('Programs/TeloPort/Outputs/' + directory + '/blast_out/blastRdna/' + directory + 'blastRdnaOut6.txt', 7)
MorepeatsSingles = dictMaker('Programs/TeloPort/Outputs/' + directory + '/blast_out/blastMorepeats/' + directory + 'blastMorepeatsOut6.txt', 7)
ran = {}
moFilter()
#moFilter(moterSingles, 'Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSinglesMoter.fasta', False)
#moFilter(RdnaSingles, 'Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSinglesRdna.fasta', False)
#moFilter(MorepeatsSingles, 'Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSinglesMorepeats.fasta', False)
#moFilter(moterSingles, 'Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSinglesGenome.fasta', True)
os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Programs/TeloPort/Outputs/' + directory + '/interrogate_out/single_reads/' + directory + 'trueSinglesGenome.fasta -subject Genomes/db/' + blastD + ' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue slen" >> Programs/TeloPort/Outputs/' + directory + '/blast_out/blastGenomes/' + directory + 'blastGenomeOut6.txt')
telContigDictMaker()
os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Genomes/db/' + blastD + ' -subject  Genomes/db/MoTeRs.gb -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qlen" >> Programs/TeloPort/Outputs/' + directory + '/blast_out/blastMoters/' + directory + 'moterContigsOut6.txt')
moterContigs = dictMaker('Programs/TeloPort/Outputs/' + directory + '/blast_out/blastMoters/' + directory + 'moterContigsOut6.txt', 7)
write = open('Programs/TeloPort/Outputs/' + directory + '/blast_out/blastMoters/' + directory + 'MotersDictionary.txt', 'a')
write.write(str(moterContigs))
write.close()
blastInterrogation()
gffBuilder()
resultsBuilder()

#clean up

write = open('Programs/TeloPort/Outputs/' + directory + '/oldClusters.txt', 'a')

for x in range(len(clusterDict)):
	write.write('cluster' + str(x) + ': ' + clusterDict['cluster ' + str(x)])


print('final time : ' + str(time.time() - startTime))
