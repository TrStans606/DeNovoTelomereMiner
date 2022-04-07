import os 


#determines if an interleaved of seperated input is needed

choiceMade = False

while not choiceMade:
	choice = input('Please enter if you are using sperated (s) genome files or an interleaved (i) file ')
	if choice == 'i' or choice == 's':
		choiceMade = True
	else:
		print('that is not a valid input. Please enter i for interleaved or s for seperated')
if choice == 's':

	genomeR1 = input('Specify the name of the R1 file you want: ')

	genomeR2= input('Specify the name of the R2 file you want: ')

elif choice == 'i':
	genomeI = input('Specify the name of the interleaved fastq file ')



directory = input('Specify the output directory name ')

blastD = input('Specify genome file to be blasted against the single reads ')

#test to make sure the file path is valid

if choice == 's':
	test= open('Genomes/' + genomeR1, 'r')
	print(test.readline())
	test.close()
	test= open('Genomes/' + genomeR2, 'r')
	print(test.readline())
elif choice == 'i':
	test = open('Genomes/' + genomeI, 'r')
	print(test.readline())

#calls all TeloPort programs up to wcd interrogate

os.system('mkdir Programs/TeloPort/' + directory)

if choice == 's':
	os.system('Programs/TeloPort/build/apps/telomereFinder -s ' + 'Genomes/' + genomeR1 + ' ' + 'Genomes/' + genomeR2 + ' -o ' + 'Programs/TeloPort/' + directory + '/find_out')
elif choice == 'i':
	os.system('Programs/TeloPort/build/apps/telomereFinder -i ' + genomeI + ' -o ' + 'Programs/TeloPort/' + directory + '/find_out')

os.system('Programs/TeloPort/build/apps/junctionFinder -i Programs/TeloPort/' + directory + '/find_out/pairReads.fastq -o Programs/TeloPort/' + directory + '/junction_out --revc true --splitJunc true --splitDir false')

os.system('Programs/TeloPort/build/apps/junctionFinder -i Programs/TeloPort/' + directory + '/find_out/pairReads.fastq -o Programs/TeloPort/' + directory + '/full_telo --revc true --splitJunc false --splitDir false')

os.system('mkdir Programs/TeloPort/' + directory + '/quality_out/') 

os.system('Programs/TeloPort/build/apps/sequenceQuality -i Programs/TeloPort/' + directory + '/junction_out/telAdjSeq.fastq -o Programs/TeloPort/' + directory + '/quality_out/' + directory + '.fasta --ofmt fasta -c 60 -l 40')


os.system('mkdir Programs/TeloPort/' + directory + '/wcd_out')

os.system('wcd Programs/TeloPort/' + directory + '/quality_out/' + directory + '.fasta -l 40 -T 5 -H 0 --show_clusters --histogram > Programs/TeloPort/' + directory + '/wcd_out/clusters.wcd')

os.system('Programs/TeloPort/build/apps/wcdInterrogate -w Programs/TeloPort/' + directory + '/wcd_out/clusters.wcd -i Programs/TeloPort/' + directory + '/quality_out/' + directory + '.fasta -s Programs/TeloPort/' + directory + '/junction_out/telAdjSeq.fastq -f fastq -o Programs/TeloPort/' + directory + '/Clusters.out -r Programs/TeloPort/' + directory + '/interrogate_out/ --indices --sort --size 1')

#This will turn all single clusters into single reads
 
intOut = open('Programs/TeloPort/' + directory + '/Clusters.out', 'r')

#this counts the last cluster in the directory

cnt = 0

for line in intOut:
	if line[0] == 'c':
		cnt = cnt + 1
last = cnt - 1

intOut.close()

#this will count the first single read

intOut = open('Programs/TeloPort/' + directory + '/Clusters.out', 'r')

cnt = 0

for line in intOut:
	if len(line) == 19:
		if line[12:18] == 'size=1':
			first = cnt
			break
	elif len(line) == 18:
		if line[11:17] == 'size=1':
			first = cnt
			break
	if line[0] == 'c':
		cnt = cnt + 1
firstEX = first

#this is code from wcdAutoSort that will rename all single clusters into single reads

os.system('mkdir Programs/TeloPort/' + directory + '/interrogate_out/single_reads/')

singleRead = 0


while first <= last:
        os.system('mv Programs/TeloPort/' + directory + '/interrogate_out/cluster' + str(first) + '.fasta Programs/TeloPort/' + directory + '/interrogate_out/single_reads/single_read' + str(singleRead) + '.fasta')
        first = first + 1
        singleRead = singleRead + 1


#code from AutoMuscle to create a consensus sequence for every cluster

start = 0

end = firstEX - 1

os.system('mkdir Programs/TeloPort/' + directory + '/muscle_out')

os.system('mkdir Programs/TeloPort/' + directory + '/cons_out')

while start <= end:
        os.system('Programs/bin/muscle -in Programs/TeloPort/' + directory + '/interrogate_out/cluster' + str(start) + '.fasta -out Programs/TeloPort/' + directory + '/muscle_out/cons' + str(start) + '.msa')
        start = start + 1

start = 0

while start <= end:
        os.system('cons -sequence Programs/TeloPort/' + directory + '/muscle_out/cons' + str(start) + '.msa -outseq Programs/TeloPort/' + directory + '/cons_out/clusterCons' + str(start) + '.fasta -name cluster' + str(start))
        start = start + 1

#This will concatonate all consensus sequences into a single file to blast

fileCnt = os.listdir('Programs/TeloPort/' + directory + '/cons_out/')

num = len(fileCnt)

start = 0

command = ''

while start < num:
	command = command + ('Programs/TeloPort/' + directory + '/cons_out/clusterCons' + str(start) + '.fasta ')
	start = start + 1

os.system('mkdir Programs/TeloPort/' +  directory + '/blast_out/') 

os.system('cat ' + command + '> Programs/TeloPort/' + directory + '/blast_out/catCons.fasta')


#this will blast each single read against the cat file

fileCnt = os.listdir('Programs/TeloPort/' + directory + '/interrogate_out/single_reads')

num = len(fileCnt)

start = 0

while start < num:
	os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Programs/TeloPort/' + directory + '/interrogate_out/single_reads/single_read' + str(start) + '.fasta -subject Programs/TeloPort/' + directory + '/blast_out/catCons.fasta -outfmt 6 > Programs/TeloPort/' + directory + '/blast_out/blastCons' + str(start) + '.txt') 
	start = start + 1


#this will read each blast outfmt 6 file and sort all those with no matches to a true single reads folder

os.system('mkdir Programs/TeloPort/' + directory + '/interrogate_out/single_reads/true_singles')

start = 0

true_single = 0

while start < num:
	read = open('Programs/TeloPort/' + directory + '/blast_out/blastCons' + str(start) + '.txt', 'r')
	one_char = read.read(1)
	read.close()
	read = open('Programs/TeloPort/' + directory + '/blast_out/blastCons' + str(start) + '.txt', 'r')
	lineCnt = read.readlines()
	if not one_char or int(lineCnt[0].split('\t')[6]) > 9 :
		os.system('mv Programs/TeloPort/' + directory + '/interrogate_out/single_reads/single_read' + str(start) + '.fasta Programs/TeloPort/' + directory + '/interrogate_out/single_reads/true_singles/true_single' + str(true_single) + '.fasta')
		true_single = true_single + 1
	else:
		if len(lineCnt) == 1:
			os.system('cat Programs/TeloPort/' + directory + '/interrogate_out/single_reads/single_read' + str(start) + '.fasta >> Programs/TeloPort/' + directory + '/interrogate_out/' + lineCnt[0].split('\t')[1] + '.fasta')
		else:
			read.close()
			lineMax = 0.0
			read = open('Programs/TeloPort/' + directory + '/blast_out/blastCons' + str(start) + '.txt', 'r')
			for line in read:
				if float(line.split('\t')[3]) > lineMax:
					lineMax = float(line.split('\t')[3])
					hold = line.split('\t')[1]
			os.system('cat Programs/TeloPort/' + directory + '/interrogate_out/single_reads/single_read' + str(start) + '.fasta >> Programs/TeloPort/' + directory + '/interrogate_out/' + hold + '.fasta')
	read.close()
	start = start + 1

#concatonate all true single reads

start = 0

command = ""

num = len(os.listdir('Programs/TeloPort/' + directory + '/interrogate_out/single_reads/true_singles/'))

while start < num:
        command = command + ' Programs/TeloPort/' + directory + '/interrogate_out/single_reads/true_singles/true_single' + str(start) + '.fasta'
        start = start + 1
os.system('cat' + command + ' >> Programs/TeloPort/' + directory + '/interrogate_out/single_reads/trueSingles.fasta')


#this will generate an updated collection of cluster sizes

fileCnt = os.listdir('Programs/TeloPort/' + directory + '/interrogate_out/')

num = len(fileCnt) - 1

start = 0


append = open('Programs/TeloPort/' + directory + '/interrogate_out/clusterHistograph.txt','a')


while start < num:
	clusterCnt = 0
	read = open('Programs/TeloPort/' + directory + '/interrogate_out/cluster' + str(start) + '.fasta', 'r')
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


#This will blast the true singles against the genome

os.system('mkdir Programs/TeloPort/' + directory + '/blast_out/blastGenomes/')


os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Programs/TeloPort/' + directory + '/interrogate_out/single_reads/trueSingles.fasta -subject Genomes/db/' + blastD + ' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue slen" >> Programs/TeloPort/' + directory + '/blast_out/blastGenomes/blastGenomeOut6.txt')

os.system('mkdir Programs/TeloPort/' + directory + '/blast_out/blastMoters')

os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Programs/TeloPort/' + directory + '/interrogate_out/single_reads/trueSingles.fasta -subject Genomes/db/MoTeRs.gb -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue slen" >> Programs/TeloPort/' + directory + '/blast_out/blastMoters/blastMoterOut6.txt')

#This will create a director of every telcontig

os.system('mkdir Programs/TeloPort/' + directory + '/blast_out/blastTel/')

os.system('Programs/ncbi-blast-2.11.0+/bin/blastn -query Genomes/db/' + blastD + ' -subject Genomes/db/telRepeats.fasta -dust no -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue slen" >> Programs/TeloPort/' + directory + '/blast_out/blastTel/telContigsOut6.txt')

read = open('Programs/TeloPort/' + directory + '/blast_out/blastTel/telContigsOut6.txt', 'r')

telContigs = {}

posistion = ""

telContigCnt = 1

for line in read:
	if int(line.split('\t')[8]) < 25:
		posistion = "front"
	elif int(line.split('\t')[8]) > int(line.split('\t')[11]) - 25:
		posistion = "end"
	else:
		posistion = "middle"

	try:
		telContigs[line.split('\t')[0]]
	except:
		telContigs[line.split('\t')[0]] = posistion
		telContigCnt = 1
	else:
		telContigCnt = telContigCnt + 1

	telContigs[line.split('\t')[0]] = (posistion, telContigCnt)

print(telContigs)

#generates results txt file

append = open('Programs/TeloPort/' + directory + '/results.txt', 'a')

append.close()

true_singles = len(os.listdir('Programs/TeloPort/' + directory + '/interrogate_out/single_reads/true_singles'))

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

read = open('Programs/TeloPort/' + directory + '/find_out/telReads.fastq', 'r')

for line in read:
	if line[0] == '@':
		teloCnt = teloCnt + 1
read.close()

command = command + ('Number of Telomere reads processed: ' + str(teloCnt) + '\n')

command = command + ('Number of de novo telomeres: ' + str(true_singles) + '\n')

clusterCnt = len(os.listdir('Programs/TeloPort/' + directory + '/interrogate_out')) - 2

command = command + ('Number of clusters: ' + str(clusterCnt) + '\n')

genomesCnt = 0

read = open('Programs/TeloPort/' + directory + '/blast_out/blastGenomes/blastGenomeOut6.txt', 'r')

for line in read:
	genomesCnt = genomesCnt + 1
read.close()

command = command + ('Number of de novo telomeres that match up with the genome: ' + str(genomesCnt) + '\n')

moterCnt = 0

read = open('Programs/TeloPort/' + directory + '/blast_out/blastMoters/blastMoterOut6.txt', 'r')

for line in read:
	moterCnt = moterCnt + 1
read.close()

command = command + ('Number of de novo telomeres that match up with the MOTER sequences: ' + str(moterCnt) + '\n')

command = command + ('Distribution of clusters \n')

append = open('Programs/TeloPort/' + directory + '/results.txt', 'a')

append.write(command)

append.close()

os.system('cat Programs/TeloPort/' + directory + '/interrogate_out/clusterHistograph.txt >> Programs/TeloPort/' + directory + '/results.txt')


#blast interrogation 

append = open('Programs/TeloPort/' + directory+ '/blast_out/blastGenomes/blastGenomeOut6.txt', 'a')

append.write('M01478:141:000000000-ACU7R:1:2109:10561:22207	Bm88324_contig02502	1451\n')

append.write('M01478:141:000000000-ACU7R:1:2109:10561:22207	Bm88324_scaffold00150	yhjhbjb')

append.close()

read = open('Programs/TeloPort/' + directory+ '/blast_out/blastGenomes/blastGenomeOut6.txt', 'r')


lineNum = 0

for line in read:
	lineNum = lineNum + 1

read.close()
cnt = 0

telContigCheck = False

while cnt < lineNum:

	read = open('Programs/TeloPort/' + directory+ '/blast_out/blastGenomes/blastGenomeOut6.txt', 'r')

	readLines = read.readlines()

	queryName = readLines[cnt].split('\t')[0]
	print(queryName)

	subjectName = readLines[cnt].split('\t')[1]
	print(subjectName)

	read.close()

	try:
		telContigs[subjectName]
	except:
		write = open('Programs/TeloPort/' + directory+ '/blast_out/blastInt.txt', 'a')
		write.write('Telomere ' + queryName + ' is not a tel contig \n \n')
	else:
		write = open('Programs/TeloPort/' + directory+ '/blast_out/blastInt.txt', 'a')
		write.write('Telomere ' + queryName + ' matched to tel contig ' + subjectName + ' \n \n')

	cnt = cnt + 1



