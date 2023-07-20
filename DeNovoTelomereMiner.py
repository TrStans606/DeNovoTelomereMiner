import subprocess
import os
import argparse
import glob
import urllib.parse


#Command line argument set up
parser = argparse.ArgumentParser(
    prog='DenovoTelomereMiner',
    argument_default='try --simple',
    description=(
        'An all-in-one pipeline and anaylsis tool for de-novo Telomeres'
    ),
    epilog=(
        'The -s and -i options are mutally exclusive. '
        'If simple mode is not used options -b is required.'
    ),
)

#command line option for seperated files 
parser.add_argument(
    '-s',
    action='store',
    nargs=2,
    help=(
        'User defined name for two seperated fastq files.' 
        'This argument exepcts two file names seperated by a space.'
    ),
    default='n/a'
)

#command line option for interleaved files
parser.add_argument(
    '-i',
    action='store',
    help=(
        'User defined name for an interleaved fastq file'
        'This argument expects one file name.'
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
        'The name of the directory where all output files will be placed.'
        'The name defaults to the name of the fasta file provided'
        'Directory names must be unique'
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
        'List all of the seqeunces you want removed from reads.'
        'Must be in fasta format.'    
    ),
    default='n/a'
)

#command line argument for determining if a cluster is a canidate de-novo 
parser.add_argument(
	'--cut', 
    action='store', 
    type=int, 
    help= ('The number of reads in a cluster before it is'
           'labeled a canidate de-novo.'
           'The default is five.'
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
    help=('List all addition seqeunces you want compared to the reads.',
          'All must be in fasta format'
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


args = parser.parse_args()

#config file maker

def config():
    if args.config:
        print('Manual Time setup')
        readDir = input('Raw Reads Directory: ')
        genomeDir = input('Assembled genome Directory: ')
        filterDir = input('Filters directory: ')
        addDir = input('Additional Blast directory: ')
        write = open('config.ini', 'w')
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
    os.remove(
        f'Outputs/{directory}/teloPortOut/pairReads.fastq'
        )
    junctionFinder()
    sequenceQuality()
    wcdestCall()
    wcdInterrogate()
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
    subprocess.run(' '.join(command),
        shell=True,
        check=True)
#tags every read with a unique id number for tracking
def idTagger():
    pairLines, pairLineCnt = fileListBuilder(
        os.path.join('Outputs',
            directory,
            'teloPortOut/pairReads.fastq')
        )
    newFile = os.path.join('Outputs',
        directory,
        'teloPortOut/subTelReads.fastq')
    with open(newFile, 'a') as write:
        x = 0
        for line in pairLines:
            if x % 4 == 0:
                newLine = f"@{directory}{x}\t{line.lstrip('@')}"
                x += 1
                write.write(newLine)
#calls junctionfinder
def junctionFinder():
    command = ['Programs/bin/junctionFinder',
                    '-i',
                    os.path.join('Outputs',
                                 directory,
                                 'teloPortOut/subTelReads.fastq'),
                    '-o',
                    os.path.join('Outputs',
                                 directory,
                                 'teloPortOut/'),
                    '--revc 0',
                    '--splitJunc 1',
                    '--splitDir 0']
    subprocess.run(' '.join(command),
                   shell=True,
                   check=True)
#calls sequenceQuality
def sequenceQuality():
    command = ['Programs/bin/sequenceQuality',
                    '-i',
                    os.path.join('Outputs',
                                 directory,
                                 'teloPortOut/telAdjSeq.fastq'),
                    '-o',
                    os.path.join('Outputs',
                                 directory,
                                 'teloPortOut/hiQualityTelAdjSeq.fastq'),
                    '--ofmt',
                    'fasta',
                    '-c',
                    '0',
                    '-l',
                    '0']
    subprocess.run(' '.join(command),
                   shell=True,
                   check=True)
#calls WCDest
def wcdestCall():
    command = ['Programs/bin/wcd',
                    os.path.join('Outputs',
                                 directory,
                                 'teloPortOut/hiQualityTelAdjSeq.fastq'),
                    '-l',
                    '-40',
                    '-T',
                    '5',
                    '-H',
                    '0',
                    '--show_clusters',
                    '--histogram',
                    '>',
                    os.path.join(
                        'Outputs',
                        directory,
                        'teloPortOut',
                        f'{directory}clusters.wcd')]
    subprocess.run(' '.join(command),
                    shell=True,
                    check=True)
#calls wcd interrogate
def wcdInterrogate():
    command = ['Programs/bin/wcdInterrogate',
                    '-w',
                    os.path.join('Outputs',
                                 directory,
                                 'teloPortOut',
                                 f'{directory}clusters.wcd'),
                    '-i',
                    os.path.join('Outputs',
                                 directory,
                                 'teloPortOut',
                                 'hiQualityTelAdjSeq.fastq'),
                    '-s',
                    os.path.join('Outputs',
                                 directory,
                                 'teloPortOut',
                                 'telAdjSeq'),
                    '-f',
                    'fastq',
                    '-o',
                    os.path.join('Outputs',
                                 directory,
                                 'teloPortOut',
                                 f'{directory}Clusters.out'),
                    '-r',
                    os.path.join('Outputs',
                                 directory,
                                'telomereClusters',
                                'telomereReads/'),
                    '--indices',
                    '--sort',
                    '--size'
                    '1']
    subprocess.run(' '.join(command),
                   shell=True,
                   check=True)
#renames all canidate denovos from clusters
def deNovoRename():
    cntCluster = len(glob.glob(os.path.join('Outputs/',
                                  directory,
                                  'telomereReads',
                                  'telomereClusters',
                                  'cluster*.fasta')))
    firstDeNovo = firstDeNovoCount()
    cluster = os.path.join('Outputs/',
                                  directory,
                                  'telomereReads',
                                  'telomereClusters',
                                  'cluster')
    for i in range(firstDeNovo,cntCluster):
        with open(f'{cluster}{i}.fasta', 'r') as read:
            deNovoLines = read.readlines()
        with open(os.path.join('Outputs/',
                  directory,
                  'telomereReads',
                  'deNovoTelomeres',
                  f'{directory}deNovos.fasta'), 'a') as write:
            for line in deNovoLines:
                if line[0]=='>':
                    write.write(f'{directory}\t{line.lstrip(">")}')
                else:
                    write.write(line)
            os.remove(f'{cluster}{i}.fasta')
        newClusterCnt = len(glob.glob(
            os.path.join('Outputs/',
                                      directory,
                                      'telomereReads',
                                      'telomereClusters',
                                      'cluster*.fasta')))
        for i in range(0,newClusterCnt):
            with open(f'{cluster}{i}.fasta', 'r') as read:
                clusterLines = read.readlines()
            os.remove(f'{cluster}{i}.fasta')
            with open(f'{cluster}{i}.fasta', 'a') as write:
                for line in clusterLines:
                    if line[0]=='>':
                        write.write(f'{directory}\t{line.lstrip(">")}')
                    else:
                        write.write(line)
                  
                                  
def firstDeNovoCount():
    with open(os.path.join('Outputs',
                           directory,
                           'teloPortOut',
                           f'{directory}Clusters.out')) as intOut:
        clusterLines = intOut.readlines()
    cnt = 0
    for line in clusterLines:
        if line[0:1]=='cl':
            if line.split('=')[1] <= str(cutOff):
                break    
            cnt += 1
        return cnt

#runs muscle and EMBOSS cons
def autoMuscle():
    clusters = os.path.join('Outputs/',
                                        directory,
                                        'telomereReads',
                                        'telomereClusters')
    cntCluster = len(glob.glob(f'{clusters}cluster*.fasta'))
    for i in range(0,cntCluster):
        command = ['Programs/bin/muscle5',
                        '-allign',
                        f'{clusters}cluster{i}.fastq',
                        '-output',
                        os.path.join('Outputs/',
                                     directory,
                                     'muscleOut',
                                     f'{directory}cluster{i}.msa')]
        subprocess.run(' '.join(command),
                       shell=True,
                       check=True)

    allignments = os.path.join('Outputs/',
                               directory, 
                               'muscleOut')
    cons = os.path.join('Outputs/',
                               directory, 
                               'consOut')
    cntAllignment = len(glob.glob(f'{allignments}/{directory}cluster*.msa'))
    for i in range(0, cntAllignment):
        command = ['Programs/bin/cons',
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
    cons = os.path.join('Outputs/',
                               directory, 
                               'consOut')
    blast = os.path.join('Outputs/',
                               directory, 
                               'blastOut')
    singles = os.path.join('Outputs/',
                               directory, 
                               'telomereReads',
                               'deNovoTelomeres')
    cntCons = len(glob.glob(f'{cons}/cons*.fasta'))
    for i in range(0,cntCons):
        with open(f'{cons}/cons{i}.fasta', 'r') as read:
            consLines = read.readlines()
        with open(f'{blast}/{directory}catCons.fasta', 'a') as write:
            write.write(consLines)
    blast(f'{singles}/deNovos.fasta',
          f'{blast}/{directory}deNovos.fasta',
          f'{blast}/{directory}blastCons.txt',
          '')
#runs blast            
def blast(query, subject, output, dust):
    outfmt ='-outfmt \
                "6 qseqid sseqid pident length mismatch \
                gapopen qstart qend sstart send evalue qlen"'
    command = ['Programs/bin/blastn',
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
    clusters = os.path.join('Outputs/',
                               directory,
                               'telomereReads',
                               'telomereClusters')
    deNovoPath = os.path.join('Outputs/',
                               directory, 
                               'telomereReads',
                               'deNovoTelomeres')
    trueDeNovos = 0
    falsePos = dictMaker(os.path.join('Outputs/',
                                      directory,
                                      'blastOut',
                                      f'{directory}blastCons.txt'), 
                         3)
    singles, cnt = fileListBuilder(
        os.path.join('Outputs/',
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
    return trueDeNovos
                                            
def dictMaker(file, secondVal):               
    with open(file, 'r') as read:
        dictionary = {}
        for line in read:
            item = line.split('\t')
            if item[0] not in dictionary:
                dictionary[item[0]] = (item[1], item[secondVal])
            elif secondVal == 3:
                if item[3] > dictionary[item[1]][1]:
                    dictionary[item[0]]=(item[1], item[3])
                    
def histogramBuilder():
    clusters = os.path.build('Outputs/',
                 directory,
                 'telomereReads',
                 'telomereClusters')
    clustersCnt = len(glob.glob(f'{clusters}/cluster*.fasta'))
    with open(os.path.join('Outputs',
                           directory,
                           '/clusterHistogram.txt', 'a')) as \
                              write:
                                  write.write('Distribution of Clusters')
    for x in range(0,clustersCnt):
        readCnt = 0
        with open(f'{clusters}/cluster{x}.fasta', 'r') as read:
            for line in read:
                if line[0] == '>':
                    readCnt += 1
        with open(os.path.join('Outputs',
                               directory,
                               '/clusterHistogram.txt'), 'a') as \
                                  write:
                                      line = f'cluster{x}\t{readCnt}\t'
                                      line = f'{line}#' * readCnt
                                      write.write(line)
                
def deNovoFilter():
    if filtersCnt > 0:
        for file in args.f:
            blast(os.path.join('Outputs',
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
                os.path.join('Outputs/',
                  directory,
                  'telomereReads',
                  'deNovoTelomeres',
                  f'{directory}trueDeNovos.fasta'))
        
        for i in range(0,lineCnt):
            if lines[i][0] == '>':
                ID = lines[i].split('\t')[0].lstrip('>')
                if ID in filterSingles:
                    with open(os.path.join('Outputs/',
                      directory,
                      'telomereReads',
                      'deNovoTelomeres',
                      f'{directory}filteredDeNovos.fasta'), 'a') as write:
                        
                        write.write(lines[i])
                else:
                    with open(os.path.join('Outputs/',
                      directory,
                      'telomereReads',
                      'deNovoTelomeres',
                      f'{directory}blastDeNovos.fasta'), 'a') as write:
                        filteredDeNovos += 1
                        write.write(lines[i])
        return filteredDeNovos
    else:
        lines, lineCnt = fileListBuilder(
                os.path.join('Outputs/',
                  directory,
                  'telomereReads',
                  'deNovoTelomeres',
                  f'{directory}trueDeNovos.fasta'))
        with open(os.path.join('Outputs/',
          directory,
          'telomereReads',
          'deNovoTelomeres',
          f'{directory}blastDeNovos.fasta'), 'a') as write:
            filteredDeNovos = lineCnt / 2
            write.write(lines)
        return filteredDeNovos

def telContigDictMaker():
    with open(os.path.join('Outputs/',
                           directory,
                           'blastOut',
                           'blastDict',
                           f'{directory}telRepeats'), 'a') as write:
        write.write(f'>{blastGenome} telomere repeats')
        write.write(telRepeat * 8)
    blast(f'{genomeDir}/{blastGenome}', 
          os.path.join('Outputs/',
                       directory,
                       'blastOut',
                       'blastDict',
                       f'{directory}telRepeats'),
          os.path.join('Outputs/',
                       directory,
                       'blastOut',
                       'blastDict',
                       f'{directory}blastDictOut6.txt'),
          '-dust no')
    telContigs = {}
    posistionStart = 0
    posistionEnd = 0
    with open(os.path.join('Outputs/',
                 directory,
                 'blastOut',
                 'blastDict',
                 f'{directory}blastDictOut6.txt'), 'r') as read:
        for line in read:
            chromosome = line.split('\t')[0]
            qEnd = line.split('\t')[7]
            qStart = line.split('\t')[6]
            contigLen = line.split('\t')[11]
            if f'{chromosome}*s' not in telContigs:
                posistionStart = 0
            if f'{chromosome}*e' not in telContigs:
                posistionEnd = 0
            #this determines if a telomere match is at the beginning or ending
            #of a contig by seeing if it is in the first or last third
            if int(qEnd) < (int(contigLen)/3):
                if posistionStart < qEnd:
                    telContigs[f'{chromosome}*s'] = \
                        ('start', qEnd)
                    posistionStart = qEnd
            if int(qEnd) > (int(contigLen) / 3):
                if posistionEnd < qEnd:
                    telContigs[f'{chromosome}*e'] = \
                        ('end', qStart)
                    posistionEnd = qStart
    return telContigs
        
def blastInterrogate():

    readLines, linesCnt = fileListBuilder(os.path.join( \
                        'Outputs/',
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
            telDisS = abs(int(readLines[x].split('\t')[8]) - posistion)
            
        telContigCheck = False
        telDisE = 0
        
        if f'{subjectName}*e' in telContigs:
            blastAn += "\tTelYe"
            telContigCheck = True
        
        if telContigCheck:
            posistion = telContigs[readLines[x].split('\t')[1] + '*e'][1]
            telDisE =  abs(int(readLines[x].split('\t')[8]) - posistion)
        
        telContigCheck = False
        
        if min(telDisS, telDisE) > 0:
            blastAn += f'\t{min(telDisS, telDisE)}\t'
            
        elif min(telDisS, telDisE) == 0:
            blastAn += f'\t{max(telDisS, telDisE)}\t'
        
        revS = False

        revQ = False
        
        if int(readLines[x].split('\t')[6]) \
            > int(readLines[x].split('\t')[7]):
                revS = True
                
        if int(readLines[x].split('\t')[8]) \
            > int(readLines[x].split('\t')[9]):
                revQ = True
                
        if not revS and revQ:
            rev = "revT"
        else:
            rev = "revF"
            
        blastAn += f'{rev}\n'
        
        blastAn = readLines[x].lstrip('\n') + blastAn
        
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
                    f'{directory}blastGenomeOut6Annotated.txt'),\
              'a') as append:
        append.write('##gff-version 3.1.26')
        cnt = '0'
        for line in lines:
            if lines.split('\t')[15] == 'revF\n':
                strand = '+'
            elif lines.split('\t')[15] == 'revT\n':
                strand = '-'
            name = lines.split("\t")[1]
            length = \
                abs(int(lines.split('\t')[9]) - int(lines.split('\t')[8]))
            
            
            gffLine = '\t.'.join(
                            [urllib.parse.quote(name),
                            'teloPortWrapper',
                            'telomere',
                            lines.split('\t')[8],
                            lines.split('\t')[9],
                            lines.split('\t')[10],
                            strand,
                            '.'])
            gffLine += f'\tID=Telomere{cnt};Name={name};Length={str(length)}'
            telDist = lines.split('\t')[14]
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
    
    with open(os.path.join('Outputs/',
              directory,
              'teloPortOut',
              'telReads.fastq'), 'r') as read:
        for line in read:
            if line[0] == '@':
                teloCnt += 1
                
    resultsLine += f'Number of Telomere reads processed: {teloCnt}\n'
    
    deNovoCnt = 0
    
    with open(os.path.join('Outputs/',
              directory,
              'telomereReads',
              'deNovoTelomeres',
              f'{directory}deNovos.fasta'), 'r') as read:
        for line in read:
            if line[0] == '>':
                deNovoCnt += 1
                
    resultsLine += f'Number of de-novo telomeres found by TeloPort: {deNovoCnt}\n'
                
    resultsLine += f'Number of true de-novo telomeres: {trueDeNovos}\n'
    
    resultsLine += f'Number of clusters: {trueClusterCnt}\n'
    
    filteredOut = trueDeNovos - filteredDeNovos
    
    resultsLine += f'Number of true de-novos filtered out: {filteredOut}\n'
    
    resultsLine += \
        f'Number of true de-novo compared to the genome: {filteredDeNovos}\n'
    
    genomesCnt = 0
    
    with open(os.path.join('Outputs/',
                  directory,
                  'blastOut',
                  'blastGenome',
                  f'{directory}blastGenomeOut6.txt'), 'r') as read:
        genomesCnt = len(read.readlines())
    
    resultsLine += \
        f'Number of de novo telomere matches to the genome: {genomesCnt}\n'
        
    resultsLine += 'Distribution of clusters \n'
    
    with open(os.path.join('Outputs',
                           directory,
                           '/clusterHistogram.txt'), 'r') as read:
        histogram = read.readLines()
        
    with open(os.path.join('Outputs',
                           directory,
                           f'{directory}blastGenomeOut6Annotated.txt'), 'r') \
        as read:
            annotation = read.readlines()
            
    with open(os.path.join('Outputs',
                           directory,
                           f'{directory}results.txt'), 'a') as append:
        append.append(resultsLine)
        
        for line in histogram:
            append.append(line)
            
        for line in annotation:
            append.append(line)

#assigns the path varibles from the config file
config()
with open('config.ini','r') as read:
    configLines = read.readlines()
    readDir = configLines[0].split('=')[1].split('\n')[0]
    genomeDir = configLines[1].split('=')[1].split('\n')[0]
    filterDir = configLines[2].split('=')[1].split('\n')[0] 
    addDir = configLines[3].split('=')[1]

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
deNovoRename()
trueClusterCnt = autoMuscle()
consBlast()
trueDeNovos = falsePosFilter()
histogramBuilder()
filteredDeNovos = deNovoFilter()
telContigs = telContigDictMaker()
blast(os.path.join('Outputs/',
                   directory,
                   'telomereReads',
                   'deNovoTelomeres',
                   f'{directory}blastDeNovos.fasta'),
      f'{genomeDir}/{blastGenome}',
      os.path.join('Outputs/',
                   directory,
                   'blastOut',
                   'blastGenome',
                   f'{directory}blastGenomeOut6.txt'),
      '')
if addBlastCnt > 0:
    for file in args.add:
        blast(f'{addDir}file',
              f'{genomeDir}/{blastGenome}',
              os.path.join('Outputs/',
                           directory,
                           'blastOut',
                           'blastAdd',
                           f'{directory}{file.split(".")[0]}blast.txt'),
              '')

blastInterrogate()
gffBuilder()
resultsBuilder()              
                   



#print(args.s)
