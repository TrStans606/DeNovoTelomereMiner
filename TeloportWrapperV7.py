import subprocess
import os
import argparse
import tempfile
import glob


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

#command line argument for determining if a cluster is a canidate de novo 
parser.add_argument(
	'--cut', 
    action='store', 
    type=int, 
    help= ('The number of reads in a cluster before it is'
           'labeled a canidate de novo.'
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

args = parser.parse_args()

#config file maker

def config():
    with open('Programs/TeloPort/config.ini', 'r') as read:
        if read.readline().split('\n')[0] == 'setup=yes':
            read.close()
            print('First Time setup')
            readDir = input('Raw Reads Directory: ')
            genomeDir = input('Assembled genome Directory: ')
            filterDir = input('Filters directory: ')
            addDir = input('Additional Blast directory: ')
            write = open('Programs/TeloPort/config.ini', 'w')
            write.write(
                f'setup=no\nreadDirectory={readDir}\ngenomeDirectory= \
                {genomeDir}\nfilterDirectory={filterDir}\naddDirectory= \
                    {addDir}')
                #'setup=no\nreadDirectory='  \
                 #   + readDir \
                  #      + '\ngenomeDirectory=' \
                   #            + '\nfilterDirectory=' \
                    #                + filterDir \
                     #                   + '\naddDirectory=' \
                      #                      + addDir)

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
        os.makedirs(os.path.join('Programs/TeloPort/Outputs', 
            directory, 
            path))
    if choice == 's':
        telomereFinderSeperated()
    elif choice =='i':
        telomereFinderInter()
    idTagger()
    os.remove(
        f'Programs/TeloPort/OutPuts/{directory}/teloPortOut/pairReads.fastq'
        )
    junctionFinder()
    sequenceQuality()
    wcdestCall()
    wcdInterrogate()
#runs telomereFinder with the -s argument
#this is used for seperated fastq files
def telomereFinderSeperated():
    subprocess.run(['Programs/TeloPort/build/apps/telomereFinder', 
        '-s',
        os.path.join(readDir, genomeR1),
        os.path.join(readDir, genomeR1),
        '-o',
        'Programs/TeloPort/Outputs/' + directory + '/teloPortOut/',
        '-t',
        telRepeat],
        shell=True,
        check=True)
#runs telomereFinder with the -i argument
#this is used for interleaved fastq files
def telomereFinderInter():
    subprocess.run(['Programs/TeloPort/build/apps/telomereFinder',
        '-i',
        os.path.join(readDir, genomeI),
        '-o',
        f'Programs/TeloPort/Outputs/{directory}/teloPortOut/',
        '-t',
        telRepeat],
        shell=True,
        check=True)
#tags every read with a unique id number for tracking
def idTagger():
    pairLines, pairLineCnt = fileListBuilder(
        os.path.join('Programs/TeloPort/Outputs/',
            directory,
            '/teloPortOut/pairReads.fastq')
        )
    newFile = os.path.join('Programs/TeloPort/Outputs/',
        directory,
        '/teloPortOut/subTelReads.fastq')
    with open(newFile, 'a') as write:
        x = 0
        for line in pairLines:
            if x % 4 == 0:
                newLine = f"@{directory}{x}\t{line.lstrip('@')}"
                x += 1
                write.write(newLine)
#calls junctionfinder
def junctionFinder():
    subprocess.run(['Programs/TeloPort/build/apps/junctionFinder',
                    '-i',
                    os.path.join('Programs/TeloPort/OutPuts',
                                 directory,
                                 'teloPortOut/subTelReads.fastq'),
                    '-o',
                    os.path.join('Programs/TeloPort/OutPuts',
                                 directory,
                                 'teloPortOut/'),
                    '--revc 0',
                    '--splitJunc 1',
                    '--splitDir 0'],
                   shell=True,
                   check=True)
#calls sequenceQuality
def sequenceQuality():
    subprocess.run(['Programs/TeloPort/build/apps/sequenceQuality',
                    '-i',
                    os.path.join('Programs/TeloPort/OutPuts',
                                 directory,
                                 'teloPortOut/telAdjSeq.fastq'),
                    '-o',
                    os.path.join('Programs/TeloPort/OutPuts',
                                 directory,
                                 'teloPortOut/hiQualityTelAdjSeq.fastq'),
                    '-ofmt',
                    'fasta',
                    '-c',
                    '0',
                    '-l',
                    '0'],
                   shell=True,
                   check=True)
#calls WCDest
def wcdestCall():
    subprocess.run(['wcd',
                    os.path.join('Programs/TeloPort/OutPuts',
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
                        'Programs/TeloPort/OutPuts',
                        directory,
                        'teloPortOut',
                        f'{directory}clusters.wcd')],
                    shell=True,
                    check=True)
#calls wcd interrogate
def wcdInterrogate():
    subprocess.run(['Programs/TeloPort/build/apps/wcdInterrogate',
                    '-w',
                    os.path.join('Programs/TeloPort/OutPuts',
                                 directory,
                                 'teloPortOut',
                                 f'{directory}clusters.wcd'),
                    '-i',
                    os.path.join('Programs/TeloPort/OutPuts',
                                 directory,
                                 'teloPortOut',
                                 'hiQualityTelAdjSeq.fastq'),
                    '-s',
                    os.path.join('Programs/TeloPort/OutPuts',
                                 directory,
                                 'teloPortOut',
                                 'telAdjSeq'),
                    '-f',
                    'fastq',
                    '-o',
                    os.path.join('Programs/TeloPort/OutPuts',
                                 directory,
                                 'teloPortOut',
                                 f'{directory}Clusters.out'),
                    '-r',
                    os.path.join('Programs/TeloPort/OutPuts',
                                 directory,
                                'telomereClusters',
                                'telomereReads/'),
                    '--indices',
                    '--sort',
                    '--size'
                    '1'],
                   shell=True,
                   check=True)
#renames all canidate denovos from clusters
def deNovoRename():
    cntCluster = glob.glob(os.path.join('Programs/TeloPort/Outputs/',
                                  directory,
                                  'telomereReads',
                                  'telomereClusters',
                                  'cluster*.fasta'))
    firstDeNovo = firstDeNovoCount()
    cluster = os.path.join('Programs/TeloPort/Outputs/',
                                  directory,
                                  'telomereReads',
                                  'telomereClusters',
                                  'cluster')
    for i in range(firstDeNovo,cntCluster):
        with open(f'{cluster}{i}.fasta', 'r') as read:
            deNovoLines = read.readlines()
        with open(os.path.join('Programs/TeloPort/Outputs/',
                  directory,
                  'telomereReads',
                  'deNovoTelomeres',
                  'deNovos.fasta'), 'a') as write:
            for line in deNovoLines:
                if line[0]=='>':
                    write.write(f'{directory}\t{line.lstrip(">")}')
                else:
                    write.write(line)
            os.remove(f'{cluster}{i}.fasta')
        newClusterCnt = glob.glob(os.path.join('Programs/TeloPort/Outputs/',
                                      directory,
                                      'telomereReads',
                                      'telomereClusters',
                                      'cluster*.fasta'))
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
    with open(os.path.join('Programs/TeloPort/Outputs',
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
 #runs muscle   
def autoMuscle():
    clusters = os.path.join('Programs/TeloPort/Outputs/',
                                        directory,
                                        'telomereReads',
                                        'telomereClusters')
    cntCluster = glob.glob(f'{clusters}cluster*.fasta')
    for i in range(0,cntCluster):
        subprocess.run(['Programs/muscle5.1.linux_intel64',
                        '-allign',
                        f'{clusters}cluster{i}.fastq',
                        '-output',
                        os.path.join('Programs/TeloPort/Outputs/',
                                     directory,
                                     'muscleOut',
                                     f'{directory}cluster{i}.msa')],
                       shell=True,
                       check=True)
                        
              
              
                           

#assigns the path varibles from the config file
config()
with open('Programs/TeloPort/config.ini','r') as read:
    configLines = read.readlines()
    readDir = configLines[1].split('=')[1].split('\n')[0]
    genomeDir = configLines[2].split('=')[1].split('\n')[0]
    filterDir = configLines[3].split('=')[1].split('\n')[0] 
    addDir = configLines[4].split('=')[1]

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
autoMuscle()


#print(args.s)
