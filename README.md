# DeNovoTelomereMiner

## Table of contents
1. [Original Paper](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#source)
2. [Dependencies](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#dependencies)
3. [Docker Installation Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#docker-instructions-for-linux-macos-and-windows)
4. [Manual Installation Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#manual-installation)
5. [DeNovoTelomereMiner.py Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#denovotelomereminer-usage)
6. [config.ini Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#configini-guide)

## Source

Stansfield, Trey, "TeloPortWrapper: A New Tool for Understanding the Dynamic World of Fungal Telomere Ends" (2023). Mahurin Honors College Capstone Experience/Thesis Projects. Paper 999.
https://digitalcommons.wku.edu/stu_hon_theses/999

# DeNovoTelomereMiner Quick Start

DeNovoTelomereMiner is built for POSIX systems and has been tested on Fedora, Ubuntu, and macOS 15. For an easy install try the new [docker](https://www.docker.com/get-started/) option which works on Linux, macOS (Intel and ARM), and Windows. For Windows users, try [Windows Subsystem Linux](https://learn.microsoft.com/en-us/windows/wsl/install) offered by Microsoft. It has been tested on Python version 13.1. It is recommended to install dependencies via [Homebrew](https://brew.sh/) or [Conda](https://docs.anaconda.com/miniconda/install/), both of which are package managers on macOS and Linux. Instructions for installation are included on the project’s homepages.

# **Dependencies**

Python 3.8 or greater

[TeloPort20XX](https://github.com/TrStans606/TeloPort20xx) a fork of [TeloPort by Seth Baunach](https://github.com/sabaunach/TeloPort) 

[Blast+ 2.11.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/) 

[MUSCLE 3.8](https://drive5.com/muscle/downloads_v3.htm)

[Mmseqs2](https://github.com/soedinglab/MMseqs2/)

[Emboss 6.6.0](http://emboss.sourceforge.net/download/)

## Docker instructions for Linux, macOS and Windows
1. Pull the image from Github
```bash
docker pull ghcr.io/trstans606/denovotelomereminer:latest
```

2. Run the container

```bash
docker run --name de_novos -it 82b6ad51a816 bash
```

It's that simple! Now the teloport environment is set up perfectly for running. To copy your data from your host machine to the docker image use this command on your host machine:

```shell
docker cp <any file or directory> de_novos:docker/
```

You can also reverse this transfer:
```shell
docker cp de_novos:docker/<insert file path> <any file or directory> 
```

The image also comes equiped with wget and the NCBI sra-toolkit to easily add data from FTP servers or the Sequence Read Archive (SRA) directly from the container itself.

## Building Docker from Source

If for whatever reason the ghcr.io service goes down you can also use Docker to manually rebuild the image. For instructions check the [wiki](https://github.com/TrStans606/DeNovoTelomereMiner/wiki/Building-Docker-from-Source) page.

## Manual Installation

For manual installation instructions check out the [Installation wiki](https://github.com/TrStans606/DeNovoTelomereMiner/wiki/Manual-Installation-Guide) guide

## DeNovoTelomereMiner Usage

### DeNovoTelomereMiner Command Line Arguments

For this command line, only include the file name. The directory of each file type will be set in the config file.

```bash
python3 DeNovoTelomereMiner.py [-h] [-s S S] [-i I] -g G [-d D] [-t T] [-f [F ...]] [--cut CUT] [--simple] [--add [ADD ...]] [--config]
```

```-h```: show the help message and exits

```-s S S```: The two separated raw read FASTQ files. This argument expects two file names separated by a space.

```-i I```: The interleaved raw read FASTQ file. This option is mutually exclusive with -s.

```-g G```: The name of the assembled genome FASTA file used for alignment. If simple mode is not used, options -g are required.

```-d D```: The name of the directory where all output files will be placed. The name defaults to the name of the fasta file provided. Directory names must be unique and must not already exist.

```-t T```: The foreward tel repeat you want to search for. Default is CCCTAA.

```-f [F ...]```: List all of the seqeunce files you want removed from reads. This actively removes some sequences from anaylsis. Must be in fasta format. You can list any number of files.

```--add [ADD ...]```: List all addition seqeunces files you want compared to the reads. This adds in additional files which the raw reads will compared to via BLAST. All must be in fasta format.You can list any number of files.

```--cut CUT```: The number of reads in a cluster before it is labeled a canidate de-novo. The default is five.

```--simple```: Activates Simple mode, a step by step interactive input mode.

```--config```: Allows the user to interactively alter the config file.

### DeNovoTelomereMiner Sample Command

This command is assuming your raw read files are in Files/ and your genome is in testData/.

```bash
python3 DeNovoTelomereMiner.py -s ERR2061621_1.fastq ERR2061621_2.fastq -g B71v2sh.fasta -d ERR2061621
```

Your results will then be found in the ERR2061621 folder in Outputs.

## Config.ini Guide

DeNovoTelomereMiner finds where your files are via the config.ini file. The file can either be directly edited in the text editor of your choice or it can be edited within DeNovoTelomereMiner itself via the --config flag.

```
readDirectory=Files/
genomeDirectory=testData/
filterDirectory=Files/
addDirectory=Files/
```

```readDirectory```= is where the raw fastq reads are (the -s or -i flags).

```genomeDirectory```= is where the assembled fasta genome is (the -g flag).

```filterDirectory```= is where the fasta files for filtering sequences are (the -f flag).

```addDirectory```= is where the additional fasta files used for comparison are (the --add flag).

These file paths can be absolute or relavtive but the path must exist and must not have any spaces before the equal sign and the file path, so:

```readDirectory=/home/usr/Downloads/Genomes``` is okay but:
```readDirectory= /home/usr/Downloads/Genomes``` isn't as there is a space between the equals sign and path.


