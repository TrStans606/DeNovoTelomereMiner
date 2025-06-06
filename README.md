# DeNovoTelomereMiner

TODO: SEED FINDING
Scan the ref genome, BV15, not the raw reads
Add more in padding before something is consiered a seed
Stop seeds from being reversed


## Table of contents
1. [Original Paper](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#source)
2. [Docker Installation Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#docker-instructions-for-linux-macos-and-windows)
3. [Manual Installation Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#manual-installation)
4. [DeNovoTelomereMiner.py Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#denovotelomereminer-usage)
5. [config.ini Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#configini-guide)
6. Outputs

## Source

Stansfield, Trey, "TeloPortWrapper: A New Tool for Understanding the Dynamic World of Fungal Telomere Ends" (2023). Mahurin Honors College Capstone Experience/Thesis Projects. Paper 999.
https://digitalcommons.wku.edu/stu_hon_theses/999

# DeNovoTelomereMiner Quick Start

DeNovoTelomereMiner is built for POSIX systems and has been tested on Fedora, Ubuntu, and macOS 15. For an easy install try the new [docker](https://www.docker.com/get-started/) option which works on Linux, macOS (Intel and ARM), and Windows. This image should also work with Podman. For Windows users, try [Windows Subsystem Linux](https://learn.microsoft.com/en-us/windows/wsl/install) offered by Microsoft. It has been tested on Python version 13.1. It is recommended to install dependencies via [Homebrew](https://brew.sh/) or [Conda](https://docs.anaconda.com/miniconda/install/), both of which are package managers on macOS and Linux. Instructions for installation are included on the project’s homepages.

## Docker instructions for Linux, macOS and Windows
1. Pull the image from Github

Windows, Intel Macs, x86 Linux
```bash
docker pull ghcr.io/trstans606/denovotelomereminer:x86@sha256:3ba7f0478c55956eb81709658243fa13ab8bfbf4af7b94eeda20d9d7890fcb5e
```

Apple Silicon Mac, arm64 Linux
```bash
docker pull ghcr.io/trstans606/denovotelomereminer:arm64@sha256:4f6b988680b3db68b64f2fa5d73dd76405d5e5ebf21093c878d1d5bd0ef785e6
```

2. Run the container

Run docker images to find the image id:
```bash
docker images
```
You will get some result like:
REPOSITORY                               TAG       IMAGE ID       CREATED          SIZE
ghcr.io/trstans606/denovotelomereminer   <none>    bcc98deeb377   43 minutes ago   36.8GB

Copy the value of IMAGE ID into the below command

```bash
docker run --name de_novos -it <IMAGE ID> bash
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

### Using Podman

This image should also work using Podman on Linux but testing is first and foremost done with Docker.

1. Pull the image from Github
```bash
podman pull ghcr.io/trstans606/denovotelomereminer:latest
```
2. Run the container

```bash
podman run --name de_novos -it 82b6ad51a816 bash
```

Everything else should be the same with podman just make sure you use the appropriate Podman command for file copying.

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

## Outputs

The main output file is the {directory}results.txt file in your specified output directory in Outputs. This 
