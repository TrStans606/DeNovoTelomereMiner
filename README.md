# DeNovoTelomereMiner

## Table of contents
1. [Original Paper](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#source)
2. [DeNovoTelomereMiner Quick Start]

## Source

Stansfield, Trey, "TeloPortWrapper: A New Tool for Understanding the Dynamic World of Fungal Telomere Ends" (2023). Mahurin Honors College Capstone Experience/Thesis Projects. Paper 999.
https://digitalcommons.wku.edu/stu_hon_theses/999

# **Dependencies**

[TeloPort20XX](https://github.com/TrStans606/TeloPort20xx) a fork of [TeloPort by Seth Baunach](https://github.com/sabaunach/TeloPort) 

[Blast+ 2.11.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/) 

[MUSCLE 3.8](https://drive5.com/muscle/downloads_v3.htm)

[Mmseqs2](https://github.com/soedinglab/MMseqs2/)

[Emboss 6.6.0](http://emboss.sourceforge.net/download/)

# DeNovoTelomereMiner Quick Start

DeNovoTelomereMiner is built for POSIX systems and has been tested on Fedora, Ubuntu, and macOS 15. For an easy install try our new docker option which works on Linux, macOS (Intel and ARM), and Windows. For Windows users, try [Windows Subsystem Linux](https://learn.microsoft.com/en-us/windows/wsl/install) offered by Microsoft. It has been tested on Python version 13.1. It is recommended to install dependencies via [Homebrew](https://brew.sh/) or [Conda](https://docs.anaconda.com/miniconda/install/), both of which are package managers on macOS and Linux. Instructions for installation are included on the project’s homepages.

## Docker instructions for Linux, macOS and Windows
1. Run this build command in terminal (macOS and Linux) or PowerShell (Windows):

```shell
docker build -t denovo_telomere_miner .
```

2. Activate the Docker image with the command:

```shell
docker run -it denovo_telomere_miner bash
```
3. Run the shell script inside the new container enviroment:

```shell
bash install.sh
```

It's that simple! Now the teloport environment is set up perfectly for running. To copy your data from your host machine to the docker image use this command on your host machine:
```shell
docker cp <any file or directory> denovo_telomere_miner:DeNovoTelomereMiner
```

## Manual Installation

For manual installation instructions check out the [Installation.md](Installation.md) guide

# Setting up TeloPort

1. Create directory called Programs in root folder
2. Unzip TeloPort into the directory Programs
3. While in the TeloPort directory type make

#Dependencies for DeNovoTelomereMiner

#Inputs for TeloPortWrapper
an interleved or two seperated fastq files
a fasta genome file
 
#How to use DeNovoTelomereMiner
