# DeNovoTelomereMiner

## Table of contents
1. [Original Paper](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#source)
2. [Dependencies](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#dependencies)
3. [Docker Installation Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#docker-instructions-for-linux-macos-and-windows)
4. [Manual Installation Guide](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#manual-installation)
5. [DeNovoTelomereMiner.py Guide]()
6. [config.ini Guide]()

## Source

Stansfield, Trey, "TeloPortWrapper: A New Tool for Understanding the Dynamic World of Fungal Telomere Ends" (2023). Mahurin Honors College Capstone Experience/Thesis Projects. Paper 999.
https://digitalcommons.wku.edu/stu_hon_theses/999

# DeNovoTelomereMiner Quick Start

DeNovoTelomereMiner is built for POSIX systems and has been tested on Fedora, Ubuntu, and macOS 15. For an easy install try the new [docker](https://www.docker.com/get-started/) option which works on Linux, macOS (Intel and ARM), and Windows. For Windows users, try [Windows Subsystem Linux](https://learn.microsoft.com/en-us/windows/wsl/install) offered by Microsoft. It has been tested on Python version 13.1. It is recommended to install dependencies via [Homebrew](https://brew.sh/) or [Conda](https://docs.anaconda.com/miniconda/install/), both of which are package managers on macOS and Linux. Instructions for installation are included on the project’s homepages.

# **Dependencies**

[TeloPort20XX](https://github.com/TrStans606/TeloPort20xx) a fork of [TeloPort by Seth Baunach](https://github.com/sabaunach/TeloPort) 

[Blast+ 2.11.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/) 

[MUSCLE 3.8](https://drive5.com/muscle/downloads_v3.htm)

[Mmseqs2](https://github.com/soedinglab/MMseqs2/)

[Emboss 6.6.0](http://emboss.sourceforge.net/download/)

## Docker instructions for Linux, macOS and Windows
1. Pull the image from Github
```bash
docker pull ghcr.io/trstans606/denovotelomereminer:latest@sha256:82b6ad51a816c2d48e46b7289eaff1ab9e62b72ca5a145b021b2a5258f8034c2
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

If for whatever reason the ghcr.io service goes down you can also use Docker to manually rebuild the image: 

1. Run this build command in terminal (macOS and Linux) or PowerShell (Windows):

```shell
docker build -t denovo_telomere_miner .
```

2. Activate the Docker image with the command:

```shell
docker run --name de_novo -it denovo_telomere_miner bash
```
3. Run the shell script inside the new container enviroment:

```shell
bash install.sh
```

It's that simple! Now the teloport environment is set up perfectly for running. To copy your data from your host machine to the docker image use this command on your host machine:
```shell
docker cp <any file or directory> de_novo:docker/
```
The image also comes equiped with wget and the NCBI sra-toolkit to easily add data from FTP servers or the Sequence Read Archive (SRA) directly from the container itself.

## Manual Installation

For manual installation instructions check out the [Installation.md](Installation.md) guide
