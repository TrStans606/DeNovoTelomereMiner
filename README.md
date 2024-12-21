# DeNovoTelomereMiner

## Table of contents
1. [Original Paper](https://github.com/TrStans606/DeNovoTelomereMiner/tree/main#source)
2. [DeNovoTelomereMiner Quick Start]

## Source

Stansfield, Trey, "TeloPortWrapper: A New Tool for Understanding the Dynamic World of Fungal Telomere Ends" (2023). Mahurin Honors College Capstone Experience/Thesis Projects. Paper 999.
https://digitalcommons.wku.edu/stu_hon_theses/999

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
### 1. Setting up the Repository
#### Linux & macOS
Make sure you have [git](https://git-scm.com/) installed:

```shell
git clone https://github.com/TrStans606/DeNovoTelomereMiner
```

### 2. Setting Up Dependencies

If you are unsure if you have an Intel or Arm-based Mac, then click on the Apple menu () in the top-left corner of your screen. If the processor information says anything with Intel in the name, then it is Intel-based. If it says Apple M, then it is ARM-based.

#### Linux & Intel/x86/OSX-64 Macs 

I would recommend creating a conda environment for this as it cleanly installs everything in one command, which can be uninstalled easily later. 

```shell
conda create -n de_novo python=13.1 bioconda::blast bioconda::emboss bioconda::muscle bioconda::mmseqs2 conda-forge::pandas
```

```shell
conda activate de_novo
```

#### Apple Silicon/ARM/OSX-arm64 Macs

Many of the dependencies here do not have an ARM build in Bioconda despite binaries existing. Because of this, you have two options: either sticking with Conda via Apple's x86 emulation or using Homebrew. Just know that using the Conda method on ARM Macs will be slower due to the emulation overhead.

##### Installing Via Homebrew

1. Enable the BrewSci tap 

```shell
brew tap brewsci/bio
```

2. Install all of the needed dependencies

```shell
brew install blast brewsci/bio/emboss brewsci/bio/muscle mmseqs2
```

3. Install Pandas via a Conda or Pip environment

```shell
python3 -m venv .venv
```

```shell
source .venv/bin/activate
```

```shell
pip install pandas
```

or 

```shell
conda create --name de_novo pandas
```

##### Installing Via Conda

This method will be easier but will run slower due to the entire environment running in x86 emulation mode. I also cannot guarantee that future updates will work with this method. Use at your own risk.
1. Install Apple's Rosetta 2 if you haven't already.

```shell
softwareupdate --install-rosetta
```

2. Set up the Conda environment.

```shell
conda create -n de_novo --platform osx-64 python=13.1 bioconda::blast bioconda::emboss bioconda::muscle bioconda::mmseqs2 conda-forge::pandas
```

```shell
conda activate de_novo
```
### 3. Compiling TeloPort
#### Linux & macOS

TeloPort20XX needs to be compiled from source. To aid in this, there is the install.sh script. This does require a C++ compiler and has been tested on GCC and Apple clang. The boost library is needed for compilation and can be installed from your native package manager, Conda, or Homebrew 

##### Installing Boost via Homebrew

```shell
brew install boost
```

##### Installing Boost via Conda
Make sure this is done in the same conda environment used in step 2.

```shell
conda install conda-forge::boost
```

###### Installing Boost via Apt (Debian/Ubuntu/Linux Mint)
```shell
sudo apt install libboost-all-dev
```

###### Installing Boost via DNF (Red Hat/Fedora/CentOS)
```shell
sudo dnf install boost-devel
```
###### Installing Boost via Pacman (Arch)
```shell
sudo pacman -S boost
```

##### Compiling TeloPort

```shell
bash install.sh
```
If you are getting errors that boost isn't found even though it has been installed via Homebrew then you may need to add it your shell profile:
For zsh users (this is the default on macOS) edit ~/.zshrc:
```shell
nano ~/.zshrc
```
Then add the lines to the file:
```shell
export CPATH=/opt/homebrew/include
export LIBRARY_PATH=/opt/homebrew/lib
```
For Bash users do the same but for ~/.bashrc. Then close out of and reopen your terminal for the changes to take effect.


# **Dependencies**

[TeloPort20XX](https://github.com/TrStans606/TeloPort20xx) a fork of [TeloPort by Seth Baunach](https://github.com/sabaunach/TeloPort) 

[Blast+ 2.11.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/) 

[MUSCLE 3.8](https://drive5.com/muscle/downloads_v3.htm)

[Mmseqs2](https://github.com/soedinglab/MMseqs2/)

[Emboss 6.6.0](http://emboss.sourceforge.net/download/)


# Setting up TeloPort

1. Create directory called Programs in root folder
2. Unzip TeloPort into the directory Programs
3. While in the TeloPort directory type make

#Dependencies for DeNovoTelomereMiner

#Inputs for TeloPortWrapper
an interleved or two seperated fastq files
a fasta genome file
 
#How to use DeNovoTelomereMiner
