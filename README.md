# DeNovoTelomereMiner

## Source

Stansfield, Trey, "TeloPortWrapper: A New Tool for Understanding the Dynamic World of Fungal Telomere Ends" (2023). Mahurin Honors College Capstone Experience/Thesis Projects. Paper 999.
https://digitalcommons.wku.edu/stu_hon_theses/999

# DeNovoTelomereMiner Quick Start

DeNovoTelomereMiner is built for POSIX systems and has been tested on Fedora, Ubuntu, and macOS 15. For Windows users, try [Windows Subsystem Linux](https://learn.microsoft.com/en-us/windows/wsl/install) offered by Microsoft. It has been tested on Python version 13.1. It is recommended to install dependencies via [Homebrew](https://brew.sh/) or [Conda](https://docs.anaconda.com/miniconda/install/), both of which are package managers on macOS and Linux. Instructions for installation are included on the project’s homepages.

## 1. Setting up the Repository
### Linux & macOS
Make sure you have [git](https://git-scm.com/) installed:

```git clone https://github.com/TrStans606/DeNovoTelomereMiner```

## 2. Setting Up Dependencies

If you are unsure if you have an Intel or Arm-based Mac, then click on the Apple menu () in the top-left corner of your screen. If the processor information says anything with Intel in the name, then it is Intel-based. If it says Apple M, then it is ARM-based.

### Linux & Intel/x86/OSX-64 Macs 

I would recommend creating a conda environment for this as it cleanly installs everything in one command, which can be uninstalled easily later. 

```conda create -n de_novo python=13.1 bioconda::blast bioconda::emboss bioconda::muscle bioconda::mmseqs2 conda-forge::pandas```

```conda activate de_novo```

### Apple Silicon/ARM/OSX-arm64 Macs

Many of the dependencies here do not have an ARM build in Bioconda despite binaries existing. Because of this, you have two options: either sticking with Conda via Apple's x86 emulation or using Homebrew. Just know that using the Conda method on ARM Macs will be slower due to the emulation overhead.

#### Installing Via Homebrew

1. Enable the BrewSci tap 

```brew tap brewsci/bio```

2. Install all of the needed dependencies

```brew install blast brewsci/bio/emboss brewsci/bio/muscle mmseqs2```

3. Install Pandas via a Conda or Pip environment

```python3 -m venv .venv```

``` source .venv/bin/activate```

```pip install pandas```

or 

```conda create --name de_novo pandas```

#### Installing Via Conda

This method will be easier but will run slower due to the entire environment running in x86 emulation mode. I also cannot guarantee that future updates will work with this method. Use at your own risk.
1. Install Apple's Rosetta 2 if you haven't already.

```softwareupdate --install-rosetta```

2. Set up the Conda environment.

```conda create -n de_novo --platform osx-64 python=13.1 bioconda::blast bioconda::emboss bioconda::muscle bioconda::mmseqs2 conda-forge::pandas```

```conda activate de_novo```

### Linux & macOS

TeloPort20XX needs to be compiled from source. To aid in this, there is the install.sh script. This does require a C++ compiler and has been tested on GCC and Apple clang.

#### Installing Boost via Homebrew

```brew install boost```

#### Installing Boost via Conda
Make sure this is done in the same conda environment used in step 2.

```conda install conda-forge::boost```

#### Compiling TeloPort

```bash install.sh```

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
