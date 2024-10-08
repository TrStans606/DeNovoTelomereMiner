# DeNovoTelomereMiner

## Source

Stansfield, Trey, "TeloPortWrapper: A New Tool for Understanding the Dynamic World of Fungal Telomere Ends" (2023). Mahurin Honors College Capstone Experience/Thesis Projects. Paper 999.
https://digitalcommons.wku.edu/stu_hon_theses/999

DeNovoTelomereMiner Quick Start

Possible new workflow:
- conda
- replace wcdest with mmseq2

For Apple Silicon/Arm:
```conda create --name myenv --platform osx-64```


# **Dependencies**

[TeloPort20XX](https://github.com/TrStans606/TeloPort20xx) a fork of [TeloPort by Seth Baunach](https://github.com/sabaunach/TeloPort) 

[Blast+ 2.11.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/) 

[MUSCLE 3.8](https://drive5.com/muscle/downloads_v3.htm)

[WCD Express 6.3](https://code.google.com/archive/p/wcdest/downloads)

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
