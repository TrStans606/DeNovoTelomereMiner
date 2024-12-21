# Manual Installation Guide
## 1. Setting up the Repository
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
