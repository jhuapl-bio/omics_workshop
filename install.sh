#!/bin/bash -i


f="nothing"

while getopts "f" opt; do
  case "$opt" in
    f)  f="overwrite"
      ;;
  esac
done



# prompt the user with y or n (defauult is n) to set f to overwrite
if [[ $f == "nothing" ]]; then
    read -t 2 -p "This script will install the omics_workshop environment and download necessary databases. Do you want to overwrite existing environments and databases? (y/n). Default is (n) OR type "Enter/Return" to use the default " f
    if [[ $f == "y" ]] || [[ $f == "Y" ]]; then
        f="overwrite"
    else
        f="nothing"
    fi
fi


if [[ $f == "overwrite" ]]; then
    echo "Overwriting existing environments and databases for a fresh (re)install..."
fi
unameOut="$(uname -s)"

case "${unameOut}" in
    Linux*)     machine="";OS="Linux";;
    Darwin*)    machine="CONDA_SUBDIR=osx-64"; OS="mac";;
    *)          machine=""; OS="Linux";;
esac
echo Your Machine:${machine}, OS:${OS}  



source $HOME/miniconda3/etc/profile.d/conda.sh
if [[ $f != "nothing" ]]; then
    rm -r ~/miniconda3
fi

if ! command -v conda; then
    echo "Conda not found. Installing Miniconda..."
    

    OS=$(uname -s)
    ARCH=$(uname -m)
    echo "OS: $OS, ARCH: $ARCH"
    if [[ "$OS" == "Linux" && "$ARCH" == "x86_64" ]]; then
        echo "Downloading Miniconda for Linux x86_64..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    elif [[ "$OS" == "Linux" && "$ARCH" == "aarch64" ]]; then
        echo "Downloading Miniconda for Linux aarch64..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O ~/miniconda.sh
    elif [[ "$OS" == "Darwin" && "$ARCH" == "x86_64" ]]; then
        echo "Downloading Miniconda for macOS x86_64..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
    elif [[ "$OS" == "Darwin" && "$ARCH" == "arm64" ]]; then
        echo "Downloading Miniconda for macOS arm64..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -O ~/miniconda.sh
    else
        echo "Unsupported OS or architecture: $OS $ARCH"
        exit 1
    fi
    # # Install Miniconda
    rm -r ~/miniconda3
    bash ~/miniconda.sh -f -b -u -p ~/miniconda3

    echo "Initializing Miniconda... for OS: $OS"
    source $HOME/miniconda3/etc/profile.d/conda.sh

    if [[ $OS == "Darwin" ]]; then
        echo "Adding darwin"
        # Check if the line already exists in the file
        grep -qxF 'bash $HOME/miniconda3/etc/profile.d/conda.sh' ~/.zshrc || echo 'bash $HOME/miniconda3/etc/profile.d/conda.sh' >> ~/.zshrc
        grep -qxF 'export PATH=$HOME/miniconda3/bin:$PATH' ~/.zshrc || echo 'export PATH=$HOME/miniconda3/bin:$PATH' >> ~/.zshrc
    elif [[ $OS == "Linux" ]]; then
        echo "Adding Linux"
        # Check if the line already exists in the file
        grep -qxF 'bash $HOME/miniconda3/etc/profile.d/conda.sh' ~/.bashrc || echo 'bash $HOME/miniconda3/etc/profile.d/conda.sh' >> ~/.bashrc
        grep -qxF 'export PATH=$HOME/miniconda3/bin:$PATH' ~/.bashrc || echo 'export PATH=$HOME/miniconda3/bin:$PATH' >> ~/.bashrc
    fi
    conda_base=$HOME/miniconda3
else
    echo "Conda is already installed."
    

fi
conda_base=$(conda info --base) 
source $conda_base/etc/profile.d/conda.sh

if ! command -v conda; then
    echo "Conda not found. Please install Miniconda and rerun this script."
    exit 1
fi
# # Wait for the user to close and reopen their shell or source their bashrc
# echo "Please close and reopen your terminal, or run 'source ~/.bashrc', then rerun this script."

# Check if the omics_workshop environment exists

if  conda env list | grep -q 'omics_workshop' && [[ $f == "nothing" ]]; then
    echo "Env found, updating omics_workshop environment..."
    # 
    eval "conda activate omics_workshop && $machine conda install -y -c bioconda bowtie2 minimap2 kraken2 krona fastqc samtools bcftools git fastp python"
else
    echo "Creating new env"
    conda env remove -y -n omics_workshop
    eval "$machine conda create -y -c bioconda -n omics_workshop bowtie2 minimap2 kraken2 krona fastqc samtools bcftools git fastp python"
fi

eval "$machine conda install  -y git"

conda activate omics_workshop


# check if git is a command and if not then exit
if ! command -v git; then
    echo "Git not found. Please install Git and rerun this script."
    exit 1
fi

if [[ $machine -eq "CONDA_SUBDIR=osx-64" ]]; then    
    conda activate omics_workshop && \
        python -c "import platform;print(platform.machine())" && \
        conda config --env --set subdir osx-64
fi



# Check if ~/omics_workshop exists and if it does then cd and git pull force changes otherwise clone new one
REPO_URL="https://github.com/jhuapl-bio/omics_workshop"
TARGET_DIR="$HOME/omics_workshop"

# Check if the directory exists and is a Git repository
if [ -d "$TARGET_DIR/.git" ] || [[ $f != "nothing" ]]; then
    # Try to pull the repository
    echo "Attempting to pull the latest changes for $TARGET_DIR..."
    git -C "$TARGET_DIR" pull || {
        echo "Failed to pull changes. Removing directory and attempting to re-clone..."
        rm -rf "$TARGET_DIR"
        git clone "$REPO_URL" "$TARGET_DIR" || {
            echo "Failed to clone repository."
            exit 1
        }
    }
else
    # If not a valid Git repository, remove directory and clone fresh
    echo "Directory does not exist as a Git repository. Cloning afresh..."
    rm -rf "$TARGET_DIR"  # Make sure to remove it if it exists but is not a repo
    git clone "$REPO_URL" "$TARGET_DIR" || {
        echo "Failed to clone repository."
        exit 1
    }
fi

# 

# if [ -d ~/omics_workshop ]; then
#     echo "Updating 'omics_workshop' GitHub repository..."
#     cd ~/omics_workshop && git pull
# else
#     echo "Cloning 'omics_workshop' GitHub repository..."
#     git clone https://github.com/jhuapl-bio/omics_workshop ~/omics_workshop
# fi

## make sure we update Taxonomy for krona to work! Requires internet
echo "Downloading taxonomy information. Requires internet connection. This will take some time...."
# bash ktUpdateTaxonomy.sh
if [ -s $(dirname $(which ktImportTaxonomy))/../opt/krona/taxonomy/taxonomy.tab ] && [[ $f == "nothing" ]]; then
    echo "Krona taxonomy.tab exists, skipping...."
else
    mkdir -p  ~/omics_workshop/downloads/ && conda activate omics_workshop \
        && wget https://github.com/jhuapl-bio/datasets/raw/main/databases/ncbi/taxonomy.tab.gz \
        -O ~/omics_workshop/downloads/taxonomy.tab.gz \
        && gzip -f -d ~/omics_workshop/downloads/taxonomy.tab.gz  && \
        mv ~/omics_workshop/downloads/taxonomy.tab $(dirname $(which ktImportTaxonomy))/../opt/krona/taxonomy/
fi


tar  -xvzf ~/omics_workshop/databases/test_metagenome.tar.gz -C ~/omics_workshop/databases/


gzip -f -d --keep ~/omics_workshop/references/test.fasta.gz

# check if ~/omics_workshop/databases/k2_viral/hash.k2d and is size >0 then skip else download
if [ -s ~/omics_workshop/databases/k2_viral/hash.k2d ] && [[ $f == "nothing" ]]; then
    echo "Kraken2 viral database already exists."
else
    echo "Downloading Kraken2 viral database..."
    wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz \
        -O ~/omics_workshop/downloads/k2_viral.tar.gz && \
        mkdir -p ~/omics_workshop/downloads/k2_viral && \
        tar -xvzf ~/omics_workshop/downloads/k2_viral.tar.gz -C ~/omics_workshop/downloads/k2_viral  && \
        mv ~/omics_workshop/downloads/k2_viral ~/omics_workshop/databases/k2_viral
fi

echo "Setup complete."
