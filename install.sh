#!/bin/bash -i


while getopts "f" opt; do
  case "$opt" in
    f)  f="overwrite"
      ;;
  esac
done

if [[ $f == "overwrite" ]]; then
    echo "Overwriting existing environments and databases for a fresh (re)install..."
fi
# Check for Conda installation
source ~/.bashrc 
if ! command -v conda; then
    echo "Conda not found. Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    bash ~/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda.sh
    echo "Initializing Miniconda..."
    ~/miniconda3/bin/conda init bash
else
    echo "Conda is already installed."
fi

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine="";;
    Darwin*)    machine="CONDA_SUBDIR=osx-64";;
    *)          machine=""
esac
echo ${machine}

conda_base=$(conda info --base) 
echo $conda_base
source $conda_base/etc/profile.d/conda.sh

# # Wait for the user to close and reopen their shell or source their bashrc
# echo "Please close and reopen your terminal, or run 'source ~/.bashrc', then rerun this script."

# Check if the omics_workshop environment exists

if  conda env list | grep -q 'omics_workshop' && [[ $f == "" ]]; then
    echo "Env found, updating omics_workshop environment..."
    # 
    eval "conda activate omics_workshop && $machine conda install -y -c bioconda bowtie2 minimap2 kraken2 krona fastqc samtools bcftools git fastp python"
else
    if [[ $f != "" ]]; then 
        conda env remove -y -n omics_workshop
        eval "$machine conda create -y -c bioconda -n omics_workshop bowtie2 minimap2 kraken2 krona fastqc samtools bcftools git fastp python"
    fi 
fi



conda activate omics_workshop


if [[ $machine -eq "CONDA_SUBDIR=osx-64" ]]; then    
    conda activate omics_workshop && \
        python -c "import platform;print(platform.machine())" && \
        conda config --env --set subdir osx-64
fi
if [[ $f != "" ]]; then 
    if [ -d ~/omics_workshop ]; then
        echo "Removing existing 'omics_workshop' directory..."
        rm -rf ~/omics_workshop
    fi 
fi
# Check if ~/omics_workshop exists and if it does then cd and git pull force changes otherwise clone new one
if [ -d ~/omics_workshop ]; then
    echo "Updating 'omics_workshop' GitHub repository..."
    cd ~/omics_workshop && git pull
else
    echo "Cloning 'omics_workshop' GitHub repository..."
    git clone https://github.com/jhuapl-bio/omics_workshop ~/omics_workshop
fi

## make sure we update Taxonomy for krona to work! Requires internet
echo "Downloading taxonomy information. Requires internet connection. This will take some time...."
# bash ktUpdateTaxonomy.sh
if [ -s $(dirname $(which ktImportTaxonomy))/../opt/krona/taxonomy/taxonomy.tab ] && [[ $f == "" ]]; then
    echo "Krona taxonomy.tab exists, skipping...."
else
    mkdir -p  ~/omics_workshop/downloads/ && conda activate omics_workshop \
        && wget https://github.com/jhuapl-bio/datasets/raw/main/databases/ncbi/taxonomy.tab.gz \
        -O ~/omics_workshop/downloads/taxonomy.tab.gz \
        && gzip -f -d ~/omics_workshop/downloads/taxonomy.tab.gz  && \
        mv ~/omics_workshop/downloads/taxonomy.tab $(dirname $(which ktImportTaxonomy))/../opt/krona/taxonomy/
fi
# check if ~/omics_workshop/databases/k2_viral/hash.k2d and is size >0 then skip else download
if [ -s ~/omics_workshop/databases/k2_viral/hash.k2d ] && [[ $f == "" ]]; then
    echo "Kraken2 viral database already exists."
else
    echo "Downloading Kraken2 viral database..."
    wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz \
        -O ~/omics_workshop/downloads/k2_viral.tar.gz && \
        mkdir -p ~/omics_workshop/downloads/k2_viral && \
        tar -xvzf ~/omics_workshop/downloads/k2_viral.tar.gz -C ~/omics_workshop/downloads/k2_viral  && \
        mv ~/omics_workshop/downloads/k2_viral ~/omics_workshop/databases/k2_viral
fi

tar -k -xvzf ~/omics_workshop/databases/test_metagenome.tar.gz -C ~/omics_workshop/databases/

echo "Setup complete."
