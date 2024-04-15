#!/bin/bash -i

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
    Darwin*)    machine="CONDA_SUBDIR=osx064";;
    *)          machine=""
esac
echo ${machine}

$HOME/miniconda3/bin/conda init bash
source ~/miniconda3/etc/profile.d/conda.sh

# Wait for the user to close and reopen their shell or source their bashrc
echo "Please close and reopen your terminal, or run 'source ~/.bashrc', then rerun this script."
conda env list
# Check if the omics_workshop environment exists
if  conda env list | grep -q 'omics_workshop'; then
    echo "Removing existing 'omics_workshop' environment..."
    conda env remove -y -n omics_workshop
fi

echo "Creating 'omics_workshop' environment..."



$machine conda create -y -c bioconda -n omics_workshop bowtie2 minimap2 kraken2 krona fastqc samtools bcftools git python

if [ -d ~/omics_workshop ]; then
    echo "Removing existing 'omics_workshop' directory..."
    rm -rf ~/omics_workshop
fi

echo "Cloning 'omics_workshop' GitHub repository..."
git clone https://github.com/jhuapl-bio/omics_workshop ~/omics_workshop

## make sure we update Taxonomy for krona to work! Requires internet
echo "Downloading taxonomy information. Requires internet connection. This will take some time...."
# bash ktUpdateTaxonomy.sh
mkdir -p  ~/omics_workshop/downloads/ && conda activate omics_workshop \
    && wget https://github.com/jhuapl-bio/datasets/raw/main/databases/ncbi/taxonomy.tab.gz \
    -O ~/omics_workshop/downloads/taxonomy.tab.gz \
    && gzip -f -d ~/omics_workshop/downloads/taxonomy.tab.gz  && \
    mv ~/omics_workshop/downloads/taxonomy.tab $(dirname $(which ktImportTaxonomy))/../opt/krona/taxonomy/

wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz \
    -O ~/omics_workshop/downloads/k2_viral.tar.gz && \
    mkdir -p ~/omics_workshop/downloads/k2_viral && \
    tar -xvzf ~/omics_workshop/downloads/k2_viral.tar.gz -C ~/omics_workshop/downloads/k2_viral  && \
    mv ~/omics_workshop/downloads/k2_viral ~/omics_workshop/databases/k2_viral

tar -k -xvzf ~/omics_workshop/databases/test_metagenome.tar.gz -C ~/omics_workshop/databases/

echo "Setup complete."