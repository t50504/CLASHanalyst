CLASH Analyst README

## Pre-requisite
- workable OS: Linux ubuntu 16.04 (LTS)
- All the toolkit need to be installed
---
- SRA-toolkit
```bash=
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-ubuntu64.tar.gz
tar -zxvf sratoolkit.2.9.2-ubuntu64.tar.gz
export PATH=$PATH:/{path to program}/sratoolkit.2.9.2-ubuntu64/bin
`````
- Bowtie
```bash=
sudo apt-get update
sudo apt-get install bowtie
sudo apt-get install bowtie2
`````
- ViennaRNA
```bash=
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.4.tar.gz
tar -zxvf ViennaRNA-2.4.4.tar.gz
cd ViennaRNA-2.4.4
./configure 
make
sudo make install
export PATH=$PATH:/{path to program}/ViennaRNA-2.4.4/bin
`````
- FastQC
```bash=
sudo apt-get install fastqc
`````
- cutadapt
```bash=
pip install cutadapt
`````
- TrimGalore
```bash=
# Check that cutadapt is installed
cutadapt --version
# Check that FastQC is installed
fastqc -v
# Install Trim Galore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
export PATH=$PATH:/{path to program}/TrimGalore-0.6.6
`````
- Flexbar
```bash=
sudo apt install aptitude
sudo aptitude install flexbar
`````
- Libgtextutils
```bash=
wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz 
tar -zxvf libgtextutils-0.7.tar.gz 
cd libgtextutils-0.7 
./configure 
make 
sudo make install
`````
- Fastx
```bash=
wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2 
tar -xjf fastx_toolkit-0.0.14.tar.bz2 
cd fastx_toolkit-0.0.14/ 
./configure 
make #if fail, run sudo apt-get install gcc g++ pkg-config 
sudo make install
`````
- pblat
```bash=
wget https://github.com/icebert/pblat/archive/master.zip
unzip master.zip
cd pblat-master/
make
export PATH=$PATH:/{path to program}/pblat-master
`````
- blat
```bash=
wget -nc http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip
unzip blatSrc35.zip
uname -a# check linux system type
export MACHTYPE=x86_64 #according to your system type
mkdir -p ~/bin/x86_64
cd blastrc
make #if fail, run sudo apt-get install libpng-dev libssl-dev
export PATH=${HOME}/bin/x86_64:$PATH
`````
- blastn
```bash=
sudo apt-get install ncbi-blast+
`````
- parallel
```bash=
sudo apt-get install parallel
`````
## how to quick start?
- Download
```bash=
git clone git@github.com:t50504/CLASHanalyst.git
`````
- Modify the run.sh
```bash=
cd ./CLASHanalyst/run_sample
vim run.sh
`````
- In the run.sh
```bash=
#revise first line path
DIR="/{absolute path to CLASHanlyst folder}/"
`````
- Run the sample with supplied file.
```bash=
bash run.sh
`````
- Result will generated in /run_sample/{hyb, clan, pir} folder, with RNA interactions information in hyb_file_step5.csv file. 
