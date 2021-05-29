DIR="{absolute path to CLASHanalyst folder}/CLASHanalyst" #etc. /home/t50504

bash un_zip.sh $DIR/run_sample/ #directory to inputfile including: target.zip, regulator.zip, read.zip 

cd $DIR/run_sample/

#### preprocess and build ref
/usr/bin/time -f "	%E real,	%U user,	%S sys" -a -o $DIR/run_sample/time_log make -f $DIR/suite/bin/makefile preprocess qc=trim_galore trim=30 link=AGATCGGAAGAG len=17 slen=70 rc=None fd=None in=hyb_file.fastq

/usr/bin/time -f "	%E real,	%U user,	%S sys" -a -o $DIR/run_sample/time_log make -f $DIR/suite/bin/makefile build reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq

### run three pipeline

mkdir $DIR/run_sample/pir $DIR/run_sample/hyb $DIR/run_sample/clan
# run pirtarbase algorithm
/usr/bin/time -f "	%E real,	%U user,	%S sys" -a -o $DIR/run_sample/time_log make -f $DIR/suite/bin/makefile detect way=pir llen=17 reg_mis=0 tran_mis=0 hmax=10 reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq

/usr/bin/time -f "	%E real,	%U user,	%S sys" -a -o $DIR/run_sample/time_log make -f $DIR/suite/bin/makefile analyse reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq way=pir

mv hyb_file_step5.csv pir/

rm hyb_file_step4.csv
# run hyb algorithm
/usr/bin/time -f "	%E real,	%U user,	%S sys" -a -o $DIR/run_sample/time_log make -f $DIR/suite/bin/makefile detect way=hyb hval=0.1 hmax=10 gmax=4 reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq

/usr/bin/time -f "	%E real,	%U user,	%S sys" -a -o $DIR/run_sample/time_log make -f $DIR/suite/bin/makefile analyse reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq way=hyb

mv hyb_file_step5.csv hyb/

rm hyb_file_step4.csv
# run clan algorithm
/usr/bin/time -f "	%E real,	%U user,	%S sys" -a -o $DIR/run_sample/time_log make -f $DIR/suite/bin/makefile detect way=clan llen=17 hmax=10 gmax=4 reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq

/usr/bin/time -f "	%E real,	%U user,	%S sys" -a -o $DIR/run_sample/time_log make -f $DIR/suite/bin/makefile analyse reg=reg_file.fasta tran=tran_file.fasta in=hyb_file.fastq way=clan

mv hyb_file_step5.csv clan/

rm hyb_file_step4.csv


