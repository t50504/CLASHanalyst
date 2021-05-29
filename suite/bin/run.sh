
#for i in $(ls -d *|grep "^[G|S]")
#for i in $(ls -d SRR6512655_up*)
#for i in $(ls -d *|grep "^[G|S].........$")
#do
    #echo $i
    #name=$i

    #cd $name
    #bash ../fun_comp/compare.sh $name
    #bash ../fun_comp_step6/compare.sh $name
    #cd ../
#done


#name=GSM1219487
#mkdir $name

#cd $name
#bash ../fun_comp/compare.sh $name
#bash ../fun_comp_step6/compare.sh $name
#cd ../


for i in $(ls -d *|grep "^[G|S].........$")
do
    echo $i
    name=$i

    cd $name
    python ../fun_comp/compare.py boxplot $name

    cd ../
done


#printf "" > all_count.txt
#for i in $(ls -d *|grep "^[G|S].........$")
#do
    #echo $i
    #name=$i

    #cd $name
    #pir=$(cat pir.tab|wc -l)
    #hyb=$(cat hyb.tab|wc -l)
    #clan=$(cat clan.tab|wc -l)
    #cd ../
    #echo ${name:8:2},$pir,$hyb,$clan >> all_count.txt
#done

################## hits ################
#printf "" > hit_fold_count.csv
#printf "" > hit_readCount_count.csv
#for i in $(ls -d *|grep "^[G|S].........$")
#do
    #echo $i
    #name=$i

    #cd $name
    #cat hits_RNAfold >> ../hit_fold_count.csv
    #cat hits_read_count >> ../hit_readCount_count.csv
    #cd ../
#done


