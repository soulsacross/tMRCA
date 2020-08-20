#!/bin/bash
: '
The following script calculates
the tMRCA for selected alleles
'

pop=''
file=''
out=''
verbose='false'

print_usage() {
  printf "Usage: -p specifies population \n -f specifies input .csv file \n -o specifies output dir \n";
}

while getopts 'p:f:o:v' flag; 
do
  case "${flag}" in
    p) pop="${OPTARG}" ;;
    f) file="${OPTARG}" ;;
    o) out="${OPTARG}" ;;
    v) verbose='true' ;;
    *) print_usage
       exit 1 ;;
    esac
done

echo "===================================================================================="
echo "                         Getting EHH for markers."
echo "                                                                                   "

#i=0
#while  IFS=$',\t\r\n'  read -r l   
#    do chr=$(echo $l | awk -F, '{print $1}')
#    pos=$(echo $l | awk -F, '{print $2}')
#    master_array[i]=${ARRAY[i]}
#    ((i++))
#    [ -f $out/$chr.$pop.ehh.$pos.out] || selscan --ehh $pos --tped ~/Dropbox/input_chr16/TPED/$chr.$pop.scans.input.tped.gz --out $out/$chr.$pop --cutoff 0.005 # running selscan (output is completly inaccurate)
#     [ -f $out/$pop.$chr.$pos.ehh.gz ] || ~/Dropbox/input_chr16/PopLDdecay/bin/PopLDdecay -InVCF /Volumes/PEDRO\ 70GB/Genetics/VCFphased/VCFphased/ALL."${chr:3:2}".PHASED.vcf.gz -OutStat $out/$pop.$chr.$pos  -SubPop ~/Dropbox/input_chr16/Pops/$pop.txt -EHH "${chr:3:2}":$pos
#done < $file

echo "                                                                                   "
echo "There are`ls $out/*.ehh.gz | wc -l` putatively selected variant(s) in your dataset."
echo "===================================================================================="
echo "                         Beggining TMRCA calculation."
echo "                                                                                   "


#for f in $(ls $out | grep $pop.*ehh.gz | sort -V)
#    do zcat < $out/$f |
#        awk '{
#        if($4>0 && $4<0.25 && $3<0){ 
#        upstream=$0;                   # < this will keep last record
##        print left;
#        }else if($4>0 && $4<0.25 && $3>0){
#        downstream=$0;
##        print right;
#        exit                           # < this will keep first record
#    }
#    }END{print upstream "\n" downstream}' | tee -a $out/$pop.output.txt
#done
#
#chmod 777 $out/$pop.output.txt

echo "                                                                                   "
echo "===================================================================================="
echo "                       Getting genetic distances in cM."

while read l 
    do query_chr=$(echo $l | awk '{print "chr"$1}') 
    query_pos=$(echo $l | awk '{print $2}') 
    cM=$(zcat < ~/Dropbox/input_chr16/TPED/$query_chr.$pop.scans.input.tped.gz | grep -w $query_pos | awk '{print $3}') 
    echo $l " "$cM
done < $out/output.txt | 

sed -e 's/  */ /g' | # space to tab substitution not currently working in Mac os X, should be replaced by \t in Linux
awk '{print $0 } NR%2{p=$7;next}{print($7-p)" "((-25*100*log(0.25))/($7-p)) " " "<<<- tMRCA" }' | # awk -v n=1 '1; NR % n == 0 {print ""}' | #commented cmd would skip evey other line in print statement (I would like than to to this in an entirely new column, which I couldn't get)
tee $out/GeneticDistances_tMRCA.txt # 

#[ -f $out/output.txt] && rm $out/output.txt

#bc -l <<< "-25*100*l(0.25)/0.236672"

