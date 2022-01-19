#!/bin/bash
: '
The following script calculates the tMRCA for selected alleles.
Provide coordinates in (-f), population name (1kG, e.g. "-p CLM"), output dir (-o). Needs also dir to PopLDdecay binary, vcf file and population list (l40) and TPED (l76) edits in-script.
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

[ -f $out/$pop.output.txt ] && rm $out/$pop.output.txt

# Reading coordinates for selected variants
{
read
while  IFS=$',\t\r\n'  read -r l   
   do chr=$(echo $l | awk -F, '{print $1}')
   pos=$(echo $l | awk -F, '{print $2}')
   [ -f $out/$pop.$chr.$pos.ehh.gz ] || ~/Dropbox/input_chr16/PopLDdecay/bin/PopLDdecay -InVCF /Volumes/PEDRO\ 70GB/Genetics/VCFphased/VCFphased/ALL."${chr:3:2}".PHASED.vcf.gz -OutStat $out/$pop.$chr.$pos  -SubPop ~/Dropbox/input_chr16/Pops/$pop.txt -MaxDist 1000 -EHH "${chr:3:2}":$pos
done 
}< $file

echo "                                                                                   "
echo "There are`ls $out/*.ehh.gz | wc -l` putatively selected variant(s) in your dataset."
echo "===================================================================================="
echo "                       Getting haplotype breakpoints."
echo "                                                                                   "

# Get coordinates of flanking variants that have EHH < 0.25, that is, Pr[Homoz] drops below 25% and thus >75% of chromosomes have recombined off of the selected haplotype, following that we'll then go get the recombination distance (r) 
# These correspond to putative recombination breaks, get physical positions on them (bp)
for f in $(ls $out | grep $pop.*ehh.gz | sort -V)
   do zcat < $out/$f |
       awk '{
       if($4>0 && $4<0.25 && $3<0){ 
       upstream=$0;                   # < this will keep last record
#        print upstream;
       }else if($4>0 && $4<0.25 && $3>0){
       downstream=$0;
#        print downstream;
       exit                           # < this will keep first record
   }
   }END{print upstream "\n" downstream}' | tee -a $out/$pop.output.txt
done
#
chmod 777 $out/$pop.output.txt

echo "                                                                                   "
echo "===================================================================================="
echo "             Getting genetic distances in cM and beginning TMRCA calculation."

#Get genetic distances for such breakpoints of selected haplotypes and calculate tMRCA.
while read l 
    do query_chr=$(echo $l | awk '{print "chr"$1}') 
    query_pos=$(echo $l | awk '{print $2}') 
    cM=$(cat < ~/Dropbox/input_chr16/GRCh37/plink.$query_chr.GRCh37.map | grep -w $query_pos | awk '{print $3}') 
    echo $l " "$cM
done < $out/$pop.output.txt | 

#cat $out/$pop.output.txt |
#sed -e 's/  */ /g' | # space to tab substitution not currently working in Mac OSX, should be replaced by \t in Linux
awk '{print $0} 
    NR%2{ p=$7;next}{if(($7-p) != 0){
        print ($7-p)" "(-25*100*log(0.25))/($7-p) " " "<<<- tMRCA"
        }
}' | # awk -v n=1 '1; NR % n == 0 {print ""}' | #commented cmd would skip every other line in print statement (I would like then to print this in an entirely new column, which I couldn't get)
tee $out/$pop.GeneticDistances_tMRCA.txt # 


echo "                       DONE."
