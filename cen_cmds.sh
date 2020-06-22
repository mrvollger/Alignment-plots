set -euxo pipefail
#~mvollger/software/Winnowmap/bin/winnowmap -W bad_k15_mers.txt -t 4 --no-long-join -DP -c --eqx -x map-ont $1 $2 > $3 

samtools faidx $2

bedtools makewindows -g $2.fai -w 5000 | bedtools getfasta -fi $2 -bed - > tmp.fasta

minimap2 -f 0.00002 -t 160 --no-long-join -DP -c --eqx -x map-ont $1 tmp.fasta > $3 

../convert_paf.py $3 $3.tbl


