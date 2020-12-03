set -euxo pipefail

samtools faidx $2
bedtools makewindows -g $2.fai -w 5000 -s 1000 | bedtools getfasta -fi $2 -bed - > tmp.fasta
#minimap2 -f 0.00001 -t 160 --no-long-join -DP -c --eqx -x map-ont $1 tmp.fasta > $3 
minimap2 -t 160 -r 500000 --secondary=no --eqx -c -x map-ont $1 tmp.fasta > $3 
./convert_paf.py -i 97 --frac $3 $3.tbl


