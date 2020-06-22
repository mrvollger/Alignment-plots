set -euxo pipefail
~mvollger/software/Winnowmap/bin/meryl count k=15 output merylDB_k15 $1
~mvollger/software/Winnowmap/bin/meryl print greater-than distinct=0.9998 merylDB_k15 > bad_k15_mers.txt

#~mvollger/software/Winnowmap/bin/winnowmap -W bad_k15_mers.txt -t 4 --no-long-join -DP -c --eqx -x map-ont $1 $2 > $3 
minimap2 -f 100000 -t 4 --no-long-join -DP -c --eqx -x map-ont $1 $2 > $3 

