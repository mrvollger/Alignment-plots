set -euxo pipefail
minimap2 -f 0.00002 -t 160 --no-long-join -DP -c --eqx -x map-ont $1 $2 > $3 
convert_paf.py $3 $3.tbl


