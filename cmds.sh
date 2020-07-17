set -euxo pipefail
minimap2 -f 50000 -t 160 -s 1000 --no-long-join -DP -c --eqx -x asm20 $1 $2 > $3 
./convert_paf.py $3 $3.tbl
./make_html.py $3.tbl $3.html


