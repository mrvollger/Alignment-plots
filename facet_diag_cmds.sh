set -euxo pipefail

cat cens.fofn | parallel --will-cite -j 1 -n 2 \
    "samtools faidx {2}; bedtools makewindows -g {2}.fai -w 5000 -s 1000 | bedtools getfasta -fi {2} -bed - | minimap2 -f 0.00001 -t 160 -r 500000 --secondary=no --eqx -c -x map-ont {1} /dev/stdin" > all.cen.paf
./convert_paf.py --frac all.cen.paf all.cen.paf.tbl


