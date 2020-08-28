r=/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/t2t_rel3_glchr8/v9_hybrid/data/glchr8v9_cen.fa
q=/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/hg00733_chr8/cen/mat/data/hg00733_mat_cen8v1.trimmed.fa
sm=13vs733


./cen_cmds.sh $r $q nobackups/$sm.paf


q2=/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/hg00733_chr8/cen/pat/v3/data/hg00733_pat_cen8v3.trimmed.fa
./cen_cmds.sh $r $q2 nobackups/13vs733pat.paf



