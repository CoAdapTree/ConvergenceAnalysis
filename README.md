# ConvergenceAnalysis


**At this stage this analysis isvery preliminary and subject to change** 


A repo containing the scripts to perform the convergence analysis on GEA results. What you need is a set of WZA results. I used the files "DFI_all_wza_tc.txt, WL_all_wza_tc.txt" etc.

The first thing I do is to add the orthogroup info to each file. I do that using the ```bin/annotateWZAfiles.py``` script.

For example, for Jack pine: 

```
python ../bin/annotateWZAfiles.py --anno Orthogroups.tsv --wza  JP_all_wza_tc.txt --species Pi --output JP_all_wza_tc_OG.txt

```
Where, Orthogroups.tsv - is the output from Orthofinder; JP_all_wza_tc.txt - is the WZA output file; Pi - is the species-code for the pine pan-genome; P_all_wza_tc_OG.txt - is the name of the output script

A file like that is generated for each species. The species files are then analysed using the ```bin/WZA_pMaxAnalysis.R``` script.
