# ConvergenceAnalysis

A repo containing the scripts to perform the convergence analysis on GEA results. What you need is a set of WZA results. These files are assumed to be in 

The first to do is to add the orthogroup info to each file. I do that using the ```bin/annotateWZAfiles.py``` script. This is a vitally important step.

This script can be used like this, for example:

```
python ../bin/annotateWZAfiles.py --anno Orthogroups.tsv --wza  JP_all_wza_tc.txt --species Pi --output JP_all_wza_tc_OG.txt

```
Where, Orthogroups.tsv - is the output from Orthofinder; JP_all_wza_tc.txt - is the WZA output file; Pi - is the species-code for the pine pan-genome; P_all_wza_tc_OG.txt - is the name of the output script

Once you've run that, you'll have a set of files with both the Zw scores from the WZA as well as orthogroup information for each species. We'll put that into the script analysing convergent evolution.

The script ```bin/convergenceAnalysis.py``` implements both multi-species tests for convergent evolution as well as pairwise tests. You'll have to take care which species you are providing to the program.

The script has help info accessed via ```python bin/convergenceAnalysis.py --help``` if you want to take a look at the list of arguments it takes. 

Here's an example of how you might use the convergence script:

```
python bin/convergenceAnalysis.py -i species_1.ENV.csv species_2.ENV.csv species_3.ENV.csv species_4.ENV.csv -o ENV.convergence --alpha_1 0.01 --alpha_2 0.05 -v --allPairs
```

That command will run the script on data for species 1, 2, 3, 4 using alpha_1 = 0.01 and alpha_2 = 0.05. The ```-v``` is for verbose mode, so the script will tell you what it's doing. ```--allPairs``` runs the analysis using the multi-species test as well as on all pair wie comparisons.