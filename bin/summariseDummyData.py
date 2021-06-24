## A script to summarise the convergence analysis results

import argparse, glob
import pandas as pd
import numpy as np
from scipy.stats import rankdata

## A function grabbed from SO to do perform the BH procedure on the p-Values you're working with
def p_adjust_bh( p ):

	p = np.asfarray(p)

	by_descend = p.argsort()[::-1]

	by_orig = by_descend.argsort()

	steps = float(len(p)) / np.arange(len(p), 0, -1 )

	q = np.minimum( 1, np.minimum.accumulate( steps * p[by_descend] ) )

	return q[by_orig]

def main():

## Define command line args
	parser = argparse.ArgumentParser(description="This takes resultas of the WZA sims and ")

	parser.add_argument("--input", "-i",
			required = True,
			dest = "input",
			type = str, 
			help = "Directory containing CSV files for the WZA results.")

	parser.add_argument("--output", "-o",
			required = True,
			dest = "output",
			type = str, 
			help = "The name you want to give to the output file")

	args = parser.parse_args()

	print( args.input + "random_set_*_multi_species.csv" )

	multi_species_output_list = []

## Do the stuff with multi-species analyses

	for g in glob.glob(args.input + "/random_set_*.multi_species.csv"):
		
		multi = pd.read_csv(g) 
		
		true_converges = (multi.converge > 2).sum()
## Get the genes that were tested for convergence
## i.e. those that had at least one hit
		multi_hits_strim = multi[ ( ~np.isnan( multi.converge_pVal ) ) ].copy()

## Do the multiple-tests correction
		multi_hits = multi[ ( ~np.isnan( multi.converge_pVal ) ) ]

		multi_hits_strim["converge_pVal_adjust"] = p_adjust_bh( np.array( multi_hits.converge_pVal ) )

		identified = np.argmin( multi_hits_strim.iloc[:,1:5].values, axis = 1) + 2

		multi_hits_strim["ident"] = identified

		multi_hits_strim_sig = multi_hits_strim[ multi_hits_strim.converge_pVal_adjust < 0.05 ]

		false_positives = (multi_hits_strim_sig.converge < 2).sum()

		good_guesses = (multi_hits_strim_sig.converge == multi_hits_strim_sig.ident).sum()

		multi_species_output_list.append( [ true_converges, good_guesses,  false_positives, multi_hits_strim_sig.shape[0] ] )
		
	output_header = ["TruePositives",
						"NumberOfGoodGuesses",
						"FalsePositives",
						"NumberOfHits"]
	
	result = pd.DataFrame( np.array(multi_species_output_list), columns = output_header)
	
	result.to_csv( args.output, index = False )
	

if __name__ == "__main__":
	main()
