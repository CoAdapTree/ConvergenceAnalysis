import argparse, itertools
import pandas as pd
import numpy as np
from scipy.stats import poisson, binom, chi2

## A simple function for doing the Dunn-Sidak correction for multiple hits on the minimum p-value
def p_value_adjust( p_vec ):
	p_adj = 1 - (1 - p_vec.min()) ** len(p_vec)
	return p_adj
	
def FishersCombinedProb( p_vec ):
## Get the length of the p-value vector
	k = p_vec.shape[0]
	
## Calculate the test statistic
	test_stat = -2. * np.log( p_vec ).sum()
	
## Convert the test_stat to a pValue
	p_val = 1 - chi2.cdf( test_stat, 2*k )
	
	return( p_val )

## This function will read a dataframe and calculate empirical p-values from the summary statistic named
def parse_result(input_name, statistic, descending = True, DunnSidak = False):

	df = pd.read_csv(input_name) 
	
	name = input_name.split("/")[-1].split(".csv")[0]

# calculate the empirical p-values for each gene

	df["empirical_p"] = df[statistic].rank(pct = True, method = "average", ascending =  not descending) 

	if DunnSidak:
# Group by orthogroup and apply the p_value_adjust function (i.e. the Dunn-Sidak correction)
		new_df = df.groupby("orthogroup")["empirical_p"].apply( p_value_adjust ).to_frame(name).reset_index()
	else:
# Group by orthogroup and apply the p_value_adjust function (i.e. the Dunn-Sidak correction)
		new_df = df.groupby("orthogroup")["empirical_p"].apply( FishersCombinedProb ).to_frame(name).reset_index()
	
	return new_df



# This function will parse the results stored in each dataframe and combine them into a single dataframe
def prep_results(input_data, summary_stat):
	
	my_DFs = [ parse_result(dataframe, summary_stat) for dataframe in input_data ] 
	
	all_DFs = [dataF.set_index(["orthogroup"]) for dataF in my_DFs]

	merged = pd.concat(all_DFs, axis = 1 , sort = True).rename_axis("orthogroup").reset_index()

	return(merged)



# This function will parse the results stored in each dataframe and combine them into a single dataframe
def prep_results_tests( input_data, descending = True ):
	
# Read the input dataframe
	my_DF = pd.read_csv(input_data) 

	genes =  ["ortho" + str(ll) for ll in range(my_DF.shape[0]) ]
	
# Make the data 
	p_val_dataframe = pd.DataFrame( [ my_DF["Z_1"].rank(pct = True, method = "average", ascending =  not descending),
					my_DF["Z_2"].rank(pct = True, method = "average", ascending =  not descending),
					my_DF["Z_3"].rank(pct = True, method = "average", ascending =  not descending),
					my_DF["Z_4"].rank(pct = True, method = "average", ascending =  not descending),
					my_DF["Z_5"].rank(pct = True, method = "average", ascending =  not descending)] ).transpose()
	p_val_dataframe.insert( loc=0, value= genes, column="orthogroup" )

	num_LA = my_DF.converge.copy()
	

	return( p_val_dataframe, num_LA )



## Generate permutations for 2-way test

def permute_overlaps_for_nullZ( pValArray, alpha_1, alpha_2, perms = 1000 ):

# Make a container for the results
	perms_array = np.zeros(perms)

	for i in range( perms ):

## Permute the genes - would only have to do one, but I'll do both just for kicks 
		pValArray[:,0] = np.random.permutation( pValArray[:,0] )
		pValArray[:,1] = np.random.permutation( pValArray[:,1] )

## Remove rows with missing data
		perm_p_value_array_na_omit = pValArray[~np.isnan(pValArray).any(axis = 1),:]

## Sort the rows of the dataframe
		perm_sorted_array = np.sort(perm_p_value_array_na_omit)

## Get the rows that have the first value less than alpha_1 and the second value less than alpha_2
		perm_hits =  perm_sorted_array[(perm_sorted_array[:,0]<alpha_1)&(perm_sorted_array[:,1]<alpha_2 )]

## Store the results in an array
		perms_array[i] = perm_hits.shape[0]

## Calculate outlier thresholds from the permutations 
	thresholds = np.quantile( perms_array, 
								np.array([ 0.95, 0.99, 0.999] ) ) 

	return( thresholds )

def null_Z( DF, alpha_1, alpha_2 , species = [ 0 , 1 ]):
## Convert the DF of p-values to a numpy array
	p_value_array = DF[ list(species)].values 

## Get the outlier  thresholds for 
	null_Z_thresholds = permute_overlaps_for_nullZ( p_value_array.copy(), alpha_1, alpha_2 )

## Get a count of the orthogroups before removing missing data
	n_orthogroups_before_naomit = p_value_array.shape[0]

## Remove rows with missing data
	p_value_array_na_omit = p_value_array[~np.isnan(p_value_array).any(axis = 1),:]

## Get a count of the orthogroups after removing missing data
	n_orthogroups_after_naomit = p_value_array_na_omit.shape[0]

## Sort the rows of the dataframe
	sorted_array = np.sort(p_value_array_na_omit)

## Get the rows that have the first value less than alpha_1 and the second value less than alpha_2
	hits =  sorted_array[(sorted_array[:,0]<alpha_1)&(sorted_array[:,1]<alpha_2 )]

## calculate the expected number of overlaps
#	expected_no_missing = alpha_1 * alpha_2 * 2 * n_orthogroups_before_naomit
	expected_missing = alpha_1 * alpha_2 * 2 * n_orthogroups_after_naomit

## calculate the expected number of overlaps
#	expected_no_missing = alpha_1 * alpha_2 * 2 * n_orthogroups_before_naomit
#	expected_missing = alpha_1 * alpha_2 * 2 * n_orthogroups_after_naomit

## calculate the expected number of overlaps at 95% confidence
#	expected_no_missing_95 =  poisson.ppf( 0.95, expected_no_missing) 
#	expected_no_missing_99 =  poisson.ppf( 0.99, expected_no_missing) 

## calculate the expected number of overlaps at 99% confidence
#	expected_missing_95 =  poisson.ppf( 0.95, expected_missing) 
#	expected_missing_99 =  poisson.ppf( 0.99, expected_missing) 

## Return result
	return({ "number_of_overlaps":hits.shape[0], 
#			"expected_no_missing":expected_no_missing, 
			"expected_missing":expected_missing,
			"expected_missing_95":null_Z_thresholds[0],
			"expected_missing_99":null_Z_thresholds[1],
			"expected_missing_999":null_Z_thresholds[2],
#			"expected_no_missing_95":expected_no_missing_95,
#			"expected_no_missing_99":expected_no_missing_99,
#			"expected_missing_95":expected_missing_95,
#			"expected_missing_99":expected_missing_99,
			"alpha_1":alpha_1,
			"alpha_2":alpha_2})



def p_jth(pVec, j):

	j = j
	
	k = pVec.shape[0]

	jth = pVec[j-1]

	if j > k:
		return np.nan

	if j < k:
		prob = jth/pVec[j]  ## This index is really annoying as it is making use of 0-indexing instead of adding 1
							## Here's the R Code to demonstrate: prob = pjth / sort(plist)[j+1]
	else:
		prob = pVec.max()

	return binom.pmf( j, j, prob)


## Here's the function with the actual pMax calculations:
def pMax_guts( pVal_array ):

#	print()
## Things we'll want to return:
# 1. The pValue of convergence
# 2. The vector of pValues tested - to see where the signal came from
# So the output will have k+1 values where k is the length of the input array

	num_na = np.isnan(pVal_array).sum()
## If there are only nans in the array, return NaN
## Or, if there is only one NaN in the array, return NaN
	if num_na >= pVal_array.shape[0]-1:
		out_array = np.empty(pVal_array.shape[0] + 1)
		out_array[:] = np.nan

	else:
		p_array = pVal_array[~np.isnan(pVal_array)]
		pMax_results = [ p_jth(p_array, i+1) for i in range(len(p_array))]

		out_array = np.array( [p_value_adjust(np.array(pMax_results))] +pMax_results + [np.nan for q in range( num_na )] )

	return out_array

def MultiSpeciesConvergence(big_DF, alpha):

## Convert the pandas dataframe into an array of p-values
	p_value_array = big_DF.iloc[:,1:].values  

## Sort the rows of the dataframe
	sorted_array = np.sort(p_value_array)

## Masks out those genes that have no hit in the first column of the sorted array
	masking_array = sorted_array[:,0] >= alpha ## Need to do multi-hit correction here?
	sorted_array[masking_array,:] = np.nan

## Only retain information about the k-1th highest p-values	
	pMax_input = sorted_array[:,1:]

## Apply the pMax function to the genes in the top set
	pMax_results = []

	for ar in pMax_input:
		pMax_results.append( pMax_guts( ar ) )

## Make a DF of the output
	pMax_DF = pd.DataFrame( pMax_results )

## Rename the header of the output
	pMax_DF.columns = ["converge_pVal"] + ["p"+str(r) for r in range(2, sorted_array.shape[1]+1)]

## Add orthogroup info to the file
	pMax_DF["orthogroup"] = big_DF["orthogroup"]
	
	return pMax_DF




def main():

## Define command line args
	parser = argparse.ArgumentParser(description="This script implements the statistical tests for convergence that we (Singh, Booker et al. have developed for the analysis of the data coming out of the CoAdapTree project. If you provide two input files only, the script will run the two-way analysis. If you provide more than 2 files the script runs the multi-spcies pMax test ONLY as a default behaviour. Use the [ --all ] flag if you want to do all pairwise comparisons as well as the multi-species test.")

	parser.add_argument("--input", "-i",
			required = True,
			dest = "input",
			type = str, 
			nargs = "+",
			help = "CSV files containing the results from the multiple tests. Give the CSV files nice names like 'DF_AHM.csv'. That way the final dataframe will have nice formatting")

	parser.add_argument("--alpha_1", 
			required = True,
			dest = "alpha_1",
			type = float, 
			help = "The initial threshold to identify outliers")

	parser.add_argument("--alpha_2", 
			required = False,
			dest = "alpha_2",
			type = float, 
			help = "The second threshold to identify outliers (only used for the 2-way test)")

	parser.add_argument("--outputPrefix", "-o",
			required = True,
			dest = "outputPrefix",
			type = str, 
			help = "The name you want to give to the output file[s] (give a prefix)")

	parser.add_argument("--stat", "-s",
			required = False,
			dest = "stat",
			type = str, 
			help = "The of the statistic you want to use for testing convergent evolution. This HAS to have the same name in all of the files you provide",
			default = "Z")

	parser.add_argument("--allPairs", 
			required = False,
			dest = "allPairs",
			action = "store_true",
			help = "Do you want to examine convergence between all pairs of species in the analysis? - NOTE This only has an effect if you supply >2 input files")
			
	parser.add_argument("--verbose", "-v",
			required = False,
			dest = "verbose",
			action = "store_true",
			help = "Run the program in verbose mode?")

	parser.add_argument("--test", "-t",
			required = False,
			dest = "test",
			action = "store_true",
			help = "Run the program in test mode? This assumes the results are all in one file")
	args = parser.parse_args()


	
	## Quit if the user only gave one input file
	if len(args.input) <= 1 and not args.test:
		
		if args.verbose: print("You only provided 1 file. You need to provide at least 2, bozo")
		return
		
	elif len( args.input ) == 2:
		if args.verbose: print("\nTwo species mode\n")
	
	elif len( args.input ) > 2:
		if args.verbose: print("\nMulti-species mode\n")
		if args.allPairs:
			if args.verbose: print("Running all comparisons")
		
	if args.verbose: print("\nCombining results from each file into a single dataframe for all species \n ")

	if args.test:
		p_value_dataframe, converge = prep_results_tests( args.input[0] ) 
		print( p_value_dataframe )
	else:
		p_value_dataframe = prep_results( args.input , args.stat) 
	
	if args.verbose: print("Here's a peek at the combined dataframe:", p_value_dataframe)

	species = list(p_value_dataframe)[1:]

## If more than 2 species AND you want to do the pMax test
	if len(args.input) > 2 or args.test:
		if args.verbose: print("\nPerforming the pMax test on the following data:", species)

## Perform the pMax test on the data that you provided
		multi_species_result = MultiSpeciesConvergence(p_value_dataframe, args.alpha_1)

## Save the result to disc using the name that was given
		if args.test:
			multi_species_result["converge"] = converge

		multi_species_result.to_csv(args.outputPrefix + ".multi_species.csv", index = False)

## Do the null_Z test on pairs of species...
	if args.allPairs or len(args.input) == 2:

		if args.verbose: print("\nThis is how the species you are analysing will appear in the output:", species)
		null_Z_results_for_pairs = []

		for pair in itertools.combinations( species, 2):
			if args.verbose: print("\nAnalysing:", pair,  " using the Null-Z test\n")
## Run the two-way test for the current pair 
			zull_Z_result_for_pair = null_Z( p_value_dataframe, 
											alpha_1 = args.alpha_1, 
											alpha_2 = args.alpha_2,
											species = pair)
			zull_Z_result_for_pair["species_1"] = pair[0]
			zull_Z_result_for_pair["species_2"] = pair[1]
			null_Z_results_for_pairs.append(zull_Z_result_for_pair)
## Make a DF for the two-way results and write to file
		
		pd.DataFrame(null_Z_results_for_pairs).to_csv(args.outputPrefix + ".nullZ.csv", index = False)


if __name__ == "__main__":
	main()
