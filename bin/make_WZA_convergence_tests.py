## A script to make datasets that contain information from multiple datasets
import argparse
import pandas as pd
import numpy as np


def shuffle_DF( input_DF ):
	
	return input_DF.set_index("rep").loc[
		np.random.permutation( input_DF.rep.unique() )
		].reset_index()

def main():

## Define command line args
	parser = argparse.ArgumentParser(description="This takes resultas of the WZA sims and ")

	parser.add_argument("--input", "-i",
			required = True,
			dest = "input",
			type = str, 
			nargs = "+",
			help = "CSV files containing the WZA results.")


	parser.add_argument("--output", "-o",
			required = True,
			dest = "output",
			type = str, 
			help = "The name you want to give to the output file")



	args = parser.parse_args()

## Use the three WZA files for the three maps
	WZA_1 = pd.read_csv( args.input[0] )
	WZA_2 = pd.read_csv( args.input[1] )
	WZA_3 = pd.read_csv( args.input[2] )

	reps_not_chimed = True

	while reps_not_chimed:
		shuffled_1 = shuffle_DF( WZA_1.copy() ) 
		shuffled_2 = shuffle_DF( WZA_1.copy() ) 
		shuffled_3 = shuffle_DF( WZA_2.copy() ) 
		shuffled_4 = shuffle_DF( WZA_2.copy() ) 
		shuffled_5 = shuffle_DF( WZA_3.copy() ) 

		BC_chimes = (shuffled_1.rep.unique() == shuffled_2.rep.unique()).sum()
		cline_chimes = (shuffled_3.rep.unique() == shuffled_4.rep.unique()).sum() 

#		print( BC_chimes, cline_chimes )
		
		if BC_chimes == 0 and cline_chimes ==0:
			reps_not_chimed = False

	thresh = 0.005
	 
	temp_DF_list =  [ shuffled_1.LA > thresh,
					shuffled_2.LA > thresh,
					shuffled_3.LA > thresh,
					shuffled_4.LA > thresh,
					shuffled_5.LA > thresh ]

	temp_header = ["Z_1",
				"Z_2",
				"Z_3",
				"Z_4",
				"Z_5"]

	temp_DF = pd.concat( temp_DF_list, axis = 1, keys = temp_header )

	LA_indicator = temp_DF.sum( axis = 1 )

	print( temp_DF )

	out_DF_list =  [ shuffled_1.empR_Z,
					shuffled_2.empR_Z,
					shuffled_3.empR_Z,
					shuffled_4.empR_Z,
					shuffled_5.empR_Z,
					LA_indicator ]

	out_header = ["Z_1",
				"Z_2",
				"Z_3",
				"Z_4",
				"Z_5",
				"converge"]

	out_DF = pd.concat( out_DF_list, axis = 1, keys = out_header )

	out_DF.to_csv(args.output, index = False)






if __name__ == "__main__":
	main()
