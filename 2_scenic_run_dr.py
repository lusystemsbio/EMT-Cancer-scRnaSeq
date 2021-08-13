if __name__ == '__main__':
	import os
	import glob
	import pickle
	import pandas as pd
	import numpy as np
	import json
	import seaborn as sns
	import sys

	from dask.diagnostics import ProgressBar
	from arboreto.utils import load_tf_names
	from arboreto.algo import grnboost2
	from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
	from pyscenic.utils import modules_from_adjacencies, load_motifs
	from pyscenic.prune import prune2df, df2regulons
	from pyscenic.aucell import aucell


	inputFilename = sys.argv[1] # Name of the specific data file provided as input through the shell script file
	RESULT_FOLDER = "results"
	DATA_FOLDER = "data"
	RESOURCES_FOLDER = "resources"
	DATABASE_FOLDER = "databases"
	DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg19*.mc9nr.feather")
	MOTIF_ANNOTATIONS_FNAME = os.path.join(DATABASE_FOLDER, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
	MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'hs_hgnc_curated_tfs.txt')
	SC_EXP_FNAME = os.path.join(DATA_FOLDER, (inputFilename + ".tsv"))
	MODULES_FNAME = os.path.join(RESULT_FOLDER, ("modules_" + inputFilename + ".p"))
	REGULONS_FNAME = os.path.join(RESULT_FOLDER, "regulons_" + inputFilename + ".p")
	MOTIFS_FNAME = os.path.join(RESULT_FOLDER, "motifs_" + inputFilename + ".csv")

	ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
	print(ex_matrix.shape)
	tf_names = load_tf_names(MM_TFS_FNAME)
	print("tf names loaded")
	out_file = os.path.join(RESULT_FOLDER, "grn_output_" + inputFilename + ".tsv")
	db_fnames = glob.glob(DATABASES_GLOB)


	def name(fname):
		return os.path.splitext(os.path.basename(fname))[0]


	dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
	print(dbs)
	print("running grnboost")
	print("tf_names head")
	print(tf_names[1:5])
	#print("gene names head")
	#print(ex_matrix.iloc[1:5,1:5])
	adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
	adjacencies.head()
	print("identify modules")
	adjacencies.to_csv(out_file, sep='\t', index=False, header=False)
	print("grnboost done")
	modules = list(modules_from_adjacencies(adjacencies, ex_matrix, rho_mask_dropouts=True))

	#print("writing modules")
	#with open(MODULES_FNAME, 'wb') as f:
	#	pickle.dump(modules, f)

	print("Finding Enriched modules")
	# Calculate a list of enriched motifs and the corresponding target genes for all modules.

	with ProgressBar():
		df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
		df.head()

	# Create regulons from this table of enriched motifs.
	print("creating regulons")
	regulons = df2regulons(df)

	print("writing regulons")
	# Save the enriched motifs and the discovered regulons to disk.
	#df.to_csv(MOTIFS_FNAME)
	#with open(REGULONS_FNAME, "wb") as f:
	#	pickle.dump(regulons, f)

	print("Finding AUC of cells")
	auc_mtx = aucell(ex_matrix, regulons, num_workers=1)
	auc_file = os.path.join(RESULT_FOLDER, "AUC_" + inputFilename + ".csv")
	auc_mtx.to_csv(auc_file)

	# Write to JSON for import to R
	rdict = {} 
	for reg in regulons:
		targets = [ target for target in reg.gene2weight ]
		rdict[reg.name] = targets


	rjson = json.dumps(rdict)
	jsonFile = os.path.join(RESULT_FOLDER, "regulons_" + inputFilename + ".json")
	f = open(jsonFile,"w")
	f.write(rjson)
	f.close()

	# Plot the results for visualization
	#sns_plot = sns.clustermap(auc_mtx, figsize=(8,8))
	#plotFile = os.path.join(RESULT_FOLDER, "AUC_" + inputFilename + ".png")
	#sns_plot.savefig(plotFile)
