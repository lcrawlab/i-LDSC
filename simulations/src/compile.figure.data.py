import pandas as pd
import numpy as np 
import pickle
import sys


EXPPATH = str(sys.argv[1])


RESOUT=f'/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/sim_results/{EXPPATH}/'
SIMCACHE = f'/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/sim_cache/{EXPPATH}/'

job2filename = pickle.load(open(f'{RESOUT}job2filename.{EXPPATH}.pkl','rb'))

firstpheno =job2filename[0][0]
metadata = pickle.load(open(f'{SIMCACHE}{firstpheno}.pkl', 'rb'))
simdf = pd.DataFrame.from_dict(metadata, orient = 'index').T
simdf = simdf.rename(index = {0:firstpheno})

estimates = pd.read_csv(f'{RESOUT}{firstpheno}.estimates.txt', sep = '\t')
estimates = estimates.rename(index = {0:firstpheno})

simdf.loc[firstpheno,'total_meld'] = estimates.loc[firstpheno,'total_meld']
simdf.loc[firstpheno,'ldsc_add'] = estimates.loc[firstpheno,'ldsc_add']
simdf.loc[firstpheno,'meld_score'] = estimates.loc[firstpheno,'meld_score']

outdf = simdf

for filename in job2filename:
	if filename == 0:
		phenos = [i for i in job2filename[filename][1:]]
	else:
		phenos = [i for i in job2filename[filename]]
	for pheno in phenos:
		metadata = pickle.load(open(f'{SIMCACHE}{pheno}.pkl', 'rb'))
		simdf = pd.DataFrame.from_dict(metadata, orient = 'index').T
		simdf = simdf.rename(index = {0:pheno})

		estimates = pd.read_csv(f'{RESOUT}{pheno}.estimates.txt', sep = '\t')
		estimates = estimates.rename(index = {0:pheno})

		simdf.loc[pheno,'total_meld'] = estimates.loc[pheno,'total_meld']
		simdf.loc[pheno,'ldsc_add'] = estimates.loc[pheno,'ldsc_add']
		simdf.loc[pheno,'meld_score'] = estimates.loc[pheno,'meld_score']
		
		outdf = pd.concat([outdf,simdf])

for column in outdf.columns:
	if '_pve' in column:
		outdf[column] = outdf[column].astype(float)/100
print(outdf)
outdf.to_pickle(f'/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld-annotations-run/figures/figure_cache/{EXPPATH}.figure.data.pkl')


if EXPPATH == 'GxAncestry' or EXPPATH == 'GxAncestry_noCASS':

	job2filename = pickle.load(open(f'{RESOUT}job2filename.{EXPPATH}.pkl','rb'))

	firstpheno =job2filename[0][0]
	metadata = pickle.load(open(f'{SIMCACHE}{firstpheno}.pkl', 'rb'))
	simdf = pd.DataFrame.from_dict(metadata, orient = 'index').T
	simdf = simdf.rename(index = {0:firstpheno})

	estimates = pd.read_csv(f'{RESOUT}{firstpheno}.pc.corrected.estimates.txt', sep = '\t')
	estimates = estimates.rename(index = {0:firstpheno})

	simdf.loc[firstpheno,'total_meld'] = estimates.loc[firstpheno,'total_meld']
	simdf.loc[firstpheno,'ldsc_add'] = estimates.loc[firstpheno,'ldsc_add']
	simdf.loc[firstpheno,'meld_score'] = estimates.loc[firstpheno,'meld_score']

	outdf = simdf

	for filename in job2filename:
		if filename == 0:
			phenos = [i for i in job2filename[filename][1:]]
		else:
			phenos = [i for i in job2filename[filename]]
		for pheno in phenos:
			metadata = pickle.load(open(f'{SIMCACHE}{pheno}.pkl', 'rb'))
			simdf = pd.DataFrame.from_dict(metadata, orient = 'index').T
			simdf = simdf.rename(index = {0:pheno})

			estimates = pd.read_csv(f'{RESOUT}{pheno}.pc.corrected.estimates.txt', sep = '\t')
			estimates = estimates.rename(index = {0:pheno})

			simdf.loc[pheno,'total_meld'] = estimates.loc[pheno,'total_meld']
			simdf.loc[pheno,'ldsc_add'] = estimates.loc[pheno,'ldsc_add']
			simdf.loc[pheno,'meld_score'] = estimates.loc[pheno,'meld_score']
			
			outdf = pd.concat([outdf,simdf])

	for column in outdf.columns:
		if '_pve' in column:
			outdf[column] = outdf[column].astype(float)/100
	print(outdf)
	outdf.to_pickle(f'/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld-annotations-run/figures/figure_cache/{EXPPATH}.pc.corrected.figure.data.pkl')
