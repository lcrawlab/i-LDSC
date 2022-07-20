"""
Implements fast calculation of MELD scores
Gregory Darnell and Samuel Pattillo Smith
Extended from initial script written by Shadi Zabad, April 2020
"""

import numpy as np
import pandas as pd
import time
import sys
import os
import errno
import argparse
from pandas_plink import read_plink1_bin
from subprocess import check_call
import csv
from numba import njit, prange
from multiprocessing import Pool
from itertools import product


def makedir(cdir):
    try:
        os.makedirs(cdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def read_plink_files(input_fname, chr_num):

    # Read/transform genotype matrices:
    try:
        gt_ac = read_plink1_bin(input_fname % chr_num + ".bed")
    except Exception as e:
        raise e

    gt_ac = np.abs(gt_ac.values - 2).astype(np.int64)
    ngt_ac = (gt_ac - gt_ac.mean(axis=0)) / gt_ac.std(axis=0)

    # Read the .bim file:
    try:
        gt_meta = pd.read_csv(input_fname % chr_num + ".bim",
                              names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], sep='\t')
    except Exception as e:
        raise e

    maf = gt_ac.sum(axis=0) / (2. * gt_ac.shape[0])
    maf = np.round(np.where(maf > .5, 1. - maf, maf), float_precision)
    gt_meta['MAF'] = maf

    gt_meta = gt_meta[['CHR', 'SNP', 'CM', 'BP', 'MAF']]

    if weights:
        sel_snp_idx = np.where(gt_meta['SNP'].isin(snp_list.squeeze()))[0]

        fgt_meta = gt_meta.iloc[sel_snp_idx, ].reset_index(drop=True)
        fgt_ac = gt_ac[:, sel_snp_idx]
        fngt_ac = ngt_ac[:, sel_snp_idx]

        return fgt_ac, fngt_ac, fgt_meta
    else:
        return gt_ac, ngt_ac, gt_meta


# --------------- Auxiliary Functions ---------------


@njit(parallel=True)
def numba_count(a, out, m, n):
    for i in prange(m):
        for j in prange(n):
            out[a[i, j], i] += 1


@njit(parallel=True)
def bincount2D_numba(a, bin_num=9):

    m, n = a.shape
    out = np.zeros((bin_num, m), dtype=np.int_)

    numba_count(a, out, m, n)

    return out

@njit
def d_squared_unphased(counts, n):
    """
    Implementation by Aaron Ragsdale
    """

    n1 = counts[0, :]
    n2 = counts[1, :]
    n3 = counts[2, :]
    n4 = counts[3, :]
    n5 = counts[4, :]
    n6 = counts[5, :]
    n7 = counts[6, :]
    n8 = counts[7, :]
    n9 = counts[8, :]

    numer = ((n2 * n4 - n2 ** 2 * n4 + 4 * n3 * n4 - 4 * n2 * n3 * n4 - 4 * n3 ** 2 * n4 - n2 * n4 ** 2 -
              4 * n3 * n4 ** 2 + n1 * n5 - n1 ** 2 * n5 + n3 * n5 + 2 * n1 * n3 * n5 - n3 ** 2 * n5 -
              4 * n3 * n4 * n5 - n1 * n5 ** 2 - n3 * n5 ** 2 + 4 * n1 * n6 - 4 * n1 ** 2 * n6 + n2 * n6 -
              4 * n1 * n2 * n6 - n2 ** 2 * n6 + 2 * n2 * n4 * n6 - 4 * n1 * n5 * n6 - 4 * n1 * n6 ** 2 - n2 * n6 ** 2 +
              4 * n2 * n7 - 4 * n2 ** 2 * n7 + 16 * n3 * n7 - 16 * n2 * n3 * n7 - 16 * n3 ** 2 * n7 -
              4 * n2 * n4 * n7 - 16 * n3 * n4 * n7 + n5 * n7 + 2 * n1 * n5 * n7 -
              4 * n2 * n5 * n7 - 18 * n3 * n5 * n7 - n5 ** 2 * n7 + 4 * n6 * n7 + 8 * n1 * n6 * n7 - 16 * n3 * n6 * n7 -
              4 * n5 * n6 * n7 - 4 * n6 ** 2 * n7 - 4 * n2 * n7 ** 2 - 16 * n3 * n7 ** 2 - n5 * n7 ** 2 -
              4 * n6 * n7 ** 2 + 4 * n1 * n8 - 4 * n1 ** 2 * n8 + 4 * n3 * n8 + 8 * n1 * n3 * n8 -
              4 * n3 ** 2 * n8 + n4 * n8 - 4 * n1 * n4 * n8 + 2 * n2 * n4 * n8 - n4 ** 2 * n8 -
              4 * n1 * n5 * n8 - 4 * n3 * n5 * n8 + n6 * n8 + 2 * n2 * n6 * n8 - 4 * n3 * n6 * n8 +
              2 * n4 * n6 * n8 - n6 ** 2 * n8 - 16 * n3 * n7 * n8 - 4 * n6 * n7 * n8 - 4 * n1 * n8 ** 2 -
              4 * n3 * n8 ** 2 - n4 * n8 ** 2 - n6 * n8 ** 2 + 16 * n1 * n9 - 16 * n1 ** 2 * n9 +
              4 * n2 * n9 - 16 * n1 * n2 * n9 - 4 * n2 ** 2 * n9 + 4 * n4 * n9 - 16 * n1 * n4 * n9 + 8 * n3 * n4 * n9 -
              4 * n4 ** 2 * n9 + n5 * n9 - 18 * n1 * n5 * n9 - 4 * n2 * n5 * n9 + 2 * n3 * n5 * n9 -
              4 * n4 * n5 * n9 - n5 ** 2 * n9 - 16 * n1 * n6 * n9 -
              4 * n2 * n6 * n9 + 8 * n2 * n7 * n9 + 2 * n5 * n7 * n9 - 16 * n1 * n8 * n9 - 4 * n4 * n8 * n9 -
              16 * n1 * n9 ** 2 - 4 * n2 * n9 ** 2 -
              4 * n4 * n9 ** 2 - n5 * n9 ** 2) / 16. +
             (-((n2 / 2. + n3 + n5 / 4. + n6 / 2.) * (n4 / 2. + n5 / 4. + n7 + n8 / 2.)) +
             (n1 + n2 / 2. + n4 / 2. + n5 / 4.) * (n5 / 4. + n6 / 2. + n8 / 2. + n9)) ** 2)

    return 4. * numer / (n * (n - 1) * (n - 2) * (n - 3))


# --------------------------------------------------
# --------------- LD Score Functions ---------------
# --------------------------------------------------

def compute_modified_ld_score(j, max_cm_dist=1.):

    # Obtain neighboring SNPs information:
    # --------------------------------------------
    # Condition to exclude focal snp: (gt_meta.index != gt_meta.iloc[j, ].name) &
    neighb_snps = gt_meta.loc[(np.abs(gt_meta['CM'] - gt_meta.iloc[j, ]['CM']) <= max_cm_dist), ]

    neighb_snps_annot = neighb_snps.iloc[:, annot_start_idx:].values
    neighb_snps_idx = neighb_snps.index.values
    var_xk = neighb_snps['VAR'].values
    var_xj = gt_meta.iloc[j, ]['VAR']

    # --------------------------------------------
    # Compute D^2
    #gt_counts = gt_ac[:, j, np.newaxis] * 3 + gt_ac[:, neighb_snps_idx]
    #count_mat = bincount2D_numba(gt_counts.T)

    ## D^2 vector with all neighboring SNPs:
    #D2 = d_squared_unphased(count_mat[::-1, :], N)
    #D2 = (4. / var_xj) * D2

    # --------------------------------------------
    # Compute r^2
    uncr_r2 = (np.dot(ngt_ac[:, j], ngt_ac[:, neighb_snps_idx]) / N)**2
    r2 = uncr_r2 - (1. - uncr_r2)/(N - 2)

    # --------------------------------------------
    # Compute scores based on different estimators/assumptions:

    # = = = = = = D^2 based estimators = = = = = =

    scores = [] # list of numpy arrays of shape (1 x n_annot)

    for lds in scores_to_compute.values():

        if lds['estimator'] == 'D2':
            scores.append(
                np.dot((neighb_snps_annot * (var_xk.reshape(-1, 1)**(-lds['alpha']))).T,
                       D2)
            )
        elif lds['estimator'] == 'R2':
            scores.append(
                np.dot((neighb_snps_annot * (var_xk.reshape(-1, 1) ** (1. - lds['alpha']))).T,
                       r2)
            )
        elif lds['estimator'] == 'NR2':
            scores.append(
                np.dot((neighb_snps_annot * (var_xk.reshape(-1, 1) ** (1. - lds['alpha']))).T,
                       uncr_r2)
            )
        else:
            raise Exception(f"LD estimator {lds['estimator']} not implemented!")

    return j, scores

def compute_meld_score(j, win_size, max_cm_dist=1.):

    # Obtain neighboring SNPs information:
    # --------------------------------------------
    # Condition to exclude focal snp: (gt_meta.index != gt_meta.iloc[j, ].name) &
    neighb_snps = gt_meta.loc[(np.abs(gt_meta['CM'] - gt_meta.iloc[j, ]['CM']) <= max_cm_dist), ]
    neighb_snps = neighb_snps.drop_duplicates(subset = ['SNP'],keep = 'first')
    neighb_snps = neighb_snps.reset_index(drop = True)
    # condition to exclude focal snp below
    #neighb_snps = gt_meta.loc[(np.abs(gt_meta['CM'] - gt_meta.iloc[j, ]['CM']) <= max_cm_dist) \
    #                          & (gt_meta.index != gt_meta.iloc[j,].name), ]

    neighb_snps_annot = neighb_snps.iloc[:, annot_start_idx:].values
    # reduce annotation by 1 column for MELD
    neighb_snps_annot = neighb_snps_annot[:,0][:,None]

    neighb_snps_idx = neighb_snps.index.values

    var_xk = neighb_snps['VAR'].values
    var_xj = gt_meta.iloc[j, ]['VAR']

    # --------------------------------------------
    # Compute D^2
    #gt_counts = gt_ac[:, j, np.newaxis] * 3 + gt_ac[:, neighb_snps_idx]
    #count_mat = bincount2D_numba(gt_counts.T)

    ## D^2 vector with all neighboring SNPs:
    #D2 = d_squared_unphased(count_mat[::-1, :], N)
    #D2 = (4. / var_xj) * D2

    # --------------------------------------------
    # Compute r^2
    uncr_r2 = (np.dot(ngt_ac[:, j], ngt_ac[:, neighb_snps_idx]) / N)**2
    r2 = uncr_r2 - (1. - uncr_r2)/(N - 2)

    # --------------------------------------------
    # Compute scores based on different estimators/assumptions:

    # = = = = = = D^2 based estimators = = = = = =

    scores = [] # list of numpy arrays of shape (1 x n_annot)

    # COMPUTE MELD SCORES
    meld_half_win = int(win_size/2)

    min_win = max(0,j - meld_half_win)
    max_win = min(len(gt_meta)-1,j + meld_half_win)

    neighb_snps_meld = gt_meta.iloc[min_win:max_win+1, ]
    # condition to exclude focal snp below
    #neighb_snps_meld = neighb_snps_meld.loc[neighb_snps_meld.index != gt_meta.iloc[j,].name]

    neighb_snps_idx_meld = neighb_snps_meld.index.values
    var_xk_meld = neighb_snps_meld['VAR'].values

    # create W matrix
    uncr_r2_meld = (np.dot(ngt_ac[:, j]**2,ngt_ac[:, neighb_snps_idx_meld]) / N)**2
    r2_meld = uncr_r2_meld - (1. - uncr_r2_meld)/(N - 2)
    r2_meld *= var_xk_meld

    meld_score = sum(r2_meld)
        
    for lds in scores_to_compute.values():

        if lds['estimator'] == 'D2':
            scores.append(
                np.dot((neighb_snps_annot * (var_xk.reshape(-1, 1)**(-lds['alpha']))).T,
                       D2)
            )
        elif lds['estimator'] == 'R2':
            scores.append(
                np.dot((neighb_snps_annot * (var_xk.reshape(-1, 1) ** (1. - lds['alpha']))).T,
                       r2)
            )
        elif lds['estimator'] == 'NR2':
            scores.append(
                np.dot((neighb_snps_annot * (var_xk.reshape(-1, 1) ** (1. - lds['alpha']))).T,
                       uncr_r2)
            )
        else:
            raise Exception(f"LD estimator {lds['estimator']} not implemented!")
        
        scores[-1] = np.append(scores[-1],meld_score)

    return j, scores

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='LD Score Regression Using 1000 Genomes Project Data')
    parser.add_argument('--pop', dest='pop', type=str, default='EUR',
                        help='The population name')
    parser.add_argument('--weights', dest='weights', action='store_true',
                        help='Calculate the weights for the LDSC')
    parser.add_argument('--chrom', dest='chrom', type=int)
    parser.add_argument('--win', dest='win', type=int)
    parser.add_argument('--alpha',dest = 'alpha',type = float)
    parser.add_argument('--plink_dir',dest = 'plink_dir',type = str)
    parser.add_argument('--out',dest = 'out',type = str)
    parser.add_argument('--filter_list',dest = 'list',type = str)
    parser.add_argument('--annotations',dest = 'annotations',type = str)

    args = parser.parse_args()
    chrom = args.chrom
    win = args.win

    # Global parameters
    # ---------------------------------------------------------
    dist_measure = "cM"
    annot_start_idx = 6
    weights = args.weights

    ld_estimator = ['R2']
    alpha = [args.alpha]

    scores_to_compute = {
        lde + '_' + str(a): {
            'estimator': lde,
            'alpha': a
        }
        for lde in ld_estimator for a in alpha
    }

    population = args.pop

    # = = = = = = = = = =
    # Computational configurations:
    num_proc = 4
    float_precision = 15
    os.environ["OMP_NUM_THREADS"] = "2"
    os.environ["NUMBA_NUM_THREADS"] = "2"

    # = = = = = = = = = =
    # Input:
    plink_dir = args.plink_dir + '%s'
    w_snp_filter = args.list
    annotations = args.annotations + '.%d.annot'

    output_dir = args.out
    output_dirs = [os.path.join(output_dir, sn) for sn in scores_to_compute]
    [makedir(od) for od in output_dirs]

    # = = = = = = = = = =
    # Read the snp filter file:
    snp_filter = False
    if not w_snp_filter == None:
        snp_list = pd.read_csv(w_snp_filter, sep="\t")['SNP'].values
        snp_filter = True

    # ---------------------------------------------------------

    for chr_num in range(chrom, chrom+1):

        output_files = [os.path.join(od, f"MELD.win_{win}.{str(chr_num)}.l2.ldscore") for od in output_dirs]

        print("Processing chromosome %s..." % str(chr_num))

        # Read the genotype file:
        try:
            gt_ac, ngt_ac, gt_meta = read_plink_files(plink_dir, str(chr_num))
        except Exception as e:
            continue

        N, M = gt_ac.shape

        gt_meta['VAR'] = 2.*gt_meta['MAF']*(1. - gt_meta['MAF'])

        # Read the annotations file:
        if weights:
            gt_meta['base'] = 1.
        else:
            try:
                annot_df = pd.read_csv(annotations % chr_num, sep="\s+").drop(['CHR', 'BP', 'CM'], axis=1)
                gt_meta = pd.merge(gt_meta, annot_df, on='SNP')
            except Exception as exp:
                gt_meta['base'] = 1.

        output_colnames = [[cn + sn for cn in gt_meta.columns[annot_start_idx:]]
                           for sn in scores_to_compute]

        print("Computing LD Scores...")

        if not weights:
            alphas = [scores_to_compute[sn]['alpha'] for sn in scores_to_compute]
            for idx,of in enumerate(output_files):
                alpha = 1.0 - alphas[idx]
                # when using an \alpha parameter, need to correct .M files
                # with total number of SNPs corrected by variance^alpha
                # in order for LDSC heritabilities to be computed correctly
                M_tot = gt_meta.iloc[:,annot_start_idx:].multiply(gt_meta["VAR"],axis="index")
                M_tot = M_tot**alpha
                M_tot = M_tot.sum(axis=0).values
                M_5_50 = gt_meta.loc[gt_meta['MAF'] >= .05, ].iloc[:, annot_start_idx:].multiply(gt_meta["VAR"],axis="index")
                M_5_50 = M_5_50**alpha
                M_5_50 = M_5_50.sum(axis=0).values

                np.savetxt(of.replace('.ldscore', '.M'), M_tot.reshape(1, -1), delimiter="\t", fmt='%.1f')
                np.savetxt(of.replace('.ldscore', '.M_5_50'), M_5_50.reshape(1, -1), delimiter="\t", fmt='%.1f')

        start = time.time()
        pool = Pool(num_proc)

        open_files = [open(outf, 'w') for outf in output_files]
        csv_writers = [csv.writer(outf, delimiter='\t') for outf in open_files]

        # Write the column names:
        for cw, col in zip(csv_writers, output_colnames):
            cw.writerow(list(gt_meta.columns[:annot_start_idx - 1]) + col)

        # Select the subset of snps to compute the the LD scores for:
        if snp_filter:
            snps_to_process = list(np.where(gt_meta['SNP'].isin(snp_list.squeeze()))[0])
        else:
            snps_to_process = gt_meta.index

        # Compute the LD Scores:
        # last parameter (1) is the distance for ld score neighbors (in cM)
        for idx, (snp_idx, ld_scores) \
                in enumerate(pool.starmap(compute_meld_score, product(snps_to_process, [win], [1])), 1):

            for cw, lds in zip(csv_writers, ld_scores):
                cw.writerow(list(gt_meta.iloc[snp_idx, :annot_start_idx - 1]) +
                            list(np.round(lds, float_precision)))

            if idx % 1000 == 0:
                print("Computed LD Score for %d variants" % idx)
                sys.stdout.flush()

        [outf.close() for outf in open_files]

        pool.close()
        pool.join()

        end = time.time()
        print("Processing Time:", end - start)

        # Gzip the output file
        [check_call(['gzip', '-f', of]) for of in output_files]
