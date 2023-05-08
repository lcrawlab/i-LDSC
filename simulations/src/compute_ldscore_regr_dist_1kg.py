import os
import sys
import pickle
import pandas as pd

def generate_annotations(chromosomes=[1]):
    for chrom in chromosomes:
        # replace last two columns with 1
        cmd = 'sed "s/[A-Z]\\t[A-Z]$/1\\t1/"'
        cmd += f' {DIR}{chrom}.bim'
        cmd += f' > {DIR}{chrom}.annot'
        os.system(cmd)
        # add column headers
        cmd = 'sed -i \'1iCHR\\tSNP\\tCM\\tBP\\tLD\\tML\''
        cmd += f' {DIR}{chrom}.annot'
        os.system(cmd)

def generate_ldscores(chromosomes=[1]):
    for chrom in chromosomes:
        # dev for MELD
        cmd = './ldsc-master-dev/ldsc.py --l2'
        cmd += f' --bfile {DIR}{chrom}'
        cmd += f' --ld-wind-cm {LDWIND}'
        cmd += f' --out {LDOUT}chr_{chrom}'
        cmd += f' --annot {DIR}{chrom}.annot'
        os.system("source /users/gdarnel1/py2/bin/activate && " + cmd)

        # with no annotation for ldsc weights
        cmd = './ldsc-master/ldsc.py --l2' 
        cmd += f' --bfile {DIR}{chrom}'
        cmd += f' --ld-wind-cm {LDWIND}'
        cmd += f' --out {LDOUT}w_chr_{chrom}' 
        # source python 2 environment for ldsc
        os.system("source /users/gdarnel1/py2/bin/activate && " + cmd)

def get_pcs(filename):

    simids = pd.read_csv(f'/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/sim_cache/{EXPPATH}/{filename}.txt', sep = '\t')
    simids = simids[['fid']].rename(columns = {'fid':'FID'})

    allindivs = pd.read_csv('/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/data/UK_ALL_REMAINING.pcs', delim_whitespace = True)

    mergepcs = allindivs.merge(simids, on = 'FID', how = 'inner')
    mergepcs = mergepcs[['FID', 'IID'] + ['PC' + str(i) for i in range(1,11)]]

    mergepcs = mergepcs.rename(columns = {'FID':'fid', 'IID':'iid'})
    mergepcs.to_csv(f'/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/sim_cache/{EXPPATH}/{filename}.qcovar', sep = '\t', index = False, header = False)


def compute_sumstats(filename_prefix):
    for filename in filename_prefix:
        print(filename)
        # check if sumstats have already been created
        if not os.path.exists(f'{RESOUT}{filename}.fastGWA'):
        #if True:
            cmd = '( cd /users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/gcta_1.93.2beta/ && ./gcta64'
            cmd += f' --mbfile {UKBData}chr_list.txt'
            cmd += f' --pheno {PHENOCACHE}{filename}.txt'
            cmd += f' --out {RESOUT}{filename}'
            cmd += ' --fastGWA-lr)'

            os.system(cmd)

        if EXPPATH == 'GxAncestry' or EXPPATH == 'GxAncestry_noCASS':
            get_pcs(filename)
            if not os.path.exists(f'{RESOUT}{filename}.pc.corrected.fastGWA'):
                cmd = '( cd /users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/gcta_1.93.2beta/ && ./gcta64'
                cmd += f' --mbfile {UKBData}chr_list.txt'
                cmd += f' --pheno {PHENOCACHE}{filename}.txt'
                cmd += f' --qcovar {PHENOCACHE}{filename}.qcovar'
                cmd += f' --out {RESOUT}{filename}.pc.corrected'
                cmd += ' --fastGWA-lr)'

                os.system(cmd)

            


def munge_sumstats(filename_prefix):
    for filename in filename_prefix:
        print(filename)
        # check if munging has already happened
        
    #if True:
        cmd = 'python /users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/ldsc-master/munge_sumstats.py'
        cmd += f' --sumstats {RESOUT}{filename}.fastGWA'
        #cmd += f' --frq AF1'
        #cmd += f' --maf-min 0.05'
        cmd += f' --out {RESOUT}{filename}'
        os.system("source activate ldsc && " + cmd)

        if EXPPATH == 'GxAncestry' or EXPPATH == 'GxAncestry_noCASS':
            if not os.path.exists(f'{RESOUT}{filename}.pc.corrected.sumstats.gz'):
                cmd = 'python /users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/ldsc-master/munge_sumstats.py'
                cmd += f' --sumstats {RESOUT}{filename}.pc.corrected.fastGWA'
                #cmd += f' --frq AF1'
                #cmd += f' --maf-min 0.05'
                cmd += f' --out {RESOUT}{filename}.pc.corrected'
                os.system("source activate ldsc && " + cmd)

def run_ldsr(filename_prefix, exp_alpha):
    """
    win_size specified in kb
    """
    # run with MELD and LD scores
    print(EXPPATH)
    print(len(filename_prefix))
    
    for filename in filename_prefix:

        cmd = 'python /users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/ldsc-master/ldsc.py'
        cmd += f' --h2 {RESOUT}{filename}.sumstats.gz'
        #cmd += f' --h2 {RESOUT}{filename}_maf.sumstats.gz'
        
        # window size will determine MELD scores
        if exp_alpha == 1.0:
            cmd += f' --ref-ld-chr {LDOUT}/unild_meld_scores/merged.'
            cmd += f' --w-ld-chr {LDOUT}/weights_hm3_no_hla/weights.'

        elif exp_alpha == 0.5:
            cmd += f' --ref-ld-chr {LDOUT}/unild_meld_scores/merged.'
            cmd += f' --w-ld-chr {LDOUT}/weights_hm3_no_hla/weights.'

        else:
            cmd += f' --ref-ld-chr {LDOUT}/unild_meld_scores/merged.'
            cmd += f' --w-ld-chr {LDOUT}/weights_hm3_no_hla/weights.'


        cmd += f' --out {RESOUT}{filename}.meld.ldsr'
        cmd += f' --frqfile-chr {LDOUT}1000G_Phase3_frq/1000G.EUR.QC.'
        cmd += ' --chisq-max 200'
        cmd += ' --print-coefficients'
        cmd += ' --overlap-annot'
        cmd += ' --not-M-5-50'
        print(cmd)
        os.system("source activate ldsc && " + cmd)

        if EXPPATH == 'GxAncestry' or EXPPATH == 'GxAncestry_noCASS':
            cmd = 'python /users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/ldsc-master/ldsc.py'
            cmd += f' --h2 {RESOUT}{filename}.pc.corrected.sumstats.gz'
            #cmd += f' --h2 {RESOUT}{filename}_maf.sumstats.gz'
            
            # window size will determine MELD scores
            if exp_alpha == 1.0:
                cmd += f' --ref-ld-chr {LDOUT}/unild_meld_scores/merged.'
                cmd += f' --w-ld-chr {LDOUT}/weights_hm3_no_hla/weights.'

            elif exp_alpha == 0.5:
                cmd += f' --ref-ld-chr {LDOUT}/unild_meld_scores/merged.'
                cmd += f' --w-ld-chr {LDOUT}/weights_hm3_no_hla/weights.'

            else:
                cmd += f' --ref-ld-chr {LDOUT}/unild_meld_scores/merged.'
                cmd += f' --w-ld-chr {LDOUT}/weights_hm3_no_hla/weights.'


            cmd += f' --out {RESOUT}{filename}.pc.corrected.meld.ldsr'
            cmd += f' --frqfile-chr {LDOUT}1000G_Phase3_frq/1000G.EUR.QC.'
            cmd += ' --chisq-max 200'
            cmd += ' --print-coefficients'
            cmd += ' --overlap-annot'
            cmd += ' --not-M-5-50'
            print(cmd)
            os.system("source activate ldsc && " + cmd)

    # run again for regular LDScore (additive only)
    for filename in filename_prefix:

        cmd = 'python /users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/ldsc-master/ldsc.py'
        cmd += f' --h2 {RESOUT}{filename}.sumstats.gz'
        
        cmd += f' --ref-ld-chr {LDOUT}/unild_trick_sldsc/merged.'
        cmd += f' --w-ld-chr {LDOUT}/weights_hm3_no_hla/weights.'
        cmd += ' --not-M-5-50'
        cmd += f' --frqfile-chr {LDOUT}1000G_Phase3_frq/1000G.EUR.QC.'
        cmd += f' --out {RESOUT}{filename}.add.ldsr'
        cmd += ' --chisq-max 200'
        cmd += ' --print-coefficients'
        os.system("source activate ldsc && " + cmd)

        if EXPPATH == 'GxAncestry' or EXPPATH == 'GxAncestry_noCASS':

            cmd = 'python /users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/ldsc-master/ldsc.py'
            cmd += f' --h2 {RESOUT}{filename}.pc.corrected.sumstats.gz'
            
            cmd += f' --ref-ld-chr {LDOUT}/unild_trick_sldsc/merged.'
            cmd += f' --w-ld-chr {LDOUT}/weights_hm3_no_hla/weights.'
            cmd += ' --not-M-5-50'
            cmd += f' --frqfile-chr {LDOUT}1000G_Phase3_frq/1000G.EUR.QC.'
            cmd += f' --out {RESOUT}{filename}.pc.corrected.add.ldsr'
            cmd += ' --chisq-max 200'
            cmd += ' --print-coefficients'
            os.system("source activate ldsc && " + cmd)


def process_ldsr_results(filename_prefix,):


    for filename in filename_prefix:
      print(filename)
      meldlog = [i.strip() for i in open(f'{RESOUT}{filename}.meld.ldsr.log', 'r')]
      print(meldlog)
      meldestimatesline = [i for i in meldlog if 'Total Observed scale h2:' in i][0].split(' ')
      
      totalmeldh2 = float(meldestimatesline[-2])
      

      # meld_partitioned = 
      ldsclog = [i.strip()for i in open(f'{RESOUT}{filename}.add.ldsr.log', 'r')]
      ldscestimatesline = [i for i in ldsclog if 'Total Observed scale h2:' in i][0].split(' ')
      print(ldscestimatesline)
      ldscaddh2 = float(ldscestimatesline[-2])
      

      meldh2 = totalmeldh2 - ldscaddh2
      outfile = open(f'{RESOUT}{filename}.estimates.txt', 'w')
      outfile.write('\t'.join(['total_meld', 'ldsc_add', 'meld_score']) + '\n')
      outfile.write('\t'.join([str(k) for k in [totalmeldh2,ldscaddh2, meldh2]]))

      if EXPPATH == 'GxAncestry' or EXPPATH == 'GxAncestry_noCASS':
          print(filename)
          meldlog = [i.strip() for i in open(f'{RESOUT}{filename}.pc.corrected.meld.ldsr.log', 'r')]
          print(meldlog)
          meldestimatesline = [i for i in meldlog if 'Total Observed scale h2:' in i][0].split(' ')
          
          totalmeldh2 = float(meldestimatesline[-2])
          

          # meld_partitioned = 
          ldsclog = [i.strip()for i in open(f'{RESOUT}{filename}.pc.corrected.add.ldsr.log', 'r')]
          ldscestimatesline = [i for i in ldsclog if 'Total Observed scale h2:' in i][0].split(' ')
          print(ldscestimatesline)
          ldscaddh2 = float(ldscestimatesline[-2])
          

          meldh2 = totalmeldh2 - ldscaddh2
          outfile = open(f'{RESOUT}{filename}.pc.corrected.estimates.txt', 'w')
          outfile.write('\t'.join(['total_meld', 'ldsc_add', 'meld_score']) + '\n')
          outfile.write('\t'.join([str(k) for k in [totalmeldh2,ldscaddh2, meldh2]]))



def run_pipeline(filename_prefix, exp_alpha):
    chr = list(range(1,23))

    # generate_annotations(chr)
    # generate_ldscores(chr)

    compute_sumstats(filename_prefix)
    munge_sumstats(filename_prefix)

    run_ldsr(filename_prefix, exp_alpha)
    process_ldsr_results(filename_prefix)

if __name__ == "__main__":
    global EXPPATH
    global DIR
    global UKBData
    global LDOUT
    global PHENOCACHE
    global RESOUT
    global LDWIND

    # will throw indexerror if no args given
    job_num = int(sys.argv[1])
    EXPPATH = str(sys.argv[2])
    exp_alpha = float(sys.argv[3])

    # global variable definitions
    DIR="/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/1kg_data/1000G_EUR_Phase3_plink/1000G.EUR.QC."
    UKBData="/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/data/"
    LDOUT="/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/1kg_data/ldscores/"
    PHENOCACHE=f'/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/sim_cache/{EXPPATH}/'
    RESOUT=f'/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/sim_results/{EXPPATH}/'
    LDWIND = 1000

    job2filename = pickle.load(open(f'{RESOUT}job2filename.{EXPPATH}.pkl','rb'))
    filename_prefix = job2filename[job_num]

    if len(filename_prefix) < 1:
        raise Exception('no files for this job')

    run_pipeline(filename_prefix, exp_alpha)

