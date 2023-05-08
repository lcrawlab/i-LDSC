from pysnptools.snpreader import Bed
from pysnptools.standardizer import Identity

class SNPData(object):

    def __init__(self, chromosomes, n_samples, block_size=5000):
        self.chromosomes = chromosomes
        self.n_chromosomes = len(self.chromosomes)
        self.n_samples = n_samples
        self.n_snps = 0
        self.curr_loc = -1
        self.block_size = block_size
        self.chrom_pos_end = []

        # UKB data for generating synthetic phenotypes
        self.snp_dir = "/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/data/"

        # 1000G data for estimating LD scores using reference panel
        #self.snp_dir = "/users/gdarnel1/data/gdarnel1/1kg_data/1000G_EUR_Phase3_plink/"

        for chrom_num in chromosomes:

            # actual UKB individual-level data for generating simulation data
            snpreader = Bed(f"{self.snp_dir}/UK_ALL_REMAINING_{chrom_num}.bed",count_A1=True)
            #snpreader = Bed(f"{self.snp_dir}/ukb_cal_chr{chrom_num}_v2.bed",count_A1=True)
            #snpreader = Bed(f"{self.snp_dir}/chr{chrom_num}.bed",count_A1=True)
            # 1000G data subsetted to UKB SNPs for estimating MELD in simulations
            #snpreader = Bed(f"{self.snp_dir}/1000G.UKB.SNPS.{chrom_num}.bed",count_A1=True)
            # full 1000G data for estimating MELD in real data
            #snpreader = Bed(f"{self.snp_dir}/1000G.EUR.QC.{chrom_num}.bed",count_A1=True)

            n_snps_chrom = snpreader.shape[1]
            self.n_snps += n_snps_chrom
            self.chrom_pos_end.append(n_snps_chrom) # account for 0-based indexing

        self.curr_chrom_idx = 0
        self.get_next_chrom()

    def get_next_chrom(self):
        """
        side-effects only
        """

        # UKB data
        self.snpreader = Bed(f"{self.snp_dir}/UK_ALL_REMAINING_{self.chromosomes[self.curr_chrom_idx]}.bed",count_A1=True)
        #self.snpreader = Bed(f"{self.snp_dir}/ukb_cal_chr{self.chromosomes[self.curr_chrom_idx]}_v2.bed",count_A1=True)
        #self.snpreader = Bed(f"{self.snp_dir}/chr{self.chromosomes[self.curr_chrom_idx]}.bed",count_A1=True)
        # 1000G data subsetted to UKB SNPs
        #self.snpreader = Bed(f"{self.snp_dir}/1000G.UKB.SNPS.{self.chromosomes[self.curr_chrom_idx]}.bed",count_A1=True)
        # full 1000G data
        #self.snpreader = Bed(f"{self.snp_dir}/1000G.EUR.QC.{self.chromosomes[self.curr_chrom_idx]}.bed",count_A1=True)

    def get_next_snps(self):
        # no more chromosomes to load
        if self.curr_chrom_idx >= self.n_chromosomes:
            return None, None, None

        self.curr_loc += 1
        block_end = self.curr_loc + self.block_size
        block_spillover = False
        # next block would end beyond chromosome boundary
        if block_end > self.chrom_pos_end[self.curr_chrom_idx]:
            block_end = self.chrom_pos_end[self.curr_chrom_idx]
            self.curr_chrom_idx += 1
            block_spillover = True

        snpreader = self.snpreader[:self.n_samples,self.curr_loc:block_end].read().standardize()#standardize(Identity())
        self.curr_loc = block_end-1

        if block_spillover and self.curr_chrom_idx < self.n_chromosomes:
            self.curr_loc = -1
            self.get_next_chrom()

        #pos = snpreader.pos[:,1] # pos information in cm
        pos = snpreader.pos[:,2] # pos information in kb
        sid = snpreader.sid[:]
        return snpreader.read().val, pos, sid

