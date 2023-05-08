import numpy as np
from sklearn.preprocessing import scale
from SNPData import SNPData
import pandas as pd
import sys

class GxEArchitecture(object):

    def __init__(self, n_samples, block_size, pve, rho, phi, win_size, sparsity, phi_sparsity, scalar, exp_alpha, null_sim = False):
        print("pve {}, rho {}, phi {}, simtype: {}".format(pve, rho, phi, 'GxE'))
        ## broad sense heritability
        self.pve = pve
        
        # proportion of linear effects contributing to broad-sense heritability
        self.rho = rho
        chrom_range = np.arange(1,23)
        
        self.n_samples = n_samples
        self.block_size = block_size
        
        self.snpdata = SNPData(chrom_range, self.n_samples, self.block_size)

        self.null_sim = null_sim
        self.scalar = scalar
        #size of windows that the CASS score is constructed with
        self.window_size = win_size
        self.sparsity = sparsity
        
        #Initialize effect size vectors for each type of generative architectures
        self.y_marginal = np.zeros((self.n_samples))
        self.y_err = np.zeros((self.n_samples))
        self.y_epi = np.zeros((self.n_samples))
        self.y_gxe = np.zeros((self.n_samples))


        self.y = np.zeros((self.n_samples))
        # relationship between MAF-ES
        self.exp_alpha = exp_alpha
        
        # for testing how many SNPs make up the cis window
        self.w_shapes = []
        idx = 0

        next_snps = np.zeros((0))
        # counter = 1
        while type(next_snps) == np.ndarray:
            # print(counter)
            if not self.snpdata.curr_chrom_idx >= self.snpdata.n_chromosomes:
                curr_chrom = self.snpdata.chromosomes[self.snpdata.curr_chrom_idx]
                maf = pd.read_csv(f"/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/data/UK_ALL_REMAINING_{curr_chrom}.maf", delim_whitespace = True)

            next_snps, snp_pos, snp_id = self.snpdata.get_next_snps()
            # print(next_snps)
            if type(next_snps) == np.ndarray:
                self.n_snps = next_snps.shape[1]
                self.n_epistatic_hubs = int(self.n_snps*(sparsity/100))
                self.n_gxe_hubs = int(self.n_snps*(phi_sparsity/100))
 
                maf_vals = maf.set_index('SNP').loc[snp_id].MAF.values

                # nan to zero
                next_snps[np.isnan(next_snps)] = 0

                y_marginal, y_epi, y_gxe, y_err = self.generate_effects(rho, phi, pve, next_snps, snp_pos, maf_vals)
                
                self.y_marginal += y_marginal
                self.y_epi += y_epi
                self.y_err += y_err
                self.y_gxe += y_gxe
                idx += 1
            # counter += 1

        if 1-self.rho-self.phi > 0 and self.phi != 0:
            self.y_marginal = self.y_marginal*np.sqrt((self.pve*self.rho)/np.var(self.y_marginal))
            self.y_epi = self.y_epi*np.sqrt((self.pve*(1-self.rho-self.phi))/np.var(self.y_epi))
            self.y_gxe = self.y_gxe*np.sqrt((self.pve*(self.phi))/np.var(self.y_gxe))
            self.y_err=self.y_err*np.sqrt((1-self.pve)/np.var(self.y_err))
        
        elif self.phi == 0:
            self.y_marginal = self.y_marginal*np.sqrt((self.pve*self.rho)/np.var(self.y_marginal))
            self.y_epi = self.y_epi*np.sqrt((self.pve*(1-self.rho-self.phi))/np.var(self.y_epi))
            self.y_gxe = np.zeros((self.n_samples))
            self.y_err=self.y_err*np.sqrt((1-self.pve)/np.var(self.y_err))
        
        else:
            self.y_marginal = self.y_marginal*np.sqrt((self.pve*self.rho)/np.var(self.y_marginal))
            self.y_epi = np.zeros((self.n_samples))
            self.y_gxe = self.y_gxe*np.sqrt((self.pve*(self.phi))/np.var(self.y_gxe))
            self.y_err=self.y_err*np.sqrt((1-self.pve)/np.var(self.y_err))
        
        self.y = self.y_marginal + self.y_epi + self.y_gxe + self.y_err
        self.y = scale(self.y)
        
        self.linear_pve = np.var(self.y_marginal.flatten())/np.var(self.y.flatten())
        self.epistatic_pve = np.var(self.y_epi.flatten())/np.var(self.y.flatten())
        self.gxe_pve = np.var(self.y_gxe.flatten())/np.var(self.y.flatten())
        print("linear pve: ",self.linear_pve)
        print("epistatic pve: ",self.epistatic_pve)
        print("GxE pve: ", self.gxe_pve)
        print("error pve: ",np.var(self.y_err.flatten())/np.var(self.y.flatten()))

    def generate_effects(self, rho, phi, pve, X, snp_pos, maf): 
        #self.pve_per_hub = []
        self.rho = float(rho)
        self.phi = float(phi)
        self.pve = float(pve)

        # randomly choose causal epistatic snps
        n_snps_draw = [i for i in range(self.n_snps)]
        
        self.causal_snp_idx_epi = np.random.choice(n_snps_draw,size=self.n_epistatic_hubs,replace=False)
        self.causal_snp_idx_epi = self.causal_snp_idx_epi.tolist()
        
        self.res = [i for i in n_snps_draw if i not in self.causal_snp_idx_epi]
        
        # self.tempsnps = self.n_snps.remove[self.causal_snp_idx_epi]
        # randomly choose causal snps for gxe sources of architecture
        self.causal_snp_idx_gxe = np.random.choice(self.res ,size=self.n_gxe_hubs,replace=False)
        self.causal_snp_idx_gxe = self.causal_snp_idx_gxe.tolist()

        # not used now
        self.null_alt = np.array([int(i not in self.causal_snp_idx_epi) for i in range(self.n_snps)])


        # genetic architecture change -- effect sizes related to MAF
        self.beta = np.random.normal(0,np.sqrt((2*maf*(1.0-maf))**self.exp_alpha),size=self.n_snps)

        # print("Generating effects")
        y_marginal = np.dot(X,self.beta)
        # print("Linear effects generated")
       
        # renormalize additive effects
        y_marginal = y_marginal * np.sqrt(self.pve*(self.rho)/(np.var(y_marginal)))
        # renormalize betas
        #self.beta = self.beta * np.sqrt(self.pve*self.rho/(np.var(self.y_marginal)))
        #self.y_marginal = np.dot(X,self.beta)

        # null effect draw and re-"normalize" epistatic y to have the correct PVE
        self.theta = []
        y_epi = np.zeros((self.n_samples))

        if (1-self.rho-self.phi) != 0.0:
            # print("get marginal epistatic")
            y_epi = self.get_marginal_epistatic(X, snp_pos, maf)

            # print("get marginal epistatic renorm")
            correction_factor = np.sqrt(self.pve*(1.0-self.rho-self.phi)/np.var(y_epi))# / self.n_epistatic_hubs
            y_epi *= correction_factor

        # null effect draw and re-"normalize" gxe effects to have the correct PVE
        if self.phi > 0:
        
            # print("get GxE effects")
            y_gxe = self.get_marginal_GxE(X, snp_pos, maf, self.scalar)
            # print("get GxE effects renorm")
            correction_factor = np.sqrt(self.pve*(self.phi)/np.var(y_gxe))# / self.n_epistatic_hubs
            y_gxe *= correction_factor
        else:
            y_gxe = 0
        
        #sets the error to 1-pve
        # if self.pve > 0:
        #     eps = (1.0-self.pve)*np.var(y_marginal+y_epi+y_gxe)/self.pve
        # else:
        #     eps = 1
        # y_err = np.random.normal(0,np.sqrt(eps),size=self.n_samples)

        y_err = np.random.normal(0,size=self.n_samples)
        y_err=y_err*np.sqrt((1-self.pve)/np.var(y_err))
    
        y = y_marginal + y_epi + y_gxe + y_err
        # print(np.mean(y_marginal),np.mean(y_epi),np.mean(y_gxe),np.mean(y_err))
        #print("rbf pve: ",np.var(self.y_rbf)/np.var(self.y))
        return y_marginal, y_epi, y_gxe, y_err

    def output_data(self):
        return self.y, self.null_alt, self.linear_pve, self.epistatic_pve, self.gxe_pve
    
    def get_epistatic_range(self,snp_idx,snp_pos):
        """
        snp_pos: numpy array of snp positions (in kb)
        """
        epi_snp_pos = snp_pos[snp_idx]
        epi_left_pos = np.searchsorted(snp_pos, epi_snp_pos - self.window_size)
        epi_right_pos = np.searchsorted(snp_pos, epi_snp_pos + self.window_size)
        return epi_left_pos, epi_right_pos
    # get epistatic SNP matrix
    
    def get_W(self, X, snp_idx, snp_pos, maf):
        (min_window,max_window) = self.get_epistatic_range(snp_idx, snp_pos)

        epi_snp = np.expand_dims(X[:,snp_idx],1)

        left_window = X[:,min_window:snp_idx]
        right_window = X[:,snp_idx+1:max_window]

        W = np.multiply(epi_snp, np.hstack((left_window, right_window)))
        self.w_shapes.append(W.shape[1])
        maf_W = np.concatenate((maf[min_window:snp_idx],maf[snp_idx+1:max_window]))

        return W, maf_W

    def get_marginal_epistatic(self, X, snp_pos, maf, correction_factor=1):
        theta = []
        y_epi = np.zeros((self.n_samples))
        for i in range(self.n_epistatic_hubs):
            snp_idx = self.causal_snp_idx_epi[i]
            if correction_factor == 1: 
                W_i, maf_W = self.get_W(X, snp_idx, snp_pos, maf)
                theta_i = np.random.normal(0,np.sqrt((2*maf_W*(1.0-maf_W))**self.exp_alpha),size=W_i.shape[1])
                
            else:
                
                theta_i = self.theta[i]*correction_factor
                W_i = self.get_W(X, snp_idx)
                
            theta.append(theta_i)
            
            y_epi += np.dot(W_i,theta_i)
        
        self.theta = theta
        return y_epi


    def get_marginal_GxE(self, X, snp_pos, maf, scalar,correction_factor=1):
        #This simulation is for GxE, set rho = 1 in parameters to have no epistatic effect
        kappa = []
        E = X[:,self.causal_snp_idx_gxe]

        splitter = int(self.n_samples/2)
        
        y_env1 = np.zeros((splitter))
        y_env2 = np.zeros((splitter))
        # variance = np.sqrt((2*maf*(1.0-maf)))
        # print(variance)
        maf = np.array(maf)[self.causal_snp_idx_gxe]
        kappa_1 = np.random.normal(0,np.sqrt((2*maf*(1.0-maf))**self.exp_alpha),size=len(self.causal_snp_idx_gxe))

        kappa_2 = kappa_1*float(scalar)
            
        kappa.append([kappa_1,kappa_2])
        
                    
        #first half of the individuals are from one environment with a kappa
        y_env1 += np.dot(E[:splitter,:],kappa_1)
        #second half come from a different environment where the effects are perturbed in the same direction
        y_env2 += np.dot(E[splitter:,:],kappa_2)

        y_gxe = np.concatenate((y_env1,y_env2))

        return y_gxe

#n_samples, block_size, pve, rho, psi, win_size, sparsity, phi_sparsity, simtype, exp_alpha = 0.0, null_sim = False
# x = AdditionalArchitecture(50000, 5000, 0.6, 0.5, 0.2, 10000, 5, 5, 1.1, exp_alpha = 0.0, null_sim = False)
# x