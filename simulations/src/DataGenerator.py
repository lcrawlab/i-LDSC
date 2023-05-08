import numpy as np
from sklearn.preprocessing import scale
from SNPData import SNPData
import pandas as pd

class DataGenerator(object):

    def __init__(self, n_samples, block_size, pve, rho, win_size, sparsity, exp_alpha = 0.0, null_sim = False):
        print("pve {}, rho {}".format(pve, rho))
        ## broad sense heritability
        self.pve = pve
        # proportion of linear effects contributing to broad-sense heritability
        self.rho = rho
        chrom_range = np.arange(1,23)
        #chrom_range = np.arange(1,2)
        self.n_samples = n_samples
        self.block_size = block_size
        
        self.snpdata = SNPData(chrom_range, self.n_samples, self.block_size)

        self.null_sim = null_sim
        # half window size
        #self.window_size = 10 #50 # in terms of number of snps
        #self.window_size = 10000 #100000 # in terms of kb
        self.window_size = win_size
        self.sparsity = sparsity
        self.y_marginal = np.zeros((self.n_samples))
        self.y_epi = np.zeros((self.n_samples))
        self.y_err = np.zeros((self.n_samples))
        self.y = np.zeros((self.n_samples))
        # relationship between MAF-ES
        self.exp_alpha = exp_alpha
        # for testing how many SNPs make up the cis window
        self.w_shapes = []
        idx = 0

        next_snps = np.zeros((0))
        while type(next_snps) == np.ndarray:
            # print(self.snpdata.curr_chrom_idx)
            if not self.snpdata.curr_chrom_idx >= self.snpdata.n_chromosomes:
                curr_chrom = self.snpdata.chromosomes[self.snpdata.curr_chrom_idx]
                maf = pd.read_csv(f"/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/data/UK_ALL_REMAINING_{curr_chrom}.maf",delim_whitespace = True)

            next_snps, snp_pos, snp_id = self.snpdata.get_next_snps()

            if type(next_snps) == np.ndarray:
                # print(next_snps.shape)
                self.n_snps = next_snps.shape[1]
                self.n_epistatic_hubs = int(self.n_snps*(sparsity/100)) 
                # print(maf)
                maf_vals = maf.set_index('SNP').loc[snp_id].MAF.values

                # nan to zero
                next_snps[np.isnan(next_snps)] = 0

                y_marginal, y_epi, y_err = self.generate_effects(rho, pve, next_snps, snp_pos, maf_vals)
                self.y_marginal += y_marginal
                self.y_epi += y_epi
                self.y_err += y_err
                idx += 1
            #if idx > 10:
            #    import pdb; pdb.set_trace() 

        self.y_marginal = self.y_marginal*np.sqrt((self.pve*self.rho)/np.var(self.y_marginal))
        self.y_epi = self.y_epi*np.sqrt((self.pve*(1-self.rho))/np.var(self.y_epi))
        self.y_err=self.y_err*np.sqrt((1-self.pve)/np.var(self.y_err))
        
        self.y = self.y_marginal + self.y_epi + self.y_err
        self.y = scale(self.y)
        
        self.linear_pve = np.var(self.y_marginal.flatten())/np.var(self.y.flatten())
        self.epistatic_pve = np.var(self.y_epi.flatten())/np.var(self.y.flatten())
        print("linear pve: ",self.linear_pve)
        print("epistatic pve: ",self.epistatic_pve)
        print("error pve: ",np.var(self.y_err.flatten())/np.var(self.y.flatten()))
        # self.y = scale(self.y)
        #total_var = np.var(self.y)
        #self.y /= np.sqrt(total_var)

        #self.pve_per_hub = []

    def generate_effects(self, rho, pve, X, snp_pos, maf): 
        #self.pve_per_hub = []
        self.rho = rho
        self.pve = pve

        # randomly choose causal epistatic snps
        self.causal_snp_idx = np.random.choice(self.n_snps,size=self.n_epistatic_hubs,replace=False)
        self.causal_snp_idx = self.causal_snp_idx.tolist()

        # not used now
        self.null_alt = np.array([int(i not in self.causal_snp_idx) for i in range(self.n_snps)])
        # print("causal idx")
        # print(self.causal_snp_idx)
        #self.beta = np.random.normal(0,1,size=self.n_snps)

        # genetic architecture change -- effect sizes related to MAF
        # N(mean,stddev)
        #self.beta = np.random.normal(0,(2*maf*(1.0-maf))**0.5,size=self.n_snps) # this implies alpha=1
        #self.beta = np.random.normal(0,(2*maf*(1.0-maf))**0.25,size=self.n_snps) # this implies alpha=0.5

        self.beta = np.random.normal(0,np.sqrt((2*maf*(1.0-maf))**self.exp_alpha),size=self.n_snps)

        #self.y_marginal = np.zeros((self.n_samples))
        # print("Generating effects")
        y_marginal = np.dot(X,self.beta)
        # print("effects generated")
        # renormalize additive effects
        y_marginal = y_marginal * np.sqrt(self.pve*self.rho/(np.var(y_marginal)))
        # renormalize betas
        #self.beta = self.beta * np.sqrt(self.pve*self.rho/(np.var(self.y_marginal)))
        #self.y_marginal = np.dot(X,self.beta)

        # null
        self.alpha = []
        y_epi = np.zeros((self.n_samples))
        if not self.null_sim and not self.rho == 1.0:
            # print("get marginal epistatic")
            y_epi = self.get_marginal_epistatic(X, snp_pos, maf)

        #self.y_marginal = np.dot(self.X_additive,self.beta)
        #self.beta = self.beta * np.sqrt(self.pve*self.rho/(np.var(self.y_marginal)))
        #self.y_marginal = np.dot(self.X_additive,self.beta)

        # re-"normalize" epistatic y to have the correct PVE
        if not self.null_sim and not self.rho == 1.0:
            # print("get marginal epistatic renorm")
            correction_factor = np.sqrt(self.pve*(1.0-self.rho)/np.var(y_epi))# / self.n_epistatic_hubs
            y_epi *= correction_factor
            #self.y_epi = self.get_marginal_epistatic(X, correction_factor)
            #self.alpha = self.alpha * np.sqrt(self.pve*(1.0-self.rho)/np.var(self.y_epi))

        # average pve per snp
        #self.pve_per_hub = np.mean(self.pve_per_hub)

        # testing rbf kernel
        #from sklearn.metrics.pairwise import rbf_kernel
        #print("rbf")
        ##K = rbf_kernel(self.X)
        #K = np.dot(self.X,self.X.T)/self.X.shape[1] 
        #K = np.multiply(K,K)
        #self.y_rbf = np.dot(K,np.expand_dims(np.random.normal(0,1,self.X.shape[0]),1)).flatten()
        #print("rbf")
        
        # if self.pve > 0:
        #     eps = (1.0-self.pve)*np.var(y_marginal+y_epi)/self.pve
        # else:
        #     eps = 1
        # y_err = np.random.normal(0,np.sqrt(eps),size=self.n_samples)

        y_err = np.random.normal(0,size=self.n_samples)
        y_err = y_err*np.sqrt((1-self.pve)/np.var(y_err))
        

        # y = y_marginal + y_epi + y_err
        #self.y += self.y_rbfs
        #self.y = scale(self.y,with_mean=False,with_std=True)

        #self.linear_pve = np.var(np.dot(self.X,self.beta))/np.var(self.y)


        return y_marginal, y_epi, y_err

    def output_data(self):
        return self.y, self.null_alt, self.linear_pve, self.epistatic_pve

    def get_epistatic_range(self,snp_idx,snp_pos):
        """
        snp_pos: numpy array of snp positions (in kb)
        """
        epi_snp_pos = snp_pos[snp_idx]
        epi_left_pos = np.searchsorted(snp_pos, epi_snp_pos - self.window_size)
        epi_right_pos = np.searchsorted(snp_pos, epi_snp_pos + self.window_size)
        return epi_left_pos, epi_right_pos
        #min_window = max(0,snp_idx-self.window_size)
        #max_window = min(self.n_snps,snp_idx+self.window_size)
        #return (min_window,max_window)

    # get epistatic SNP matrix
    def get_W(self, X, snp_idx, snp_pos, maf):
        (min_window,max_window) = self.get_epistatic_range(snp_idx, snp_pos)
        #X_subset = self.X[:5000,:]
        #epi_snp = np.expand_dims(X_subset[:,snp_idx],1)
        ## exclude the hub SNP
        #W = np.multiply(epi_snp,np.hstack((X_subset[:,min_window:min_window+self.window_size],X_subset[:,min_window+self.window_size+1:max_window])))

        epi_snp = np.expand_dims(X[:,snp_idx],1)
        #epi_snp = X[:,snp_idx]
        # exclude the hub SNP
        #W = np.multiply(epi_snp,np.hstack((self.X[:,min_window:min_window+self.window_size],self.X[:,min_window+self.window_size+1:max_window])))

        left_window = X[:,min_window:snp_idx]
        right_window = X[:,snp_idx+1:max_window]
        #left_window = X[:,min_window:min_window+self.window_size]
        #right_window = X[:,min_window+self.window_size+1:max_window]

        W = np.multiply(epi_snp, np.hstack((left_window, right_window)))
        self.w_shapes.append(W.shape[1])
        maf_W = np.concatenate((maf[min_window:snp_idx],maf[snp_idx+1:max_window]))

        #W = np.multiply(epi_snp,np.hstack((self.X[:,min_window:min_window+self.window_size],self.X[:,min_window+self.window_size+1:max_window])))
        # TODO: should we scale W?? for the correct genetic architecture of epistatic effects
        #if W.shape[1] > 0:
        #    W = scale(W)
        #W = np.multiply(epi_snp,self.X[:,min_window:max_window])
        return W, maf_W

    # get y epistatic with potential correction factor to ensure the correct
    # overall PVE
    # modifies this classes alpha
    # return y_epi
    def get_marginal_epistatic(self, X, snp_pos, maf, correction_factor=1):
        alpha = []
        y_epi = np.zeros((self.n_samples))
        for i in range(self.n_epistatic_hubs):
            snp_idx = self.causal_snp_idx[i]
            if correction_factor == 1: 
                #(min_window,max_window) = self.get_epistatic_range(snp_idx)
                #realized_window_size = max_window-min_window
                W_i, maf_W = self.get_W(X, snp_idx, snp_pos, maf)
                #self.W.append(W_i)
                #W_i = self.W[snp_idx]
                #alpha_i = np.random.normal(0,1,size=W_i.shape[1])
                alpha_i = np.random.normal(0,np.sqrt((2*maf_W*(1.0-maf_W))**self.exp_alpha),size=W_i.shape[1])
                #alpha_i = np.random.normal(0,np.sqrt((2*maf_W*(1.0-maf_W))**self.exp_alpha*\
                #                             (2*maf[snp_idx]*(1.0-maf[snp_idx]))**self.exp_alpha)\
                #                             ,size=W_i.shape[1])
                #print(np.var(np.dot(W_i,alpha_i))/np.var(self.y_marginal))
            else:
                #(min_window,max_window) = self.get_epistatic_range(snp_idx)
                #alpha_i = self.alpha[i*self.window_size*2:(i+1)*self.window_size*2]*correction_factor
                #alpha_i = self.alpha[min_window:max_window]*correction_factor
                alpha_i = self.alpha[i]*correction_factor
                W_i = self.get_W(X, snp_idx)
                #W_i = self.W[i]
                #pve_per_hub = np.var(np.dot(W_i,alpha_i))/np.var(self.y_marginal)
                #print(pve_per_hub)
                #self.pve_per_hub.append(pve_per_hub)
            alpha.append(alpha_i)
            #W_i = self.get_W(snp_idx)
            y_epi += np.dot(W_i,alpha_i)
        #self.alpha = np.hstack(alpha)
        self.alpha = alpha
        return y_epi
