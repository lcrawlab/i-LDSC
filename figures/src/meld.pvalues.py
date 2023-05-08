import pandas as pd 
import numpy as np 
from scipy.stats import t, norm

continuous_ukb = ['Height_norm','BMI_norm','MCV_norm','Platelet_norm','DBP_norm','SBP_norm','WBC_norm','RBC_norm','Hemoglobin_norm','Hematocrit_norm','MCH_norm','MCHC_norm','Lymphocyte_norm','Monocyte_norm','Neutrophil_norm','Eosinophil_norm','Basophil_norm','Urate_norm','HBA1C_norm','EGFR_norm','CRP_norm','Triglyceride_norm','HDL_norm','LDL_norm','Cholesterol_norm']
continuous_bbj = ['Height','BMI','MCV','Platelet','DBP','SBP','WBC','RBC','Hemoglobin','Hematocrit','MCH','MCHC','Lymphocyte','Monocyte','Neutrophil','Eosinophil','Basophil','Urate','HBA1C','EGFR','CRP','Triglyceride','HDL','LDL','Cholesterol']


def qqplot(trait,tstat):
    df = []   
    p = 2.0*(norm.sf(np.abs(float(tstat))))
    
    return [trait,p]
   
# df = []
# for trait in continuous_ukb:
#     new = [i.strip().split() for i in open('../../UKBB/output_bivariateldsc/' + trait + '.results', 'r')]
#     df.append(qqplot(trait,new[-2][-1]))
# pd.DataFrame(np.array(df)).to_csv('UKB.meld.pvalues.txt',sep = '\t', index = None, header = None)

# df = []
# for trait in continuous_bbj:
#     new = [i.strip().split() for i in open('../../BBJ/output_bivariateldsc/' + trait + '.results', 'r')]
#     df.append(qqplot(trait,new[-2][-1]))
# pd.DataFrame(np.array(df)).to_csv('BBJ.meld.pvalues.txt',sep = '\t', index = None, header = None)

newdf = []
df = pd.read_csv('../figure_cache/UKB.sLDSC.zscores.txt', sep = '\t').set_index('Trait')

for trait in continuous_bbj:
	
	newdf.append(qqplot(trait,df.loc[trait]['Z']))

pd.DataFrame(np.array(newdf)).to_csv('UKB.sLDSC.meld.pvalues.txt',sep = '\t', index = None, header = None)



