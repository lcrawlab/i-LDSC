
import matplotlib
from matplotlib import pyplot as plt
import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import math
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import t, norm, pearsonr
import argparse
import statsmodels.api as sm
import sys
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

mypalette = ['#1b9e77','#d95f02','#7570b3','#e7298a']
windowpalette = ['#1b9e77','#d95f02','#7570b3','#e7298a']
windowpalette_dict = {1:'#1b9e77',2:'#d95f02',3:'#7570b3',4:'#e7298a'}
alpha = 0.6
rbgscaler = float(255)
windowpalette_rgba = {1:(27/rbgscaler,158/rbgscaler,119/rbgscaler,alpha),2:(217/rbgscaler,95/rbgscaler,2/rbgscaler,alpha),3:(117/rbgscaler,112/rbgscaler,179/rbgscaler,alpha),4:(231/rbgscaler,41/rbgscaler,138/rbgscaler,alpha)}

alphaparams = {'sim1':0,'sim2':0,'sim3':0.5,'sim4':1,'sim5':0,'sim6':0.5,'sim7':1}
windowparams = {'sim1':10,'sim2':100,'sim3':10,'sim4':10,'sim5':10,'sim6':100,'sim7':100}

def MAE_calculator(df,rho):

    outdf = pd.DataFrame(np.zeros((5,4)),columns = pd.unique(df['rho']), index = pd.unique(df['est_win']))
    # df = df.groupby(['est_win'])
    # print(df)
    df['h2add'] = df['pve'] * df['rho']
    df['mae'] = np.abs(df['h2add'] - df["Observed_scale_h2"])
    returndf = df

    averagedf = df[(df['est_win'] == 10) | (df['est_win'] == 20) | (df['est_win'] == 50) | (df['est_win'] == 100)]
    print(averagedf.groupby(['rho']).sem())

    temp1 = df.groupby(['est_win','rho']).mean()
    temp = df.groupby(['est_win','rho']).sem()
    temp1['sem_mae'] = temp['mae']
    temp1['h2_sem'] = temp['Observed_scale_h2']
    df = temp1

    # print(df)
    df.to_csv('../figures/mae.table.' + str(rho) + '.txt',sep = '\t')
    return returndf

def qqplot(stats,ax,color):

    stats["t"] = stats.Coefficients/stats.Coefficient_SE
    stats["P"] = 2.0*(1.0 - norm.cdf(np.abs(stats.t)))

    stats = stats.P
    n = len(stats) # number of tests
    obs_p = -np.log10(np.sort(stats))
    th_p = np.arange(1/float(n),1+(1/float(n)),1/float(n))

    th_p = -np.log10(th_p)
    ax.scatter(th_p,obs_p,color = color,alpha = 0.9)
    x = np.linspace(*ax.get_xlim())

    ax.plot(x, x,color = 'black')
    if color == '#1b9e77':
        cix = norm.ppf(np.arange(1/float(n),1+(1/float(n)),1/float(n))[:-1])
        a, b = np.polyfit(th_p, th_p, deg=1)
        y_est = a * x + b
        y_err = x.std() * np.sqrt(1/len(x) + (x - x.mean())**2 / np.sum((x - x.mean())**2))
        ax.fill_between(x, y_est - y_err, y_est + y_err, alpha=0.5,color = 'grey')


    return fig

results_dir = "/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld-annotations-run/figures/figure_cache/"
plot_dir = "/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld-annotations-run/figures/figures/"

parser = argparse.ArgumentParser()
parser.add_argument('-f','--figs',type=str,help='figures to make')
parser.add_argument('-d','--dir',type=str,help='figure directory')
args = parser.parse_args()
figs = [int(x) for x in args.figs.split(",")]
print(figs)
fig_dir = args.dir
# stats = pd.read_pickle(results_dir + fig_dir + ".pkl")
# try:
#     stats_add = pd.read_pickle(results_dir + fig_dir + "_add.pkl")
# except: # if no additive results
#    pass 



"""
FIG 1
"""
if 1 in figs:
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16,9))
    #fig.suptitle('MELD main text fig 1')

    facetA = stats.loc[(stats.Categories.str.contains("ML")) & (stats.rho == 0.5) & \
                       (stats.pve == 0.6)]
    sns.barplot(ax=ax3,x="sparsity",y="sig",hue="est_win",\
                data=facetA,
                capsize=.15,
                palette=windowpalette,
                edgecolor = 'white',
                errwidth=1.5,
                saturation = 0.8,)
    ax3.set_title(r"$H^2=0.6, \rho=0.5$")
    ax3.set_ylabel("Power")
    ax3.set_xlabel("Percentage sparsity")
    ax3.legend().set_visible(False)
    rectanglecount = 1

    # for j in ax1.get_children():
    #   if type(j) == matplotlib.patches.Rectangle and rectanglecount <= 12:
    #       j.set_facecolor(windowpallette_rgba[math.ceil(rectanglecount/3)])
    #       rectanglecount +=1

    facetB = stats.loc[(stats.Categories.str.contains("ML")) & (stats.rho == 0.8) & \
                       (stats.pve == 0.6)]
    sns.barplot(ax=ax4,x="sparsity",y="sig",hue="est_win",\
                data=facetB,
                capsize=.15,
                palette=windowpalette,
                edgecolor = 'white',
                errwidth=1.5,
                saturation = 0.8,)
    ax4.set_title(r"$H^2=0.6, \rho=0.8$")
    ax4.set_ylabel("Power")
    ax4.set_xlabel("Percentage sparsity")
    ax4.legend().set_visible(False)

    facetC = stats.loc[(stats.Categories.str.contains("ML")) & (stats.rho == 0.5) & \
                       (stats.pve == 0.3)]
    sns.barplot(ax=ax1,x="sparsity",y="sig",hue="est_win",\
                data=facetC,
                capsize=.15,
                palette=windowpalette,
                edgecolor = 'white',
                errwidth=1.5,
                saturation = 0.8,)
    ax1.set_title(r"$H^2=0.3, \rho=0.5$")
    ax1.set_ylabel("Power")
    ax1.set_xlabel("Percentage sparsity")
    ax1.legend().set_visible(False)

    facetD = stats.loc[(stats.Categories.str.contains("ML")) & (stats.rho == 0.8) & \
                       (stats.pve == 0.3)]
    sns.barplot(ax=ax2,x="sparsity",y="sig",hue="est_win",\
                data=facetD,
                capsize=.15,
                palette=windowpalette,
                edgecolor = 'white',
                errwidth=1.5,
                saturation = 0.8,)
    ax2.set_title(r"$H^2=0.3, \rho=0.8$")
    ax2.set_ylabel("Power")
    ax2.set_xlabel("Percentage sparsity")
    ax2.legend(loc='upper right',title="window size")

    for ax in fig.get_axes():
        ax.label_outer()
    sns.despine()
    fig.savefig(os.path.join(plot_dir,"meld_fig_1_" + fig_dir + "_alpha" + str(alphaparams[fig_dir]) + "_window" + str(windowparams[fig_dir]) + ".pdf"))


"""
FIG 2
"""
if 2 in figs:


    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(16,9))

    facetA = stats.loc[(stats.Categories.str.contains("ML")) & (stats.rho == 1.0) &\
                       (stats.pve == 0.6)]

    sns.barplot(ax=ax1,x="est_win",y="sig",\
                data=facetA,
                capsize=.45,
                edgecolor = 'black',
                errwidth = 1.5,
                linewidth = 2,
                palette=windowpalette,
                alpha = 0.8)

    rectanglecount = 1
    for j in ax1.get_children():
        if type(j) == matplotlib.patches.Rectangle and rectanglecount <= 4:
            j.set_facecolor(windowpalette_rgba[rectanglecount])
            j.set_edgecolor(windowpalette_rgba[rectanglecount])

            rectanglecount +=1


    ax1.axhline(0.05, color='black', linestyle='--')
    ax1.set_ylabel("False discovery rate")


    h2_col = "Observed_scale_h2"
    sns.boxplot(x='est_win',y=h2_col,ax = ax2,
                data=facetA,fliersize=0,linewidth = 2)

    counter = 0
    subcount=1
    for j in ax2.get_children():
        if type(j) == matplotlib.patches.PathPatch and subcount <= 4:

            j.set_facecolor(windowpalette_rgba[subcount])
            j.set_edgecolor(windowpalette_dict[subcount])

            subcount += 1
        if type(j) == matplotlib.lines.Line2D:
            counter+=1
            
            j.set_color(windowpalette_dict[math.ceil(counter/float(6))])



    ax2.axhline(0.0, color='black', linestyle='--')
    ax2.set_ylabel(r'$\sigma$')

    facetC = stats.loc[(stats.Categories.str.contains("LD")) & (stats.rho == 1.0) & \
                       (stats.pve == 0.6)]
    h2_col = "Observed_scale_h2"
    sns.boxplot(ax=ax3,x="est_win",y=h2_col,
                data=facetC,\
                showfliers=False,linewidth = 2)


    counter = 0
    subcount=1
    for j in ax3.get_children():
        if type(j) == matplotlib.patches.PathPatch and subcount <= 4:
            j.set_facecolor(windowpalette_rgba[subcount])
            j.set_edgecolor(windowpalette_dict[subcount])
            subcount += 1

        if type(j) == matplotlib.lines.Line2D:
            counter+=1

            j.set_color(windowpalette_dict[math.ceil(counter/float(5))])

    ax3.axhline(0.6, color='black', linestyle='--')
    ax3.set_ylabel(r'$h^2$')

    ax1.set_xlabel("Estimation Window Size")
    ax2.set_xlabel("Estimation Window Size")
    ax3.set_xlabel("Estimation Window Size")

    qqplot(facetA.loc[(facetA.est_win==10)],ax=ax4,color = windowpalette[0])
    qqplot(facetA.loc[(facetA.est_win==20)],ax=ax4,color = windowpalette[1])
    qqplot(facetA.loc[(facetA.est_win==50)],ax=ax4,color = windowpalette[2])
    qqplot(facetA.loc[(facetA.est_win==100)],ax=ax4,color = windowpalette[3])

    ax4.set_xlim(left = 0)
    ax4.set_ylim(bottom = 0)
    ax4.set_xlabel(r'Expected $-log_{10}(p)$')
    ax4.set_ylabel(r'Observed $-log_{10}(p)$')

    sns.despine()
    fig.savefig(os.path.join(plot_dir,"meld_fig_2_" + fig_dir + "_alpha" + str(alphaparams[fig_dir]) + "_window" + str(windowparams[fig_dir]) + "_37_H2_0.6.pdf"))

methodpalette = ['#e41a1c','#377eb8','#4daf4a']
methodpalette_dict = {'LDSC':'#e41a1c','MELD':'#377eb8','LDER':'#984ea3'} 
methodpalette_rgba = {'LDSC':(228/rbgscaler,26/rbgscaler,28/rbgscaler,alpha),'MELD':(55/rbgscaler,126/rbgscaler,184/rbgscaler,alpha),'LDER':(152/rbgscaler,78/rbgscaler,163/rbgscaler,alpha)}
methodpalette3_dict = {1:'#e41a1c',0:'#377eb8',2:'#984ea3'} 
methodpalette3_rgba = {1:(228/rbgscaler,26/rbgscaler,28/rbgscaler,alpha),0:(55/rbgscaler,126/rbgscaler,184/rbgscaler,alpha),2:(152/rbgscaler,78/rbgscaler,163/rbgscaler,alpha)}
methodpalette3 = ['#377eb8','#e41a1c','#4daf4a']

windowpalette = ['#e6ab02','#1b9e77','#d95f02','#7570b3','#e7298a']
windowpalette_dict = {0:'#e41a1c',1:'#1b9e77',2:'#d95f02',3:'#7570b3',4:'#e7298a'}
windowpalette_rgba = {0:(228/rbgscaler,26/rbgscaler,28/rbgscaler,alpha),1:(27/rbgscaler,158/rbgscaler,119/rbgscaler,alpha),2:(217/rbgscaler,95/rbgscaler,2/rbgscaler,alpha),3:(117/rbgscaler,112/rbgscaler,179/rbgscaler,alpha),4:(231/rbgscaler,41/rbgscaler,138/rbgscaler,alpha)}

"""
FIG 3
"""
if 3 in figs:

    fig, ((ax1, ax2)) = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(8,6))

    stats = pd.read_pickle(results_dir + fig_dir + '.figure.data.pkl')

    LD = stats.loc[(stats.Categories == "LDR2_1.0") ]
    ML = stats.loc[(stats.Categories == "MLR2_1.0") ]
    MLLD = LD.merge(ML,on="sim_idx")

    # print(MLLD.columns)

    # print(MLLD[['Observed_scale_h2_x','Observed_scale_h2_y']])
    MLLD['total_meld'] = MLLD['Observed_scale_h2_x'] + MLLD['Observed_scale_h2_y']
    MLLD['ldsc_add'] = MLLD['Observed_scale_h2_x']

    MLLD_plot = MLLD[['pve_x','rho_x','total_meld','ldsc_add','est_win_x']]
    # print(MLLD_plot)

    h2_06 = MLLD_plot[(MLLD_plot['pve_x'] == 0.6) & (MLLD_plot['est_win_x'] == 100)]
    h2_06_total_meld = h2_06[['total_meld','rho_x']].rename(columns = {'total_meld':'estimate', 'rho_x':'rho'})
    h2_06_total_meld['label'] = 'total_meld'

    h2_06_total_ldsc = h2_06[['ldsc_add','rho_x']].rename(columns = {'ldsc_add':'estimate','rho_x':'rho'})
    h2_06_total_ldsc['label'] = 'total_ldsc'
    h2_06 = pd.concat([h2_06_total_meld,h2_06_total_ldsc])

    sns.boxplot(ax=ax1,x="rho",y='estimate', hue="label",\
                data=h2_06,
                
                fliersize = 0)
    ax1.set_title(r"$h^2=0.6$")
    ax1.set_ylabel(r"$\hat{h}^2$")
    ax1.set_xlabel(r"$\rho$")
    [ax1.axhline(y=i*0.6,xmin=idx*0.25,xmax=(idx+1)*0.25,linestyle='--',color = 'grey') for idx,i in enumerate(np.array([0.2,0.4,0.6,0.8]))]
    ax1.axhline(y=0.6,linestyle='--',color = 'black')


    h2_03 = MLLD_plot[(MLLD_plot['pve_x'] == 0.3) & (MLLD_plot['est_win_x'] == 100)]
    h2_03_total_meld = h2_03[['total_meld','rho_x']].rename(columns = {'total_meld':'estimate', 'rho_x':'rho'})
    h2_03_total_meld['label'] = 'total_meld'

    h2_03_total_ldsc = h2_03[['ldsc_add','rho_x']].rename(columns = {'ldsc_add':'estimate','rho_x':'rho'})
    h2_03_total_ldsc['label'] = 'total_ldsc'
    h2_03 = pd.concat([h2_03_total_meld,h2_03_total_ldsc])

    sns.boxplot(ax=ax2,x="rho",y='estimate', hue="label",\
                data=h2_03,
                fliersize = 0)
    ax2.set_title(r"$h^2=0.3$")
    ax2.set_ylabel(r"$\hat{h}^2$")
    ax2.set_xlabel(r"$\rho$")
    [ax2.axhline(y=i*0.3,xmin=idx*0.25,xmax=(idx+1)*0.25,linestyle='--',color = 'grey') for idx,i in enumerate(np.array([0.2,0.4,0.6,0.8]))]
    ax2.axhline(y=0.3,linestyle='--',color = 'black')

    counter = 1
    subcount = 0
    for j in ax1.get_children():
        if type(j) == matplotlib.patches.PathPatch:
            j.set_facecolor(methodpalette3_rgba[subcount%2])
            j.set_edgecolor(methodpalette3_dict[subcount%2]) 
            subcount += 1

        if type(j) == matplotlib.lines.Line2D and counter < 48:
            j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%2])
            counter+=1

    
    counter = 1
    subcount = 0
    for j in ax2.get_children():
        # print(j,counter)
        if type(j) == matplotlib.patches.PathPatch:
            # print(subcount)
            j.set_facecolor(methodpalette3_rgba[subcount%2])
            j.set_edgecolor(methodpalette3_dict[subcount%2])
            subcount += 1

        if type(j) == matplotlib.lines.Line2D and counter < 48:
             # and counter <= 48
            j.set_color(methodpalette3_rgba[(subcount+1)%2])
            counter+=1

    for ax in fig.get_axes():
        ax.label_outer()

    total_line = mlines.Line2D([], [], color='grey', linestyle = '--',label=r'Generative $h^2$')
    additive_line = mlines.Line2D([], [], color='black', linestyle = '--',label='Generative\nadditive pve')
    ildsc_patch = Patch(edgecolor=(55/rbgscaler,126/rbgscaler,184/rbgscaler,1), facecolor = (55/rbgscaler,126/rbgscaler,184/rbgscaler,0.5), label='iLDSC')
    ldsc_patch = Patch(edgecolor=(228/rbgscaler,26/rbgscaler,28/rbgscaler,1), facecolor = (228/rbgscaler,26/rbgscaler,28/rbgscaler,0.5), label='LDSC')

    ax1.legend(handles=[total_line, additive_line, ildsc_patch, ldsc_patch], ncol = 1,loc = 'upper right',bbox_to_anchor=(1.3, 1.0))
    ax2.legend().set_visible(False)

    fig.subplots_adjust(right=0.8)

    sns.despine()
    fig.savefig(os.path.join(plot_dir,"fig_3_" + fig_dir + ".pdf"))
    
    ####stratify the plot based on estimation window size
    plt.clf()

    fig, ((ax1, ax2)) = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(8,6))
    
    MLLD_plot = MLLD[['pve_x','rho_x','est_win_x','total_meld','ldsc_add']]
    

    h2_06 = MLLD_plot[MLLD_plot['pve_x'] == 0.6]
    h2_06_total_meld = h2_06[['total_meld','rho_x','est_win_x']].rename(columns = {'total_meld':'estimate', 'rho_x':'rho'})
    h2_06_total_meld['label'] = h2_06_total_meld['est_win_x']

    h2_06_total_ldsc = h2_06[['ldsc_add','rho_x','est_win_x']].rename(columns = {'ldsc_add':'estimate','rho_x':'rho'})
    h2_06_total_ldsc['label'] = 'total_ldsc'
    h2_06 = pd.concat([h2_06_total_meld,h2_06_total_ldsc])

    # print(h2_06)
    sns.boxplot(ax=ax1,x="rho",y='estimate', hue="label",\
                data=h2_06,
                fliersize = 0, 
                hue_order = ['total_ldsc',10,20,50,100])

    ax1.set_title(r"$h^2=0.6$")
    ax1.set_ylabel(r"$\hat{h}^2$")
    ax1.set_xlabel(r"$\rho$")
    [ax1.axhline(y=i*0.6,xmin=idx*0.25,xmax=(idx+1)*0.25,linestyle='--',color = 'grey') for idx,i in enumerate(np.array([0.2,0.4,0.6,0.8]))]
    ax1.axhline(y=0.6,linestyle='--',color = 'black')

    h2_03 = MLLD_plot[MLLD_plot['pve_x'] == 0.3]
    h2_03_total_meld = h2_03[['total_meld','rho_x','est_win_x']].rename(columns = {'total_meld':'estimate', 'rho_x':'rho'})
    h2_03_total_meld['label'] = h2_03_total_meld['est_win_x']

    h2_03_total_ldsc = h2_03[['ldsc_add','rho_x','est_win_x']].rename(columns = {'ldsc_add':'estimate','rho_x':'rho'})
    h2_03_total_ldsc['label'] = 'total_ldsc'
    h2_03 = pd.concat([h2_03_total_meld,h2_03_total_ldsc])
    print(h2_03['label'].unique())
    sns.boxplot(ax=ax2,x="rho",y='estimate', hue="label",\
                data=h2_03,
                fliersize = 0,
                hue_order = ['total_ldsc',10,20,50,100])

    ax2.set_title(r"$h^2=0.3$")
    ax2.set_ylabel(r"$\hat{h}^2$")
    ax2.set_xlabel(r"$\rho$")
    [ax2.axhline(y=i*0.3,xmin=idx*0.25,xmax=(idx+1)*0.25,linestyle='--',color = 'grey') for idx,i in enumerate(np.array([0.2,0.4,0.6,0.8]))]
    ax2.axhline(y=0.3,linestyle='--',color = 'black')


    counter = 0
    subcount=0
    for j in ax1.get_children():

        if type(j) == matplotlib.lines.Line2D and counter < 144:
            j.set_color(windowpalette_dict[((subcount-1)%5)])
            
        if type(j) == matplotlib.patches.PathPatch:
            j.set_facecolor(windowpalette_rgba[subcount%5])
            j.set_edgecolor(windowpalette_dict[subcount%5])
            subcount += 1

        counter+=1
    
    counter = 0
    subcount = 0
    for j in ax2.get_children():
        if type(j) == matplotlib.lines.Line2D and counter < 144:
            j.set_color(windowpalette_dict[((subcount-1)%5)])
            
        if type(j) == matplotlib.patches.PathPatch:
            j.set_facecolor(windowpalette_rgba[subcount%5])
            j.set_edgecolor(windowpalette_dict[subcount%5])
            subcount += 1

        counter+=1

    total_line = mlines.Line2D([], [], color='black', linestyle = '--',label=r'Generative $h^2$')
    additive_line = mlines.Line2D([], [], color='grey', linestyle = '--',label='Generative\nadditive pve')
    
    ldsc_patch = Patch(edgecolor=(228/rbgscaler,26/rbgscaler,28/rbgscaler,1), facecolor = (228/rbgscaler,26/rbgscaler,28/rbgscaler,0.5), label='LDSC')
    
    patch10 = Patch(edgecolor=(27/rbgscaler,158/rbgscaler,119/rbgscaler,1), facecolor = (27/rbgscaler,158/rbgscaler,119/rbgscaler,0.5), label='+/- 10 SNPs')
    patch20 = Patch(edgecolor=(217/rbgscaler,95/rbgscaler,2/rbgscaler,1), facecolor = (217/rbgscaler,95/rbgscaler,2/rbgscaler,0.5), label='+/- 20 SNPs')
    patch50 = Patch(edgecolor=(117/rbgscaler,112/rbgscaler,179/rbgscaler,1), facecolor = (117/rbgscaler,112/rbgscaler,179/rbgscaler,0.5), label='+/- 50 SNPs')
    patch100 = Patch(edgecolor=(231/rbgscaler,41/rbgscaler,138/rbgscaler,1), facecolor = (231/rbgscaler,41/rbgscaler,138/rbgscaler,0.5), label='+/- 100 SNPs')
    
    ax1.legend(handles=[total_line, additive_line, ldsc_patch, patch10, patch20, patch50, patch100], ncol = 1,loc = 'upper right',bbox_to_anchor=(1.3, 1.0))
    ax2.legend().set_visible(False)

    fig.subplots_adjust(right=0.8)

    sns.despine()
    fig.savefig(os.path.join(plot_dir,"fig_3_" + fig_dir + "_windows.pdf"))


'''Figure 4'''

if 4 in figs:
    # tof = ['Height','Triglyceride']

    ukbmeld = pd.read_csv(results_dir + 'UKB_MELD_estimates.txt',sep = '\t')
    ukbmeld = ukbmeld[['Trait','h2','intercept','lambdaGC']]
    ukbmeld.columns = ['Trait','h2_meld','intercept','lambda']
    ukbmeld['Trait'] = ukbmeld['Trait'].str.replace('_norm','')
    
    ukbldsc = pd.read_csv(results_dir + 'UKB_LDSC_estimates.txt',sep = '\t')
    ukbldsc = ukbldsc[['Trait','h2','intercept','lambdaGC']]
    ukbldsc.columns = ['Trait','h2_ldsc','intercept','lambda']
    ukbldsc['Trait'] = ukbldsc['Trait'].str.replace('_norm','')
    # ukbmerged = ukblder.merge(ukbmeld,on = 'Trait',how = 'inner').dropna()

    bbjmeld = pd.read_csv(results_dir + 'BBJ_MELD_estimates.txt',sep = '\t')
    bbjmeld = bbjmeld[['Trait','h2','intercept','lambdaGC']] 
    bbjmeld.columns = ['Trait','h2_meld','intercept','lambda']

    bbjldsc = pd.read_csv(results_dir + 'BBJ_LDSC_estimates.txt',sep = '\t')
    bbjldsc = bbjldsc[['Trait','h2','intercept','lambdaGC']]
    bbjldsc.columns = ['Trait','h2_ldsc','intercept','lambda']

    ukbmeld['meldsigma'] = ukbmeld['h2_meld'] - ukbldsc['h2_ldsc']
    bbjmeld['meldsigma'] = bbjmeld['h2_meld'] - bbjldsc['h2_ldsc']
    
    finaldata = ukbmeld[['Trait','meldsigma']].merge(bbjmeld[['Trait','meldsigma']], on = 'Trait')
    finaldata.columns = ['Trait','ukb','bbj']
    # print(finaldata)

    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2, figsize = (6,6))

    #begin panel A
    ukbfitdata = ukbldsc[['Trait','h2_ldsc']].merge(ukbmeld[['Trait','h2_meld','meldsigma']], on = 'Trait')
    bbjfitdata = bbjldsc[['Trait','h2_ldsc']].merge(bbjmeld[['Trait','h2_meld','meldsigma']], on = 'Trait')


    ax1.scatter(bbjfitdata['h2_ldsc'],bbjfitdata['h2_meld'], color = '#984ea3',s = 1.5)
    ax1.scatter(ukbfitdata['h2_ldsc'],ukbfitdata['h2_meld'], color = '#4daf4a',s = 1.5)
    xlimtemp = ax1.get_xlim()
    ylimtemp = ax1.get_ylim()
    

    b, a = np.polyfit(bbjfitdata['h2_ldsc'], bbjfitdata['h2_meld'], deg=1)
    xseq = np.linspace(0,1, num=100)
    ax1.plot(xseq, a + b * xseq, color="#984ea3", lw=1)

    print('panel A bbj ldsc-meld',bbjfitdata[['h2_ldsc','h2_meld']].corr()**2)
    print('panel A bbj ldsc-meld',pearsonr(bbjfitdata['h2_ldsc'],bbjfitdata['h2_meld']))

    print('panel A ukb ldsc-meld',ukbfitdata[['h2_ldsc','h2_meld']].corr()**2)
    print('panel A ukb ldsc-meld',pearsonr(ukbfitdata['h2_ldsc'],ukbfitdata['h2_meld']))

    b, a = np.polyfit(ukbfitdata['h2_ldsc'], ukbfitdata['h2_meld'], deg=1)
    ax1.plot(xseq, a + b * xseq, color="#4daf4a", lw=1)

    ax1.set_xlim(xlimtemp)
    ax1.set_ylim(ylimtemp)
    xseq = np.linspace(0,0.8, num=100)
    ax1.plot(xseq,xseq,color = 'black',linestyle='--')

    ax1.set_xlabel(r"LDSC $\hat{h}^2$")
    ax1.set_ylabel(r"i-LDSC $\hat{h}^2$")
    
    ###Begin panel b - comparing h2_add estimates
    additives = ukbmeld[['Trait','h2_meld']].merge(bbjmeld[['Trait','h2_meld']], on = 'Trait')
    additives.columns = ['Trait','ukb','bbj']

    ldsc_additives = bbjldsc[['Trait','h2_ldsc']].merge(ukbldsc[['Trait','h2_ldsc']], on = 'Trait')
    ldsc_additives.columns = ['Trait','bbj','ukb']

    ax2.scatter(additives['ukb'],additives['bbj'],color = '#377eb8', s = 1.5)
    ax2.scatter(ldsc_additives['ukb'],ldsc_additives['bbj'],color = '#e41a1c', s = 1.5)
    
    xlimtemp = ax2.get_xlim()
    ylimtemp = ax2.get_ylim()

    b, a = np.polyfit(additives['ukb'], additives['bbj'], deg=1)
    ax2.plot(xseq, a + b * xseq, color="#377eb8", lw=1)

    b, a = np.polyfit(ldsc_additives['ukb'], ldsc_additives['bbj'], deg=1)
    ax2.plot(xseq, a + b * xseq, color="#e41a1c", lw=1)
    ax2.set_xlabel(r"UK Biobank $\hat{h}^2$")
    ax2.set_ylabel(r"Biobank Japan $\hat{h}^2$")

    ax2.set_xlim(xlimtemp)
    ax2.set_ylim(ylimtemp)
    print('panel B meld',additives[['ukb','bbj']].corr()**2)
    print('panel B meld',pearsonr(additives['ukb'],additives['bbj']))

    print('panel B ldsc',ldsc_additives[['ukb','bbj']].corr()**2)
    print('panel B ldsc',pearsonr(ldsc_additives['ukb'],ldsc_additives['bbj']))


    ###Begin panel c - comparison of meld estimates
    xseq = np.linspace(-1,1, num=201)

    bias = ukbmeld[['Trait','meldsigma']].merge(bbjmeld[['Trait','meldsigma']], on = 'Trait')
    bias.columns = ['Trait','ukb','bbj']
    print(bias)
    ax3.scatter(bias['ukb'],bias['bbj'],color = '#377eb8', s = 1.5)
    xlimtemp = ax3.get_xlim()
    ylimtemp = ax3.get_ylim()

    print('meld sigma',bias.corr()**2)
    print('meld sigma',pearsonr(bias['ukb'],bias['bbj']))

    b, a = np.polyfit(bias['ukb'], bias['bbj'], deg=1)
    ax3.plot(xseq, a + b * xseq, color="#377eb8", lw=1)
    ax3.set_xlabel(r"UK Biobank $\hat{\sigma}$")
    ax3.set_ylabel(r"Biobank Japan $\hat{\sigma}$")

    ax3.set_xlim(xlimtemp)
    ax3.set_ylim(ylimtemp)

    #Make intercept plot

    ukbfitdata = ukbldsc[['Trait','intercept']].merge(ukbmeld[['Trait','intercept']], on = 'Trait')
    ukbfitdata.columns = ['Trait','LDSC','MELD']

    bbjfitdata = bbjldsc[['Trait','intercept']].merge(bbjmeld[['Trait','intercept']], on = 'Trait')
    bbjfitdata.columns = ['Trait','LDSC','MELD']

    print('ukb intercept',ukbfitdata[['LDSC','MELD']].corr()**2)
    print('ukb intercept', pearsonr(ukbfitdata['LDSC'],ukbfitdata['MELD']))

    print('bbj intercept',bbjfitdata[['LDSC','MELD']].corr()**2)
    print('bbj intercept', pearsonr(bbjfitdata['LDSC'],bbjfitdata['MELD']))
    
    ax4.scatter(bbjfitdata['LDSC'],bbjfitdata['MELD'], color = '#984ea3',s = 1.5)
    ax4.scatter(ukbfitdata['LDSC'],ukbfitdata['MELD'], color = '#4daf4a',s = 1.5)

    b, a = np.polyfit(bbjfitdata['LDSC'], bbjfitdata['MELD'], deg=1)
    xseq = np.linspace(1,1.8, num=100)
    ax4.plot(xseq, a + b * xseq, color="#984ea3", lw=1)

    b, a = np.polyfit(ukbfitdata['LDSC'], ukbfitdata['MELD'], deg=1)
    ax4.plot(xseq, a + b * xseq, color="#4daf4a", lw=1)

    ax4.plot(xseq,xseq,color = 'black',linestyle='--')
    xlimtemp = ax4.get_xlim()
    ylimtemp = ax4.get_ylim()

    ax4.set_xlabel(r"LDSC intercept")
    ax4.set_ylabel(r"i-LDSC intercept")
    
    ax1.spines[['right','top']].set_visible(False)
    ax2.spines[['right','top']].set_visible(False)
    ax3.spines[['right','top']].set_visible(False)
    ax4.spines[['right','top']].set_visible(False)    

    plt.tight_layout()
    fig.savefig(os.path.join(plot_dir,"figure4.pdf"))

    plt.clf()
    
if 5 in figs:

    if fig_dir in ['GxE_noCASS','GxE']:
        fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(4, 1, sharex=True, figsize=(8,12))

        simdata = pd.read_pickle(results_dir + fig_dir + '.figure.data.pkl')
        
        simdata['meldh2_total_pve_diff'] = simdata['pve'] - simdata['total_meld'] 
        simdata['ldsch2_total_pve_diff'] = simdata['pve'] - simdata['ldsc_add'] 
        simdata['truemeld_meldestimate_pve_diff'] = simdata['meld_score'] - simdata['epistatic_pve']

        h2_06 = simdata[simdata['pve'] == 0.6]
        h2_06_total_meld = h2_06[['total_meld','amplifier']].rename(columns = {'total_meld':'estimate'})
        h2_06_total_meld['label'] = 'total_meld'

        h2_06_total_ldsc = h2_06[['ldsc_add','amplifier']].rename(columns = {'ldsc_add':'estimate'})
        h2_06_total_ldsc['label'] = 'total_ldsc'

        h2_06 = pd.concat([h2_06_total_meld,h2_06_total_ldsc,])
        
        methodpalette3_dict = {1:'#e41a1c',0:'#377eb8'} 
        methodpalette3_rgba = {1:(228/rbgscaler,26/rbgscaler,28/rbgscaler,alpha), 0:(55/rbgscaler,126/rbgscaler,184/rbgscaler,alpha),}

        sns.boxplot(data = h2_06, x = 'amplifier', y = 'estimate', hue = 'label', ax = ax1, fliersize = 0)

        fliersize = 0
        ax1.set_title("Generative model\n" + r"$h^2=0.6$, additive pve = 0.3, CASS pve = 0, GxE pve = 0.3")
        ax1.set_ylabel(r"$\hat{h}^2$")
        ax1.set_xlabel(r"Amplifier")
        ax1.legend().set_visible(False)

        counter = 1
        subcount = 0
        for j in ax1.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                j.set_facecolor(methodpalette3_rgba[subcount%2])
                j.set_edgecolor(methodpalette3_dict[subcount%2])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%2])
                counter+=1
        ax1.axhline(y=0.6,linestyle='--',color = 'black')

        h2_06 = simdata[simdata['pve'] == 0.6]

        h2_06_total_meld = h2_06[['meld_score','amplifier']]
        sns.boxplot(data = h2_06_total_meld, x = 'amplifier', y = 'meld_score', ax = ax2, fliersize = 0)

        counter = 1
        subcount = 0
        for j in ax2.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                j.set_facecolor(methodpalette3_rgba[subcount%1])
                j.set_edgecolor(methodpalette3_dict[subcount%1])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%1])
                counter+=1
        
        ax2.set_ylabel("iLDSC component\nestimate")
        ax2.set_xlabel(r"Amplifier")
    
        h2_03 = simdata[simdata['pve'] == 0.3]
        h2_03_total_meld = h2_03[['total_meld','amplifier']].rename(columns = {'total_meld':'estimate'})
        h2_03_total_meld['label'] = 'total_meld'

        h2_03_total_ldsc = h2_03[['ldsc_add','amplifier']].rename(columns = {'ldsc_add':'estimate'})
        h2_03_total_ldsc['label'] = 'total_ldsc'
        h2_03 = pd.concat([h2_03_total_meld,h2_03_total_ldsc])
        
        methodpalette3_dict = {1:'#e41a1c',0:'#377eb8'} 
        methodpalette3_rgba = {1:(228/rbgscaler,26/rbgscaler,28/rbgscaler,alpha), 0:(55/rbgscaler,126/rbgscaler,184/rbgscaler,alpha),}

        sns.boxplot(data = h2_03, x = 'amplifier', y = 'estimate', hue = 'label', ax = ax3, fliersize = 0)

        fliersize = 0
        ax3.set_title("Generative model\n" + r"$h^2=0.3$, additive pve = 0.15, CASS pve = 0, GxE pve = 0.15")
        ax3.set_ylabel(r"$\hat{h}^2$")
        ax3.set_xlabel(r"Amplifier")
        ax3.legend().set_visible(False)

        counter = 1
        subcount = 0
        for j in ax3.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                j.set_facecolor(methodpalette3_rgba[subcount%2])
                j.set_edgecolor(methodpalette3_dict[subcount%2])
                subcount += 1
            
            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%2])
                counter+=1

        ax3.axhline(y=0.3,linestyle='--',color = 'black')

        h2_03 = simdata[simdata['pve'] == 0.3]
        h2_03_total_meld = h2_03[['meld_score','amplifier']]        
        sns.boxplot(data = h2_03_total_meld, x = 'amplifier', y = 'meld_score', ax = ax4, fliersize = 0)

        counter = 1
        subcount = 0
        for j in ax4.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                # print(subcount)
                j.set_facecolor(methodpalette3_rgba[subcount%1])
                j.set_edgecolor(methodpalette3_dict[subcount%1])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%1])
                counter+=1
        
        ax4.set_xlabel(r"Amplifier")
        ax4.set_ylabel("iLDSC component\nestimate")

        ildsc_patch = Patch(edgecolor=(55/rbgscaler,126/rbgscaler,184/rbgscaler,1), facecolor = (55/rbgscaler,126/rbgscaler,184/rbgscaler,0.5), label='iLDSC')
        ldsc_patch = Patch(edgecolor=(228/rbgscaler,26/rbgscaler,28/rbgscaler,1), facecolor = (228/rbgscaler,26/rbgscaler,28/rbgscaler,0.5), label='LDSC')

        ax1.tick_params(axis='both', which='both', labelbottom=True)
        ax2.tick_params(axis='both', which='both', labelbottom=True)
        ax3.tick_params(axis='both', which='both', labelbottom=True)
        ax4.tick_params(axis='both', which='both', labelbottom=True)

        ax1.legend(handles=[ildsc_patch, ldsc_patch], ncol = 1, loc = 'upper right',bbox_to_anchor=(1.28, 1))
        sns.despine()
        plt.tight_layout()
        fig.subplots_adjust(right=0.8)

        plt.savefig(plot_dir + 'simulations.' + fig_dir + '.pdf')


    if fig_dir in ['GxAncestry', 'GxAncestry.pc.corrected', 'GxAncestry_noCASS', 'GxAncestry_noCASS.pc.corrected']:
        plot_params = {'GxAncestry':{'0.6':['0.3','0.15','0.15'], '0.3':['0.15','0.075','0.075']}, 'GxAncestry.pc.corrected':{'0.6':['0.3','0.15','0.15'], '0.3':['0.15','0.075','0.075']}, 'GxAncestry_noCASS':{'0.6':['0.3','0','0.3'], '0.3':['0.15','0','0.15']}, 'GxAncestry_noCASS.pc.corrected':{'0.6':['0.3','0','0.3'], '0.3':['0.15','0','0.15']}}
        fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(4, 1, sharex=True, figsize=(8,12))

        simdata = pd.read_pickle(results_dir + fig_dir + '.figure.data.pkl')

        simdata['meldh2_total_pve_diff'] = simdata['pve'] - simdata['total_meld'] 
        simdata['ldsch2_total_pve_diff'] = simdata['pve'] - simdata['ldsc_add'] 

        eigenvalues = pd.read_csv(results_dir + 'eigenvalues.txt', header = None)
        eigenvalues = eigenvalues.reset_index()
        eigenvalues['proportion'] = eigenvalues[0]/np.sum(eigenvalues[0])
        eigenvalues['index'] = eigenvalues['index'].astype(float)
        eigenvalues['index'] = eigenvalues['index'] + 1 
        # eigenvalues = np.array(eigenvalues)

        h2_06 = simdata[simdata['pve'] == 0.6]
        h2_06_total_meld = h2_06[['total_meld','pc']].rename(columns = {'total_meld':'estimate'})
        h2_06_total_meld['label'] = 'total_meld'

        h2_06_total_ldsc = h2_06[['ldsc_add','pc']].rename(columns = {'ldsc_add':'estimate'})
        h2_06_total_ldsc['label'] = 'total_ldsc'
        
        h2_06 = pd.concat([h2_06_total_meld,h2_06_total_ldsc])
        
        methodpalette3_dict = {1:'#e41a1c',0:'#377eb8'} 
        methodpalette3_rgba = {1:(228/rbgscaler,26/rbgscaler,28/rbgscaler,alpha), 0:(55/rbgscaler,126/rbgscaler,184/rbgscaler,alpha),}

        sns.boxplot(data = h2_06, x = 'pc', y = 'estimate', hue = 'label', ax = ax1, fliersize = 0)
        ax5 = ax1.twinx()

        # print(eigenvalues[:10,0])
        ploteigen = eigenvalues.head(n=10)[['index','proportion']]
        ax5.scatter(ploteigen['index']-1,ploteigen['proportion'], color = '#984ea3')
        
        fliersize = 0
        ax1.set_title("Generative model\n" + r"$h^2=0.6$, additive pve = " + plot_params[fig_dir]['0.6'][0] + ", CASS pve = " + plot_params[fig_dir]['0.6'][1] + ", GxAncestry pve = " + plot_params[fig_dir]['0.6'][2])
        ax1.set_ylabel(r"$\hat{h}^2$")
        ax1.set_xlabel(r"Principal Component")
        ax5.set_ylabel("Eigenvalue Proportion\nof variance", color = '#984ea3')

        ax1.legend().set_visible(False)
        

        counter = 1
        subcount = 0
        for j in ax1.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                j.set_facecolor(methodpalette3_rgba[subcount%2])
                j.set_edgecolor(methodpalette3_dict[subcount%2])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%2])
                counter+=1
        ax1.axhline(y=0.6,linestyle='--',color = 'black')
        h2_06 = simdata[simdata['pve'] == 0.6]
        h2_06_total_meld = h2_06[['meld_score','pc']]
        sns.boxplot(data = h2_06_total_meld, x = 'pc', y = 'meld_score', ax = ax2, fliersize = 0)

        counter = 1
        subcount = 0
        for j in ax2.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                j.set_facecolor(methodpalette3_rgba[subcount%1])
                j.set_edgecolor(methodpalette3_dict[subcount%1])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%1])
                counter+=1
        
        ax2.set_ylabel("iLDSC component\nestimate")
        ax2.set_xlabel(r"Principal Component")

        h2_03 = simdata[simdata['pve'] == 0.3]
        h2_03_total_meld = h2_03[['total_meld','pc']].rename(columns = {'total_meld':'estimate'})
        h2_03_total_meld['label'] = 'total_meld'

        h2_03_total_ldsc = h2_03[['ldsc_add','pc']].rename(columns = {'ldsc_add':'estimate'})
        h2_03_total_ldsc['label'] = 'total_ldsc'
        h2_03 = pd.concat([h2_03_total_meld,h2_03_total_ldsc,])
        
        methodpalette3_dict = {1:'#e41a1c',0:'#377eb8'} 
        methodpalette3_rgba = {1:(228/rbgscaler,26/rbgscaler,28/rbgscaler,alpha), 0:(55/rbgscaler,126/rbgscaler,184/rbgscaler,alpha),}

        sns.boxplot(data = h2_03, x = 'pc', y = 'estimate', hue = 'label', ax = ax3, fliersize = 0)
        ax6 = ax3.twinx()
        ax6.scatter(ploteigen['index']-1,ploteigen['proportion'], color = '#984ea3')
        ax6.set_ylabel("Eigenvalue Proportion\nof variance", color = '#984ea3')

        fliersize = 0
        ax1.set_title("Generative model\n" + r"$h^2=0.6$, additive pve = " + plot_params[fig_dir]['0.3'][0] + ", CASS pve = " + plot_params[fig_dir]['0.3'][1] + ", GxAncestry pve = " + plot_params[fig_dir]['0.3'][2])
        ax3.set_ylabel(r"$\hat{h}^2$")
        ax3.set_xlabel(r"Principal Component")
        ax3.legend().set_visible(False)
        

        counter = 1
        subcount = 0
        for j in ax3.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                j.set_facecolor(methodpalette3_rgba[subcount%2])
                j.set_edgecolor(methodpalette3_dict[subcount%2])
                subcount += 1
            
            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%2])
                counter+=1
        ax3.axhline(y=0.3,linestyle='--',color = 'black')

        h2_03 = simdata[simdata['pve'] == 0.3]
        h2_03_total_meld = h2_03[['meld_score','pc']]

        sns.boxplot(data = h2_03_total_meld, x = 'pc', y = 'meld_score', ax = ax4, fliersize = 0)

        counter = 1
        subcount = 0
        for j in ax4.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                
                j.set_facecolor(methodpalette3_rgba[subcount%1])
                j.set_edgecolor(methodpalette3_dict[subcount%1])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%1])
                counter+=1

        ax4.set_xlabel(r"Principal Component")
        ax4.set_ylabel("iLDSC component\nestimate")

        ildsc_patch = Patch(edgecolor=(55/rbgscaler,126/rbgscaler,184/rbgscaler,1), facecolor = (55/rbgscaler,126/rbgscaler,184/rbgscaler,0.5), label='iLDSC')
        ldsc_patch = Patch(edgecolor=(228/rbgscaler,26/rbgscaler,28/rbgscaler,1), facecolor = (228/rbgscaler,26/rbgscaler,28/rbgscaler,0.5), label='LDSC')
        eigen_dot = Line2D([0], [0], marker='o', color='#984ea3', label='Eigen pve',
                          markerfacecolor='#984ea3', markersize=15)

        ax1.legend(handles=[ildsc_patch, ldsc_patch], ncol = 1, loc = 'upper right',bbox_to_anchor=(1.3, 0.0))
        sns.despine(ax = ax1, top = True, right = False)
        sns.despine(ax = ax2, top = True, right = True)
        sns.despine(ax = ax3, top = True, right = False)
        sns.despine(ax = ax4, top = True, right = True)
        plt.tight_layout()
        fig.subplots_adjust(right=0.8)

        ax1.tick_params(axis='both', which='both', labelbottom=True)
        ax2.tick_params(axis='both', which='both', labelbottom=True)
        ax3.tick_params(axis='both', which='both', labelbottom=True)
        ax4.tick_params(axis='both', which='both', labelbottom=True)
        fig.subplots_adjust(hspace=0.2)
        plt.savefig(plot_dir + 'simulations.' + fig_dir + '.pdf')

    if fig_dir == 'ld2':

        fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(4, 1, sharex=True, figsize=(8,12))

        simdata = pd.read_pickle(results_dir + 'ld2.figure.data.pkl')
        simdata['percentage'] = np.where(simdata['l2sample'] == 'low',
                                           simdata['percentage'] * -1,
                                           simdata['percentage'])        

        h2_06 = simdata[simdata['pve'] == 0.6]
        h2_06_total_meld = h2_06[['total_meld','percentage']].rename(columns = {'total_meld':'estimate'})
        h2_06_total_meld['label'] = 'total_meld'

        h2_06_total_ldsc = h2_06[['ldsc_add','percentage']].rename(columns = {'ldsc_add':'estimate'})
        h2_06_total_ldsc['label'] = 'total_ldsc'
        
        h2_06 = pd.concat([h2_06_total_meld,h2_06_total_ldsc])
        
        methodpalette3_dict = {1:'#e41a1c',0:'#377eb8'} 
        methodpalette3_rgba = {1:(228/rbgscaler,26/rbgscaler,28/rbgscaler,alpha), 0:(55/rbgscaler,126/rbgscaler,184/rbgscaler,alpha),}

        sns.boxplot(data = h2_06, x = 'percentage', y = 'estimate', hue = 'label', ax = ax1, fliersize = 0)

        fliersize = 0
        ax1.set_title(r"$h^2=0.6$")
        ax1.set_ylabel(r"$h^2$")
        ax1.set_xlabel(r"Percentile of ld score with trait effect")
        ax1.legend().set_visible(False)
        

        counter = 1
        subcount = 0
        for j in ax1.get_children():
            if type(j) == matplotlib.patches.PathPatch and subcount >= 2:
                j.set_facecolor(methodpalette3_rgba[subcount%2])
                j.set_edgecolor(methodpalette3_dict[subcount%2])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%2])
                counter+=1

        ax1.axhline(y=0.6,linestyle='--',color = 'black')
        h2_06 = simdata[simdata['pve'] == 0.6]
        sns.boxplot(data = h2_06, x = 'percentage', y = 'meld_score', ax = ax2, fliersize = 0)

        counter = 1
        subcount = 0
        for j in ax2.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                j.set_facecolor(methodpalette3_rgba[subcount%1])
                j.set_edgecolor(methodpalette3_dict[subcount%1])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%1])
                counter+=1
        
        ax2.set_ylabel(r"iLDSC estimate")
        ax2.set_xlabel(r"Percentile of ld score with trait effect")
        
        h2_03 = simdata[simdata['pve'] == 0.3]
        h2_03_total_meld = h2_03[['total_meld','percentage']].rename(columns = {'total_meld':'estimate'})
        h2_03_total_meld['label'] = 'total_meld'

        h2_03_total_ldsc = h2_03[['ldsc_add','percentage']].rename(columns = {'ldsc_add':'estimate'})
        h2_03_total_ldsc['label'] = 'total_ldsc'
        
        h2_03 = pd.concat([h2_03_total_meld,h2_03_total_ldsc])
        
        methodpalette3_dict = {1:'#e41a1c',0:'#377eb8'} 
        methodpalette3_rgba = {1:(228/rbgscaler,26/rbgscaler,28/rbgscaler,alpha), 0:(55/rbgscaler,126/rbgscaler,184/rbgscaler,alpha),}

        sns.boxplot(data = h2_03, x = 'percentage', y = 'estimate', hue = 'label', ax = ax3, fliersize = 0)

        fliersize = 0
        ax3.set_title(r"$h^2=0.3$")
        ax3.set_ylabel(r"$h^2$")
        ax3.set_xlabel(r"Percentile of ld score with trait effect")
        ax3.legend().set_visible(False)

        counter = 1
        subcount = 0
        for j in ax3.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                j.set_facecolor(methodpalette3_rgba[subcount%2])
                j.set_edgecolor(methodpalette3_dict[subcount%2])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%2])
                counter+=1

        ax3.axhline(y=0.3,linestyle='--',color = 'black')
        h2_03 = simdata[simdata['pve'] == 0.3]
        
        sns.boxplot(data = h2_03, x = 'percentage', y = 'meld_score', ax = ax4, fliersize = 0)

        counter = 1
        subcount = 0
        for j in ax4.get_children():
            if type(j) == matplotlib.patches.PathPatch:
                j.set_facecolor(methodpalette3_rgba[subcount%1])
                j.set_edgecolor(methodpalette3_dict[subcount%1])
                subcount += 1

            if type(j) == matplotlib.lines.Line2D:
                j.set_color(methodpalette3_dict[(math.ceil(counter/float(6))-1)%1])
                counter+=1
        
        ax4.set_xlabel(r"Percentile of ld score with trait effect")
        ax4.set_ylabel(r"iLDSC estimate")

        ildsc_patch = Patch(edgecolor=(55/rbgscaler,126/rbgscaler,184/rbgscaler,1), facecolor = (55/rbgscaler,126/rbgscaler,184/rbgscaler,0.5), label='iLDSC')
        ldsc_patch = Patch(edgecolor=(228/rbgscaler,26/rbgscaler,28/rbgscaler,1), facecolor = (228/rbgscaler,26/rbgscaler,28/rbgscaler,0.5), label='LDSC')

        ax1.legend(handles=[ildsc_patch, ldsc_patch,], ncol = 1, loc = 'upper right',bbox_to_anchor=(1.28, 1))
        sns.despine()
        plt.tight_layout()
        fig.subplots_adjust(right=0.8)

        plt.savefig(plot_dir + 'simulations.ld2.pdf')








        










