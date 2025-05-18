# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 11:06:59 2023

@author: danie
"""
import os
from pathlib import Path
import numpy as np
import scipy.stats as sts
import cobra
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"
from sklearn.linear_model import LinearRegression as LR
import scipy.stats as sts
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

class GeneExpr:
    
    def __init__(self, 
                 geneFolder,
                 tpmFile,
                 groupLabels,
                 groupIDX,
                 groupComparison,
                 featuresFile,
                 sbmlModel,
                 rootID = 'ncbi',
                 species = 'bt'):
        
        self.rootID = rootID
        self.species = species
        self.geneFolder = geneFolder
        self.groupLabels = groupLabels
        self.groupIDX = groupIDX #in tpm file
        self.groupIDX_ = []
        counter = -1
        for idxs in self.groupIDX:
            a = []
            for c in idxs:
                counter+=1
                a.append(counter)
            self.groupIDX_.append(a[:])
        self.groupComparison = groupComparison
        self.genes, self.tpm = self.__parseTPM(tpmFile)
        self.__parsePatricFeatures(featuresFile)
        self.means = self.__getMeans(self.tpm)
        self.powerM, self.powerSTD = self.__getPowerT(self.tpm)
        self.pvalsDeseq = self.__getPvalsDeseq()
        self.fc = self.__getFCDeseq()
        self.genes2reactions, self.reactions2genes = self.__parseModelGenes(sbmlModel)
        
        
    def __parseTPM(self, tpmFile):
        with open(os.path.join(self.geneFolder, tpmFile)) as f:
            d = {}
            genes = []
            f.readline()
            for line in f:
                a = line.strip().split('\t')
                
                if self.species=='bt':
                    idx = 0
                    
                elif self.species == 'bh':
                    idx = 1
                
                elif self.species == 'ri':
                    idx=1
                
                genes.append(a[idx])
                d[a[idx]] = {}
                for i, group in enumerate(self.groupLabels):
                    #the patric id is assumed at position 1 of the tpm file
                    d[a[idx]][group] = np.array(list(map(self.floatTry, a[self.groupIDX[i][0]: self.groupIDX[i][-1] + 1])))
        return np.array(genes), d
    
    def __parsePatricFeatures(self, featuresFile):
        
        ncbi2patric = {}
        patric2function = {}
        
        with open(os.path.join(self.geneFolder, featuresFile)) as f:
            f.readline()
            for line in f:
                a = line.strip().split('\t')
                
                ncbi2patric[a[6].replace('"', '')] = a[5].replace('"', '')
                patric2function[a[5].replace('"', '')] = a[12]
        self.ncbi2patric = ncbi2patric
        self.patric2function = patric2function
            
    def __parseModelGenes(self, model):
        g2r = {i:[] for i in self.genes}
        reacs = {}
        
        for reaction in model.reactions:
            reacs[reaction.id] = []
            if reaction.genes != frozenset():
                for gene in reaction.genes:
                    
                    if self.species=='bt':
                        geneID = gene.id
                    else:
                        geneID = 'fig|' + gene.id
                    
                    g2r[geneID].append(reaction.id)
                    reacs[reaction.id].append(geneID)
        return g2r, reacs
        
        
        
    def __getGeneMatrix(self, geneD):
        m = []
        
        for gene in self.genes:
            m.append(np.concatenate([geneD[gene][group] for group in self.groupLabels]))
            
        return np.array(m)
    
    def __getGeneDict(self, geneM):
        d = {}
        
        for gi, gene in enumerate(self.genes):
            d[gene] = {}
            for gri, group in enumerate(self.groupLabels):
                d[gene][group] = np.array(geneM[gi][self.groupIDX_[gri][0]: self.groupIDX_[gri][-1] + 1])
                
        return d
    
                
    
    def __getMeans(self, geneD):
        
        return {group:np.array([np.mean(geneD[gene][group]) for gene in self.genes]) for group in self.groupLabels}
    
    def __getStds(self, geneD):
        
        return {group:np.array([np.std(geneD[gene][group]) for gene in self.genes]) for group in self.groupLabels}
    
    
    def __getPowerT(self, geneD):
        expr = self.__getGeneMatrix(geneD).T
        trExpr = np.array([sts.yeojohnson(i)[0] for i in expr]).T
        
        trD = self.__getGeneDict(trExpr)
        
        return self.__getMeans(trD), self.__getStds(trD)
        
    def __getPvalsDeseq(self):
        pvals = {gene:{comp:np.nan for comp in self.groupComparison} for gene in self.genes}
        self.notInPatric = []
        
        for comp in self.groupComparison:
            with open(os.path.join(self.geneFolder, self.groupComparison[comp])) as f:
                f.readline()
                for line in f:
                    a = line.strip().split('\t')
                    geneID = a[0].replace('"', '')
                    
                    if self.rootID == 'ncbi':
                        if geneID not in self.ncbi2patric:
                            self.notInPatric.append(geneID)
                        else:
                            if abs(self.floatTry(a[2]))>1:
                                pvals[self.ncbi2patric[geneID]][comp] = self.floatTry(a[-1], 1.0)
                            else:
                                pvals[self.ncbi2patric[geneID]][comp] = 1
                    else:
                        
                        pvals[geneID][comp] = self.floatTry(a[-1], 1.0)
                        
                        # else:
                        #     pvals[geneID][comp] = 1
    
    
        return pvals
    
    
    def __getFCDeseq(self):
        
        fc = {gene:{comp:np.nan for comp in self.groupComparison} for gene in self.genes}
        self.notInPatric = []
        
        for comp in self.groupComparison:
            with open(os.path.join(self.geneFolder, self.groupComparison[comp])) as f:
                f.readline()
                for line in f:
                    a = line.strip().split('\t')
                    geneID = a[0].replace('"', '')
                    
                    if self.rootID == 'ncbi':
                        if geneID not in self.ncbi2patric:
                            pass
                        else:
                            if abs(self.floatTry(a[2]))>1:
                                fc[self.ncbi2patric[geneID]][comp] = self.floatTry(a[2])
                            else:
                                fc[self.ncbi2patric[geneID]][comp] = 0
                    else:
                        
                        fc[geneID][comp] = self.floatTry(a[2])
                        
                        # else:
                        #     pvals[geneID][comp] = 1
    
    
        return fc
        
    
        
        
                    
    
    @staticmethod
    def floatTry(string, r=0.0):
        try:
            return float(string)
        
        except:
            return r






def extracReactions_fisher(exprObj, reactionList, group):
    rxn_fc, rxn_p = [], []

    for reac in reactionList:
        genes = exprObj.reactions2genes.get(reac, [])
        if genes:
            # collect gene‐level stats
            fcs = []; pvs = []
            for g in genes:
                fcs.append(exprObj.fc[g][group])
                pvs.append(exprObj.pvalsDeseq[g][group])
            fcs = np.array(fcs); pvs = np.array(pvs)

            # median FC
            rxn_fc.append(np.median(fcs))

            # Fisher’s method to combine p-values
            # χ² = –2 ∑ ln p_i ; df = 2·k
            chi2 = -2.0 * np.sum(np.log(np.clip(pvs, 1e-300, 1.0)))
            df   = 2 * len(pvs)
            rxn_p.append(sts.chi2.sf(chi2, df))
        else:
            # no genes → no change
            rxn_fc.append(0.0)
            rxn_p.append(1.0)

    return np.array(rxn_fc), np.array(rxn_p)


def extracReactions_stouffer(exprObj, reactionList, group):
    rxn_fc, rxn_p = [], []

    for reac in reactionList:
        genes = exprObj.reactions2genes.get(reac, [])
        if genes:
            fcs = []; zs = []; ws = []
            for g in genes:
                fc = exprObj.fc[g][group]
                p  = exprObj.pvalsDeseq[g][group]
                # convert p to z (two-sided)
                z = sts.norm.isf(p / 2.0)
                w = 1.0  # or 1/np.var(fc_estimate) if known
                fcs.append(fc)
                zs.append(z * np.sign(fc))  # keep direction
                ws.append(w)

            fcs = np.array(fcs); zs = np.array(zs); ws = np.array(ws)
            # weighted mean fold-change
            rxn_fc.append(np.average(fcs, weights=ws))

            # Stouffer’s combined z
            z_comb = zs.dot(ws) / np.sqrt(np.sum(ws**2))
            p_comb = 2 * sts.norm.sf(abs(z_comb))
            rxn_p.append(p_comb)
        else:
            rxn_fc.append(0.0)
            rxn_p.append(1.0)

    return np.array(rxn_fc), np.array(rxn_p)

def extracReactions_median(exprObj, reactionList, group):
    rxn_fc, rxn_p = [], []
    for reac in reactionList:
        genes = exprObj.reactions2genes.get(reac, [])
        if genes:
            fcs = [ exprObj.fc[g][group]             for g in genes ]
            pvs = [ exprObj.pvalsDeseq[g][group]     for g in genes ]
            rxn_fc.append(np.median(fcs))
            rxn_p.append(np.median(pvs))
        else:
            rxn_fc.append(0.0); rxn_p.append(1.0)
    return np.array(rxn_fc), np.array(rxn_p)


def extracReactions(exprObj, reactionList, group):
    x, p, g = [],[],[]
    
    for reac in reactionList:
        genes = np.array(exprObj.reactions2genes[reac])
        
        if len(genes)>0:
            
            
            pv = []
            fc = []
            gids = []
            
            
            for gene in genes:
                gidx = np.arange(len(exprObj.genes))[exprObj.genes==gene]
                fc.append(exprObj.fc[gene][group])
                pv.append(exprObj.pvalsDeseq[gene][group])
                gids.append(gene)
            
            fc = np.array(fc).flatten()
            pv = np.array(pv).flatten()
            gids = np.array(gids).flatten()
            
            x.append(fc[np.abs(pv)==min(np.abs(pv))][0])
            p.append(pv[np.abs(pv)==min(np.abs(pv))][0])
            g.append(gids[np.abs(pv)==min(np.abs(pv))][0])
        
        else:
            x.append(0)
            p.append(1)
            g.append('NF')

    return np.array(x), np.array(p)#, np.array(g)
    
    
    

            
#example usage
geneFolder = os.path.join(Path(os.getcwd()).parents[0], 'files', 'bh', 'genes')
modelFolder = os.path.join(Path(os.getcwd()).parents[0], 'files', 'bh', 'model') 
model = cobra.io.read_sbml_model(os.path.join(modelFolder, 'bh_final.xml'))

bhGE = GeneExpr(geneFolder = geneFolder, 
                tpmFile = 'bh_tpm.txt', 
                groupLabels = ['t14', 't32', 't72'], 
                groupIDX = [[4,5,6],
                            [7,8,9],
                            [10,11,12]
                            ], 
                groupComparison = {('t14', 't32'):'t14vst32_deseq.txt',
                                   ('t14', 't72'): 't14vst72_deseq.txt',
                                   ('t32', 't72'): 't32vst72_deseq.txt'}, 
                featuresFile =  'bh_BVBRC_features.txt', 
                sbmlModel = model,
                species='bh'
                )   




reactionList = [i.id for i in model.reactions]

wc_reactions = []

with open(os.path.join(geneFolder, 'wcReactions.txt')) as f:
    for line in f:
        wc_reactions.append(line.strip())


reactionList = np.array(reactionList)



group = ('t14', 't72')

x,p = extracReactions(bhGE, reactionList, group)



# 1) Build DataFrame
df = pd.DataFrame({
    'reaction':    reactionList,
    'fold_change': x,
    'p_value':     p
})
df['neg_log10_p'] = -np.log10(df['p_value'])

# 2) Significance & regulation
df['regulation'] = 'Not significant'
sig = df['p_value'] < 0.01
df.loc[sig & (df['fold_change'] >  0), 'regulation'] = 'Upregulated'
df.loc[sig & (df['fold_change'] <  0), 'regulation'] = 'Downregulated'

# 3) Highlight wc_reactions
df['highlight'] = df['reaction'].isin(wc_reactions)

# 4) Create a combined category for coloring
def cat(r, h):
    if h and r=='Upregulated':   return 'WC Up'
    if h and r=='Downregulated': return 'WC Down'
    return r

df['cat'] = df.apply(lambda row: cat(row['regulation'], row['highlight']), axis=1)

# 5) Pull reaction equations
eqs = []
for rid in df['reaction']:
    try:
        eq = model.reactions.get_by_id(rid).build_reaction_string(use_metabolite_names=True)
    except KeyError:
        eq = ''
    eqs.append(eq)
df['equation'] = eqs

# 6) Jitter
np.random.seed(0)
jitter_x = 0.5  # adjust as needed
jitter_y = 0.5
df['x_jit'] = df['fold_change'] + np.random.uniform(-jitter_x, jitter_x, size=len(df))
df['y_jit'] = df['neg_log10_p'] + np.random.uniform(-jitter_y, jitter_y, size=len(df))

# 7) Marker sizing
df['size'] = np.where(df['highlight'], 14, 6)

# 8) Color map
color_map = {
    'Upregulated':     '#3A68AE',    # original blue
    'Downregulated':   '#44A043',    # original green
    'Not significant': 'lightgray',
    'WC Up':           '#00CCFF',    # neon cyan/blue
    'WC Down':         '#39FF14'     # neon green
}

# 9) Plotly volcano
fig = px.scatter(
    df,
    x='x_jit',
    y='y_jit',
    color='cat',
    color_discrete_map=color_map,
    size='size',
    hover_name='reaction',
    hover_data={'equation':True, 'p_value':True, 'size':False},
    labels={
      'x_jit': 'Log₂ Fold Change',
      'y_jit': '-log₁₀(p-value)'
    },
    title=f"Volcano Plot: {group[1]} vs {group[0]}"
)

# 10) Threshold lines
fig.add_hline(y=-np.log10(0.01), line_dash='dash', line_color='gray',
              annotation_text="p = 0.01", annotation_position='top left')
fig.add_vline(x=0, line_dash='dash', line_color='black')

fig.update_traces(marker=dict(opacity=0.8), selector=dict(mode='markers'))
fig.update_layout(legend_title_text='Regulation/Highlight', template='simple_white')



# … after your px.scatter and any fig.update_layout calls …

# # 12) Overlay black stars for special reactions
# special = ['TREpts', 'GLCabc']
# for reac in special:
#     row = df[df['reaction'] == reac]
#     if not row.empty:
#         x0 = float(row['x_jit'])
#         y0 = float(row['y_jit'])
#         fig.add_trace(go.Scatter(
#             x=[x0], y=[y0],
#             mode='markers',
#             marker=dict(
#                 symbol='star',
#                 size=30,           # bigger so it’s obvious
#                 color='black',     # solid black
#                 line=dict(width=1, color='white')  # white outline for contrast
#             ),
#             showlegend=False,
#             hovertemplate=(
#                 f"<b>{reac}</b><br>"
#                 "Log₂ FC: %{x:.2f}<br>"
#                 "-log₁₀(p): %{y:.2f}<extra></extra>"
#             )
#         ))




# 11) Save
out_file = "volcano_plot_jitter_14_72.html"
fig.write_html(out_file, include_plotlyjs='cdn')
print(f"🚀 Saved interactive plot to {out_file}")

