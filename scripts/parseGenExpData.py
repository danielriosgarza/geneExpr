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
    
    
    
def analyze_group(groupA, groupB, geneExprObj, gsmm, reactionList, wc_reactions, species_name):
    group = (groupA, groupB)
    x, p = extracReactions(geneExprObj, reactionList, group)

    # 1) Build DataFrame
    df = pd.DataFrame({
        'reaction':    reactionList,
        'fold_change': x,
        'p_value':     p
    })
    df['neg_log10_p'] = -np.log10(df['p_value'])

    # 2) Significance & regulation
    df['regulation'] = 'Not significant'
    sig = df['p_value'] < 0.05
    df.loc[sig & (df['fold_change'] >  0), 'regulation'] = f'OverExpressed({groupA}_vs_{groupB})'
    df.loc[sig & (df['fold_change'] <  0), 'regulation'] = f'UnderExpressed({groupA}_vs_{groupB})'

    # 3) Highlight wc_reactions
    df['In_WC_model'] = df['reaction'].isin(wc_reactions)

    # 4) Combined category for coloring
    def cat(reg, in_wc):
        if in_wc and reg.startswith('OverExpressed'):
            return 'WC OverExpressed'
        if in_wc and reg.startswith('UnderExpressed'):
            return 'WC UnderExpressed'
        return reg
    df['cat'] = df.apply(lambda row: cat(row['regulation'], row['In_WC_model']), axis=1)

    # 5) Pull reaction equations
    eqs = []
    for rid in df['reaction']:
        try:
            eq = gsmm.reactions.get_by_id(rid) \
                      .build_reaction_string(use_metabolite_names=True)
        except Exception:
            eq = ''
        eqs.append(eq)
    df['equation'] = eqs

    # 6) Jitter
    np.random.seed(0)
    df['x_jit'] = df['fold_change'] + np.random.uniform(-0.5, 0.5, size=len(df))
    df['y_jit'] = df['neg_log10_p'] + np.random.uniform(-0.5, 0.5, size=len(df))

    # 7) Marker sizing: large for ALL wc_reactions, small otherwise
    df['size'] = np.where(df['In_WC_model'], 14, 6)

    # 8) Neon color map
    color_map = {
        f'OverExpressed({groupA}_vs_{groupB})': '#3A68AE',
        f'UnderExpressed({groupA}_vs_{groupB})': '#44A043',
        'Not significant':                     'lightgray',
        'In WC model OverExpressed({groupA}_vs_{groupB})':                    '#00BFFF',
        'In WC model UnderExpressed({groupA}_vs_{groupB})':                   '#39FF14'
    }

    # 9) Plotly volcano
    fig = px.scatter(
        df,
        x='x_jit', y='y_jit',
        color='cat',
        color_discrete_map=color_map,
        size='size',
        size_max=16,
        hover_name='reaction',
        hover_data={'equation':True, 'p_value':True},
        labels={'x_jit':'Log₂ Fold Change','y_jit':'-log₁₀(p-value)'},
        title=f"Volcano Plot {species_name}: {group[0]} vs {group[1]}"
    )
    fig.add_hline(y=-np.log10(0.05), line_dash='dash', line_color='gray',
                  annotation_text="p = 0.05", annotation_position='top left')
    fig.add_vline(x=0, line_dash='dash', line_color='black')
    fig.update_traces(marker=dict(opacity=0.8))
    fig.update_layout(legend_title_text='Regulation/Highlight', template='simple_white')
    
    fig.add_annotation(
    dict(
        xref='paper', yref='paper',
        x=1.02, y=0.95,   # position just to the right of the top‐right corner
        text="<b>Marker size</b><br>"
             "Large = in WC model<br>"
             "Small = not in WC model",
        showarrow=False,
        align='left',
        font=dict(size=12, color='black'),
        bordercolor="black",
        borderwidth=1,
        bgcolor="white",
        opacity=0.8
            )
        )

    # 10) Display inline
    fig.show()

    # 11) Return both figure and a slimmed‐down table
    df_out = df[['reaction','equation','fold_change','p_value','In_WC_model']].copy()
    return fig, df_out


