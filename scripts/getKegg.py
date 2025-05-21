#import gurobipy
#from gurobipy import Model
import cobra
import numpy as np
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'scripts'))



from parseGenExpData import *

#path with files
bh_geneFolder = os.path.join(os.path.join(Path(os.getcwd()).parents[0],'files', 'bh', 'genes'))
bh_modelFolder = os.path.join(os.path.join(Path(os.getcwd()).parents[0],'files', 'bh', 'model'))
bh_gsmm = cobra.io.read_sbml_model(os.path.join(bh_modelFolder, 'bh_final.xml'))

bhGE_obj = GeneExpr(geneFolder = bh_geneFolder,
                tpmFile = 'bh_tpm.txt',
                groupLabels = ['t14', 't32', 't72'],
                groupIDX = [[4,5,6],
                            [7,8,9],
                            [10,11,12]
                            ],
                groupComparison = {('14h', '32h'):'t14vst32_deseq.txt',
                                   ('14h', '72h'): 't14vst72_deseq.txt',
                                   ('32h', '72h'): 't32vst72_deseq.txt'},
                featuresFile =  'bh_BVBRC_features.txt',
                sbmlModel = bh_gsmm,
                species='bh'
                )