#!/bin/env python3

import sys
import pandas as pd
import numpy as np

for f in sys.argv[1:]:
    df = pd.read_csv(f,sep='\t',index_col=0)
    ratio_cols = [ h for h in df.keys() if 'Log Ratio' in h ]
    shrc = [ h for h in ratio_cols if 'Unshared' not in h ]
    unshrc = [ h for h in ratio_cols if 'Unshared' in h ]
    dfsh = df[shrc]
    dfunsh = df[unshrc]
    print("Matrix:              ",f)
    print("Samples:             ",dfsh.shape[1])
    print("Genes:               ",dfsh.shape[0])
    # print("Genes per Sample Min:",np.min(np.sum(~np.isnan(dfsh),axis=0)))
    # print("Genes per Sample Max:",np.max(np.sum(~np.isnan(dfsh),axis=0)))
    # print("50% Genes:           ",np.sum((np.sum(np.isnan(dfsh),axis=1)<dfsh.shape[1]*.5)))
    print("80% Genes:           ",np.sum((np.sum(np.isnan(dfsh),axis=1)<dfsh.shape[1]*.2)))
    print("Full Genes:          ",np.sum((np.sum(np.isnan(dfsh),axis=1)==0)))
    print("#Nan:                ",np.sum(np.sum(np.isnan(dfsh))),"(%.2f%%)"%(100*np.sum(np.sum(np.isnan(dfsh)))/(dfsh.shape[1]*dfsh.shape[0]),))
    print("Min:                 ",np.nanmin(np.nanmin(dfsh)))
    print("Max:                 ",np.nanmax(np.nanmax(dfsh)))
    print()
