import scipy.io, pandas as pd, numpy as np, geomux, os
def main():
    d='results/comprehensive_bench'; outd=d+'/geomux_pr'; os.makedirs(outd, exist_ok=True)
    m=scipy.io.mmread(d+'/guide_counts.mtx').tocsr()
    cells=pd.read_csv(d+'/cells.csv'); guides=pd.read_csv(d+'/guides.csv')
    X=np.asarray(m.T.todense())
    gx=geomux.Geomux(matrix=X, cell_names=cells['barcode'].values, guide_names=guides['guide'].values, n_jobs=1)
    gx.test()
    for q in [0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.35,0.5,0.7]:
        a=gx.assignments(pvalue_threshold=q, lor_threshold=0.0)   # lor=0: pure p-value curve
        rows=[(r['cell_id'],g) for _,r in a.iterrows() if isinstance(r['assignment'],(list,np.ndarray)) for g in r['assignment']]
        pd.DataFrame(rows,columns=['cell_barcode','guide']).to_csv(f'{outd}/q{q}.csv',index=False)
        print(f'q={q}: {len(rows)} calls',flush=True)
if __name__=='__main__': main()
