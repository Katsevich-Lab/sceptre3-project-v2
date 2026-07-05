import scipy.io, pandas as pd, numpy as np, geomux, os, sys
def main():
    base='results/barnyard_cohort_export'
    for s in ['Cropseq_mix0hr','Cropseq_mix72hr','DirectCapture_mix0hr','DirectCapture_mix72hr']:
        d=os.path.join(base,s)
        m=scipy.io.mmread(d+'/guide_counts.mtx').tocsr()
        cells=pd.read_csv(d+'/cells.csv'); guides=pd.read_csv(d+'/guides.csv')
        X=np.asarray(m.T.todense())
        gx=geomux.Geomux(matrix=X, cell_names=cells['barcode'].values, guide_names=guides['guide'].values, n_jobs=1)
        gx.test()
        for tag,lor in [('default',10.0),('nolor',0.0)]:
            a=gx.assignments(pvalue_threshold=0.05, lor_threshold=lor)
            rows=[(r['cell_id'],g) for _,r in a.iterrows() if isinstance(r['assignment'],(list,np.ndarray)) for g in r['assignment']]
            pd.DataFrame(rows,columns=['cell_barcode','guide']).to_csv(f'{d}/geomux_{tag}_calls.csv',index=False)
            print(f'{s} {tag}: {len(rows)} calls',flush=True)
if __name__=='__main__':
    main()
