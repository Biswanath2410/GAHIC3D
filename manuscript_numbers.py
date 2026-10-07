#every number quoted in the manuscript, recomputed from the result tables (python manuscript_numbers.py)
import os
import glob
import numpy as np
import pandas as pd
from scipy import stats
pd.set_option('display.width', 250)
ROOT = os.path.dirname(os.path.abspath(__file__))


def paired(a, b):
    k = a.index.intersection(b.index)
    a, b = a.loc[k], b.loc[k]
    ok = ~(a.isna() | b.isna())
    a, b = a[ok], b[ok]
    return len(a), a.median(), b.median(), int((a < b).sum()), stats.wilcoxon(a, b).pvalue if len(a) >= 5 else np.nan


def main():
    d = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(ROOT, 'results', 'evals', '*.csv'))], ignore_index=True)
    d['CHROMOSOME'] = d['CHROMOSOME'].astype(str)
    d['NB'] = d['TAG'].str.contains('nb')
    d['COV'] = d['TAG'].str.extract(r'cov([0-9.]+)')[0]
    extra = d.METHOD.str.contains(r'-(?:pop|relax|pivot|iter|a)[0-9]') | ~d.COV.isin(['1.0', '0.01'])
    d = d[~extra].copy()                      # run_extra.py results are summarised by plot_extra.py
    main = d[~d.NB]
    print('== Table 1: median RMSD (n)')
    t = main.groupby(['METHOD', 'RESOLUTION', 'COV'])['RMSD'].agg(['median', 'count']).round(3)
    print(t.unstack([1, 2]).to_string())
    print('\n== median Pearson r of distances')
    print(main.groupby(['METHOD', 'RESOLUTION', 'COV'])['PCoeff'].median().round(3).unstack([1, 2]).to_string())
    print('\n== paired tests GAHIC3D vs method (RMSD; n, med_G, med_o, G better, p)')
    for (r, c), sub in main.groupby(['RESOLUTION', 'COV']):
        g = sub[sub.METHOD == 'GAHIC3D'].set_index(['CELL', 'CHROMOSOME']).RMSD
        for m in sorted(set(sub.METHOD) - {'GAHIC3D'}):
            o = sub[sub.METHOD == m].set_index(['CELL', 'CHROMOSOME']).RMSD
            n, mg, mo, w, p = paired(g, o)
            if n:
                print(r, c, m, n, round(mg, 3), round(mo, 3), w, '{0:.2g}'.format(p))
    print('\n== multistart by chromosome (400 kb)')
    for (c, ch), sub in main[(main.RESOLUTION == 400000)].groupby(['COV', 'CHROMOSOME']):
        g = sub[sub.METHOD == 'GAHIC3D'].set_index('CELL').RMSD
        o = sub[sub.METHOD == 'MultiStart'].set_index('CELL').RMSD
        print(c, ch, paired(g, o))
    print('\n== robustness (alpha -3 + NB), 1 Mb')
    rb = d[d.NB]
    print(rb.groupby(['METHOD', 'COV'])['RMSD'].median().round(3).unstack().to_string())
    for c, sub in rb.groupby('COV'):
        g = sub[sub.METHOD == 'GAHIC3D'].set_index(['CELL', 'CHROMOSOME']).RMSD
        for m in sorted(set(sub.METHOD) - {'GAHIC3D'}):
            o = sub[sub.METHOD == m].set_index(['CELL', 'CHROMOSOME']).RMSD
            n, mg, mo, w, p = paired(g, o)
            print('robust', c, m, n, round(mg, 3), round(mo, 3), w, '{0:.2g}'.format(p))
    print('\n== mirrored fraction')
    mm = main[main.METHOD.isin(['GAHIC3D', 'pastis-pm2', 'pastis-pm1', 'pastis-mds', 'lordg', 'shrec3d', 'minimds',
                                'hsa1', 'chromsde', 'chromosome3d', 'gem'])]
    fr = (mm.RMSD_noreflect > mm.RMSD + 1e-6)
    print(fr.groupby(mm.METHOD).mean().round(2).to_string(), '\noverall', round(fr.mean(), 3))
    ratio = (mm.RMSD_noreflect / mm.RMSD)[fr]
    print('inflation when mirrored: median', round(ratio.median(), 2), 'max', round(ratio.max(), 2))
    print('\n== runtime (s) median')
    rt = pd.read_csv(os.path.join(ROOT, 'results', 'tables', 'runtime_raw.csv'))
    print(rt.groupby(['METHOD', 'RESOLUTION'])['seconds'].median().round(1).unstack().to_string())
    # real data
    r = pd.read_csv(os.path.join(ROOT, 'results', 'real', 'real_all.csv'))
    print('\n== real A (median Pearson, Spearman, n)')
    A = r[r.ANALYSIS == 'A']
    print(A.groupby(['RESOLUTION', 'METHOD'])[['PCoeff', 'SCoeff']].agg(['median', 'count']).round(3).to_string())
    print('\n== real B (median AUC, n)')
    B = r[r.ANALYSIS == 'B']
    print(B.groupby(['RESOLUTION', 'METHOD'])['AUC'].agg(['median', 'count']).round(3).to_string())
    print('\n== real C')
    C = r[r.ANALYSIS == 'C']
    print(C.groupby(['RESOLUTION', 'METHOD'])['SCoeff_mean'].median().round(3).unstack(0).to_string())
    ma = pd.read_csv(os.path.join(ROOT, 'results', 'real', 'model_agreement.csv'))
    bc = pd.read_csv(os.path.join(ROOT, 'results', 'real', 'between_cell_agreement.csv'))
    print('models agree', ma.groupby('RESOLUTION').PCoeff_models.median().round(3).to_dict(),
          'between cells', bc.groupby('RESOLUTION').PCoeff.median().round(3).to_dict())
    st = pd.read_csv(os.path.join(ROOT, 'data', 'real', 'sc_matrix_stats.tsv'), sep='\t')
    print(st.groupby(['res', 'chrom'])[['contacts', 'bins_with_contacts']].median())
    fails = [os.path.basename(f) for f in glob.glob(os.path.join(ROOT, 'results', 'real', 'evals', '*.failed'))]
    print('failures', pd.Series([f.split('_')[0] + '_' + f.split('_')[1] for f in fails]).value_counts().to_dict())


if __name__ == '__main__':
    main()


def native_pastis_common_bins():
    """Supplementary Note 1: PASTIS run natively (raw counts + biases, 4 % bins filtered),
    compared with GAHIC3D on the bins PASTIS kept"""
    from RMSD import compare
    from simulateHiC import load_groundtruth
    rows = []
    for f in glob.glob(os.path.join(ROOT, 'results', 'evals', 'pastis-pm*-native_*.csv')):
        r = pd.read_csv(f).iloc[0]
        cond = os.path.join(ROOT, 'results', 'structures', r.TAG, '{0}_{1}_{2}'.format(r.CELL, r.CHROMOSOME, r.RESOLUTION))
        nat = np.loadtxt(os.path.join(cond, r.METHOD + '_struct.txt'))
        keep = nat[:, 0].astype(int)
        gi, g = load_groundtruth(os.path.join(ROOT, 'data', 'groundtruth', 'groundtruth_{0}_{1}_{2}.txt'.format(
            r.CELL, r.CHROMOSOME, r.RESOLUTION)), str(r.CHROMOSOME), str(r.RESOLUTION))
        e = np.loadtxt(os.path.join(cond, 'GAHIC3D_struct.txt'))
        sel = np.isin(e[:, 0].astype(int), keep)
        rows.append([r.METHOD, r.RESOLUTION, r.TAG, r.RMSD, compare(e[sel, 0].astype(int), e[sel, 1:4], gi, g)['RMSD']])
    t = pd.DataFrame(rows, columns=['METHOD', 'RES', 'TAG', 'RMSD_native', 'RMSD_GAHIC3D_same_bins'])
    for k, s in t.groupby(['METHOD', 'RES', 'TAG']):
        p = stats.wilcoxon(s.RMSD_GAHIC3D_same_bins, s.RMSD_native).pvalue
        print(k, len(s), 'native', round(s.RMSD_native.median(), 3), 'GAHIC3D', round(s.RMSD_GAHIC3D_same_bins.median(), 3),
              'GAHIC3D better', int((s.RMSD_GAHIC3D_same_bins < s.RMSD_native).sum()), '{0:.2g}'.format(p))
    t.to_csv(os.path.join(ROOT, 'results', 'tables', 'tableS3_pastis_native_common_bins.csv'), index=False)


if __name__ == '__main__':
    print('\n== native PASTIS, common bins')
    native_pastis_common_bins()
