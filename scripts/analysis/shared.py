from os.path import join
from numpy import log10, sqrt
from pandas import read_csv, read_pickle, cut

plotsdir = join(
    '..',
    '..',
    'plots'
)

datadir = join(
    '..',
    '..',
    'data'
)

fnconus = join(
    datadir,
    'exp_pro',
    'conus_angles.csv'
)

fnmars = join(
    datadir,
    'exp_pro',
    'mars_angles.csv'
)

def load_conus():
    try:
        df = read_pickle(fnconus.replace('.csv', '.pickle'))
    except FileNotFoundError:
        df = read_csv(fnconus)
        df.to_pickle(fnconus.replace('.csv', '.pickle'))
    return df

def load_mars():
    return read_csv(fnmars)

def add_logslope(df):
    a = df['slope A'].values
    b = df['slope B'].values
    a[a < 1e-5] = 1e-5
    b[b < 1e-5] = 1e-5
    df['logslope'] = log10(a)/2 + log10(b)/2
    return None

def add_maxorder(df):
    df['maxorder'] = [max(df.at[i,'order A'], df.at[i,'order B']) for i in df.index]
    return None

def add_minorder(df):
    df['minorder'] = [min(df.at[i,'order A'], df.at[i,'order B']) for i in df.index]
    return None

def add_derived_cols(df):
    add_logslope(df)
    add_maxorder(df)
    add_minorder(df)

def standardize_col(df, col):
    df[col] = (df[col] - df[col].mean())/df[col].std()
    return None

def standardize_non_angle(df, cols=None):
    if cols is None:
        cols = list(df.columns)
        if 'angle' in cols:
            cols.remove('angle')
    for col in cols:
        standardize_col(df, col)
    return None

def rename_col(df, fromcol, tocol):
    cols = list(df.columns)
    if fromcol in cols:
        cols[cols.index(fromcol)] = tocol
    df.columns = cols
    return None

def rename_cols(df):
    rename_col(df, 'ppt_annual', 'P')
    rename_col(df, 'tmean_annual', 'T')
    for col in df.columns:
        if 'ppt' in col:
            rename_col(df, col, col.replace('ppt', 'P'))
        if 'tmean' in col:
            rename_col(df, col, col.replace('tmean', 'T'))
    return None

def binmeanplot(ax, df, cols, x, bins=10):
    df = df.copy()
    if x == 'angle':
        df['angle_bin'] = cut(df.angle, bins=bins, right=True)
        df['bin_center'] = [(I.left + I.right)/2 for I in df.angle_bin]
        g = df.groupby('angle_bin')
        m = g.mean()
        s = g.std()/sqrt(g.count())
        for col in cols:
            ax.errorbar(
                m.angle,
                m[col],
                s[col],
                fmt='.-',
                label=col
            )
    elif x == 'cols':
        for col in cols:
            sl = df[['angle',col]].copy()
            sl['bin'] = cut(sl[col], bins=bins)
            sl['bin_center'] = [(I.left + I.right)/2 for I in sl['bin']]
            g = sl.groupby('bin_center')
            m = g.mean()
            s = g.std()/sqrt(g.count())
            ax.errorbar(
                m[col],
                m.angle,
                s.angle,
                fmt='.-',
                label=col
            )
    else:
        raise(StandardError('incorrect x argument'))

    return None
# %%
