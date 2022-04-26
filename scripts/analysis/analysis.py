# %%

import matplotlib.pyplot as plt

from shared import *

# %%

df = load_conus()
rename_cols(df)
add_derived_cols(df)
cols = [
    'P',
    'AI',
    #'NDVI',
    'EVI',
    'SSM',
    'SUSM',
    #'SMP',
    'logslope',
    'maxorder',
    'minorder'
]
df = df[['angle'] + cols]
standardize_non_angle(df)
# %%
