# A simple script to generate the 'predischarged', 'discharged', and (dsrg) 'accreted' dataframes.
#
# __________________________
# Last revised: 1 Dec. 2021

import tqdm

from compiler import *



keys = ['h148_13','h148_28','h148_37','h148_45','h148_68','h148_80','h148_283',
        'h148_278','h148_329','h229_20','h229_22','h229_23','h229_27','h229_55',
        'h242_24','h242_41','h242_80','h329_33','h329_137']

print('Compiling sim gas into sets (predischarged, discharged, accreted) for the following keys:', keys)

for key in keys:
    sim = str(key[:4])
    haloid = int(key[5:])
    # note that predischarged, discharged, etc. are automatically concatenated.
    predischarged, discharged, dsrg_accreted = calc_discharged(sim, haloid, save=True, verbose=False)