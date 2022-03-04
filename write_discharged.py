# A simple script to generate the 'predischarged', 'discharged', and (dsrg) 'accreted' dataframes.
#
# __________________________
# Last revised: 3 Mar. 2022

import tqdm

from compiler import *



# keys = ['h148_13','h148_28','h148_37','h148_45','h148_68','h148_80','h148_283',
#         'h148_278','h148_329','h229_20','h229_22','h229_23','h229_27','h229_55',
#         'h242_24','h242_41','h242_80','h329_33','h329_137','h242_12','h242_30',
#         'h242_44','h148_3','h148_14','h148_36']

keys = ['h148_12','h148_27','h148_34','h148_38','h148_55','h148_65','h148_249',
        'h148_251','h148_282','h229_14','h229_18','h229_20','h229_22',
        'h229_49','h242_21','h242_38','h242_69','h329_29','h329_117']

print('Compiling sim gas into sets (predischarged, discharged, accreted) for the following keys:', keys)

for key in keys:
    sim = str(key[:4])
    haloid = int(key[5:])
    # note that predischarged, discharged, etc. are automatically concatenated.
    predischarged, discharged, accreted = calc_discharged(sim, haloid, save=True, verbose=False)