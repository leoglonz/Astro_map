# A simple script to generate the 'predischarged', 'discharged', and (dsrg) 'accreted' dataframes.
#
# __________________________
# Last revised: 3 Mar. 2022

import tqdm

from compiler import *



keys = ['h148_12','h148_27','h148_34','h148_38','h148_55','h148_65','h148_249',
        'h148_251','h148_282','h229_14','h229_18','h229_20','h229_22',
        'h229_49','h242_21','h242_38','h242_69','h329_29','h329_117']

print('Compiling sim gas into set (SNgas) for the following keys:', keys)

for key in keys:
    sim = str(key[:4])
    haloid = int(key[5:])
    sngas = calc_snGas(sim, haloid, save=True, verbose=False)