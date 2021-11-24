# Archived script for compiling set of accreted gas particles; no longer used because it does not specifically screen for
# particles accreted following a discharge event.

# from analysis import *
from accretion import *

keys = ['h148_13','h148_28','h148_37','h148_45','h148_68','h148_80','h148_283',
        'h148_278','h148_329','h229_20','h229_22','h229_23','h229_27','h229_55',
        'h242_24','h242_41','h242_80','h329_33','h329_137']

print('Compiling sim gas into sets (advanced accreted) for the following keys:', keys)

for key in keys:
    sim = str(key[:4])
    haloid = int(key[5:])
    # note that heated is automatically concatenated by 'calc_ejected_expelled'.
    adv_accreted = calc_adv_accreted(sim, haloid, save=True, verbose=False)