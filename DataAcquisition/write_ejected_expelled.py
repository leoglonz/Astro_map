# A simple script to generate the 'ejected', 'cooled', 'expelled', and 'accreted' datasets.
#
# __________________________________________________________________________________________________________
# Code credit to Hollis Akins 2021;
# Github permalink: https://github.com/hollisakins/Justice_League_Code/blob/
#                    e049137edcfdc9838ebb3cf0fcaa4ee46e977cec/Analysis/RamPressure/write_ejected_expelled.py
# __________________________________________________________________________________________________________
# Last revised: 28 Nov. 2021

from compiler import *



keys = ['h148_13','h148_28','h148_37','h148_45','h148_68','h148_80','h148_283',
        'h148_278','h148_329','h229_20','h229_22','h229_23','h229_27','h229_55',
        'h242_24','h242_41','h242_80','h329_33','h329_137']

print('Compiling sim gas into sets (ejected, cooled, expelled, accreted) for the following keys:', keys)

for key in keys:
    sim = str(key[:4])
    haloid = int(key[5:])
    # note that ejected, cooled, etc. are automatically concatenated by 'calc_ejected_expelled'.
    ejected, cooled, expelled, accreted = calc_ejected_expelled(sim, haloid, save=True, verbose=False)
