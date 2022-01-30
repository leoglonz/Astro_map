# include `wait` in between commands to do them in batches (each command uses at most ~3.5 Gb of memory).
# run this in terminal via 'bash runall.sh'.

#!/bin/bash

cd ~/astro_research/Stellar_Feedback_Code/
pwd
date

python particletracking.py h242 41 & #1 
python particletracking.py h329 33 & #2
python particletracking.py h229 27 & #3
python particletracking.py h148 13 & #4
wait
python particletracking.py h229 20 & #5
python particletracking.py h242 24 & #6
python particletracking.py h229 22 & #7
wait
python particletracking.py h242 80 & #8
python particletracking.py h148 37 & #9
python particletracking.py h148 28 & #10 
python particletracking.py h148 68 & #11
wait
python particletracking.py h148 278 & #12
python particletracking.py h229 23 & #13
python particletracking.py h148 45 & #14
wait
python particletracking.py h148 283 & #15
python particletracking.py h229 55 & #16
python particletracking.py h329 137 & #17
wait
python particletracking.py h148 80 & #18
python particletracking.py h148 329 & #19
wait