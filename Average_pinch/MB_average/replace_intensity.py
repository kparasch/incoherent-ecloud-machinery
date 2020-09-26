from replaceline import replaceline_and_save
import sys

fname = sys.argv[1]
int_e11ppb = sys.argv[2]
replaceline_and_save(fname, 
                     findln='intensity =', 
                     newline=f'intensity = {int_e11ppb}e+11'
                    )
