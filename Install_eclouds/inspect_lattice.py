import pysixtrack
import numpy as np
import pickle

with open('line.pkl','rb') as fid:
    line = pysixtrack.Line.from_dict(pickle.load(fid))


ss = '81'
sector = { '12' : ['r1','l2'],
           '23' : ['r2','l3'],
           '34' : ['r3','l4'],
           '45' : ['r4','l5'],
           '56' : ['r5','l6'],
           '67' : ['r6','l7'],
           '78' : ['r7','l8'],
           '81' : ['r8','l1']
         }

mags = ['mq', 'mb']
not_mags = ['mqt','mqs']
s_els = line.get_s_elements()
slist = []
for ii in range(len(line.element_names)):
    name = line.element_names[ii]
    el = line.elements[ii]
    split_name = name.split('.')
    mag_name = split_name[0]
    if len(split_name) == 3:
        halfsect = split_name[1][-2:]
        if halfsect in sector[ss]:
            if (mag_name[0:2] in mags) and (mag_name not in not_mags):
                if mag_name == 'mq':
                    if float(split_name[1][:-2]) > 10:
                        slist.append(s_els[ii])
                        print(el)
                print(ii,name, s_els[ii])

print(slist)
print(np.array(slist[1:] - np.array(slist[:-1])))
