import pysixtrack
import pickle

with open('line.pkl','rb') as fid:
    line = pysixtrack.Line.from_dict(pickle.load(fid))


mags = ['mq', 'mb']
not_mags = ['mqt','mqs']
s_els = line.get_s_elements()
for ii in range(len(line.element_names)):
    name = line.element_names[ii]
    el = line.elements[ii]
    split_name = name.split('.')
    if len(split_name) == 3:
        if split_name[1][-2:] == 'l7':
            if (split_name[0][0:2] in mags) and (split_name[0][0:3] not in not_mags):
                print(ii,name, s_els[ii])
