import pickle

fLine = 'line_with_ecloud_markers_and_collimators.pkl'
fOptics = 'optics.pkl'
fPartCO = 'part_on_CO.pkl'
fEclouds = 'eclouds_info.pkl'

simulation_input = {}
with open(fLine, 'rb') as fid:
    simulation_input['line'] = pickle.load(fid)
with open(fOptics, 'rb') as fid:
    simulation_input['optics'] = pickle.load(fid)
with open(fPartCO, 'rb') as fid:
    simulation_input['partCO'] = pickle.load(fid)
with open(fEclouds, 'rb') as fid:
    simulation_input['eclouds_info'] = pickle.load(fid)

pickle.dump(simulation_input, open('simulation_input.pkl', 'wb'))
