import os
import numpy as np
import replaceline as rl

from scipy.constants import c as clight

tag_prefix = 'Qinj'
	
tobecopied = 'job.job launch_single beam.beam machine_parameters.input secondary_emission_parameters.input simulation_parameters.input LHC_chm_ver.mat'

target_storage='/eos/project/e/ecloud-simulations/kparasch/Quadrupole_initial_states/'
current_dir = os.getcwd()
study_folder =  current_dir.split('/config')[0]
scan_folder = study_folder+'/simulations'
#os.mkdir(scan_folder)
#os.mkdir(scan_folder+'/progress')

launch_file_lines = []
launch_file_lines +=['#!/bin/bash\n']


# scan parameters
blenf = lambda V_RF: 1.33 + (1.12-1.33)/(6.5-3.5)*(V_RF-3.5)
V_RF_vect_MeV = np.array([8.])
#V_RF_vect_MeV = np.arange(3., 8.1, 1)
V_RF_vect = V_RF_vect_MeV * 1e6
sigmat_vect_ns = blenf(V_RF_vect_MeV)/4
sigmat_vect = sigmat_vect_ns * 1e-9
sigmaz_vect = sigmat_vect * clight

sey_vect = [1.2]

#fact_beam_vect = np.arange(0.30e11,1.26e11,0.05e11)
fact_beam_vect=np.array([1.20e11,0.70e11])

prog_num = 0 
for iters in range(10):
    for fact_beam in fact_beam_vect:
    
        for sey in sey_vect:
    
            for sigmaz,V_RF in zip(sigmaz_vect,V_RF_vect):
                            
                prog_num +=1
                current_sim_ident='injection_ArcQuad_MQF_intensity_%.2fe11ppb_sey_%.2f_VRF_%dMV_%d'%(fact_beam/1e11,sey,V_RF/1e6,iters) 
                sim_tag = tag_prefix+'%03d'%prog_num
                print sim_tag, current_sim_ident
                current_sim_folder = scan_folder+'/'+current_sim_ident
                os.mkdir(current_sim_folder)
                
                
                rl.replaceline_and_save(fname = 'secondary_emission_parameters.input',
                 findln = 'del_max = ', newline = 'del_max = %.2f\n'%sey)            
                
                rl.replaceline_and_save(fname = 'beam.beam',
                 findln = 'sigmaz = ', newline = 'sigmaz = %e\n'%sigmaz)            
                
                rl.replaceline_and_save(fname = 'beam.beam',
                 findln = 'fact_beam = ', newline = 'fact_beam = %e\n'%fact_beam)                          
                        
                            
                rl.replaceline_and_save(fname = 'simulation_parameters.input',
                 findln = 'logfile_path =', newline = 'logfile_path = '+'\''+ current_sim_folder+'/logfile.txt'+'\''+'\n')
                 
                rl.replaceline_and_save(fname = 'simulation_parameters.input',
                 findln = 'progress_path =', newline = 'progress_path = '+'\'' + scan_folder+'/progress/' +sim_tag+'\''+ '\n')
            
                rl.replaceline_and_save(fname = 'simulation_parameters.input',
                 findln = 'stopfile =', newline = 'stopfile = \''+scan_folder+'/progress/stop\'\n')
                    
                rl.replaceline_and_save(fname = 'job.job',
                                     findln = 'CURRDIR=/',
                                     newline = 'CURRDIR='+current_sim_folder)
                                     
                rl.replaceline_and_save(fname = 'job.job',
                                     findln = 'TARGETDIR=/',
                                     newline = 'TARGETDIR='+target_storage+current_sim_ident)

                launch_lines = ['bsub -L /bin/bash -J '+ sim_tag + 
                                ' -o '+ current_sim_folder+'/STDOUT',
                                ' -e '+ current_sim_folder+'/STDERR',
                                ' -q %s < '%'1nw'+current_sim_folder+'/job.job\n', 'bjobs\n']
                
                
                with open('launch_single','w') as fid:
                    fid.writelines(launch_lines) 
    
                os.chmod('job.job',0755)
                os.chmod('launch_single',0755)
                                     
                os.system('cp -r %s %s'%(tobecopied, current_sim_folder))
                
                launch_file_lines += launch_lines

            
                                    
with open(study_folder+'/run_PyECLOUD', 'w') as fid:
    fid.writelines(launch_file_lines)
os.chmod(study_folder+'/run_PyECLOUD',0755)

import htcondor_config as htcc
htcc.htcondor_config(scan_folder, time_requirement_days=2.)

