import glob
import numpy as np
import datetime as dt

inAfiles = glob.glob('S2*L2A')

out_list_name = f'command_list_{dt.datetime.now().strftime("%Y%m%dT%H%M%S")}.txt'
min_delt = 5
max_delt = 25
max_v = 30 # m/d
half_source_chip = 5
modval = 5 # for S2 could be 5 (include AB and BA as well) or 10 (only AA and BB)

command_list = []
inAdts = [(dt.datetime.strptime(x.split('_')[2],"%Y%m%d"),x.split('_')[2],x) for x in inAfiles]
inAdts.sort()

for ind1 in range(len(inAdts)-1):
    dt1,dt1str,dir1 = inAdts[ind1]

    dts = [(x[0] - inAdts[ind1][0]).days for x in inAdts[ind1+1:]]

    dts_arr = np.array([(lambda x,y: x>= min_delt and x<=max_delt and y==0)(x,np.mod(x,modval)) for x in dts])

    for delt,(dt2,dt2str,dir2) in zip(np.array(dts)[dts_arr], np.array(inAdts[ind1+1:])[dts_arr] ):
        cmdstr = f'python /Users/mark/repos/pycorr_iceflow/pycorr_iceflow_v1.py -img1dir {dir1} -img2dir {dir2} B08.tif B08.tif -img1datestr {dt1str} -img2datestr {dt2str} -datestrfmt "%Y%m%d" -inc 10 -half_source_chip {half_source_chip}  -half_target_chip {np.ceil((delt*max_v/10)+half_source_chip).astype(np.int)} -plotvmax 20 -log10 -out_name_base S2{dir1[2]}{dir2[2]}_{dir1.split("_")[1]}_{delt:02d}_{dt1str}_{dt2str} -no_speed_ref -progupdates -offset_correction_lgo_mask -lgo_mask_filename AK_ocean_land_glacier_210_mask_100m_Kienholz_v3.tif -lgo_mask_file_dir ~/Dropbox'
        print(f'{cmdstr}')
        command_list.append(cmdstr)
        

with open(out_list_name,'w') as outf:
    for line in command_list:
        outf.write(f'{line}\n')
        
