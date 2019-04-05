from __future__ import print_function
import netCDF4
import numpy as np
import gdal
import os
import sys
from matplotlib.pyplot import get_cmap
import re
import argparse
import glob

# copy and modify netCDF file to change offset correction to dc x,y offset over land pixels
# and add velocity field divergence layer

############################################################################################
# from accepted answer here: https://stackoverflow.com/questions/15141563/python-netcdf-making-a-copy-of-all-variables-and-attributes-but-one
# import netCDF4 as nc  #(this appears to be wrong - "netCDF4" used in code below)
# toexclude = ['ExcludeVar1', 'ExcludeVar2']
# 
# with netCDF4.Dataset("in.nc") as src, netCDF4.Dataset("out.nc", "w") as dst:
#     # copy global attributes all at once via dictionary
#     dst.setncatts(src.__dict__)
#     # copy dimensions
#     for name, dimension in src.dimensions.items():
#         dst.createDimension(
#             name, (len(dimension) if not dimension.isunlimited() else None))
#     # copy all file data except for the excluded
#     for name, variable in src.variables.items():
#         if name not in toexclude:
#             x = dst.createVariable(name, variable.datatype, variable.dimensions)
#             dst[name][:] = src[name][:]
#             # copy variable attributes all at once via dictionary
#             dst[name].setncatts(src[name].__dict__)
########
# NOTE: had to rearrange the above code to deal with case where variable had no shape (no values)
#       this occurs in cf1.6 files where a variable's attributes are used to store all the information
#       also chose to modify (back) to python 2 version of iteritems... 
#
#       See below for modified code
############################################################################################
default_indir = '.'
default_outdir = '.'
default_speed_max = 15.0
default_cam = 1.0
default_dcam = 0.1
default_cam1 = 0.5

make_output_prdir = False
make_log10 = True

parser = argparse.ArgumentParser( \
    description="""make_new_nc_dcoff_and_divergence.py""",
    epilog='>>  <<',
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-input_dir', 
                    action='store', 
                    type=str, 
                    default=default_indir,
                    help = 'input directory for nc file [%(default)s]')
parser.add_argument('-output_dir', 
                    action='store', 
                    type=str, 
                    default=default_outdir,
                    help='output directory for nc, tif and png files [%(default)s]')
parser.add_argument('-plt_speed_max', 
                    action='store', 
                    type=float, 
                    default = default_speed_max,
                    help='max speed on plots [%(default)f]')
parser.add_argument('-cam', 
                    action='store', 
                    type=float, 
                    default = default_cam,
                    help='default cam [%(default)f]')
parser.add_argument('-dcam', 
                    action='store', 
                    type=float, 
                    default = default_dcam,
                    help='default dcam [%(default)f]')
parser.add_argument('-cam1', 
                    action='store', 
                    type=float, 
                    default = default_cam1,
                    help='default cam1 [%(default)f]')
# parser.add_argument('-span_tree_make_missing_dirs', 
#                     action='store_true', 
#                     default = False,
#                     help='find all nc files in p_r or S2_tile dirs in specified input_dir, create output dirs as needed in output_dir [%(default)s]')
parser.add_argument('-make_missing_output_dir', 
                    action='store_true', 
                    default = False,
                    help='make output_dir if it is missing [%(default)s]')
parser.add_argument('-replace_existing_output', 
                    action='store_true', 
                    default = False,
                    help='replace existing output nc file if present (False skips this file) [%(default)s]')
parser.add_argument('input_nc_file', 
                    action='store', 
                    type=str,
                    default=None,
                    help='nc file to process [None]')
args = parser.parse_args()



cam = args.cam
dcam = args.dcam
cam1 = args.cam1




variables_to_exclude = [     
                u'offset_correction',
                u'applied_bilinear_x_offset_correction_in_pixels',
                u'applied_bilinear_y_offset_correction_in_pixels',
                u'vx',
                u'vy',
                u'vv',
                u'vx_masked',
                u'vy_masked',
                u'vv_masked',
                u'applied_bilinear_x_offset_correction_in_pixels',
                u'applied_bilinear_y_offset_correction_in_pixels',
                u'offset_correction'
            ]

# files_to_process = []

indir = args.input_dir
outdir = args.output_dir
if not(os.path.exists(outdir)):
    if args.make_missing_output_dir:
        os.makedirs(outdir)
    else:
        raise ValueError('output directory {} does not exist, and -make_missing_output_dir not set on command line'.format(outdir))

in_nc_file = args.input_nc_file
out_nc_file = in_nc_file.replace('_hp','_dcd') # dcd => dc offset, divergence

if not(args.replace_existing_output) and os.path.exists(outdir + '/' + out_nc_file):
    print('file: {} already present, skipping.  Use -replace_existing_output to force replace.'.format(outdir + '/' + out_nc_file))
    sys.exit(0)

# if not(args.span_tree_make_missing_dirs):  # only one file to process, make it the sole member of the list
#     in_nc_file = args.input_nc_file
#     out_nc_file = in_nc_file.replace('_hp','_dcd') # dcd => dc offset, divergence
#     files_to_process.append( (indir, in_nc_file, outdir, out_nc_file) )
# else:  # traverse tree of input and output files_to_process

src = netCDF4.Dataset(indir + '/' + in_nc_file)
dst = netCDF4.Dataset(outdir + '/' + out_nc_file, "w", clobber=True, format='NETCDF4')

##############  See above comment block - following is based on that stackoverflow code example -
##############  modified to deal with empty variables that are used only for metadata fields (projection in cf1.6 file for example)
# copy global attributes all at once via dictionary
dst.setncatts(src.__dict__)
# copy dimensions
for name, dimension in src.dimensions.items():
    dst.createDimension(
        name, (len(dimension) if not dimension.isunlimited() else None))
# copy all file data except for the excluded
for name, variable in src.variables.iteritems():
    if name not in variables_to_exclude:
        x = dst.createVariable(name, variable.datatype, variable.dimensions, **variable.filters())
        # copy variable attributes all at once via dictionary
        dst[name].setncatts(src[name].__dict__)
        if src[name].shape != ():
            dst[name][:] = src[name][:]




lgo_masked_offset_available = None
found_valid_offset = False
################################
# values from sc_pycorr_v5p10
lgo_masked_min_num_pix=500
lgo_masked_min_percent_valid_pix_available=0.05
max_allowable_pixel_offset_correction = 3.0

corr_val_for_offset=0.3
delcorr_min_for_offset = 0.15

vel_nodata_val=np.float32(-9999.0)

plotvmax = args.plt_speed_max      # 15.0
#
################################

regex = r"(-inc) (?P<incval>\d+)"
bb = re.search(regex, src.getncattr('history'))  # looking in processing history for the increment value used by sc_pycorr to get grid spacing in original image pixels
inc_val_pixels = float(bb.groupdict()['incval']) # this will be the spacing between chip centers in the original images, in pixels in those images

grid_mapping = src['vv'].getncattr('grid_mapping')  # used below when writing new velocity grids to output nc (projection)
del_t_speedunit_str = src['vv'].getncattr('units')  # used below when writing new velocity grids to output nc (m/d typically)

corr_arr = src['corr'][:]
corr_nodata_val = src['corr'].getncattr('_FillValue')
del_corr_arr = src['del_corr'][:]
del_i = src['del_i'][:]
del_j = src['del_j'][:]

if 'lgo_mask' in src.variables.keys():
    lgo_mask_image_utm = src['lgo_mask'][:]
    
    if True:  # if args.v:
        print(out_nc_file + ': attempting to find offset correction for land (lgo mask) areas')
    # mask pixels that are not land not(lgo==1)
    if '_FillValue' in src['lgo_mask'].ncattrs():  # input speed reference has a specified no_data value.
        lgo_mask_nodata = src['lgo_mask'].getncattr('_FillValue')
        
        lgo_masked_d_i_m=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset) | \
                                         (lgo_mask_image_utm!=1)|(lgo_mask_image_utm==lgo_mask_nodata),del_i)
        lgo_masked_d_j_m=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset) | \
                                         (lgo_mask_image_utm!=1)|(lgo_mask_image_utm==lgo_mask_nodata),del_j)
        lgo_masked_num_possible_pix=np.count_nonzero(np.array((lgo_mask_image_utm==1)&(lgo_mask_image_utm!=lgo_mask_nodata)&(corr_arr!=corr_nodata_val), dtype=np.bool))
        #             lgo_masked_num_valid_pix=np.count_nonzero(np.array((speedref_vel_vv<slow_area_zero_speed_md)&(speedref_vel_vv!=speedref_vv_nodata)&(del_i!=0.0)&(prelim_vv_ma2.mask==False), dtype=np.bool))
        lgo_masked_num_valid_pix=np.count_nonzero(np.array(lgo_masked_d_i_m.mask==False, dtype=np.bool))
    else:   # lgo mask did not include a nodata value (so don't try to use it...)
        lgo_masked_d_i_m=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset) | \
                                         (lgo_mask_image_utm!=1),del_i)
        lgo_masked_d_j_m=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset) | \
                                         (lgo_mask_image_utm!=1),del_j)
        lgo_masked_num_possible_pix=np.count_nonzero(np.array((lgo_mask_image_utm==1)&(corr_arr!=corr_nodata_val), dtype=np.bool))
        #             lgo_masked_num_valid_pix = np.count_nonzero(np.array((speedref_vel_vv<slow_area_zero_speed_md) & (speedref_vel_vv!=speedref_vv_nodata) & (del_i!=0.0) & (prelim_vv_ma2.mask==False), dtype=np.bool))
        lgo_masked_num_valid_pix=np.count_nonzero(np.array(lgo_masked_d_i_m.mask==False, dtype=np.bool))
        
    if lgo_masked_num_valid_pix>0:
        lgo_masked_offset_i = -(np.median(lgo_masked_d_i_m.compressed()))
        lgo_masked_offset_j = -(np.median(lgo_masked_d_j_m.compressed()))
        lgo_masked_stdev_i = np.std(lgo_masked_d_i_m.compressed())
        lgo_masked_stdev_j = np.std(lgo_masked_d_j_m.compressed())
        lgo_masked_offset_available=True

        if True: # args.v:
            print('found lgo_mask (land pixel) offset correction (del_i: {0} del_j: {1} std_i {2} std_j {3}) using {4} pixels out of {5} possible'.format(
                    lgo_masked_offset_i,lgo_masked_offset_j,lgo_masked_stdev_i,lgo_masked_stdev_j,lgo_masked_num_valid_pix,lgo_masked_num_possible_pix))
else: # no lgo_mask available - how to do offset?                    
    print('No land/ice/ocean mask available in nc file, exiting...')
    sys.exit(0)
               
if lgo_masked_offset_available  and \
        (lgo_masked_num_valid_pix>lgo_masked_min_num_pix) and \
        ((float(lgo_masked_num_valid_pix)/lgo_masked_num_possible_pix) >= \
        lgo_masked_min_percent_valid_pix_available/100.0) and \
        (np.sqrt(lgo_masked_offset_i**2.0 + lgo_masked_offset_j**2.0) < \
        max_allowable_pixel_offset_correction):
    # use lgo_masked corection
    final_offset_correction_i=lgo_masked_offset_i
    final_offset_correction_j=lgo_masked_offset_j
    found_valid_offset=True
    offset_correction_type_applied='lgo_masked_correction'
    offset_correction_type_descritption='lgo_masked_correction for land pixels, %d valid pixels out of %d possible for scene (%f %%)'%(lgo_masked_num_valid_pix,lgo_masked_num_possible_pix,100.0*(np.float(lgo_masked_num_valid_pix)/lgo_masked_num_possible_pix))
else:
    offset_correction_type_applied='None'
    offset_correction_type_descritption='None'
          
############################################################################        
# apply offset correction if there is one
############################################################################

# dcam=args.dcam
# cam=args.cam
# cam1=args.cam1

# -dcam 0.1 -cam 1.0 -cam1 0.5 typical parameters for GoLIVE processing
# greenland -dcam 0.05 -cam 1.0 -cam1 0.0
# dcam = 0.1
# cam = 1.0
# cam1 = 0.2

image_pix_x_m = src['input_image_details'].getncattr('image1_pix_x_size')
image_pix_y_m = src['input_image_details'].getncattr('image1_pix_y_size')
del_t_val = float(src['image_pair_times'].getncattr('del_t'))

# if not(args.offset_correction_speedref or args.offset_correction_lgo_mask) or not(found_valid_offset):  # no offset to apply
if not(found_valid_offset):  # no offset to apply
    print('Not enough valid pixels to apply offset correction - proceeding without one')
    # create masked and unmasked versions of the speed, vx, and vy output
    # the masked version will be written to the nc file with no_data values in corners and also where the correlation parameters suggest masking is needed
#     vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1)|(((offset_dist_ij_arr*img1.pix_x_m)/del_t_val)>plotvmax),((offset_dist_ij_arr*img1.pix_x_m)/del_t_val))
    vx=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),((del_i*image_pix_x_m)/del_t_val))
    vx[corr_arr==corr_nodata_val] = vel_nodata_val
    vy=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),((del_j*image_pix_y_m)/del_t_val))
    vy[corr_arr==corr_nodata_val] = vel_nodata_val
    vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),np.sqrt(np.square(vx) + np.square(vy)))
    vv[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run

    # nomask arrays only have nodata values where there was no velocity attempted (filled corners of landsat scenes, for example)
    vx_nomask=((del_i*image_pix_x_m)/del_t_val)
    vx_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
    vy_nomask=((del_j*image_pix_y_m)/del_t_val)
    vy_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
    vv_nomask=np.sqrt(np.square(vx_nomask) + np.square(vy_nomask))
    vv_nomask[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run

else:  # found_valid_offset (have single i and j offsets to apply)
    # apply offset and create masked and unmasked versions of the speed, vx, and vy output
    # the masked version will be written to the nc file with no_data values in corners and also where the correlation parameters suggest masking is needed
    vx=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),(((del_i + final_offset_correction_i) * image_pix_x_m)/del_t_val))
    vx[corr_arr==corr_nodata_val] = vel_nodata_val
    vy=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),(((del_j + final_offset_correction_j) * image_pix_y_m)/del_t_val))
    vy[corr_arr==corr_nodata_val] = vel_nodata_val
#     vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1)|(np.sqrt(np.square(vx) + np.square(vy))>plotvmax),np.sqrt(np.square(vx) + np.square(vy)))
    vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),np.sqrt(np.square(vx) + np.square(vy)))
    vv[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run

    # nomask arrays only have nodata values where there was no velocity attempted (filled corners of landsat scenes, for example)
    vx_nomask=(((del_i + final_offset_correction_i) * image_pix_x_m)/del_t_val)
    vx_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
    vy_nomask=(((del_j + final_offset_correction_j) * image_pix_y_m)/del_t_val)
    vy_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
    vv_nomask=np.sqrt(np.square(vx_nomask) + np.square(vy_nomask))
    vv_nomask[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run
                            




####################################
#
# calculate divergence
#
####################################

# vx_arr[vx_arr==vx.nodatavalue]=np.nan
# vy_arr=vy.img
# vy_arr[vy_arr==vy.nodatavalue]=np.nan

# vx and vy are masked arrays, but use no_data_value too - need all these pixels to be np.nan for following calculation
# so copy them to new arrays and reset the nodata vals to nans in those - won't get output to nc file that way...
nvx = vx.copy()
nvx[vx==vel_nodata_val] = np.nan
nvy = vy.copy()
nvy[vy==vel_nodata_val] = np.nan

divergence = np.ones_like(nvx) * vel_nodata_val
strain_rate_distance = 2.0 * image_pix_x_m * inc_val_pixels # distance between chip centers in meters * 2 for derivative calculated below
num_y,num_x = nvx.shape

if 'lgo_mask' in src.variables.keys():
    ice_indices = np.where(lgo_mask_image_utm == 0)
else: # no mask, calc at all pixels...  Won't get here until offset with no lgo mask part is put in above.
    ice_indices = np.where(nvx != vel_nodata_val)
for j,i in zip(*ice_indices):
    if ( (i > 0) & (i < (num_x - 2)) & (j > 0) & (j < (num_y - 2)) ):
        if ( (lgo_mask_image_utm[j,i]==0) & (lgo_mask_image_utm[j,i-1]==0)  & (lgo_mask_image_utm[j,i+1]==0)  & (lgo_mask_image_utm[j-1,i]==0)  & (lgo_mask_image_utm[j+1,i]==0) ):  # if pixel and surrounding are ice
            divergence[j,i] = ( ((nvx[j,i+1] - nvx[j,i-1])/strain_rate_distance) + (np.sign(image_pix_y_m)*(nvy[j+1,i] - nvy[j-1,i])/strain_rate_distance) ) # np.sign to get gradient in vy calculated in proper direction
# divergence[divergence == vel_nodata_val] = np.nan
# divergence[divergence == np.nan] = vel_nodata_val  # reset the nans that made it through calculation to nodata_val for nc output
divergence_ma = np.ma.masked_where(divergence == vel_nodata_val, divergence)




# done with src nc file, close it and reopen as gdal image to get geotransform for gtiff output if log10 is used
# output nc file will get its transform from input nc file. Otherwise, to get geotransform from nc file, you have to know what
# the projection is in advance and then seek the metadata block named after the projection
src.close()


mapjet=get_cmap('jet')
mmvv=vv.copy()
mmvv[(mmvv>plotvmax) & (~mmvv.mask)]=plotvmax
vv_zo=mmvv/plotvmax
vv_rgba=mapjet(vv_zo) # this produces a 4-band 0.0-1.0 image array using the jet colormap
(out_lines,out_pixels,out_bands)=vv_rgba.shape        # NEED these values below, EVEN if only output tif is log10

        
if (make_log10):    # prepare a log10 version of the speed for output below
    mmvv=vv.copy()
    mmvv[(mmvv>plotvmax) & (~mmvv.mask)]=plotvmax
    mmvv.mask[(mmvv==0) & (~mmvv.mask)]= True
    # ignore warning for invalid points (0 or <0) - will be set to nans
    with np.errstate(invalid='ignore', divide='ignore'):
        lmmvv=np.log10(mmvv)
#     min_lmmvv=np.min(lmmvv)
    min_lmmvv=np.log10(0.1)
    max_lmmvv=np.log10(plotvmax)
    range_lmmvv=max_lmmvv - min_lmmvv
    lvv_zo=(lmmvv - min_lmmvv)/(range_lmmvv)
    lvv_rgba=mapjet(lvv_zo) # this produces a 4-band 0.0-1.0 image array using the jet colormap





##########################################################################################
#
# output section - write various output files
#
##########################################################################################


# check if output pr directory exists - if -make_output_prdir flag is set, create a missing pr directory, otherwise fail...



if (make_output_prdir):
    if not(os.path.isdir(outdir)):
        if not(os.path.isdir(os.path.dirname(outdir))):  # make sure the parent directory (one holding all the pr directories) exists
            print(out_name_base + ': ' + 'Error: >>>>>>>>>>>>>>>parent directory {}} does not exist - halting'.format(os.path.dirname(outdir)))
            sys.exit(0)
        else:
            os.makedirs(outdir)
            print(out_name_base + ': ' + '>>>>>>>>>>>>>>>created directory %s '.format(outdir))
            os.chmod(outdir, 0o755)
        
# if not(args.no_gtif):   ###### This flag is not presently used, as GTiff forms basis for all output right now...
format = "GTiff"
driver = gdal.GetDriverByName( format )
metadata = driver.GetMetadata()

pngdriver = gdal.GetDriverByName( 'PNG' )


if(make_log10):
    gd = gdal.Open('NETCDF:' + indir + '/' + in_nc_file + ':vv')  # to get geotransform and other projection info from it, no data
#     dst_filename = outdir + '/' + file_name_base + '_log10.tif'
    dst_filename = outdir + '/' + out_nc_file.replace('.nc','.tif')
    dst_ds = driver.Create( dst_filename, out_pixels, out_lines, out_bands, gdal.GDT_Byte )

    # dst_ds.SetGeoTransform( [ com_min_x, inc * img1.pix_x_m, 0, com_max_y, 0, inc * img1.pix_y_m ] ) # note pix_y_m typically negative
    dst_ds.SetGeoTransform( gd.GetGeoTransform() ) # note pix_y_m typically negative
    dst_ds.SetProjection( gd.GetProjection() )
    dst_ds.GetRasterBand(1).WriteArray( (lvv_rgba[:,:,0]*255).astype('ubyte') )
    dst_ds.GetRasterBand(2).WriteArray( (lvv_rgba[:,:,1]*255).astype('ubyte') )
    dst_ds.GetRasterBand(3).WriteArray( (lvv_rgba[:,:,2]*255).astype('ubyte') )
    dst_ds.GetRasterBand(4).WriteArray( (lvv_rgba[:,:,3]*255).astype('ubyte') )
#     dst_ds = None # done, close the dataset


############
#   Now create PNG as well, using CreateCopy - gdal png driver does not support Create at this time...
    dst_filename1 = outdir + '/' + out_nc_file.replace('.nc','.png')
    dst_ds1 = pngdriver.CreateCopy( dst_filename1, dst_ds) # copy geotiff to png



    dst_ds1 = None # done, close the dataset
    dst_ds = None # done, close the dataset
#     make output files rw-r--r-- and remove unneeded png.aux.xml file that gdal creates
    os.chmod(dst_filename, 0o644)
    os.chmod(dst_filename1, 0o644)
    os.remove(dst_filename1 + '.aux.xml')                        


####################################################
#
# now add new (corrected) velocity layers to nc file
#
####################################################

varname='vx'
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = dst.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','x_velocity')
var.setncattr('long_name','x component of velocity')
var.setncattr('units',del_t_speedunit_str)
var[:] = vx_nomask.astype('float32')

varname='vy'
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = dst.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','y_velocity')
var.setncattr('long_name','y component of velocity')
var.setncattr('units',del_t_speedunit_str)
var[:] = vy_nomask.astype('float32')

varname='vv'
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = dst.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','speed')
var.setncattr('long_name','magnitude of velocity')
var.setncattr('units',del_t_speedunit_str)
var[:] = vv_nomask.astype('float32')

varname='vx_masked'
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = dst.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','x_velocity_masked')
var.setncattr('long_name','x component of velocity (masked)')
var.setncattr('units',del_t_speedunit_str)
var.setncattr('masking_info','masked_where(((del_corr_arr<%4.3f)&(corr_arr<%4.3f))|(corr_arr<%4.3f))'%(dcam,cam,cam1))
var[:] = vx.astype('float32')

varname='vy_masked'
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = dst.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','y_velocity_masked')
var.setncattr('long_name','y component of velocity (masked)')
var.setncattr('units',del_t_speedunit_str)
var.setncattr('masking_info','masked_where(((del_corr_arr<%4.3f)&(corr_arr<%4.3f))|(corr_arr<%4.3f))'%(dcam,cam,cam1))
var[:] = vy.astype('float32')

varname='vv_masked'
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = dst.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','speed_masked')
var.setncattr('long_name','magnitude of velocity (masked)')
var.setncattr('units',del_t_speedunit_str)
var.setncattr('masking_info','masked_where(((del_corr_arr<%4.3f)&(corr_arr<%4.3f))|(corr_arr<%4.3f))'%(dcam,cam,cam1))
var[:] = vv.astype('float32')

varname='divergence_masked'
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = dst.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','divergence')
var.setncattr('long_name','divergence of velocity (masked)')
var.setncattr('units',del_t_speedunit_str.replace('m','1'))
var.setncattr('masking_info','masked_where(((del_corr_arr<%4.3f)&(corr_arr<%4.3f))|(corr_arr<%4.3f))'%(dcam,cam,cam1))
var[:] = divergence_ma.astype('float32')

# set variables
# first set up variable not as dimension, but as holder for attributes for the times of the two input images, delt, etc.
varname='offset_correction'
datatype=np.dtype('S1')
dimensions=()
FillValue=None
var = dst.createVariable(varname,datatype, dimensions, fill_value=FillValue)
if found_valid_offset:
    var.setncattr('offset_correction_type_applied',offset_correction_type_applied)
    var.setncattr('offset_correction_type_descritption',offset_correction_type_descritption)
    var.setncattr('final_offset_correction_units','pixels')
    var.setncattr('final_offset_correction_i',str(final_offset_correction_i))
    var.setncattr('final_offset_correction_j',str(final_offset_correction_j))
else:
    var.setncattr('offset_correction_type_applied',offset_correction_type_applied)
    var.setncattr('offset_correction_type_descritption',offset_correction_type_descritption)
    var.setncattr('final_offset_correction_units','pixels')
    var.setncattr('final_offset_correction_i','None')
    var.setncattr('final_offset_correction_j','None')


dst.close()


