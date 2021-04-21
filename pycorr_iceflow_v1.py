import numpy as np
from numpy.polynomial import polynomial    # this is for the bilinear offset correction
from scipy.interpolate import RectBivariateSpline as RBS
import cv2
import gdal
import gdalconst as gdc  # constants for gdal - e.g. GA_ReadOnly, GA_Update ( http://www.gdal.org )
import os
import subprocess as sp
import string
import random
import sys
import time
import datetime as dt
import argparse
import osr
from scipy.ndimage.filters import gaussian_filter
from matplotlib.pyplot import get_cmap
import netCDF4
import re

# pycorr_iceflow_v1 started life as experimental_sentinel2_sc_pycorr_nc_oc_v5p8p1_py3.py on 4/21/2021
#   added -lgo_mask_limit_land_offset so that by default land and ocean pixels would not be search limited - good for fast tidewater and surging slow over oceans

# v5p0 is an update to v4 that adds masks/offset correction for Alaska and eventually Greenland using glacier/land/ocean masks
# the point0 (p0) is to indicate that this isn't complete yet - will remove that designation when it works right...

# 5p1 got most of this done - first working lgo mask update - limited search box to 5 pixels for land.

# 5p2 added save of lgo mast to nc file, limited search box to 5 pixels for oceans as well as land.

# 5p3 added bilinear offset correction of i and j offsets for land pixels, which improved warping-caused signal;
#  compression in nc file may have been added here or prototyped in 5p2-derived version earlier. Compression here has issues with gdal
#  because of bottom-up writing of datasets

# 5p4 fixed GDAL bottom-up issue by writing data and y coordinates out to nc file in top-down, with compression.

# 5p5 just frees the memory from the hp images before it does the offset correction, so that the process memory footprint doesn't get any bigger

# 5p6 adds the possibility of use of BOTH an lgo mask and velocity mask for offset correction (called full_stationary for lack of a better term) - this is envisioned now ONLY FOR GREENLAND NOW - maybe more places later - ALSO REQUIRES new -Greenland flag set.

# 5p7 the changes to implement the goal of 5p6 (adding Greenland speed ref vv,vx,vy offset correction like that used in Antarctica) got complicated - so when I started on figuring out
#		how to reproject the vx and vy components for mid-speed area offset correction, I created this new version number so I wouldn't step all over 5p6 code...
# 

# 5p8 this now has Greenland offset correction, changes to output nc file (includes del_i and del_j, as well as bilinear offset correction fields, if they were used, 
# 		and basic info in attributes about input files and original pixel/projection data)
#		also fixed ref_velocity nodata issue - didn't make it through reprojection.
#

# psutil is to trace memory footprint, but is not on sc system, so remove these references there
try:
	import psutil
	psutil_available=True
	print('psutil found')
	def memory_usage_psutil():
		# return the memory usage in MB
		process = psutil.Process(os.getpid())
		mem = process.memory_info()[0] / float(2 ** 20)
		return mem
except:
	psutil_available=False
	print('psutil not found, proceeding without psutil memory usage reports')

# try:
# 	import resource
# 	resource_available=True
# 	print('resource found')
# 	def memory_usage_resource():
# 		# return the memory usage in MB
# 		usage = resource.getrusage(resource.RUSAGE_SELF)
# 		return('%7.3f local %7.3f stack %7.3f max memory'%(float(usage.ru_idrss)/float(2 ** 20),float(usage.ru_isrss)/float(2 ** 20),float(usage.ru_maxrss)/float(2 ** 20)))
# except:
# 	resource_available=False
# 	print('resource package not found, proceeding without resource memory usage reports')

resource_available=False


##########################################################################################
# GeoImg_noload
###############
#
# modified 9/7/2015 for LO8 fname -> date    
# modified to noload for sc application - don't read image on setup - will read image data in main code to hp filter, delete...and keep memory footprint small
#
# modified 7/10/2017 for Collection 1 file names - recognize input filename as C1, parse first yyyyMMDD for acquisition date, set doy field...
#
# modified 4/4/2018 GeoImg_noload checks for successful file open, exits with status(0) but prints an error line if open fails
#
##########################################################################################
class GeoImg_noload:  
    """geocoded image input and info
        a=GeoImg(in_file_name,indir='.', datestr=None, datefmt='%m/%d/%y')
            a.img will contain image
            a.parameter etc...
            datefmt is datetime format string dt.datetime.strptime()"""
    
    # LXSS_LLLL_PPPRRR_YYYYMMDD_yyymmdd_CC_TX, in which:
    # L = Landsat
    # X = Sensor
    # SS = Satellite
    # PPP = WRS path
    # RRR = WRS row
    # YYYYMMDD = Acquisition date
    # yyyymmdd = Processing date
    # CC = Collection number
    # TX = Collection category

    def __init__(self, in_filename,in_dir='.',datestr=None,datefmt='%m/%d/%y'):
        self.filename = in_filename
        self.in_dir_path = in_dir  #in_dir can be relative...
        self.in_dir_abs_path=os.path.abspath(in_dir)  # get absolute path for later ref if needed
        self.gd=gdal.Open(self.in_dir_path + os.path.sep + self.filename)
        if not(self.gd):
            print('Error: open of %s failed: gdal_error: %s'%(self.in_dir_path + os.path.sep + self.filename, gdal.GetLastErrorMsg()))
            sys.exit(0)
        self.nodata_value=self.gd.GetRasterBand(1).GetNoDataValue()
        self.srs=osr.SpatialReference(wkt=self.gd.GetProjection())
        self.gt=self.gd.GetGeoTransform()
        self.proj=self.gd.GetProjection()
        self.intype=self.gd.GetDriver().ShortName
        self.min_x=self.gt[0]
        self.max_x=self.gt[0]+self.gd.RasterXSize*self.gt[1]
        self.min_y=self.gt[3]+self.gt[5]*self.gd.RasterYSize
        self.max_y=self.gt[3]
        self.pix_x_m=self.gt[1]
        self.pix_y_m=self.gt[5]
        self.num_pix_x=self.gd.RasterXSize
        self.num_pix_y=self.gd.RasterYSize
        self.XYtfm=np.array([self.min_x,self.max_y,self.pix_x_m,self.pix_y_m]).astype('float')
        if (datestr is not None):   # date specified in GeoImg call directly - could be any GeoTiff...
            self.imagedatetime=dt.datetime.strptime(datestr,datefmt)
        elif (self.filename.count('_')>=7 and self.filename[0]=='L'): # looks like collection 1 landsat
            b=self.filename.split('_')
            self.sensor=b[0]
            self.path=int(b[2][0:3])
            self.row=int(b[2][3:6])
            self.year=int(b[3][0:4])
            self.imagedatetime=dt.datetime.strptime(b[3],'%Y%m%d')    
            self.doy=self.imagedatetime.timetuple().tm_yday        
        elif ((self.filename.find('LC8') == 0) | (self.filename.find('LO8') == 0) | \
                (self.filename.find('LE7') == 0) | (self.filename.find('LT5') == 0) | \
                (self.filename.find('LT4') == 0)):    # looks landsat like (old filenames) - try parsing the date from filename (contains day of year)
            self.sensor=self.filename[0:3]
            self.path=int(self.filename[3:6])
            self.row=int(self.filename[6:9])
            self.year=int(self.filename[9:13])
            self.doy=int(self.filename[13:16])
            self.imagedatetime=dt.datetime.fromordinal(dt.date(self.year-1,12,31).toordinal()+self.doy)
        elif ( (self.filename.find('S2A') == 0) | (self.filename.find('S2B') == 0) | \
                ((self.filename.find('T') == 0) & (self.filename.find('_') == 6)) ):	# looks like sentinal 2 data (old or new format) - try parsing the date from filename (contains day of year)
            if self.filename.find('S2') == 0:  # old format Sentinel 2 data
                self.sensor=self.filename[0:3]
                b=re.search('_(?P<date>\d{8})T(?P<time>\d{6})_T(?P<tile>[A-Z0-9]{5})_A(?P<orbit>\d{6})_R(?P<rel_orbit>\d{3})_',self.filename)
                self.path=np.mod(int(b.group('orbit')),143)+3  # why + 3?  there is an offset between rel_orbit and absolute orbit numbers for S2A
                self.tile=b.group('tile')
                self.imagedatetime=dt.datetime.strptime(b.group('date'),'%Y%m%d')
            else:
                self.sensor='S2'  # would have to get S2A or when it flies S2B from full file path, which I may not maintain
                b=re.search('T(?P<tile>[A-Z0-9]{5})_(?P<date>\d{8})T(?P<time>\d{6})',self.filename)
                self.tile=b.group('tile')
                self.imagedatetime=dt.datetime.strptime(b.group('date'),'%Y%m%d')
        else:
            self.imagedatetime=None  # need to throw error in this case...or get it from metadata
        #         self.img=self.gd.ReadAsArray().astype(np.float32)   # works for L8 and earlier - and openCV correlation routine needs float or byte so just use float...
    def imageij2XY(self,ai,aj,outx=None,outy=None):
        it = np.nditer([ai,aj,outx,outy],
                        flags = ['external_loop', 'buffered'],
                        op_flags = [['readonly'],['readonly'],
                                    ['writeonly', 'allocate', 'no_broadcast'],
                                    ['writeonly', 'allocate', 'no_broadcast']])
        for ii,jj,ox,oy in it:
            ox[...]=(self.XYtfm[0]+((ii+0.5)*self.XYtfm[2]));
            oy[...]=(self.XYtfm[1]+((jj+0.5)*self.XYtfm[3]));
        return np.array(it.operands[2:4])
    def XY2imageij(self,ax,ay,outi=None,outj=None):
        it = np.nditer([ax,ay,outi,outj],
                        flags = ['external_loop', 'buffered'],
                        op_flags = [['readonly'],['readonly'],
                                    ['writeonly', 'allocate', 'no_broadcast'],
                                    ['writeonly', 'allocate', 'no_broadcast']])
        for xx,yy,oi,oj in it:
            oi[...]=((xx-self.XYtfm[0])/self.XYtfm[2])-0.5;  # if python arrays started at 1, + 0.5
            oj[...]=((yy-self.XYtfm[1])/self.XYtfm[3])-0.5;  # " " " " "
        return np.array(it.operands[2:4])
#         self.img=self.gd.ReadAsArray().astype(np.uint8)        # L7 and earlier - doesn't work with plt.imshow...
#         self.img_ov2=self.img[0::2,0::2]
#         self.img_ov10=self.img[0::10,0::10]






startproctime=time.perf_counter()
startwalltime=time.time()
t_line = "="*70
# output_log_file_fid=None	# will be fid of output log file if the -nlf flag not thrown - then log to file as well as terminal.

log_output_lines=[]
def t_log(s, outlogdisabled=False,elapsed=None):
    print(t_line)
    print(time.strftime('UTC: %x_%X', time.gmtime()),'  Wall:%8.3f'%(time.time()-startwalltime),'  Proc:%8.3f'%(time.perf_counter()-startproctime), '-', s)
    if elapsed:
        print("Elapsed time: %s" % elapsed)
    print(t_line + '\n')
    if not(outlogdisabled):
    	log_output_lines.append(time.strftime('UTC: %x_%X', time.gmtime()) + '  Wall:%8.3f'%(time.time()-startwalltime) + '  Proc:%8.3f'%(time.perf_counter()-startproctime) + ' - ' + s + '\n')
    

# from accepted answer and second answer here: http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
def polyfit2d(x, y, f, deg):
        x = np.asarray(x)
        y = np.asarray(y)
        f = np.asarray(f)
        deg = np.asarray(deg)
        vander = polynomial.polyvander2d(x, y, deg)
        vander = vander.reshape((-1,vander.shape[-1]))
        f = f.reshape((vander.shape[0],))
        # print(np.linalg.lstsq(vander, f))
        c,resid = np.linalg.lstsq(vander, f)[0:2]   # [0] is coefs, [1] is residuals
        return c.reshape(deg+1)
        

# set up command line arguments
parser = argparse.ArgumentParser( \
    description="""uses image to image correlation to detect offsets in surface features; 
    				produces map of offsets in units of pixels and velocity in m/day or m/year""",
    epilog='>>  <<',
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-imgdir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='single source dir for both images [.]')
parser.add_argument('-img1dir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='source dir for image 1 [.]')
parser.add_argument('-img2dir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='source dir for image 2 [.]')
parser.add_argument('-output_dir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='output dir [.]')
parser.add_argument('-img1datestr', 
                    action='store', 
                    type=str, 
                    help='date string for image 1 [None - set from L8 filename]')
parser.add_argument('-img2datestr', 
                    action='store', 
                    type=str, 
                    help='date string for image 2 [None - set from L8 filename]')
parser.add_argument('-datestrfmt', 
                    action='store', 
                    type=str,
                    default='%m/%d/%Y',  
                    help='date string format for img1datestr and img2datestr [None - set from L8 filename] eg. %%m/%%d/%%Y - SEE: https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior')                    
parser.add_argument('img1_name', 
                    action='store', 
                    type=str, 
                    help='image 1 filename')
parser.add_argument('img2_name', 
                    action='store', 
                    type=str, 
                    help='image 2 filename')
parser.add_argument('-out_name_base', 
                    action='store', 
                    type=str,
                    default='temp_out',
                    help='output filename base')
parser.add_argument('-bbox', 
                    action='store', 
                    type=float, 
                    metavar=('min_x','min_y','max_x','max_y'), 
                    nargs=4, 
                    help='bbox for feature tracking area in projection m - minx miny maxx maxy - defaults to entire common area between images')
# parser.add_argument('-plotvmin', 
#                     action='store', 
#                     type=float, 
#                     default=0.0, 
#                     help='min vel for colormap [%(default)d]')
parser.add_argument('-plotvmax', 
                    action='store', 
                    type=float, 
                    default=15.0, 
                    help='max vel for colormap [%(default)d]')
parser.add_argument('-trackvmax', 
                    action='store', 
                    type=float, 
                    help='max vel that can be tracked (half_target_chip will be made larger if necessary, but if speed_ref is slow, that limit will be used...) [not set]')
parser.add_argument('-nodatavmax', 
                    action='store', 
                    type=float,
                    default=0.333,
                    help='max vel (m/d) that can be tracked if speedref file is used and has nodata at location (using trackvmax for this for large time separations impractical - think ocean pixels...) [%(default)f]')
parser.add_argument('-half_source_chip', 
                    action='store',
                    default=10,  
                    type=int, 
                    help='half size of source (small) chip [%(default)d]')
parser.add_argument('-half_target_chip', 
                    action='store',
                    default=20,  
                    type=int, 
                    help='half size of target (large or search) chip [%(default)d]')
parser.add_argument('-inc', 
                    action='store', 
                    type=int, 
                    default=20, 
                    help='inc(rement) or chip center grid spacing in pixels (must be even integer) [%(default)d]')
parser.add_argument('-gfilt_sigma', 
                    action='store', 
                    type=float, 
                    default=3.0, 
                    help='gaussian filter sigma (standard deviation in pixels for Gaussian kernel) [%(default)f]')
parser.add_argument('-dcam', 
                    action='store', 
                    type=float, 
                    default=0.05, 
                    help='min difference between corr peaks 1 and 2 for mask ( full mask statement: masked_where(((del_corr_arr<dcam)&(corr_arr<cam)) | (corr_arr<cam1)) ) [%(default)f]')
parser.add_argument('-cam', 
                    action='store', 
                    type=float, 
                    default=1.0, 
                    help='corr peak max value for dcam & cam mask (see full mask statement in -dcam help) [%(default)f]')
parser.add_argument('-cam1', 
                    action='store', 
                    type=float, 
                    default=0.0, 
                    help='corr peak value max for or mask (see full mask statement in -dcam help) [%(default)f]')
# parser.add_argument('-of', 
#                     action='store', 
#                     type=str, 
#                     default='GTiff', 
#                     help='output format string (GDAL short format name - presently only GTiff and ENVI supported) [GTiff]')
parser.add_argument('-no_speed_ref', 
                    action='store_true',  
                    default=False, 
                    help='Do NOT use speed reference to set search chip sizes (and optionally to correct offsets) (use default max speed instead) [False if not raised]')
parser.add_argument('-Greenland', 
                    action='store_true',  
                    default=False, 
                    help='this pair is in Greenland - needed so that vref vv mosaic is reprojected when input (unlike Antarctica) [False if not raised]')
parser.add_argument('-speed_ref_dir', 
                    action='store', 
                    type=str, 
                    default='.', 
                    help='path for directory where speed_ref mosaic (speed in m/yr) is kept [%(default)s]')
parser.add_argument('-speed_ref_filename', 
                    action='store', 
                    type=str, 
                    default='LISA_750m_vv.tif', 
                    help='name of speed_ref mosaic (speed in m/yr) [%(default)s]')
parser.add_argument('-lgo_mask_filename', 
                    action='store', 
                    type=str, 
                    default=None, 
                    help='file name for land_glacier_ocean mask file (used for lgo-mask-based offset correction) [%(default)s]')
parser.add_argument('-lgo_mask_file_dir', 
                    action='store', 
                    type=str, 
                    default=None, 
                    help='path to diretory containing land_glacier_ocean mask file (used for lgo-mask-based offset correction) [%(default)s]')
parser.add_argument('-VRT_dir', 
                    action='store', 
                    type=str, 
                    default='.', 
                    help='path to directory where temporary .vrt files will be stored, then deleted, for access to part of speed_ref or lgo-mask mosaic [%(default)s]')
parser.add_argument('-max_allowable_pixel_offset_correction', 
                    action='store', 
                    type=float, 
                    default=3.0, 
                    help='if calculated offset correction of any type is larger than this many pixels, it will NOT be applied [%(default)f pixels]')
######################################
#
# the rest of the parameters are flags - if raised, set to true, otherwise false
#
######################################
parser.add_argument('-do_not_highpass_input_images', 
                    action='store_true',
                    default=False, 
                    help='DO NOT highpass input images before correlation - [highpass images before correlation]')
parser.add_argument('-v', 
                    action='store_true',  
                    default=False, 
                    help='verbose - extra diagnostic and image info put in log file [False if not raised]')
parser.add_argument('-mpy', 
                    action='store_true',
                    default=False, 
                    help='meters per year - [meters per day if not raised]')
parser.add_argument('-log10', 
                    action='store_true',
                    default=False, 
                    help='output second rgba color image that is log10(v) (and second kmz with colorbar in addition to GTiff if requested) - [only linear color image produced]')
# parser.add_argument('-no_gtif', 
#                     action='store_true',
#                     default=False, 
#                     help='do not output default GTiff - [geoTiff always produced if not raised]')
parser.add_argument('-nlf', 
                    action='store_true',
                    default=False, 
                    help='do not output log entries in nc file that contain command line etc - [output log entries in nc file]')
parser.add_argument('-progupdates', 
                    action='store_true',
                    default=False, 
                    help='do not provide progress updates - [output progress updates]')
parser.add_argument('-debug_tifs', 
                    action='store_true',
                    default=False, 
                    help='output int offsets and other diagnostic tifs - [no]')
parser.add_argument('-offset_correction_speedref', 
                    action='store_true',
                    default=False, 
                    help='estimate offset at slow moving areas and correct output vels with constant vx and vy shifts (requires -speed_ref files (vv,vx,vy - vv specified)) - [False]')
parser.add_argument('-offset_correction_lgo_mask', 
                    action='store_true',
                    default=False, 
                    help='estimate offset from land pixels and make correct output vels with constant vx and vy shifts (requires -lgo_mask_file_fullpath - [False]')
parser.add_argument('-lgo_mask_limit_land_offset', 
                    action='store_true',
                    default=False, 
                    help='when source chip center is a land pixel, limit search distance (requires -lgo_mask_file_fullpath - [False]')
parser.add_argument('-offset_correction_bilinear_fit', 
                    action='store_true',
                    default=False, 
                    help='use planar corrections to both i and j offsets to correct misfit - currently ONLY for lgo land pixels... - [False]')
parser.add_argument('-offset_correction_bilinear_fit_min_pixels', 
                    action='store',
                    type=int,
                    default=2000, 
                    help='min valid land pixels to use planar corrections to both i and j offsets to correct misfit - [%(default)d]')
args = parser.parse_args()


outdir=args.output_dir
out_name_base=args.out_name_base
out_parent_dir='./'
out_file_dir='./'
file_name_base=out_name_base
	
if not(args.nlf):
# 	outsf=open(outdir + '/' + file_name_base+"_log.txt",'w')
# 	output_log_file_fid=outsf  # this pointer to the log file fid used in t_log...
	log_output_lines.append('# log from %s\n'% sys.argv[0])
	log_output_lines.append(" ".join(sys.argv) + "\n")
	log_output_lines.append(str(args) + "\n")
	
t_log("Start Program",outlogdisabled=args.nlf)

if(args.v):
	print('#', args)

# set up source directories
if(args.imgdir != "."): # single source directory specified
	if((args.img1dir != ".") | (args.img2dir != ".")):
		sys.exit('Error: -imgdir (specifying a single source directory for both images) cannot be used with -img1dir or -img2dir')
	else:
		image1dir=args.imgdir
		image2dir=args.imgdir
else:    # -imgdir not used - set directories to img1dir and img2dir
	image1dir=args.img1dir
	image2dir=args.img2dir


if psutil_available:
	print('psutil reports process %s using '%(args.out_name_base),memory_usage_psutil())
if resource_available:
	print('resource module reports process %s using '%(args.out_name_base),memory_usage_resource())

format = "GTiff"
driver = gdal.GetDriverByName( format )

int_limit_max = np.iinfo('int16').max  # get max (and min) value for rescaling float hp before converting back to signed int16
int_limit_minp1 = np.iinfo('int16').min + 1
int_nodata_val = np.iinfo('int16').min
img1_names=args.img1_name.split('.')

hp1filename = args.img1_name

hp2filename = args.img2_name


if psutil_available:
	print('psutil reports process %s using '%(args.out_name_base),memory_usage_psutil())
if resource_available:
	print('resource module reports process %s using '%(args.out_name_base),memory_usage_resource())

##########################################################################
# now read hp images in and generate mask for img1
##########################################################################

img1=GeoImg_noload(hp1filename,in_dir=image1dir,datestr=args.img1datestr,datefmt=args.datestrfmt)
if psutil_available:
	print('about to read hp img1 - using ',memory_usage_psutil())






if not(args.do_not_highpass_input_images):
	t_log("making low pass and high pass images (set -do_not_highpass_input_images to disable this step)")
	inimg1=img1.gd.ReadAsArray().astype(np.float32)
	lp_arr_1=inimg1.copy()
	gaussian_filter(lp_arr_1, args.gfilt_sigma, output=lp_arr_1)
	hp_arr_1=inimg1 - lp_arr_1
	lp_arr_1=None
	t_log("done with low pass and high pass image 1 - removed lp")
	img1_data_mask=np.ones_like(hp_arr_1, dtype=bool)
	# need to test for presence of nodata value - fixed in 5p2
	if img1.nodata_value:
		img1_data_mask[inimg1==np.float32(img1.nodata_value)]=False
	else:
		img1_data_mask[inimg1==0]=False
	if psutil_available:
		print('made mask - psutil reports process %s using '%(args.out_name_base),memory_usage_psutil())
	if resource_available:
		print('made mask - resource module reports process %s using '%(args.out_name_base),memory_usage_resource())
	inimg1=None		
else:
    t_log("not filtering input images because -do_not_highpass_input_images is set")
    hp_arr_1=img1.gd.ReadAsArray().astype(np.float32)
    img1_data_mask=np.ones_like(hp_arr_1, dtype=bool)
    # need to test for presence of nodata value - fixed in 5p2
    if img1.nodata_value:
        img1_data_mask[hp_arr_1==np.float32(img1.nodata_value)]=False
    else:
        img1_data_mask[hp_arr_1==0]=False

        if psutil_available:
            print('made mask - psutil reports process %s using '%(args.out_name_base),memory_usage_psutil())
        if resource_available:
            print('made mask - resource module reports process %s using '%(args.out_name_base),memory_usage_resource())

if psutil_available:
	print('read hp img1 - psutil reports process %s using '%(args.out_name_base),memory_usage_psutil())
if resource_available:
	print('read hp img1 - resource module reports process %s using '%(args.out_name_base),memory_usage_resource())



img2=GeoImg_noload(hp2filename,in_dir=image2dir,datestr=args.img2datestr,datefmt=args.datestrfmt)
if not(args.do_not_highpass_input_images):
# 	t_log("sentinel 2 data, so making low pass and high pass images")
# 	file_name_base = file_name_base + '_hp_filt_%3.1f'%(args.gfilt_sigma)
	inimg2=img2.gd.ReadAsArray().astype(np.float32)
	lp_arr_2=inimg2.copy()
	gaussian_filter(lp_arr_2, args.gfilt_sigma, output=lp_arr_2)
	hp_arr_2=inimg2 - lp_arr_2
	lp_arr_2=None
	t_log("done with low pass and high pass image 2 - removed lp")
	inimg2=None
		
else:
	hp_arr_2=img2.gd.ReadAsArray().astype(np.float32)

if psutil_available:
	print('read hp img2 - psutil reports process %s using '%(args.out_name_base),memory_usage_psutil())
if resource_available:
	print('read hp img2 - resource module reports process %s using '%(args.out_name_base),memory_usage_resource())






del_t=img2.imagedatetime - img1.imagedatetime
del_t_days=del_t.total_seconds()/(24.0*3600.0)
del_t_years=del_t_days/(365.25)
	
# switch from m/d to m/yr if -mpy
if (args.mpy):
	del_t_unit_str='years'
	del_t_speedunit_str='m/yr'
	del_t_val=del_t_years
else:
	del_t_unit_str='days'
	del_t_speedunit_str='m/d'
	del_t_val=del_t_days


# offset pixels to subtract from ij_1 to get ij_2 coordinates (how far the two image arrays are offset based only on their geolocation)
offset_2_1=img1.XY2imageij(*img2.imageij2XY(*[0,0])) 

# either work in specified bounding box, or find largest box of common pixels between the two images and use that
if (args.bbox is None):
	# 	find largest common x,y using original XY (geo cartesian coordinates) box between images
	# image coordinates are for bounding box (outer edge of pixels)
	com_max_x=np.min([img1.max_x, img2.max_x])
	com_max_y=np.min([img1.max_y, img2.max_y])
	com_min_x=np.max([img1.min_x, img2.min_x])
	com_min_y=np.max([img1.min_y, img2.min_y])
else: # read area from command line specified bounding box (still XY etc.)
	print('specifying -bbox is not implemented in this version of sc_pycorr - must allow it to use common overlap area of the two images - exiting')
	sys.exit(-1)
	com_max_x=args.bbox[2]
	com_max_y=args.bbox[3]
	com_min_x=args.bbox[0]
	com_min_y=args.bbox[1]
#############
####  should check that a bbox, if specified on the command line, fits within the common footprint of the two images - but so far we never specify a bbox on the command line
####  EITHER remove this option OR put in this check...
#############


bbox=[com_min_x, com_min_y, com_max_x, com_max_y]
# plotvmin and plotvmax - speeds for colorbar in color_based geotiffs... - ?maybe set to None for autoscaling?
# plotvmin=args.plotvmin
plotvmax=args.plotvmax
# plot_2_mult=args.plot_2_mult # vmax for second plot will be plot_2_mult * plotvmax - to show detail on low end...

half_source_chip=args.half_source_chip # must be < half_target_chip - not that anything would happen, if trackvmax is specified.  This is not enforced but could be to catch errors.
half_target_chip=args.half_target_chip  


# if trackvmax specified, make sure half_target_chip is large enough to support measuring that large an offset. This is set before the buffer around the grid is worked out, so that target chips remain in the images.
if args.trackvmax:
	if args.trackvmax > ((half_target_chip - half_source_chip) * img1.pix_x_m / del_t_days):
		half_target_chip = np.int(np.ceil(((args.trackvmax * del_t_days)/img1.pix_x_m) + half_source_chip))
		print('increased half_target_chip to %d to allow tracking up to speed -trackvmax %f'%(half_target_chip,args.trackvmax))
maxtrackablevel=((half_target_chip - half_source_chip) * img1.pix_x_m / del_t_days)

print('>>>> half_target_chip %d pixels - maximum trackable velocity %f m/d (unless speed_ref at that location is lower, then it is 1.5 * (speed_ref_vel/pix_size_x_m) * del_t_days + 5 pixels) <<<<'%(half_target_chip,maxtrackablevel))
print("time separation in days: %f years: %f - using %s"%(del_t_days,del_t_years,del_t_speedunit_str))
if not(args.nlf):
	log_output_lines.append('>>>> half_target_chip %d pixels - maximum trackable velocity %f m/d (unless speed_ref at that location is lower, then it is 1.5 * (speed_ref_vel/pix_size_x_m) * del_t_days + 5 pixels) <<<<'%(half_target_chip,maxtrackablevel))
	log_output_lines.append("time separation in days: %f years: %f - using %s\n"%(del_t_days,del_t_years,del_t_speedunit_str))

inc=args.inc
if ((inc/2.0) != int(inc/2.0)):
	print('-inc must be an even integer number of pixels for output geolocation to register properly; currently inc=%f'%(inc))
	sys.exit(-1)
	
###############################
# sub_pixel interpolation setup
###############################
rbs_halfsize=3 # size of peak area used for spline for subpixel peak loc
rbs_order=3    # polynomial order for subpixel rbs interpolation of peak location

d2xdx2_spacing=0.05 # pixels +- for calculating curvature of spline surface at peak correlation...


# use this box to find overlap range in i and j (array coordinates) for both images
com_ul_i_j_img1=img1.XY2imageij(*(com_min_x,com_max_y))+0.5
com_ul_i_j_img2=img2.XY2imageij(*(com_min_x,com_max_y))+0.5
com_lr_i_j_img1=img1.XY2imageij(*(com_max_x,com_min_y))-0.5
com_lr_i_j_img2=img2.XY2imageij(*(com_max_x,com_min_y))-0.5


# image pixels in bbox for img1 as a slice are img1.img[c1_minj:c1_maxjp1,c1_mini:c1_maxip1]
# remember (Mark) that a python slice [3:5] returns elements 3 and 4, indexed from 0

c1_minj=com_ul_i_j_img1[1]
c1_maxjp1=(com_lr_i_j_img1[1]+1) # may be one more than the max index, for use in slices (p1 means plus 1...)
c1_mini=com_ul_i_j_img1[0]
c1_maxip1=(com_lr_i_j_img1[0]+1) # same

c2_minj=com_ul_i_j_img2[1]
c2_maxjp1=(com_lr_i_j_img2[1]+1) # same
c2_mini=com_ul_i_j_img2[0]
c2_maxip1=(com_lr_i_j_img2[0]+1) # same

com_num_pix_i=c1_maxip1-c1_mini  # number of pixels per line in common (overlap) area
com_num_pix_j=c1_maxjp1-c1_minj  # number of lines in common (overlap) area

if(args.v):
	print('numpix %f numlines %f in box'%(com_num_pix_i,com_num_pix_j))
	print('ul X %f Y %f of box'%(com_min_x,com_max_y))
	print('ul image1 i %f j%f  image2 i %f j%f'%(c1_mini,c1_minj,c2_mini,c2_minj))
	print('lr X %f Y %f of box'%(com_max_x,com_min_y))
	print('lr image1 i %f j%f  image2 i %f j%f'%(c1_maxip1,c1_maxjp1,c2_maxip1,c2_maxjp1))


half_inc_rim=int(inc/2.0)  # inc was checked to be even (above) so int is redundant but included for clarity - want integer-width rim...

num_grid_i=int(com_num_pix_i/inc)		# these are integer divisions...
num_grid_j=int(com_num_pix_j/inc)

ul_i_1_chip_grid=int(c1_mini + half_inc_rim)
ul_j_1_chip_grid=int(c1_minj + half_inc_rim)
lr_i_1_chip_grid=int(c1_mini + half_inc_rim + ((num_grid_i-1) * inc))
lr_j_1_chip_grid=int(c1_minj + half_inc_rim + ((num_grid_j-1) * inc))

ul_i_2_chip_grid=int(c2_mini + half_inc_rim)
ul_j_2_chip_grid=int(c2_minj + half_inc_rim)
lr_i_2_chip_grid=int(c2_mini + half_inc_rim + ((num_grid_i-1) * inc))
lr_j_2_chip_grid=int(c2_minj + half_inc_rim + ((num_grid_j-1) * inc))

r_i_1=range(ul_i_1_chip_grid,lr_i_1_chip_grid+1,inc) # range over common area
r_j_1=range(ul_j_1_chip_grid,lr_j_1_chip_grid+1,inc)

r_i_2=range(ul_i_2_chip_grid,lr_i_2_chip_grid+1,inc) # range over common area
r_j_2=range(ul_j_2_chip_grid,lr_j_2_chip_grid+1,inc)

# find required size of buffer around chip center positons so that target chip centered there remains in bbox
##### Note that this does not take advantage of a rim of pixels in image 2 that are outside the common pixels in the two images - this could be fixed
# not sure that last comment is true any more - should be ok now.
min_ind_i_1=np.min(np.nonzero(np.array(r_i_1)>half_target_chip))# this is width of the rim of zeros on left or top in the correlation arrays
min_ind_j_1=np.min(np.nonzero(np.array(r_j_1)>half_target_chip))# this is width of the rim of zeros on left or top in the correlation arrays
min_ind_i_2=np.min(np.nonzero(np.array(r_i_2)>half_target_chip))# this is width of the rim of zeros on left or top in the correlation arrays
min_ind_j_2=np.min(np.nonzero(np.array(r_j_2)>half_target_chip))# this is width of the rim of zeros on left or top in the correlation arrays
min_ind_i = np.max([min_ind_i_1,min_ind_i_2])
min_ind_j = np.max([min_ind_j_1,min_ind_j_2])

max_ind_i_1=np.max(np.nonzero(half_target_chip<=(img1.num_pix_x - np.array(r_i_1))))  # maximum allowed i index (buffer right edge of grid to accomodate target chip size)
max_ind_j_1=np.max(np.nonzero(half_target_chip<=(img1.num_pix_y - np.array(r_j_1))))  # maximum allowed i index (buffer bottom edge of grid to accomodate target chip size)
max_ind_i_2=np.max(np.nonzero(half_target_chip<=(img2.num_pix_x - np.array(r_i_2))))  # maximum allowed i index (buffer right edge of grid to accomodate target chip size)
max_ind_j_2=np.max(np.nonzero(half_target_chip<=(img2.num_pix_y - np.array(r_j_2))))  # maximum allowed i index (buffer bottom edge of grid to accomodate target chip size)
max_ind_i=np.min([max_ind_i_1,max_ind_i_2])  # maximum allowed i index (buffer right edge of grid to accomodate target chip size)
max_ind_j=np.min([max_ind_j_1,max_ind_j_2])  # maximum allowed i index (buffer bottom edge of grid to accomodate target chip size)

peak_block_halfsize=3 # this is the size of the area zeroed out before looking for the second highest peak
# the following assumes square chips - which we use here so far...
cent_loc=half_target_chip-half_source_chip  # Note cent_loc,cent_loc is array i,j of chip center because indexing is zero based

# the -0.5 in what follows shifts from the ij of the element to it's upper left corner (chip center)
r_grid=np.meshgrid(np.array(r_i_1)-0.5,np.array(r_j_1)-0.5,indexing='xy')
chip_center_grid_xy=img1.imageij2XY(*r_grid)
# note output_array_ul_corner will be the same as [com_min_x, com_max_y] but is derived here just to check coordinates
output_array_ul_corner=chip_center_grid_xy[:,0,0]-((inc/2.0)*img1.pix_x_m,(inc/2.0)*img1.pix_y_m)
output_array_pix_x_m=inc*img1.pix_x_m
output_array_pix_y_m=inc*img1.pix_y_m	# note this will nearly always be negative, for image stored top line first, just as in world file
output_array_num_pix_x=r_i_2.__len__()
output_array_num_pix_y=r_j_2.__len__()


vel_nodata_val=np.float32(-9999.0)
corr_nodata_val=np.float32(-2.0)
curvature_nodata_value=np.float32(-999.0)

#################################
# set up (allocate) output arrays
#################################
corr_arr=np.ones((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)*corr_nodata_val		# NOTE: chips not used (image nodata areas) will have this value...
del_corr_arr=np.ones((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)*corr_nodata_val
offset_dist_ij_arr=np.zeros((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)
sp_dist_i=np.zeros((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)
sp_dist_j=np.zeros((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)
d2idx2=np.ones((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)*curvature_nodata_value
d2jdx2=np.ones((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)*curvature_nodata_value
d2dx2_mean=np.ones((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)*curvature_nodata_value
del_i=np.zeros((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)
del_j=np.zeros((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)
orig_int_del_i=np.zeros((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)
orig_int_del_j=np.zeros((output_array_num_pix_y,output_array_num_pix_x), dtype=np.float32)
count_tenths=np.zeros((output_array_num_pix_y,output_array_num_pix_x), dtype=np.int16)
count_hundredths=np.zeros((output_array_num_pix_y,output_array_num_pix_x), dtype=np.int16)

if psutil_available:
	print('set up output arrays - psutil reports process %s using '%(args.out_name_base),memory_usage_psutil())
if resource_available:
	print('set up output arrays - resource module reports process %s using '%(args.out_name_base),memory_usage_resource())

if not(args.no_speed_ref):
	# set up vrt file to sample reference speed mosaic at chip centers - this will be used to set search chip sizes for each point
	# if args.offset_correction_speedref is True, then read vx and vy reference data as well, for offset correction
	ulx=output_array_ul_corner[0]
	uly=output_array_ul_corner[1]
	lrx=output_array_ul_corner[0] + (output_array_num_pix_x * output_array_pix_x_m)
	lry=output_array_ul_corner[1] + (output_array_num_pix_y * output_array_pix_y_m)
	
	path_to_speed_ref_vv = args.speed_ref_dir+ '/' + args.speed_ref_filename
	output_vv_vrt_name= 'tmp' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) + 'vv.vrt'
	
	
	
	if not(args.Greenland):  # ANTARCTICA - this is in Antarctica, where Landsat images and vref mosaic are in the SAME projection!
		# -te needs xmin ymin xmax ymax - that is what is below, using the corners
		sp.call(["gdalbuildvrt", "-te", "%f"%ulx, "%f"%lry, "%f"%lrx, "%f"%uly, "-tr", "%9.4f"%(output_array_pix_x_m), "%9.4f"%(output_array_pix_x_m), "%s"%(args.VRT_dir + '/' + output_vv_vrt_name), path_to_speed_ref_vv])  # import subprocess as sp
		speedref_geoimg_vv=GeoImg_noload(output_vv_vrt_name,in_dir=args.VRT_dir)
		speedref_vel_vv=speedref_geoimg_vv.gd.ReadAsArray().astype(np.float32)
		if speedref_geoimg_vv.nodata_value:
			speedref_vv_nodata=speedref_geoimg_vv.nodata_value
		else:
			speedref_vv_nodata=None
		os.remove(args.VRT_dir + '/' + output_vv_vrt_name)  # get rid of the temporary vrt file now
		print('got speedref vel data')
		if not(args.nlf):
			log_output_lines.append('# using speed reference (speedref file: %s) from speedref vrt file: %s\n'%(args.speed_ref_filename,output_vv_vrt_name))

		if args.offset_correction_speedref:  # need to read vx and vy ref data as well...
			path_to_speed_ref_vx = args.speed_ref_dir+ '/' + args.speed_ref_filename.replace('vv','vx')   # MARIN need to fix this...
			path_to_speed_ref_vy = args.speed_ref_dir+ '/' + args.speed_ref_filename.replace('vv','vy')   # MARIN need to fix this...
			output_vx_vrt_name= 'tmp' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) + 'vx.vrt'
			output_vy_vrt_name= 'tmp' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) + 'vy.vrt'
			sp.call(["gdalbuildvrt", "-te", "%f"%ulx, "%f"%lry, "%f"%lrx, "%f"%uly, "-tr", "%9.4f"%(output_array_pix_x_m), "%9.4f"%(output_array_pix_x_m), "%s"%(args.VRT_dir + '/' + output_vx_vrt_name), path_to_speed_ref_vx])  # import subprocess as sp
			sp.call(["gdalbuildvrt", "-te", "%f"%ulx, "%f"%lry, "%f"%lrx, "%f"%uly, "-tr", "%9.4f"%(output_array_pix_x_m), "%9.4f"%(output_array_pix_x_m), "%s"%(args.VRT_dir + '/' + output_vy_vrt_name), path_to_speed_ref_vy])  # import subprocess as sp
			speedref_geoimg_vx=GeoImg_noload(output_vx_vrt_name,in_dir=args.VRT_dir)
			speedref_vel_vx=speedref_geoimg_vx.gd.ReadAsArray().astype(np.float32)
			if speedref_geoimg_vx.nodata_value:
				speedref_vx_nodata=speedref_geoimg_vx.nodata_value
			else:
				speedref_vx_nodata=None
		
			speedref_geoimg_vy=GeoImg_noload(output_vy_vrt_name,in_dir=args.VRT_dir)
			speedref_vel_vy=speedref_geoimg_vy.gd.ReadAsArray().astype(np.float32)
			if speedref_geoimg_vy.nodata_value:
				speedref_vy_nodata=speedref_geoimg_vy.nodata_value
			else:
				speedref_vy_nodata=None		
			os.remove(args.VRT_dir + '/' + output_vx_vrt_name)  # get rid of the temporary vrt file now
			os.remove(args.VRT_dir + '/' + output_vy_vrt_name)  # get rid of the temporary vrt file now
			if not(args.nlf):
				log_output_lines.append('# using speed reference vx and vy data as well\n')
		
	else:  # this is in GREENLAND, where Landsat images and vref mosaic are in DIFFERENT PROJECTIONS and so vref components must be reprojected!
		
		# first - open the vv mosaic but don't read it in - only want the srs from this file so we know it's projection - will read subregion at output resolution with .vrt
		tmp_speed_ref_orig=GeoImg_noload(args.speed_ref_filename,args.speed_ref_dir)
		source = img1.srs		# the local UTM of the input image
		target = tmp_speed_ref_orig.srs
		# set up transform from local utm to PS of vel ref image, to calculate the corners of the vref mosaic covering this area
		transform_utm_to_PS = osr.CoordinateTransformation(source, target)
		# set up inv transform to that above to project the endpoints of each PS vx,vy vector to find local vref vy and vy in utm reference frame
		inv_transform_PS_to_utm = osr.CoordinateTransformation(target, source)

		# set up to sample speed_ref mz at chip centers - this will be used to find mask value at each point
		# in L8 image local UTM, these are:
		utm_ulx=output_array_ul_corner[0]
		utm_uly=output_array_ul_corner[1]
	# 	lrx=output_array_ul_corner[0] + (r_i_2.__len__() * inc * img1.pix_x_m)
	# 	lry=output_array_ul_corner[1] + (r_j_2.__len__() * inc * img1.pix_y_m)
		utm_lrx=output_array_ul_corner[0] + (output_array_num_pix_x * output_array_pix_x_m)
		utm_lry=output_array_ul_corner[1] + (output_array_num_pix_y * output_array_pix_y_m)
	
		#########################
		# now we project output array corners from local UTM back to projection of speed_ref mz
		(tulx, tuly, tulz ) = transform_utm_to_PS.TransformPoint( utm_ulx, utm_uly, 0.0 )
		(turx, tury, turz ) = transform_utm_to_PS.TransformPoint( utm_lrx, utm_uly, 0.0 )
		(tlrx, tlry, tlrz ) = transform_utm_to_PS.TransformPoint( utm_lrx, utm_lry, 0.0 )
		(tllx, tlly, tllz ) = transform_utm_to_PS.TransformPoint( utm_ulx, utm_lry, 0.0 )
		#####################
		# Make .vrt for speed_ref image that is in local UTM to match Band 8 images 
		# - will reproject speed_ref img in memory to local utm after reading in just this part of data - .vrt limits the read size...
		#####################
		ulx=np.min([tulx,tllx])
		uly=np.max([tuly,tury])
		lrx=np.max([turx,tlrx])
		lry=np.min([tlly,tlry])
	
		sp.call(["gdalbuildvrt", "-te", "%f"%ulx, "%f"%lry, "%f"%lrx, "%f"%uly, "-tr", "%9.4f"%(output_array_pix_x_m), "%9.4f"%(output_array_pix_x_m), "%s"%(args.VRT_dir + '/' + output_vv_vrt_name), path_to_speed_ref_vv])  # import subprocess as sp
		speed_ref_vv_cropped_geoimg=GeoImg_noload(output_vv_vrt_name,in_dir=args.VRT_dir)
# AM I CORRECT (11/3/16)? that this next read served no purpose? commenting it out here to see...
# 		speed_ref_vv_cropped_mask=speed_ref_vv_cropped_geoimg.gd.ReadAsArray().astype(np.float32)
		if speed_ref_vv_cropped_geoimg.nodata_value:
			speedref_vv_nodata=speed_ref_vv_cropped_geoimg.nodata_value
		else:
			speedref_vv_nodata=None
		tmp_speed_ref_orig=None  # close GeoImg_noload object used on full original lgo_mask file, used to get srs above, but no input of data
	
		#####################
		# now reproject cropped input mask into local UTM so speed_ref_vv can be used to identify
		# pixels that shouldn't move
		#####################
		mem_drv = gdal.GetDriverByName( 'MEM' )
		# 	inmem_lgo_mask_utm=mem_drv.CreateCopy('',in_vx.gd)   # actually don't need the data, but the rest is good
		inmem_speed_ref_vv_utm = mem_drv.Create('', output_array_num_pix_x, output_array_num_pix_y, 1, gdal.GDT_Float32)
		inmem_speed_ref_vv_utm.SetGeoTransform( [ output_array_ul_corner[0], output_array_pix_x_m, 0, output_array_ul_corner[1], 0, output_array_pix_y_m ] ) #
		inmem_speed_ref_vv_utm.SetProjection ( img1.proj )
		#############################################
		# following is fix for mem driver not copying nodata pixels - set nodata value and fill band before reprojecting
		# FROM: https://trac.osgeo.org/gdal/ticket/6404 or stackexchange links found there
		#############################################
		rb1=inmem_speed_ref_vv_utm.GetRasterBand(1)
		rb1.SetNoDataValue(speedref_vv_nodata)
		rb1.Fill(speedref_vv_nodata)

		res = gdal.ReprojectImage( speed_ref_vv_cropped_geoimg.gd, inmem_speed_ref_vv_utm, speed_ref_vv_cropped_geoimg.srs.ExportToWkt(), img1.srs.ExportToWkt(), gdal.GRA_NearestNeighbour )
		speedref_vel_vv=inmem_speed_ref_vv_utm.ReadAsArray().astype(np.float32)

	
		print('got speed_ref_vv data and reprojected to original image local utm')
		if not(args.nlf):
			log_output_lines.append('# got speed_ref_vv data and reprojected to original image local utm')
			log_output_lines.append('# speed_ref_vv reference (speed_ref_vv mask file: %s) from speed_ref_vv vrt file: %s\n'%(path_to_speed_ref_vv,output_vv_vrt_name))
		
		if args.offset_correction_speedref:  # need to read and reproject vx and vy ref data as well...
			path_to_speed_ref_vx = args.speed_ref_dir+ '/' + args.speed_ref_filename.replace('vv','vx')
			path_to_speed_ref_vy = args.speed_ref_dir+ '/' + args.speed_ref_filename.replace('vv','vy')
			output_vx_vrt_name= 'tmp' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) + 'vx.vrt'
			output_vy_vrt_name= 'tmp' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) + 'vy.vrt'
			sp.call(["gdalbuildvrt", "-te", "%f"%ulx, "%f"%lry, "%f"%lrx, "%f"%uly, "-tr", "%9.4f"%(output_array_pix_x_m), "%9.4f"%(output_array_pix_x_m), "%s"%(args.VRT_dir + '/' + output_vx_vrt_name), path_to_speed_ref_vx])  # import subprocess as sp
			sp.call(["gdalbuildvrt", "-te", "%f"%ulx, "%f"%lry, "%f"%lrx, "%f"%uly, "-tr", "%9.4f"%(output_array_pix_x_m), "%9.4f"%(output_array_pix_x_m), "%s"%(args.VRT_dir + '/' + output_vy_vrt_name), path_to_speed_ref_vy])  # import subprocess as sp
			speed_ref_vx_cropped_geoimg=GeoImg_noload(output_vx_vrt_name,in_dir=args.VRT_dir)
			speed_ref_vy_cropped_geoimg=GeoImg_noload(output_vy_vrt_name,in_dir=args.VRT_dir)
			# reproject PS vx value into utm before converting to utm vx value
			inmem_speed_ref_vx_utm = mem_drv.Create('', output_array_num_pix_x, output_array_num_pix_y, 1, gdal.GDT_Float32)
			inmem_speed_ref_vx_utm.SetGeoTransform( [ output_array_ul_corner[0], output_array_pix_x_m, 0, output_array_ul_corner[1], 0, output_array_pix_y_m ] ) #
			inmem_speed_ref_vx_utm.SetProjection ( img1.proj )
			res = gdal.ReprojectImage( speed_ref_vx_cropped_geoimg.gd, inmem_speed_ref_vx_utm, speed_ref_vx_cropped_geoimg.srs.ExportToWkt(), img1.srs.ExportToWkt(), gdal.GRA_NearestNeighbour )
			speedref_vel_PS_vx_in_utm=inmem_speed_ref_vv_utm.ReadAsArray().astype(np.float32)
			# reproject PS vy value into utm before converting to utm vy value
			inmem_speed_ref_vy_utm = mem_drv.Create('', output_array_num_pix_x, output_array_num_pix_y, 1, gdal.GDT_Float32)
			inmem_speed_ref_vy_utm.SetGeoTransform( [ output_array_ul_corner[0], output_array_pix_x_m, 0, output_array_ul_corner[1], 0, output_array_pix_y_m ] ) #
			inmem_speed_ref_vy_utm.SetProjection ( img1.proj )
			res = gdal.ReprojectImage( speed_ref_vy_cropped_geoimg.gd, inmem_speed_ref_vy_utm, speed_ref_vy_cropped_geoimg.srs.ExportToWkt(), img1.srs.ExportToWkt(), gdal.GRA_NearestNeighbour )
			speedref_vel_PS_vy_in_utm=inmem_speed_ref_vv_utm.ReadAsArray().astype(np.float32)
			# will convert the PSvx,PSvy components (at local utm centers already) to local utmvx,utmvy components to compare with feature tracking of local utm Landsat...
			if speedref_vv_nodata:
				speedref_vel_vx=np.ones_like(speedref_vel_vv) * speedref_vv_nodata
				speedref_vel_vy=np.ones_like(speedref_vel_vv) * speedref_vv_nodata
				pts=np.where(speedref_vel_vv!=speedref_vv_nodata)  # indices of points to reproject vx and vy
			else:
				speedref_vel_vx=np.zeros_like(speedref_vel_vv)
				speedref_vel_vy=np.zeros_like(speedref_vel_vv)
				pts=np.where((speedref_vel_vv!=0.0)&(speedref_vel_vv!=-99.0))  # indices of points to reproject vx and vy
			for i,j in zip(*pts): # iterate over these points
				PSvx=speedref_vel_PS_vx_in_utm[i,j]
				PSvy=speedref_vel_PS_vy_in_utm[i,j]
				UTMx,UTMy=chip_center_grid_xy[:,i,j]
				(PSx,PSy,PSz) = transform_utm_to_PS.TransformPoint( UTMx, UTMy, 0.0 )
				(UTMx_vec_end_point, UTMy_vec_end_point, UTMz_vec_end_point) = inv_transform_PS_to_utm.TransformPoint( PSx + PSvx, PSy + PSvy, PSz)
				speedref_vel_vx[i,j] = UTMx_vec_end_point - UTMx
				speedref_vel_vy[i,j] = UTMy_vec_end_point - UTMy
		
		os.remove(args.VRT_dir + '/' + output_vx_vrt_name)  # get rid of the temporary vrt file now
		os.remove(args.VRT_dir + '/' + output_vy_vrt_name)  # get rid of the temporary vrt file now
		os.remove(args.VRT_dir + '/' + output_vv_vrt_name)  # get rid of the temporary vrt file now
##########################
# commented these out for debugging driver issue - uncomment for production to free memory
##########################
# 		inmem_speed_ref_vv_utm=None
# 		inmem_speed_ref_vx_utm=None
# 		inmem_speed_ref_vy_utm=None
		
		print('got speed_ref_vx and vy data and reprojected to original image local utm')
		if not(args.nlf):
			log_output_lines.append('# got speed_ref_vx and vy data and reprojected to original image local utm')
		
	if psutil_available:
		print('got speed_ref data - psutil reports process %s using '%(args.out_name_base),memory_usage_psutil())
	if resource_available:
		print('got speed_ref data - resource module reports process %s using '%(args.out_name_base),memory_usage_resource())
	
	
	
	
	
if args.offset_correction_lgo_mask:
	# use land(1)-glacier(0)-ocean(2) mask to find land pixels for offset correction - means reprojecting mask into Landsat image's local UTM
	#########################  Following code started from GRN_find_offset_reproject_to_PS_olg_mask_v3_fitWCS_ifnoland_multiprocess.py
	# now make a .vrt of the area of the land/ocean/glacier mask image (so you only 
	# have to read in a small part of the mask)
	# original mask is in it's own projection 
	# find local UTM corners of input, transform to mask's projection, set up .vrt of required area 
	# from mask, read it in, warp in memory to local UTM, apply land mask to find
	# pixels that should have had 0 vel, do offset correction
	# seems like a long way around, but would like to have consistent dataset tied to original image
	# Will need to reproject for mosaicking - may want to use a central UTM for each group of glaciers rather than ACC?
	#########################
	# first - open the lgo mask but don't read it in - only want the srs from this file so we know it's projection - will read subregion at output resolution with .vrt
	tmp_lgo_orig=GeoImg_noload(args.lgo_mask_filename,args.lgo_mask_file_dir)
	source = img1.srs		# the local UTM of the L8 image
# 	target = osr.SpatialReference()   # to use .ImportFromEPSG you have to create an object first - not needed from image that has been opened with gdal
	target = tmp_lgo_orig.srs
# 	target.ImportFromEPSG(102006) # Alaska Albers Conformal Conic
# 	target.ImportFromEPSG(3413) # PS
	transform = osr.CoordinateTransformation(source, target)

	# set up to sample lgo_mask at chip centers - this will be used to find mask value at each point
	# in L8 image local UTM, these are:
	utm_ulx=output_array_ul_corner[0]
	utm_uly=output_array_ul_corner[1]
# 	lrx=output_array_ul_corner[0] + (r_i_2.__len__() * inc * img1.pix_x_m)
# 	lry=output_array_ul_corner[1] + (r_j_2.__len__() * inc * img1.pix_y_m)
	utm_lrx=output_array_ul_corner[0] + (output_array_num_pix_x * output_array_pix_x_m)
	utm_lry=output_array_ul_corner[1] + (output_array_num_pix_y * output_array_pix_y_m)
	
	#########################
	# now we project output array corners from local UTM back to projection of lgo_mask
# 	x_size = output_array_pix_x_m
# 	y_size = output_array_pix_y_m
# 	(tulx, tuly, tulz ) = transform.TransformPoint( in_vx.gt[0], in_vx.gt[3])
# 	(turx, tury, turz ) = transform.TransformPoint( in_vx.gt[0] + in_vx.gt[1]*x_size, in_vx.gt[3])
# 	(tlrx, tlry, tlrz ) = transform.TransformPoint( in_vx.gt[0] + in_vx.gt[1]*x_size, in_vx.gt[3] + in_vx.gt[5]*y_size )
# 	(tllx, tlly, tllz ) = transform.TransformPoint( in_vx.gt[0], in_vx.gt[3] + in_vx.gt[5]*y_size )
	(tulx, tuly, tulz ) = transform.TransformPoint( utm_ulx, utm_uly, 0.0 )
	(turx, tury, turz ) = transform.TransformPoint( utm_lrx, utm_uly, 0.0 )
	(tlrx, tlry, tlrz ) = transform.TransformPoint( utm_lrx, utm_lry, 0.0 )
	(tllx, tlly, tllz ) = transform.TransformPoint( utm_ulx, utm_lry, 0.0 )
	#####################
	# Make .vrt for lgo image that is in local UTM to match Band 8 images 
	# - will reproject lgo mask in memory to local utm after reading in just this part of mask - .vrt limits the read size...
	#####################
	ulx=np.min([tulx,tllx])
	uly=np.max([tuly,tury])
	lrx=np.max([turx,tlrx])
	lry=np.min([tlly,tlry])
	
	path_to_lgo_mask = args.lgo_mask_file_dir+ '/' + args.lgo_mask_filename    # yes we opened this above, but just for srs - this full path is for vrt file
	output_lgo_vrt_name= 'tmp' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) + 'lgo.vrt'
	# -te needs xmin ymin xmax ymax - that is what is below, using the corners
	sp.call(["gdalbuildvrt", "-te", "%f"%ulx, "%f"%lry, "%f"%lrx, "%f"%uly, "-tr", "%9.4f"%(output_array_pix_x_m), "%9.4f"%(output_array_pix_x_m), "%s"%(args.VRT_dir + '/' + output_lgo_vrt_name), path_to_lgo_mask])  # import subprocess as sp
	lgo_mask_cropped_geoimg=GeoImg_noload(output_lgo_vrt_name,in_dir=args.VRT_dir)
# AM I CORRECT (11/3/16)? that this next read served no purpose? commenting it out here to see...
# 	lgo_mask_cropped_mask=lgo_mask_cropped_geoimg.gd.ReadAsArray().astype(np.float32)
	if lgo_mask_cropped_geoimg.nodata_value:
		lgo_mask_nodata=lgo_mask_cropped_geoimg.nodata_value
	else:
		lgo_mask_nodata=None
	os.remove(args.VRT_dir + '/' + output_lgo_vrt_name)  # get rid of the temporary vrt file now
	tmp_lgo_orig=None  # close GeoImg_noload object used on full original lgo_mask file, used to get srs above, but no input of data
	
	#####################
	# now reproject cropped input mask into local UTM so land mask can be used to identify
	# pixels that shouldn't move
	#####################
	mem_drv = gdal.GetDriverByName( 'MEM' )
	# 	inmem_lgo_mask_utm=mem_drv.CreateCopy('',in_vx.gd)   # actually don't need the data, but the rest is good
	inmem_lgo_mask_utm = mem_drv.Create('', output_array_num_pix_x, output_array_num_pix_y, 1, gdal.GDT_Float32)
	inmem_lgo_mask_utm.SetGeoTransform( [ output_array_ul_corner[0], output_array_pix_x_m, 0, output_array_ul_corner[1], 0, output_array_pix_y_m ] ) #
	inmem_lgo_mask_utm.SetProjection ( img1.proj )
	res = gdal.ReprojectImage( lgo_mask_cropped_geoimg.gd, inmem_lgo_mask_utm, lgo_mask_cropped_geoimg.srs.ExportToWkt(), img1.srs.ExportToWkt(), gdal.GRA_NearestNeighbour )
	lgo_mask_image_utm=inmem_lgo_mask_utm.ReadAsArray().astype(np.byte)
	
	
	print('got lgo mask data and reprojected to original image local utm')
	if not(args.nlf):
		log_output_lines.append('# got lgo mask data and reprojected to original image local utm')
		log_output_lines.append('# lgo mask reference (lgo mask file: %s) from lgo mask vrt file: %s\n'%(path_to_lgo_mask,output_lgo_vrt_name))



numcorr=0
##############################
# sc version of pycorr only uses hp - no choice!
##############################
img1_src_arr=hp_arr_1
img2_src_arr=hp_arr_2

offset_one_tenth_i=np.zeros([3,3],dtype="float")
offset_one_tenth_j=np.zeros([3,3],dtype="float")
offset_one_tenth_i[:,0]=-0.1
offset_one_tenth_j[0,:]=-0.1
offset_one_tenth_i[:,2]=0.1
offset_one_tenth_j[2,:]=0.1

offset_one_hundredth_i=0.1*offset_one_tenth_i
offset_one_hundredth_j=0.1*offset_one_tenth_j

np99p0=np.float32(99.0) # used in a comparison below - don't know if a variable reference is more efficient than calling np.float over and over, but try that - does interpreter only do this once???
new_half_target_chip=half_target_chip    # will only be done once unless use_wcs == True (now also no_speed_ref is False and if offset_correction_lgo_mask is used)...then smaller target chips are possible
new_cent_loc=cent_loc					# same here...  (in other words, need this second variable name if using wcs - this allows either)
if args.progupdates:
	updatestep = int((max_ind_i-min_ind_i)/10)
	progupdateintervals=range(min_ind_i+updatestep,max_ind_i,updatestep)
for i in range(min_ind_i,max_ind_i+1):
	if args.progupdates and (i in progupdateintervals):
		t_log('progress i at %4.1f%%'%(100.0*(float(i-min_ind_i)/float(max_ind_i))),outlogdisabled=args.nlf)
	for j in range(min_ind_j,max_ind_j+1):
		chip_src=img1_src_arr[(r_j_1[j]-half_source_chip):(r_j_1[j]+half_source_chip),(r_i_1[i]-half_source_chip):(r_i_1[i]+half_source_chip)].astype(np.float32)
# 		if (args.use_hp):	# hp images do have pix < 0, so need to look at original image for no data flag here...
# now test img1_data_mask - if np.min of chip is True, there are no 'nodata' pixels in chip_center_grid_xy
		test_for_masked_pix=np.min(img1_data_mask[(r_j_1[j]-half_source_chip):(r_j_1[j]+half_source_chip),(r_i_1[i]-half_source_chip):(r_i_1[i]+half_source_chip)])
# 		else:
# 			test_for_masked_pix=np.min(chip_src)
		if(test_for_masked_pix):  # all valid pixels in source chip - proceed
			numcorr+=1
			if not(args.no_speed_ref):
				speedrefvel=speedref_vel_vv[j][i]
				if speedref_vv_nodata:
					if (speedrefvel==speedref_vv_nodata):
# 						speedrefvel=0.0
						speedrefvel=args.nodatavmax # (will be reset to give half_target_chip below)
				elif ((speedrefvel < 0.0) or (speedrefvel == np99p0)):  
# 					speedrefvel = 0.0  # less than zero would be no data - but do it anyway right now.
					speedrefvel = args.nodatavmax  # less than zero would be no data and something (LISA?) at some point used 99.0 for no data - will be reset to give half_target_chip below)
# 				vmax=1.5 * speedrefvel / 365.25  # wcs speeds are in m/yr... we need m/day  And, allow for 50% greater speed... plus a few pixels
				vmax=1.5 * speedrefvel #  LISA in m/day  And, allow for 50% greater speed...
				maxpixoffset = np.int(np.ceil(vmax * del_t_days / img1.pix_x_m))  # changed to make result integer because of future subscript warning on line 895
				new_half_target_chip = half_source_chip + maxpixoffset + 5 # set target chip size to a few more than (1.5 * max v * interval)/pixsize
				if new_half_target_chip > half_target_chip:
					new_half_target_chip = half_target_chip   # max target chip size is specified half_target_chip - this ensures we stay in the image... (half_target_chip value was used in setting bounding box)
				new_cent_loc=new_half_target_chip-half_source_chip  # offset found in cv2 is referenced to returned array chip center, which has size full_source_chip - full_target_chip - so center is dif of half sizes.
# 				print wcs_vel.img[j][i],vmax,maxpixoffset,new_half_target_chip
				chip_tar=img2_src_arr[(r_j_2[j]-new_half_target_chip):(r_j_2[j]+new_half_target_chip),(r_i_2[i]-new_half_target_chip):(r_i_2[i]+new_half_target_chip)].astype(np.float32)
			elif args.offset_correction_lgo_mask and args.lgo_mask_limit_land_offset and (lgo_mask_image_utm[j][i]>=1):		# if there is a land mask and this particular source chip location is "land" (or ocean (2)), then make small search size (how far can land move?)
				new_half_target_chip = half_source_chip + 5 # set target chip size to 5 more than source chip size - image to image offset with no motion should not be larger than this...
				if new_half_target_chip > half_target_chip:
					new_half_target_chip = half_target_chip   # max target chip size is specified half_target_chip - this ensures we stay in the image... (half_target_chip value was used in setting bounding box)
				new_cent_loc=new_half_target_chip-half_source_chip  # offset found in cv2 is referenced to returned array chip center, which has size full_source_chip - full_target_chip - so center is dif of half sizes.
				chip_tar=img2_src_arr[(r_j_2[j]-new_half_target_chip):(r_j_2[j]+new_half_target_chip),(r_i_2[i]-new_half_target_chip):(r_i_2[i]+new_half_target_chip)].astype(np.float32)	
			else: # no speed reference, use -half_target_chip value from command line or default
				chip_tar=img2_src_arr[(r_j_2[j]-half_target_chip):(r_j_2[j]+half_target_chip),(r_i_2[i]-half_target_chip):(r_i_2[i]+half_target_chip)].astype(np.float32)
				new_cent_loc=cent_loc	# new_cent_loc is used below to calculate offset - need to reset it here in case it was changed before.  Could have been a bug in prior versions.
			resultCv=cv2.matchTemplate(chip_src, chip_tar, cv2.TM_CCOEFF_NORMED)
			mml_1=cv2.minMaxLoc(resultCv)
			peak1_corr=mml_1[1]
			peak1_loc=np.array(mml_1[3])
# 			peak1_corr=cv2.minMaxLoc(resultCv)[1]
# 			peak1_loc=cv2.minMaxLoc(resultCv)[3]
			tempres_peak1_zeroed=resultCv.copy()
			tempres_peak1_zeroed[(peak1_loc[1]-peak_block_halfsize):(peak1_loc[1]+peak_block_halfsize+1),(peak1_loc[0]-peak_block_halfsize):(peak1_loc[0]+peak_block_halfsize+1)]=0
			
			orig_int_del_i[j,i]=peak1_loc[0]-new_cent_loc
			orig_int_del_j[j,i]=peak1_loc[1]-new_cent_loc
			
			corr_arr[j,i]=peak1_corr
			del_corr_arr[j,i]=peak1_corr - cv2.minMaxLoc(tempres_peak1_zeroed)[1]
			# now if peak is far enough from the edges of the correlation surface, fit locale of peak with a spline and find sub-pixel offset at 1/10th of a pixel level
			if((np.array([n-rbs_halfsize for n in peak1_loc])>=np.array([0,0])).all() & 
			   (np.array([(n+rbs_halfsize) for n in peak1_loc ])<np.array(list(resultCv.shape))).all()):# 
				# at this point we have the i,j loc of integer pixel offset, need to chase peak at sub-pixel level
				rbs_p=RBS(range(-rbs_halfsize,rbs_halfsize+1),
						  range(-rbs_halfsize,rbs_halfsize+1),
						  resultCv[(peak1_loc[1]-rbs_halfsize):(peak1_loc[1]+rbs_halfsize+1),
								   (peak1_loc[0]-rbs_halfsize):(peak1_loc[0]+rbs_halfsize+1)],
						  kx=rbs_order,
						  ky=rbs_order)
				spoffi=0.0
				spoffj=0.0
				sphist=[]
				sphist.append((0,0,0.0,0.0,peak1_corr,0.0))
				count_ten=0
				peak_in_center=False
				while (not(peak_in_center) and (count_ten<10)):
					count_ten+=1
					# evaluate the eight pixels around present loc (spoffi,spoffj)
					tspsurf=rbs_p.ev((spoffj+offset_one_tenth_j).flatten(),(spoffi+offset_one_tenth_i).flatten())
					# find peak in this 3 x 3
					mml=cv2.minMaxLoc(tspsurf.reshape(offset_one_tenth_i.shape))
					if((mml[3][0]==1) and (mml[3][1]==1)):
						# peak stayed in center - done
						peak_in_center=True
						peak1_corr=mml[1]
						sarr=np.sort(tspsurf)
						sphist.append((mml[3][0],mml[3][1],spoffi,spoffj,mml[1],sarr[-1]-sarr[-2]))
					else:
						# peak not in center, move
						mmi=mml[3][0]
						mmj=mml[3][1]
						peak1_corr=mml[1]
						spoffi=spoffi+offset_one_tenth_i[mmj,mmi]
						spoffj=spoffj+offset_one_tenth_j[mmj,mmi]
						sarr=np.sort(tspsurf)
						sphist.append((mml[3][0],mml[3][1],spoffi,spoffj,mml[1],sarr[-1]-sarr[-2]))
				count_tenths[j,i]=count_ten
				###################################################
				# walk surface to find peak to hundredths of a pixel
				###################################################
# 				if (subpixel_hundredths):
				count_hund=0
				peak_in_center=False
				while (not(peak_in_center) and (count_hund<20)):
					count_hund+=1
					# evaluate the eight pixels around present loc (spoffi,spoffj)
					tspsurf=rbs_p.ev((spoffj+offset_one_hundredth_j).flatten(),(spoffi+offset_one_hundredth_i).flatten())
					# find peak is in this 3 x 3
					mml=cv2.minMaxLoc(tspsurf.reshape(offset_one_hundredth_i.shape))
					if((mml[3][0]==1) and (mml[3][1]==1)):
						# peak stayed in center - done
						peak_in_center=True
						peak1_corr=mml[1]
						sarr=np.sort(tspsurf)
						sphist.append((mml[3][0],mml[3][1],spoffi,spoffj,mml[1],sarr[-1]-sarr[-2]))
					else:
						# peak not in center, move
						mmi=mml[3][0]
						mmj=mml[3][1]
						peak1_corr=mml[1]
						spoffi=spoffi+offset_one_hundredth_i[mmj,mmi]
						spoffj=spoffj+offset_one_hundredth_j[mmj,mmi]
						sarr=np.sort(tspsurf)
						sphist.append((mml[3][0],mml[3][1],spoffi,spoffj,mml[1],sarr[-1]-sarr[-2]))
# 					if count_hund > 50:
# 						print 'count_hund {countvar:4d} exceeded 50'.format(countvar=count_hund)
				count_hundredths[j,i]=count_hund
				####################################################
				# end walk surface
				####################################################

				sp_deli=spoffi
				sp_delj=spoffj
# 				peak1_corr=mml[1]
				corr_arr[j,i]=peak1_corr
				del_corr_arr[j,i]=peak1_corr - cv2.minMaxLoc(tempres_peak1_zeroed)[1]
##########################################################
				deli=(peak1_loc[0]+sp_deli)-new_cent_loc
				delj=(peak1_loc[1]+sp_delj)-new_cent_loc
##########################################################
				d2i_pts=rbs_p.ev([spoffj,spoffj,spoffj],[spoffi-d2xdx2_spacing,spoffi,spoffi+d2xdx2_spacing])
				d2j_pts=rbs_p.ev([spoffj-d2xdx2_spacing,spoffj,spoffj+d2xdx2_spacing],[spoffi,spoffi,spoffi])
				d2idx2_t=(2.0*d2i_pts[1] - (d2i_pts[0]+d2i_pts[2]))/(np.power(d2xdx2_spacing,2.0))
				d2jdx2_t=(2.0*d2j_pts[1] - (d2j_pts[0]+d2j_pts[2]))/(np.power(d2xdx2_spacing,2.0))
######################## dump some peak-tracking histories to see how it is working...
# 				if((np.abs(j-425)<1) & (np.abs(i-640)<20) & (count_ten>=1)):
# 					print '\ni: %d j: %d count_ten: %d count_hund: %d orig_pl1 %d %d peak1_loc %d %d '%(i,j,count_ten,count_hund,orig_int_del_i[j,i],orig_int_del_j[j,i],peak1_loc[0]-new_cent_loc,peak1_loc[1]-new_cent_loc)
# 					for k in sphist: print '%d %d %4.2f %4.2f %8.6f %8.7f'%(k)
######################## #  and quit if you need too #					sys.exit(0)
			else: # peak is too close to an edge of correlation surface, don't do sub-pixel fit.
				sp_deli=0.0
				sp_delj=0.0
				deli=peak1_loc[0]-new_cent_loc
				delj=peak1_loc[1]-new_cent_loc
				d2idx2_t=curvature_nodata_value
				d2jdx2_t=curvature_nodata_value
# 		offset_dist_ij_arr[j,i]=np.sqrt(np.sum(np.power((cent_loc-peak1_loc[0],cent_loc-peak1_loc[1]),2.0)))
			offset_dist_ij_arr[j,i]=np.sqrt(np.sum(np.power((deli,delj),2.0)))
			del_i[j,i]=deli
			del_j[j,i]=delj
			sp_dist_i[j,i]=sp_deli
			sp_dist_j[j,i]=sp_delj
			d2idx2[j,i]=d2idx2_t
			d2jdx2[j,i]=d2jdx2_t
		else: # all pixels zero in chip_src - mask this one
			corr_arr[j,i]=corr_nodata_val
			del_corr_arr[j,i]=corr_nodata_val
			offset_dist_ij_arr[j,i]=0.0
			del_i[j,i]=0
			del_j[j,i]=0
			sp_dist_i[j,i]=0
			sp_dist_j[j,i]=0
			d2idx2[j,i]=curvature_nodata_value
			d2jdx2[j,i]=curvature_nodata_value
		d2dx2_mean[j,i]=(d2idx2[j,i] + d2jdx2[j,i])/2.0
#   	if not(args.npb):
#   		pbar.update(i)
# if not(args.npb):
# 	pbar.finish()

t_log('Done with image to image correlation: ' + np.str(numcorr) + ' correlations of ' + np.str(output_array_num_pix_x*output_array_num_pix_y) + ' possible',outlogdisabled=args.nlf)

# done with hp image arrays now - free all pointers that reference them and the memory
# del chip_src
# del chip_tar
del img1_src_arr
del hp_arr_1
del img2_src_arr
del hp_arr_2


###########################################################################
#
#  offset correction section
#
###########################################################################
found_valid_offset=False
found_valid_bilinear_offset=False # flag that will be used below to decide to apply bilinear surfaces as offsets...
offset_correction_type_applied=''
offset_correction_type_descritption=''

if args.offset_correction_speedref or args.offset_correction_lgo_mask:
	# place holder flags for the different types of possible offsets - defines variable to prevent errors in if X: statements
	zspeed_offset_available=None
	midspeed_offset_available=None
	lgo_masked_offset_available=None
	lgo_bilinear_masked_offset_available=None
	full_stationary_masked_offset_available=None
	full_stationary_bilinear_masked_offset_available=None
	full_vector_masked_offset_available=None
	full_vector_bilinear_masked_offset_available=None
	
	slow_area_zero_speed_myr=20.0 # m/yr for zero_speed correction
	slow_area_max_vector_speed_myr=40.0 # m/yr for mid_speed correction
	slow_area_zero_speed_md=slow_area_zero_speed_myr/365.25
	slow_area_max_vector_speed_md=slow_area_max_vector_speed_myr/365.25
	
	zero_speed_min_num_pix=500
	zero_speed_min_percent_valid_pix_available=0.5
	
	mid_speed_min_num_pix=4000
	mid_speed_min_percent_valid_pix_available=5.0

	lgo_masked_min_num_pix=500
	lgo_masked_min_percent_valid_pix_available=0.5
	
	full_stationary_min_num_pix=1000
	full_stationary_min_percent_valid_pix_available=2.0
	
	full_vector_min_num_pix=1000
	full_vector_min_percent_valid_pix_available=2.0

	corr_val_for_offset=0.3
	delcorr_min_for_offset = 0.15

	
	if args.offset_correction_speedref:
		print('attempting to find offset correction for slow areas')
		# 	slow_area_zero_speed_pixels=(del_t_days*slow_area_zero_speed_md)/img1.pix_x_m
		# 	slow_area_max_vector_speed_pixels=(del_t_days*slow_area_max_vector_speed_md)/img1.pix_x_m
		# 	if args.wcs=="":
		# 		print 'offset_correction requires -wcs use, exiting...'
		# 		sys.exit(-1)
		# 	else:
		# 		##########  find offset correction and make new output images...
		# 		# del_t_val=32.0  set from input images...
		# 		# make masked arrays of offsets, and pull wcs speeds for the unmasked points to find slow areas to correct image offsets to 0
		
		# first - 
		# Develop Haran-type filters for isolated points/sd limits if not isolated - using the speeds tracked here.  Then apply the masks to the masked i and j offsets above, and find estimated i and j offset corrections 
		# for the two (zero and slow) areas by combining the two masks and determining offsets from remaining pixels.  This approach is conservative in that it rejects a lot of data - hoping to make offsets accurate.
		#
		prelim_vv=((offset_dist_ij_arr*img1.pix_x_m)/del_t_val)
		prelim_vv_ma=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset),prelim_vv)

		ijnotmasked=np.where(prelim_vv_ma.mask==False)    # mask is true where masked, valid data where false

		std_vv_ratio=np.zeros_like(prelim_vv)

		min_std_allowed = 0.001  # prevents division by zero if all velocities in neighborhood are identical, which would be really suspicious anyway.. Haran used 0.01

		for cpi,cpj in zip(ijnotmasked[0],ijnotmasked[1]):
			i_s=[cpi-1,cpi-1,cpi-1,cpi,cpi,cpi+1,cpi+1,cpi+1]
			j_s=[cpj-1,cpj,cpj+1,cpj-1,cpj+1,cpj-1,cpj,cpj+1]
			nn=np.sum(prelim_vv_ma[(i_s,j_s)].mask==False)
			if nn>1:
				cpstd=np.std(prelim_vv_ma[(i_s,j_s)])
				if cpstd > min_std_allowed:
					delv=prelim_vv_ma[cpi,cpj]-np.mean(prelim_vv_ma[(i_s,j_s)])
					std_vv_ratio[cpi,cpj]= np.abs(delv)/cpstd
				else:
					std_vv_ratio[cpi,cpj]=0.0
			elif nn==1:
				if np.abs(prelim_vv_ma[cpi,cpj]-np.mean(prelim_vv_ma[(i_s,j_s)])) < 1.0:
					std_vv_ratio[cpi,cpj]=1.0  # set so this point survives the mask at the end (std_vv_ratio <= 3.0) - Haran used points with single neighbors that were within a m/d
				else:
					std_vv_ratio[cpi,cpj]=0.0
			else:
				std_vv_ratio[cpi,cpj]=0.0

		prelim_vv_ma2=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset) | (std_vv_ratio==0.0) | (std_vv_ratio > 3.0), prelim_vv)

		# second - find pixels where speed ref is below slow_area_zero_speed_md
		if speedref_vv_nodata:  # input speed reference has a specified no_data value.
			zspeed_d_i_m=np.ma.masked_where((speedref_vel_vv>=slow_area_zero_speed_md)|(speedref_vel_vv==speedref_vv_nodata)|(del_i==0.0)|(prelim_vv_ma2.mask==True),del_i)
			zspeed_d_j_m=np.ma.masked_where((speedref_vel_vv>=slow_area_zero_speed_md)|(speedref_vel_vv==speedref_vv_nodata)|(del_j==0.0)|(prelim_vv_ma2.mask==True),del_j)
			zspeed_num_possible_pix=np.count_nonzero(np.array((speedref_vel_vv<slow_area_zero_speed_md)&(speedref_vel_vv!=speedref_vv_nodata)&(del_i!=0.0)))
			zspeed_num_valid_pix=np.count_nonzero(np.array((speedref_vel_vv<slow_area_zero_speed_md)&(speedref_vel_vv!=speedref_vv_nodata)&(del_i!=0.0)&(prelim_vv_ma2.mask==False)))
		else:
			zspeed_d_i_m=np.ma.masked_where((speedref_vel_vv>=slow_area_zero_speed_md)|(speedref_vel_vv==0.0)|(speedref_vel_vv==-99.0)|(del_i==0.0)|(prelim_vv_ma2.mask==True),del_i)
			zspeed_d_j_m=np.ma.masked_where((speedref_vel_vv>=slow_area_zero_speed_md)|(speedref_vel_vv==0.0)|(speedref_vel_vv==-99.0)|(del_j==0.0)|(prelim_vv_ma2.mask==True),del_j)
			zspeed_num_possible_pix=np.count_nonzero(np.array((speedref_vel_vv<slow_area_zero_speed_md)&(speedref_vel_vv!=0.0)&(speedref_vel_vv!=-99.0)&(del_i!=0.0)))
			zspeed_num_valid_pix=np.count_nonzero(np.array((speedref_vel_vv<slow_area_zero_speed_md)&(speedref_vel_vv!=0.0)&(speedref_vel_vv!=-99.0)&(del_i!=0.0)&(prelim_vv_ma2.mask==False)))
		# third - find pixels where speed ref is between slow_area_zero_speed_md and slow_area_max_vector_speed_md
		
		#########################################################################################################
		#### THIS WAS NOT CORRECT FOR RUNS BEFORE 11/6/2016 - midspeed top cutoff was misapplied - masked speeds less than top cutoff, rather than greater. 
		#### 
		#########################################################################################################
		
############################################		
#	This was the code prior to v5p8 - now midspeed includes all speeds down to zero, not down to the zero_speed cutoff.  Vector matching can include the slow speeds.
#
# 		if speedref_vv_nodata:
# 			midspeed_d_i_m=np.ma.masked_where((speedref_vel_vv<slow_area_zero_speed_md)|(slow_area_max_vector_speed_md<speedref_vel_vv)| \
# 											  (speedref_vel_vv==speedref_vv_nodata)|((del_i==0.0))|(prelim_vv_ma2.mask==True),del_i)
# 			midspeed_d_j_m=np.ma.masked_where((speedref_vel_vv<slow_area_zero_speed_md)|(slow_area_max_vector_speed_md<speedref_vel_vv)| \
# 											  (speedref_vel_vv==speedref_vv_nodata)|((del_j==0.0))|(prelim_vv_ma2.mask==True),del_j)
# 			midspeed_num_possible_pix=np.count_nonzero(np.array((speedref_vel_vv>=slow_area_zero_speed_md)&(slow_area_max_vector_speed_md>=speedref_vel_vv)&(speedref_vel_vv!=speedref_vv_nodata)&(del_i!=0.0), dtype=np.bool))
# 			midspeed_num_valid_pix=np.count_nonzero(np.array((speedref_vel_vv>=slow_area_zero_speed_md)&(slow_area_max_vector_speed_md>=speedref_vel_vv)&(speedref_vel_vv!=speedref_vv_nodata)&(del_i!=0.0)&(prelim_vv_ma2.mask==False), dtype=np.bool))
# 		else:
# 			midspeed_d_i_m=np.ma.masked_where((speedref_vel_vv<slow_area_zero_speed_md)|(slow_area_max_vector_speed_md<speedref_vel_vv)| \
# 											  (speedref_vel_vv==0.0)|(speedref_vel_vv==-99.0)|((del_i==0.0))|(prelim_vv_ma2.mask==True),del_i)
# 			midspeed_d_j_m=np.ma.masked_where((speedref_vel_vv<slow_area_zero_speed_md)|(slow_area_max_vector_speed_md<speedref_vel_vv)| \
# 											  (speedref_vel_vv==0.0)|(speedref_vel_vv==-99.0)|((del_j==0.0))|(prelim_vv_ma2.mask==True),del_j)
# 			midspeed_num_possible_pix=np.count_nonzero(np.array((speedref_vel_vv>=slow_area_zero_speed_md)&(slow_area_max_vector_speed_md>=speedref_vel_vv)&(speedref_vel_vv!=0.0)&(speedref_vel_vv!=-99.0)&(del_i!=0.0), dtype=np.bool))
# 			midspeed_num_valid_pix=np.count_nonzero(np.array((speedref_vel_vv>=slow_area_zero_speed_md)&(slow_area_max_vector_speed_md>=speedref_vel_vv)&(speedref_vel_vv!=0.0)&(speedref_vel_vv!=-99.0)&(del_i!=0.0)&(prelim_vv_ma2.mask==False), dtype=np.bool))
############################################

		if speedref_vv_nodata:
			midspeed_d_i_m=np.ma.masked_where((slow_area_max_vector_speed_md<speedref_vel_vv)| \
											  (speedref_vel_vv==speedref_vv_nodata)|((del_i==0.0))|(prelim_vv_ma2.mask==True),del_i)
			midspeed_d_j_m=np.ma.masked_where((slow_area_max_vector_speed_md<speedref_vel_vv)| \
											  (speedref_vel_vv==speedref_vv_nodata)|((del_j==0.0))|(prelim_vv_ma2.mask==True),del_j)
			midspeed_num_possible_pix=np.count_nonzero(np.array((slow_area_max_vector_speed_md>=speedref_vel_vv)&(speedref_vel_vv!=speedref_vv_nodata)&(del_i!=0.0)))
			midspeed_num_valid_pix=np.count_nonzero(np.array((slow_area_max_vector_speed_md>=speedref_vel_vv)&(speedref_vel_vv!=speedref_vv_nodata)&(del_i!=0.0)&(prelim_vv_ma2.mask==False)))
		else:
			midspeed_d_i_m=np.ma.masked_where((slow_area_max_vector_speed_md<speedref_vel_vv)| \
											  (speedref_vel_vv==0.0)|(speedref_vel_vv==-99.0)|((del_i==0.0))|(prelim_vv_ma2.mask==True),del_i)
			midspeed_d_j_m=np.ma.masked_where((slow_area_max_vector_speed_md<speedref_vel_vv)| \
											  (speedref_vel_vv==0.0)|(speedref_vel_vv==-99.0)|((del_j==0.0))|(prelim_vv_ma2.mask==True),del_j)
			midspeed_num_possible_pix=np.count_nonzero(np.array((slow_area_max_vector_speed_md>=speedref_vel_vv)&(speedref_vel_vv!=0.0)&(speedref_vel_vv!=-99.0)&(del_i!=0.0)))
			midspeed_num_valid_pix=np.count_nonzero(np.array((slow_area_max_vector_speed_md>=speedref_vel_vv)&(speedref_vel_vv!=0.0)&(speedref_vel_vv!=-99.0)&(del_i!=0.0)&(prelim_vv_ma2.mask==False)))



		# fourth - generate offsets for both speed categories
		# 	zspeed_num_valid_pix=len(zspeed_d_i_m.compressed())
		if zspeed_num_valid_pix>0:
			zspeed_offset_i = -(np.median(zspeed_d_i_m.compressed()))
			zspeed_offset_j = -(np.median(zspeed_d_j_m.compressed()))
			zspeed_offset_available=True

		# 	midspeed_num_valid_pix=len(midspeed_d_i_m.compressed())
		#########################################################################################################
		#### THIS WAS NOT CORRECT FOR RUNS BEFORE 8/26/2016 - should have scaled the ref vx and vy by del_t_days to get total offset in i or j pixels, but omitted the del_t_days term, so just one day's motion subtracted... oops! 
		#### this should have been an issue with Antarctic data processed with sc_pycorr before the reprocessing of October 2016, which was the whole dataset and should replace the earlier data with any issues from this omission.
		#########################################################################################################
		if midspeed_num_valid_pix>0:
			midspeed_offset_i = -(np.median((midspeed_d_i_m - (del_t_days*speedref_vel_vx/img1.pix_x_m)).compressed()))
# 			midspeed_offset_j = -(np.median((midspeed_d_j_m - (del_t_days*speedref_vel_vy/img1.pix_x_m)).compressed()))
			midspeed_offset_j = -(np.median((midspeed_d_j_m - (del_t_days*speedref_vel_vy/img1.pix_y_m)).compressed()))
			midspeed_offset_available=True

		if zspeed_num_valid_pix>0:
			print('found zero_speed offset correction (del_i: %f del_j: %f) in slow moving areas (less than %f m/yr) using %d pixels out of %d possible'%\
					(zspeed_offset_i,zspeed_offset_j,slow_area_zero_speed_myr,zspeed_num_valid_pix,zspeed_num_possible_pix))
		if midspeed_num_valid_pix>0:
			print('found mid_speed offset correction (del_i: %f del_j: %f) in slow moving areas (more than %f m/yr and less than %f m/yr) using %d pixels out of %d possible'%\
					(midspeed_offset_i,midspeed_offset_j,slow_area_zero_speed_myr,slow_area_max_vector_speed_myr,midspeed_num_valid_pix,midspeed_num_possible_pix))
		if not(args.nlf):
			if zspeed_num_valid_pix>0:
				t_log('found zero_speed offset correction (del_i: %f del_j: %f) in slow moving areas (less than %f m/yr) using %d pixels out of %d possible (%f %%)'%\
					(zspeed_offset_i,zspeed_offset_j,slow_area_zero_speed_myr,zspeed_num_valid_pix,zspeed_num_possible_pix,100.0*(zspeed_num_valid_pix/float(zspeed_num_possible_pix))),outlogdisabled=args.nlf)
			if midspeed_num_valid_pix>0:
				t_log('found mid_speed offset correction (del_i: %f del_j: %f) in slow moving areas (more than %f m/yr and less than %f m/yr) using %d pixels out of %d possible (%f %%)'%\
				(midspeed_offset_i,midspeed_offset_j,slow_area_zero_speed_myr,slow_area_max_vector_speed_myr,midspeed_num_valid_pix,midspeed_num_possible_pix,100.0*(midspeed_num_valid_pix/float(midspeed_num_possible_pix))),outlogdisabled=args.nlf)


	if args.offset_correction_lgo_mask:   # before 5p6 this was an elif - but now do this section even if you did the last, so both can be combined below if successful
		print('attempting to find offset correction for land (lgo mask) areas')
		# mask pixels that are not land not(lgo==1)
		if lgo_mask_nodata:  # input speed reference has a specified no_data value.
			lgo_masked_d_i_m=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset) | \
			 								(lgo_mask_image_utm!=1)|(lgo_mask_image_utm==lgo_mask_nodata),del_i)
			lgo_masked_d_j_m=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset) | \
			 								(lgo_mask_image_utm!=1)|(lgo_mask_image_utm==lgo_mask_nodata),del_j)
			lgo_masked_num_possible_pix=np.count_nonzero(np.array((lgo_mask_image_utm==1)&(lgo_mask_image_utm!=lgo_mask_nodata)&(corr_arr!=corr_nodata_val)))
			# 			lgo_masked_num_valid_pix=np.count_nonzero(np.array((speedref_vel_vv<slow_area_zero_speed_md)&(speedref_vel_vv!=speedref_vv_nodata)&(del_i!=0.0)&(prelim_vv_ma2.mask==False), dtype=np.bool))
			lgo_masked_num_valid_pix=np.count_nonzero(np.array(lgo_masked_d_i_m.mask==False))
		else:   # lgo mask did not include a nodata value (so don't try to use it...)
			lgo_masked_d_i_m=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset) | \
			 								(lgo_mask_image_utm!=1),del_i)
			lgo_masked_d_j_m=np.ma.masked_where((corr_arr==corr_nodata_val) | (del_corr_arr < delcorr_min_for_offset) | (corr_arr < corr_val_for_offset) | \
			 								(lgo_mask_image_utm!=1),del_j)
			lgo_masked_num_possible_pix=np.count_nonzero(np.array((lgo_mask_image_utm==1)&(corr_arr!=corr_nodata_val)))
			# 			lgo_masked_num_valid_pix = np.count_nonzero(np.array((speedref_vel_vv<slow_area_zero_speed_md) & (speedref_vel_vv!=speedref_vv_nodata) & (del_i!=0.0) & (prelim_vv_ma2.mask==False), dtype=np.bool))
			lgo_masked_num_valid_pix=np.count_nonzero(np.array(lgo_masked_d_i_m.mask==False))
			
# 		lgo_masked_offset_available=False
		if lgo_masked_num_valid_pix>0:
			lgo_masked_offset_i = -(np.median(lgo_masked_d_i_m.compressed()))
			lgo_masked_offset_j = -(np.median(lgo_masked_d_j_m.compressed()))
			lgo_masked_stdev_i = np.std(lgo_masked_d_i_m.compressed())
			lgo_masked_stdev_j = np.std(lgo_masked_d_j_m.compressed())
			lgo_masked_offset_available=True

			print('found lgo_mask (land pixel) offset correction (del_i: %f del_j: %f std_i %f std_j %f) using %d pixels out of %d possible'%\
					(lgo_masked_offset_i,lgo_masked_offset_j,lgo_masked_stdev_i,lgo_masked_stdev_j,lgo_masked_num_valid_pix,lgo_masked_num_possible_pix))
			if not(args.nlf):
				if lgo_masked_num_valid_pix>0:
					t_log('found lgo_mask (land pixel) offset correction (del_i: %f del_j: %f std_i %f std_j %f) using %d pixels out of %d possible'%\
					(lgo_masked_offset_i,lgo_masked_offset_j,lgo_masked_stdev_i,lgo_masked_stdev_j,lgo_masked_num_valid_pix,lgo_masked_num_possible_pix),outlogdisabled=args.nlf)
			
			##############################################################
			# bilinear offset correction - if you have enough valid pixels
			##############################################################
			# lgo_bilinear_masked_offset_available=False
			if args.offset_correction_bilinear_fit and (lgo_masked_num_valid_pix>=args.offset_correction_bilinear_fit_min_pixels):
				lgo_bilinear_masked_offset_available=True
				ix,iy=lgo_masked_d_i_m.nonzero()
				coefs_i=polyfit2d(ix, iy, lgo_masked_d_i_m[ix,iy], [1,1])
				coefs_j=polyfit2d(ix, iy, lgo_masked_d_j_m[ix,iy], [1,1])
				Ind_arr = np.array(np.meshgrid(range(del_i.shape[1]),range(del_i.shape[0])))
				lgo_i_bilinear_offset_corr_arr = -1.0 * polynomial.polyval2d(Ind_arr[1], Ind_arr[0], coefs_i)
				lgo_j_bilinear_offset_corr_arr = -1.0 * polynomial.polyval2d(Ind_arr[1], Ind_arr[0], coefs_j)


	if args.Greenland:
		######################################################################################################################################################
		# special case for GREENLAND where both velocity field and lgo mask can be used together to calculate offset - otherwise use only one or none below...
		######################################################################################################################################################
		# if args.offset_correction_speedref and args.offset_correction_lgo_mask:  could do this, but actually at this point have enough info to do this:
		if zspeed_offset_available and lgo_masked_offset_available:
			print('attempting to find offset correction for combined slow and land masked areas')
			# did mask these above, but now need to do the combined mask because some points will be in both masks...
			full_stationary_point_mask = lgo_masked_d_i_m.mask & zspeed_d_i_m.mask
			full_stationary_point_mask_masked_d_i_m=np.ma.masked_where(full_stationary_point_mask,del_i)
			full_stationary_point_mask_masked_d_j_m=np.ma.masked_where(full_stationary_point_mask,del_j)

			if speedref_vv_nodata and lgo_mask_nodata:  # input speed reference has a specified no_data value, and lgo mask does as well...
				full_stationary_masked_num_possible_pix=np.count_nonzero(np.array((((speedref_vel_vv<slow_area_zero_speed_md)&(speedref_vel_vv!=speedref_vv_nodata))|((lgo_mask_image_utm==1)&(lgo_mask_image_utm!=lgo_mask_nodata)))&(del_i!=0.0)& \
																				  (corr_arr!=corr_nodata_val)))
			elif speedref_vv_nodata:
				full_stationary_masked_num_possible_pix=np.count_nonzero(np.array((((speedref_vel_vv<slow_area_zero_speed_md)&(speedref_vel_vv!=speedref_vv_nodata))|(lgo_mask_image_utm==1))&(del_i!=0.0)& \
																				  (corr_arr!=corr_nodata_val)))
			else:
				full_stationary_masked_num_possible_pix=np.count_nonzero(np.array((((speedref_vel_vv<slow_area_zero_speed_md)&(speedref_vel_vv!=0.0)&(speedref_vel_vv!=-99.0))|(lgo_mask_image_utm==1))&(del_i!=0.0)& \
																				  (corr_arr!=corr_nodata_val)))
																			  
			full_stationary_masked_num_valid_pix=np.count_nonzero(np.array(full_stationary_point_mask==False))

	# 		full_stationary_masked_offset_available=False
			if full_stationary_masked_num_valid_pix>0:
				full_stationary_masked_offset_i = -(np.median(full_stationary_point_mask_masked_d_i_m.compressed()))
				full_stationary_masked_offset_j = -(np.median(full_stationary_point_mask_masked_d_j_m.compressed()))
				full_stationary_masked_stdev_i = np.std(full_stationary_point_mask_masked_d_i_m.compressed())
				full_stationary_masked_stdev_j = np.std(full_stationary_point_mask_masked_d_j_m.compressed())
				full_stationary_masked_offset_available=True

				print('found full_stationary_mask (land and zero_speed) offset correction (del_i: %f del_j: %f std_i %f std_j %f) using %d pixels out of %d possible'%\
						(full_stationary_masked_offset_i,full_stationary_masked_offset_j,full_stationary_masked_stdev_i,full_stationary_masked_stdev_j,full_stationary_masked_num_valid_pix,full_stationary_masked_num_possible_pix))
				if not(args.nlf):
					if full_stationary_masked_num_valid_pix>0:
						t_log('found full_stationary_mask (land and zero_speed) offset correction (del_i: %f del_j: %f std_i %f std_j %f) using %d pixels out of %d possible'%\
						(full_stationary_masked_offset_i,full_stationary_masked_offset_j,full_stationary_masked_stdev_i,full_stationary_masked_stdev_j,full_stationary_masked_num_valid_pix,full_stationary_masked_num_possible_pix),outlogdisabled=args.nlf)
			
				##############################################################
				# bilinear offset correction - if you have enough valid pixels
				##############################################################
				full_stationary_bilinear_masked_offset_available=False
				if args.offset_correction_bilinear_fit and (full_stationary_masked_num_valid_pix>=args.offset_correction_bilinear_fit_min_pixels):
					full_stationary_bilinear_masked_offset_available=True
					ix,iy=full_stationary_point_mask_masked_d_i_m.nonzero()
					coefs_i=polyfit2d(ix, iy, full_stationary_point_mask_masked_d_i_m[ix,iy], [1,1])
					coefs_j=polyfit2d(ix, iy, full_stationary_point_mask_masked_d_j_m[ix,iy], [1,1])
					Ind_arr = np.array(np.meshgrid(range(del_i.shape[1]),range(del_i.shape[0])))
					full_stationary_i_bilinear_offset_corr_arr = -1.0 * polynomial.polyval2d(Ind_arr[1], Ind_arr[0], coefs_i)
					full_stationary_j_bilinear_offset_corr_arr = -1.0 * polynomial.polyval2d(Ind_arr[1], Ind_arr[0], coefs_j)



		######################################################################################################################################################
		# SECOND part of special case for GREENLAND where both velocity field and lgo mask can be used together to calculate offset - otherwise use only one or none below...
		#  this combines land and medium and low speed pixels through vector fits to find offsets
		######################################################################################################################################################
		# if args.offset_correction_speedref and args.offset_correction_lgo_mask:  could do this, but actually at this point have enough info to do this:
		if (midspeed_offset_available or (midspeed_offset_available and zspeed_offset_available)) and lgo_masked_offset_available:
			print('attempting to find vector-tied offset correction for combined medium, low and land masked areas')
			# did mask these above, but now need to do the combined mask because some points will be in both masks...
		
			# make sure vv ref mz speeds over land have zero i and j components, then find vector field offsets
			# ok to do this here because it is the last use of speedref_ vector components - only speedref_vel_vv used elsewhere below?
			if speedref_vv_nodata and (speedref_vv_nodata!=0.0):
				speedref_vel_vx[lgo_mask_image_utm==1]=0.0
				speedref_vel_vy[lgo_mask_image_utm==1]=0.0
				speedref_vel_vv[lgo_mask_image_utm==1]=0.0
			else:
				speedref_vel_vx[lgo_mask_image_utm==1]=0.000000001
				speedref_vel_vy[lgo_mask_image_utm==1]=0.000000001   # zero might be nodata - is tested that way below if speedref_vv_nodata not set?
				speedref_vel_vv[lgo_mask_image_utm==1]=0.000000001

			full_vector_point_mask = lgo_masked_d_i_m.mask & zspeed_d_i_m.mask & midspeed_d_i_m.mask
		
		
	# 		if speedref_vv_nodata:
	# 			full_vector_point_mask_masked_d_i_m=np.ma.masked_where((full_vector_point_mask | (speedref_vel_vv==speedref_vv_nodata)),(midspeed_d_i_m - (del_t_days*speedref_vel_vx/img1.pix_x_m)))
	# 			full_vector_point_mask_masked_d_j_m=np.ma.masked_where((full_vector_point_mask | (speedref_vel_vv==speedref_vv_nodata)),(midspeed_d_j_m - (del_t_days*speedref_vel_vy/img1.pix_x_m)))
	# 		else:
	# 			full_vector_point_mask_masked_d_i_m=np.ma.masked_where((full_vector_point_mask | (speedref_vel_vv!=0.0) | (speedref_vel_vv==-99.0)),(midspeed_d_i_m - (del_t_days*speedref_vel_vx/img1.pix_x_m)))
	# 			full_vector_point_mask_masked_d_j_m=np.ma.masked_where((full_vector_point_mask | (speedref_vel_vv!=0.0) | (speedref_vel_vv==-99.0)),(midspeed_d_j_m - (del_t_days*speedref_vel_vy/img1.pix_x_m)))
			if speedref_vv_nodata:
				full_vector_point_mask_masked_d_i_m=np.ma.masked_where((full_vector_point_mask | (speedref_vel_vv==speedref_vv_nodata)),(del_i - (del_t_days*speedref_vel_vx/img1.pix_x_m)))
				full_vector_point_mask_masked_d_j_m=np.ma.masked_where((full_vector_point_mask | (speedref_vel_vv==speedref_vv_nodata)),(del_j - (del_t_days*speedref_vel_vy/img1.pix_y_m)))
			else:
				full_vector_point_mask_masked_d_i_m=np.ma.masked_where((full_vector_point_mask | (speedref_vel_vv!=0.0) | (speedref_vel_vv==-99.0)),(del_i - (del_t_days*speedref_vel_vx/img1.pix_x_m)))
				full_vector_point_mask_masked_d_j_m=np.ma.masked_where((full_vector_point_mask | (speedref_vel_vv!=0.0) | (speedref_vel_vv==-99.0)),(del_j - (del_t_days*speedref_vel_vy/img1.pix_y_m)))

			if speedref_vv_nodata and lgo_mask_nodata:  # input speed reference has a specified no_data value, and lgo mask does as well...
				full_vector_masked_num_possible_pix=np.count_nonzero(np.array((((speedref_vel_vv<slow_area_max_vector_speed_md)&(speedref_vel_vv!=speedref_vv_nodata))|((lgo_mask_image_utm==1)&(del_i!=0.0)& \
																				  (corr_arr!=corr_nodata_val)))))
			elif speedref_vv_nodata:
				full_vector_masked_num_possible_pix=np.count_nonzero(np.array((((speedref_vel_vv<slow_area_max_vector_speed_md)&(speedref_vel_vv!=speedref_vv_nodata))|((lgo_mask_image_utm==1)&(del_i!=0.0)& \
																				  (corr_arr!=corr_nodata_val)))))
			else:
				full_vector_masked_num_possible_pix=np.count_nonzero(np.array((((speedref_vel_vv<slow_area_max_vector_speed_md)&(speedref_vel_vv!=0.0)&(speedref_vel_vv!=-99.0))|((lgo_mask_image_utm==1)&(del_i!=0.0)& \
																				  (corr_arr!=corr_nodata_val)))))
																			  
			full_vector_masked_num_valid_pix=np.count_nonzero(np.array(full_vector_point_mask==False))

	# 		full_vector_masked_offset_available=False
			if full_vector_masked_num_valid_pix>0:
				full_vector_masked_offset_i = -(np.median(full_vector_point_mask_masked_d_i_m.compressed()))
				full_vector_masked_offset_j = -(np.median(full_vector_point_mask_masked_d_j_m.compressed()))
				full_vector_masked_stdev_i = np.std(full_vector_point_mask_masked_d_i_m.compressed())
				full_vector_masked_stdev_j = np.std(full_vector_point_mask_masked_d_j_m.compressed())
				full_vector_masked_offset_available=True

				print('found full_vector_mask (land and zero_speed) offset correction (del_i: %f del_j: %f std_i %f std_j %f) using %d pixels out of %d possible'%\
						(full_vector_masked_offset_i,full_vector_masked_offset_j,full_vector_masked_stdev_i,full_vector_masked_stdev_j,full_vector_masked_num_valid_pix,full_vector_masked_num_possible_pix))
				if not(args.nlf):
					if full_vector_masked_num_valid_pix>0:
						t_log('found full_vector_mask (land and zero_speed) offset correction (del_i: %f del_j: %f std_i %f std_j %f) using %d pixels out of %d possible'%\
						(full_vector_masked_offset_i,full_vector_masked_offset_j,full_vector_masked_stdev_i,full_vector_masked_stdev_j,full_vector_masked_num_valid_pix,full_vector_masked_num_possible_pix),outlogdisabled=args.nlf)
			
				##############################################################
				# bilinear offset correction - if you have enough valid pixels
				##############################################################
				full_vector_bilinear_masked_offset_available=False
				if args.offset_correction_bilinear_fit and (full_vector_masked_num_valid_pix>=args.offset_correction_bilinear_fit_min_pixels):
					full_vector_bilinear_masked_offset_available=True
					ix,iy=full_vector_point_mask_masked_d_i_m.nonzero()
					coefs_i=polyfit2d(ix, iy, full_vector_point_mask_masked_d_i_m[ix,iy], [1,1])
					coefs_j=polyfit2d(ix, iy, full_vector_point_mask_masked_d_j_m[ix,iy], [1,1])
					Ind_arr = np.array(np.meshgrid(range(del_i.shape[1]),range(del_i.shape[0])))
					full_vector_i_bilinear_offset_corr_arr = -1.0 * polynomial.polyval2d(Ind_arr[1], Ind_arr[0], coefs_i)
					full_vector_j_bilinear_offset_corr_arr = -1.0 * polynomial.polyval2d(Ind_arr[1], Ind_arr[0], coefs_j)



	# finally - decide which offset to use...  Haran used 0.5 % for zspeed and 2.0 % for slow...
	found_valid_offset=False
	found_valid_bilinear_offset=False # flag that will be used below to decide to apply bilinear surfaces as offsets...


	if full_vector_bilinear_masked_offset_available and (full_vector_masked_num_valid_pix>args.offset_correction_bilinear_fit_min_pixels) and full_vector_masked_offset_available and (full_vector_masked_num_valid_pix>full_vector_min_num_pix) and ((float(full_vector_masked_num_valid_pix)/full_vector_masked_num_possible_pix)>=full_vector_min_percent_valid_pix_available/100.0) and (np.sqrt(full_vector_masked_offset_i**2.0 + full_vector_masked_offset_j**2.0)<args.max_allowable_pixel_offset_correction):
		# use bilinear lgo_masked corection - but fields applied below, where i and j dc offsets are applied for the other corrections...
		final_offset_correction_i=full_vector_masked_offset_i	# keep these anyway, to compare to bilinear solution offsets 
		final_offset_correction_j=full_vector_masked_offset_j
		i_bilinear_offset_corr_arr=full_vector_i_bilinear_offset_corr_arr
		j_bilinear_offset_corr_arr=full_vector_j_bilinear_offset_corr_arr
		found_valid_offset=True
		found_valid_bilinear_offset=True
		offset_correction_type_applied='full_vector_bilinear_masked_correction'
		offset_correction_type_descritption='full_vector_bilinear_masked_correction fields for i and j for land pixels,' + ' - using bilinear coefs_i ' + str(coefs_i.ravel()) + ' coefs_j ' + str(coefs_j.ravel()) +' %d valid pixels out of %d possible for scene (%f %%)'%(full_vector_masked_num_valid_pix,full_vector_masked_num_possible_pix,100.0*(np.float(full_vector_masked_num_valid_pix)/full_vector_masked_num_possible_pix))

	elif full_stationary_bilinear_masked_offset_available and (full_stationary_masked_num_valid_pix>args.offset_correction_bilinear_fit_min_pixels) and full_stationary_masked_offset_available and (full_stationary_masked_num_valid_pix>full_stationary_min_num_pix) and ((float(full_stationary_masked_num_valid_pix)/full_stationary_masked_num_possible_pix)>=full_stationary_min_percent_valid_pix_available/100.0) and (np.sqrt(full_stationary_masked_offset_i**2.0 + full_stationary_masked_offset_j**2.0)<args.max_allowable_pixel_offset_correction):
		# use bilinear lgo_masked corection - but fields applied below, where i and j dc offsets are applied for the other corrections...
		final_offset_correction_i=full_stationary_masked_offset_i	# keep these anyway, to compare to bilinear solution offsets 
		final_offset_correction_j=full_stationary_masked_offset_j
		i_bilinear_offset_corr_arr=full_stationary_i_bilinear_offset_corr_arr
		j_bilinear_offset_corr_arr=full_stationary_j_bilinear_offset_corr_arr
		found_valid_offset=True
		found_valid_bilinear_offset=True
		offset_correction_type_applied='full_stationary_bilinear_masked_correction'
		offset_correction_type_descritption='full_stationary_bilinear_masked_correction fields for i and j for land pixels,' + ' - using bilinear coefs_i ' + str(coefs_i.ravel()) + ' coefs_j ' + str(coefs_j.ravel()) +' %d valid pixels out of %d possible for scene (%f %%)'%(full_stationary_masked_num_valid_pix,full_stationary_masked_num_possible_pix,100.0*(np.float(full_stationary_masked_num_valid_pix)/full_stationary_masked_num_possible_pix))

# this should be moved down the stack?
	elif full_vector_masked_offset_available and (full_vector_masked_num_valid_pix>full_vector_min_num_pix) and ((float(full_vector_masked_num_valid_pix)/full_vector_masked_num_possible_pix)>=full_vector_min_percent_valid_pix_available/100.0) and (np.sqrt(full_vector_masked_offset_i**2.0 + full_vector_masked_offset_j**2.0)<args.max_allowable_pixel_offset_correction):
		final_offset_correction_i=full_vector_masked_offset_i
		final_offset_correction_j=full_vector_masked_offset_j
		found_valid_offset=True
		offset_correction_type_applied='full_vector_masked_correction'
		offset_correction_type_descritption='full_vector_masked_correction for vel mosaic speeds < %f m/yr and lgo land pixels, %d valid pixels out of %d possible for scene (%f %%)'%(slow_area_zero_speed_myr,full_vector_masked_num_valid_pix,full_vector_masked_num_possible_pix,100.0*(np.float(full_vector_masked_num_valid_pix)/full_vector_masked_num_possible_pix))

# this should be moved down the stack?
	elif full_stationary_masked_offset_available and (full_stationary_masked_num_valid_pix>full_stationary_min_num_pix) and ((float(full_stationary_masked_num_valid_pix)/full_stationary_masked_num_possible_pix)>=full_stationary_min_percent_valid_pix_available/100.0) and (np.sqrt(full_stationary_masked_offset_i**2.0 + full_stationary_masked_offset_j**2.0)<args.max_allowable_pixel_offset_correction):
		final_offset_correction_i=full_stationary_masked_offset_i
		final_offset_correction_j=full_stationary_masked_offset_j
		found_valid_offset=True
		offset_correction_type_applied='full_stationary_masked_correction'
		offset_correction_type_descritption='full_stationary_masked_correction for vel mosaic speeds < %f m/yr and lgo land pixels, %d valid pixels out of %d possible for scene (%f %%)'%(slow_area_zero_speed_myr,full_stationary_masked_num_valid_pix,full_stationary_masked_num_possible_pix,100.0*(np.float(full_stationary_masked_num_valid_pix)/full_stationary_masked_num_possible_pix))

# from here on - not a case where BOTH zero_speed and lgo_offsets are to be combined - so everywhere but GREENLAND, and greenland where both are not valid at the same time... (like all pre 5p6 runs of sc_pycorr)	
	elif zspeed_offset_available and (zspeed_num_valid_pix>zero_speed_min_num_pix) and ((float(zspeed_num_valid_pix)/zspeed_num_possible_pix)>=zero_speed_min_percent_valid_pix_available/100.0) and (np.sqrt(zspeed_offset_i**2.0 + zspeed_offset_j**2.0)<args.max_allowable_pixel_offset_correction):
		#use zspeed correction
		final_offset_correction_i=zspeed_offset_i
		final_offset_correction_j=zspeed_offset_j
		found_valid_offset=True
		offset_correction_type_applied='zero_speed_correction'
		offset_correction_type_descritption='zero_speed_correction for vel mosaic speeds < %f m/yr, %d valid pixels out of %d possible for scene (%f %%)'%(slow_area_zero_speed_myr,zspeed_num_valid_pix,zspeed_num_possible_pix,100.0*(np.float(zspeed_num_valid_pix)/zspeed_num_possible_pix))
	elif midspeed_offset_available and (midspeed_num_valid_pix>mid_speed_min_num_pix) and ((float(midspeed_num_valid_pix)/midspeed_num_possible_pix)>=mid_speed_min_percent_valid_pix_available/100.0) and (np.sqrt(midspeed_offset_i**2.0 + midspeed_offset_j**2.0)<args.max_allowable_pixel_offset_correction):
		# use midspeed corection
		final_offset_correction_i=midspeed_offset_i
		final_offset_correction_j=midspeed_offset_j
		found_valid_offset=True
		offset_correction_type_applied='mid_speed_correction'
		offset_correction_type_descritption='mid_speed_correction for vel mosaic speeds > %f m/yr and < %f m/yr, %d valid pixels out of %d possible for scene (%f %%)'%(slow_area_zero_speed_myr,slow_area_max_vector_speed_myr,midspeed_num_valid_pix,midspeed_num_possible_pix,100.0*(np.float(midspeed_num_valid_pix)/midspeed_num_possible_pix))
	elif lgo_bilinear_masked_offset_available and (lgo_masked_num_valid_pix>=args.offset_correction_bilinear_fit_min_pixels) and lgo_masked_offset_available  and (lgo_masked_num_valid_pix>lgo_masked_min_num_pix) and ((float(lgo_masked_num_valid_pix)/lgo_masked_num_possible_pix)>=lgo_masked_min_percent_valid_pix_available/100.0) and (np.sqrt(lgo_masked_offset_i**2.0 + lgo_masked_offset_j**2.0)<args.max_allowable_pixel_offset_correction):
		# use bilinear lgo_masked corection - but fields applied below, where i and j dc offsets are applied for the other corrections...
		final_offset_correction_i=lgo_masked_offset_i	# keep these anyway, to compare to bilinear solution offsets 
		final_offset_correction_j=lgo_masked_offset_j
		i_bilinear_offset_corr_arr=lgo_i_bilinear_offset_corr_arr
		j_bilinear_offset_corr_arr=lgo_j_bilinear_offset_corr_arr
		found_valid_offset=True
		found_valid_bilinear_offset=True
		offset_correction_type_applied='lgo_bilinear_masked_correction'
		offset_correction_type_descritption='lgo_bilinear_masked_correction fields for i and j for land pixels,' + ' - using bilinear coefs_i ' + str(coefs_i.ravel()) + ' coefs_j ' + str(coefs_j.ravel()) +' %d valid pixels out of %d possible for scene (%f %%)'%(lgo_masked_num_valid_pix,lgo_masked_num_possible_pix,100.0*(np.float(lgo_masked_num_valid_pix)/lgo_masked_num_possible_pix))
	elif lgo_masked_offset_available  and (lgo_masked_num_valid_pix>lgo_masked_min_num_pix) and ((float(lgo_masked_num_valid_pix)/lgo_masked_num_possible_pix)>=lgo_masked_min_percent_valid_pix_available/100.0) and (np.sqrt(lgo_masked_offset_i**2.0 + lgo_masked_offset_j**2.0)<args.max_allowable_pixel_offset_correction):
		# use lgo_masked corection
		final_offset_correction_i=lgo_masked_offset_i
		final_offset_correction_j=lgo_masked_offset_j
		found_valid_offset=True
		offset_correction_type_applied='lgo_masked_correction'
		offset_correction_type_descritption='lgo_masked_correction for land pixels, %d valid pixels out of %d possible for scene (%f %%)'%(lgo_masked_num_valid_pix,lgo_masked_num_possible_pix,100.0*(np.float(lgo_masked_num_valid_pix)/lgo_masked_num_possible_pix))
	else:
		offset_correction_type_applied='None'
		offset_correction_type_descritption='None'
		
	if found_valid_offset and not(found_valid_bilinear_offset):
		t_log('using %s offset correction %f pixels in i, %f pixels in j '%(offset_correction_type_applied,final_offset_correction_i,final_offset_correction_j) + offset_correction_type_descritption,outlogdisabled=args.nlf)
	elif found_valid_offset and found_valid_bilinear_offset:
		t_log('using %s bilinear offset correction - constant offsets would have been %f pixels in i, %f pixels in j '%(offset_correction_type_applied,final_offset_correction_i,final_offset_correction_j) + offset_correction_type_descritption,outlogdisabled=args.nlf)
	else:
		t_log('offset correction not used - not enought valid pixels',outlogdisabled=args.nlf)

############################################################################		

dcam=args.dcam
cam=args.cam
cam1=args.cam1

# if args.offset_correction:
# 	offset_dist_ij_arr=np.power((np.power(del_i_corrected,2.0) + np.power(del_j_corrected,2.0)),0.5)
# 	vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1)|(((offset_dist_ij_arr*img1.pix_x_m)/del_t_val)>plotvmax),((offset_dist_ij_arr*img1.pix_x_m)/del_t_val))
# 	vx=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),((del_i_corrected*img1.pix_x_m)/del_t_val))
# 	vy=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),((del_j_corrected*img1.pix_y_m)/del_t_val))
# else:
# 	vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1)|(((offset_dist_ij_arr*img1.pix_x_m)/del_t_val)>plotvmax),((offset_dist_ij_arr*img1.pix_x_m)/del_t_val))
# 	vx=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),((del_i*img1.pix_x_m)/del_t_val))
# 	vy=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),((del_j*img1.pix_y_m)/del_t_val))

if not(args.offset_correction_speedref or args.offset_correction_lgo_mask) or not(found_valid_offset):  # no offset to apply
	# create masked and unmasked versions of the speed, vx, and vy output
	# the masked version will be written to the nc file with no_data values in corners and also where the correlation parameters suggest masking is needed
# 	vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1)|(((offset_dist_ij_arr*img1.pix_x_m)/del_t_val)>plotvmax),((offset_dist_ij_arr*img1.pix_x_m)/del_t_val))
	vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),((offset_dist_ij_arr*img1.pix_x_m)/del_t_val))
	vv[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run
	vx=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),((del_i*img1.pix_x_m)/del_t_val))
	vx[corr_arr==corr_nodata_val] = vel_nodata_val
	vy=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),((del_j*img1.pix_y_m)/del_t_val))
	vy[corr_arr==corr_nodata_val] = vel_nodata_val

	# nomask arrays only have nodata values where there was no velocity attempted (filled corners of landsat scenes, for example)
	vv_nomask=((offset_dist_ij_arr*img1.pix_x_m)/del_t_val)
	vv_nomask[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run
	vx_nomask=((del_i*img1.pix_x_m)/del_t_val)
	vx_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
	vy_nomask=((del_j*img1.pix_y_m)/del_t_val)
	vy_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
	
elif not(found_valid_bilinear_offset):  # have single i and j offsets to apply
	# apply offset and create masked and unmasked versions of the speed, vx, and vy output
	# the masked version will be written to the nc file with no_data values in corners and also where the correlation parameters suggest masking is needed
	vx=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),(((del_i + final_offset_correction_i) * img1.pix_x_m)/del_t_val))
	vx[corr_arr==corr_nodata_val] = vel_nodata_val
	vy=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),(((del_j + final_offset_correction_j) * img1.pix_y_m)/del_t_val))
	vy[corr_arr==corr_nodata_val] = vel_nodata_val
# 	vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1)|(np.sqrt(np.square(vx) + np.square(vy))>plotvmax),np.sqrt(np.square(vx) + np.square(vy)))
	vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),np.sqrt(np.square(vx) + np.square(vy)))
	vv[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run

	# nomask arrays only have nodata values where there was no velocity attempted (filled corners of landsat scenes, for example)
	vx_nomask=(((del_i + final_offset_correction_i) * img1.pix_x_m)/del_t_val)
	vx_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
	vy_nomask=(((del_j + final_offset_correction_j) * img1.pix_y_m)/del_t_val)
	vy_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
	vv_nomask=np.sqrt(np.square(vx_nomask) + np.square(vy_nomask))
	vv_nomask[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run
	
elif found_valid_bilinear_offset:  # have bilinear i and j offset surfaces to apply
	# apply offset and create masked and unmasked versions of the speed, vx, and vy output
	# the masked version will be written to the nc file with no_data values in corners and also where the correlation parameters suggest masking is needed
	vx=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),(((del_i + i_bilinear_offset_corr_arr) * img1.pix_x_m)/del_t_val))
	vx[corr_arr==corr_nodata_val] = vel_nodata_val
	vy=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),(((del_j + j_bilinear_offset_corr_arr) * img1.pix_y_m)/del_t_val))
	vy[corr_arr==corr_nodata_val] = vel_nodata_val
# 	vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1)|(np.sqrt(np.square(vx) + np.square(vy))>plotvmax),np.sqrt(np.square(vx) + np.square(vy)))
	vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1),np.sqrt(np.square(vx) + np.square(vy)))
	vv[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run

	# nomask arrays only have nodata values where there was no velocity attempted (filled corners of landsat scenes, for example)
	vx_nomask=(((del_i + i_bilinear_offset_corr_arr) * img1.pix_x_m)/del_t_val)
	vx_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
	vy_nomask=(((del_j + j_bilinear_offset_corr_arr) * img1.pix_y_m)/del_t_val)
	vy_nomask[corr_arr==corr_nodata_val] = vel_nodata_val
	vv_nomask=np.sqrt(np.square(vx_nomask) + np.square(vy_nomask))
	vv_nomask[corr_arr==corr_nodata_val] = vel_nodata_val  # note this only eliminates pixels that had black (no data0 values from corners of input image 1 - so correlation wasn't run
	
else:
	print('should not be able to get here, but wanted to put in that last test above...what gives?')

# mapjet=plt.get_cmap('jet')
mapjet=get_cmap('jet')
mmvv=vv.copy()
mmvv[(mmvv>plotvmax) & (~mmvv.mask)]=plotvmax
vv_zo=mmvv/plotvmax
vv_rgba=mapjet(vv_zo) # this produces a 4-band 0.0-1.0 image array using the jet colormap
(out_lines,out_pixels,out_bands)=vv_rgba.shape		# NEED these values below, EVEN if only output tif is log10

        
if (args.log10):	# prepare a log10 version of the speed for output below
	mmvv=vv.copy()
	mmvv[(mmvv>plotvmax) & (~mmvv.mask)]=plotvmax
	mmvv.mask[(mmvv==0) & (~mmvv.mask)]= True
	lmmvv=np.log10(mmvv)
# 	min_lmmvv=np.min(lmmvv)
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

		
# if not(args.no_gtif):   ###### This flag is not presently used, as GTiff forms basis for all output right now...
format = "GTiff"
driver = gdal.GetDriverByName( format )
metadata = driver.GetMetadata()
# dst_filename = outdir + '/' + file_name_base + '.tif'
# dst_ds = driver.Create( dst_filename, out_pixels, out_lines, out_bands, gdal.GDT_Byte )
# if (args.v):
# 	if metadata.has_key(gdal.DCAP_CREATE) and metadata[gdal.DCAP_CREATE] == 'YES':
# 		print 'Driver %s supports Create() method.' % format
# 	if metadata.has_key(gdal.DCAP_CREATECOPY) and metadata[gdal.DCAP_CREATECOPY] == 'YES':
# 		print 'Driver %s supports CreateCopy() method.' % format
# 	img1.proj
# 	dst_ds
# 	print [com_min_x,com_max_x,com_min_y,com_max_y]
# print 'out image %s %s %d, %d, %d'%(dst_filename, format, out_pixels, out_lines, out_bands)
# if not(args.nlf):
# 	log_output_lines.append('out image %s %s %d, %d, %d\n'%(dst_filename, format, out_pixels, out_lines, out_bands))
# # dst_ds.SetGeoTransform( [ com_min_x, inc * img1.pix_x_m, 0, com_max_y, 0, inc * img1.pix_y_m ] ) # note pix_y_m typically negative
# dst_ds.SetGeoTransform( [ output_array_ul_corner[0], output_array_pix_x_m, 0, output_array_ul_corner[1], 0, output_array_pix_y_m ] ) # note pix_y_m typically negative
# dst_ds.SetProjection( img1.proj )
# dst_ds.GetRasterBand(1).WriteArray( (vv_rgba[:,:,0]*255).astype('ubyte') )
# dst_ds.GetRasterBand(2).WriteArray( (vv_rgba[:,:,1]*255).astype('ubyte') )
# dst_ds.GetRasterBand(3).WriteArray( (vv_rgba[:,:,2]*255).astype('ubyte') )
# dst_ds.GetRasterBand(4).WriteArray( (vv_rgba[:,:,3]*255).astype('ubyte') )
# dst_ds = None # done, close the dataset

if(args.log10):
	dst_filename = outdir + '/' + file_name_base + '_log10.tif'
	dst_ds = driver.Create( dst_filename, out_pixels, out_lines, out_bands, gdal.GDT_Byte )
	print('out image %s %s %d, %d, %d'%(dst_filename, format, out_pixels, out_lines, out_bands))
	if not(args.nlf):
			log_output_lines.append('out image %s %s %d, %d, %d\n'%(dst_filename, format, out_pixels, out_lines, out_bands))
	# dst_ds.SetGeoTransform( [ com_min_x, inc * img1.pix_x_m, 0, com_max_y, 0, inc * img1.pix_y_m ] ) # note pix_y_m typically negative
	dst_ds.SetGeoTransform( [ output_array_ul_corner[0], output_array_pix_x_m, 0, output_array_ul_corner[1], 0, output_array_pix_y_m ] ) # note pix_y_m typically negative
	dst_ds.SetProjection( img1.proj )
	dst_ds.GetRasterBand(1).WriteArray( (lvv_rgba[:,:,0]*255).astype('ubyte') )
	dst_ds.GetRasterBand(2).WriteArray( (lvv_rgba[:,:,1]*255).astype('ubyte') )
	dst_ds.GetRasterBand(3).WriteArray( (lvv_rgba[:,:,2]*255).astype('ubyte') )
	dst_ds.GetRasterBand(4).WriteArray( (lvv_rgba[:,:,3]*255).astype('ubyte') )
	dst_ds = None # done, close the dataset

##################################
#
# now for the netCDF data stack...
#
##################################

quiet = False
clobber = True     # overwrite existing output nc file

out_nc_filename=outdir + '/' + file_name_base + '_hp.nc'

nc_outfile = netCDF4.Dataset(out_nc_filename,'w',clobber=clobber,format='NETCDF4')


# First set global attributes that GDAL uses when it reads netCFDF files
nc_outfile.setncattr('GDAL_AREA_OR_POINT','Area')
nc_outfile.setncattr('Conventions','CF-1.6')
nc_outfile.setncattr('history','%s : %s  ====> sc_pycorr run'%(dt.datetime.now().isoformat(),' '.join(sys.argv)))
nc_outfile.setncattr('GDAL',gdal.VersionInfo("VERSION_DATE"))

# if we were copying from an nc file, we would...
# for attr in vv_nc_basefile.ncattrs():
#     nc_outfile.setncattr(attr,vv_nc_basefile.getncattr(attr))

# set dimensions
nc_outfile.createDimension('x',out_pixels)
nc_outfile.createDimension('y',out_lines)


# set variables
# first set up image_pair_times variable not as a dimension, but as holder for attributes for the times of the two input images, delt, etc.
varname='image_pair_times'
delt=str(del_t_val) # need to set this as an attribute
delt_units=del_t_unit_str
delt_speed_units=del_t_speedunit_str
datatype=np.dtype('S1')
dimensions=()
FillValue=None
lsd = None

dtstart=img1.imagedatetime
start_decimal_year=dtstart.year + (dtstart - dt.datetime(year=dtstart.year,month=1,day=1)).total_seconds()/(dt.datetime(year=dtstart.year+1,month=1,day=1)-dt.datetime(year=dtstart.year,month=1,day=1)).total_seconds()
dtmid=img1.imagedatetime + dt.timedelta(seconds=((img2.imagedatetime-img1.imagedatetime).total_seconds()/2.0))
mid_decimal_year=dtmid.year + (dtmid - dt.datetime(year=dtmid.year,month=1,day=1)).total_seconds()/(dt.datetime(year=dtmid.year+1,month=1,day=1)-dt.datetime(year=dtmid.year,month=1,day=1)).total_seconds()
dtstop=img2.imagedatetime
stop_decimal_year=dtstop.year + (dtstop - dt.datetime(year=dtstop.year,month=1,day=1)).total_seconds()/(dt.datetime(year=dtstop.year+1,month=1,day=1)-dt.datetime(year=dtstop.year,month=1,day=1)).total_seconds()

# 	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, least_significant_digit=lsd,zlib=zlib,complevel=complevel,shuffle=shuffle,fletcher32=fletcher32)
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
var.setncattr('del_t',delt)
var.setncattr('del_t_units',delt_units)
var.setncattr('del_t_speed_units',delt_speed_units)
var.setncattr('start_time_decimal_year',str(start_decimal_year))
var.setncattr('mid_time_decimal_year',str(mid_decimal_year))
var.setncattr('end_time_decimal_year',str(stop_decimal_year))
var.setncattr('start_date',dtstart.isoformat())
var.setncattr('mid_date',dtmid.isoformat())
var.setncattr('end_date',dtstop.isoformat())


# next set up input_image_details variable not as a dimension, but as holder for attributes for the times of the two input images, delt, etc.
varname='input_image_details'

datatype=np.dtype('S1')
dimensions=()
FillValue=None
lsd = None

# 	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, least_significant_digit=lsd,zlib=zlib,complevel=complevel,shuffle=shuffle,fletcher32=fletcher32)
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
var.setncattr('image1_filename',img1.filename)
var.setncattr('image1_path',img1.in_dir_abs_path)
var.setncattr('image1_pix_x_size',img1.pix_x_m)
var.setncattr('image1_pix_y_size',img1.pix_y_m)
var.setncattr('image1_proj_WKT',img1.proj)
var.setncattr('image1_proj_GeoTransform',img1.gt)
var.setncattr('image1_decimal_year',str(start_decimal_year))

var.setncattr('image2_filename',img2.filename)
var.setncattr('image2_path',img2.in_dir_abs_path)
var.setncattr('image2_pix_x_size',img2.pix_x_m)
var.setncattr('image2_pix_y_size',img2.pix_y_m)
var.setncattr('image2_proj_WKT',img2.proj)
var.setncattr('image2_proj_GeoTransform',img2.gt)
var.setncattr('image2_decimal_year',str(stop_decimal_year))


# projection variable is tricky - for all L8 not in Antarctica, it is transverse mercator; for Antarctica, polar stereo
if (img1.srs.GetAttrValue('PROJECTION') == 'Polar_Stereographic'):
	# Antarctic input file - OR this could be a non-landsat PS image from the north and it might work? because of the sign of the pole used below??? maybe...
	print('antarctic image')

	varname='polar_stereographic'
	grid_mapping='polar_stereographic'  # need to set this as an attribute for the image variables
	datatype=np.dtype('S1')
	dimensions=()
	FillValue=None
	lsd = None
# 	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, least_significant_digit=lsd,zlib=zlib,complevel=complevel,shuffle=shuffle,fletcher32=fletcher32)
	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
	# variable made, now add attributes
	out_pix_size_x=output_array_pix_x_m
	out_pix_size_y=output_array_pix_y_m
	out_geotransform = [ output_array_ul_corner[0], output_array_pix_x_m, 0, output_array_ul_corner[1], 0, output_array_pix_y_m ]
	
	var.setncattr('grid_mapping_name','polar_stereographic')
	var.setncattr('straight_vertical_longitude_from_pole',img1.srs.GetProjParm('central_meridian'))
	var.setncattr('false_easting',img1.srs.GetProjParm('false_easting'))
	var.setncattr('false_northing',img1.srs.GetProjParm('false_northing'))
	var.setncattr('latitude_of_projection_origin',np.sign(img1.srs.GetProjParm('latitude_of_origin'))*90.0)  # could hardcode this to be -90 for landsat - just making it more general, maybe...
	var.setncattr('standard_parallel',img1.srs.GetProjParm('latitude_of_origin'))
	var.setncattr('longitude_of_prime_meridian',float(img1.srs.GetAttrValue('GEOGCS|PRIMEM',1)))
	var.setncattr('semi_major_axis',float(img1.srs.GetAttrValue('GEOGCS|SPHEROID',1)))
	var.setncattr('inverse_flattening',float(img1.srs.GetAttrValue('GEOGCS|SPHEROID',2)))
	var.setncattr('spatial_ref',img1.srs.ExportToWkt())
	var.setncattr('GeoTransform',' '.join(str(x) for x in out_geotransform))  # note this has pixel size in it - set  explicitly above

elif (img1.srs.GetAttrValue('PROJECTION') == 'Transverse_Mercator'):
	# rest of the world for landsat - UTM
	# make 'transverse_mercator' variable with special gadl netCDF projection attributes
	varname='transverse_mercator'
	grid_mapping='transverse_mercator'  # need to set this as an attribute for the image variables
	datatype=np.dtype('S1')
	dimensions=()
	FillValue=None
	lsd = None
# 	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, least_significant_digit=lsd,zlib=zlib,complevel=complevel,shuffle=shuffle,fletcher32=fletcher32)
	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
	# variable made, now add attributes
	out_pix_size_x=output_array_pix_x_m
	out_pix_size_y=output_array_pix_y_m
	out_geotransform = [ output_array_ul_corner[0], output_array_pix_x_m, 0, output_array_ul_corner[1], 0, output_array_pix_y_m ]
	
	var.setncattr('grid_mapping_name','transverse_mercator')
	var.setncattr('longitude_of_central_meridian',img1.srs.GetProjParm('central_meridian'))
	var.setncattr('false_easting',img1.srs.GetProjParm('false_easting'))
	var.setncattr('false_northing',img1.srs.GetProjParm('false_northing'))
	var.setncattr('latitude_of_projection_origin',img1.srs.GetProjParm('latitude_of_origin'))
	var.setncattr('scale_factor_at_central_meridian',img1.srs.GetProjParm('scale_factor'))
	var.setncattr('longitude_of_prime_meridian',float(img1.srs.GetAttrValue('GEOGCS|PRIMEM',1)))
	var.setncattr('semi_major_axis',float(img1.srs.GetAttrValue('GEOGCS|SPHEROID',1)))
	var.setncattr('inverse_flattening',float(img1.srs.GetAttrValue('GEOGCS|SPHEROID',2)))
	var.setncattr('spatial_ref',img1.srs.ExportToWkt())
	var.setncattr('GeoTransform',' '.join(str(x) for x in out_geotransform))  # note this has pixel size in it - set  explicitly above
else:
	print('projection %s not recognized for this program'%(img1.srs.GetAttrValue('PROJECTION')))
	exit(-1)

x_vals=np.arange((out_geotransform[0] + (out_pix_size_x/2.0)),(out_geotransform[0] + (out_pix_size_x * out_pixels)),out_pix_size_x)
# y_vals=np.arange((out_geotransform[3] + (out_pix_size_y/2.0)),(out_geotransform[3] + (out_pix_size_y * out_lines)),out_pix_size_y)[::-1]  # flip the order to increasing
y_vals=np.arange((out_geotransform[3] + (out_pix_size_y/2.0)),(out_geotransform[3] + (out_pix_size_y * out_lines)),out_pix_size_y)  # not any more - top down so gdal will be happy


varname='x'
datatype=np.dtype('float64')
dimensions=('x')
FillValue=False
lsd = None
# var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, least_significant_digit=lsd,zlib=zlib,complevel=complevel,shuffle=shuffle,fletcher32=fletcher32)
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
var.setncattr('standard_name','projection_x_coordinate')
var.setncattr('long_name','x coordinate of projection')
var.setncattr('units','m')
var[:] = x_vals


varname='y'
datatype=np.dtype('float64')
dimensions=('y')
FillValue=False
lsd = None
# var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, least_significant_digit=lsd,zlib=zlib,complevel=complevel,shuffle=shuffle,fletcher32=fletcher32)
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
var.setncattr('standard_name','projection_y_coordinate')
var.setncattr('long_name','y coordinate of projection')
var.setncattr('units','m')
var[:] = y_vals

# set variables
# first set up image_pair_times variable not as dimension, but as holder for attributes for the times of the two input images, delt, etc.
varname='offset_correction'
# delt=str(del_t_val) # need to set this as an attribute
# delt_units=del_t_unit_str
# delt_speed_units=del_t_speedunit_str
datatype=np.dtype('S1')
dimensions=()
FillValue=None
var = nc_outfile.createVariable(varname,datatype, dimensions, fill_value=FillValue)
if found_valid_offset and not(found_valid_bilinear_offset):
	var.setncattr('offset_correction_type_applied',offset_correction_type_applied)
	var.setncattr('offset_correction_type_descritption',offset_correction_type_descritption)
	var.setncattr('final_offset_correction_units','pixels')
	var.setncattr('final_offset_correction_i',str(final_offset_correction_i))
	var.setncattr('final_offset_correction_j',str(final_offset_correction_j))
elif found_valid_offset and found_valid_bilinear_offset:
	var.setncattr('offset_correction_type_applied',offset_correction_type_applied)
	var.setncattr('offset_correction_type_descritption',offset_correction_type_descritption)
	var.setncattr('final_offset_correction_units','pixels')
	var.setncattr('final_offset_correction_i','coefs_i ' + str(coefs_i.ravel()))
	var.setncattr('final_offset_correction_j','coefs_j ' + str(coefs_j.ravel()))
else:
	var.setncattr('offset_correction_type_applied',offset_correction_type_applied)
	var.setncattr('offset_correction_type_descritption',offset_correction_type_descritption)
	var.setncattr('final_offset_correction_units','pixels')
	var.setncattr('final_offset_correction_i','None')
	var.setncattr('final_offset_correction_j','None')

if found_valid_offset and found_valid_bilinear_offset:  # write out bilinear correction arrays so that one can get back to original del_i,del_j if needed
	varname='applied_bilinear_x_offset_correction_in_pixels'
	print('writing %s'%(varname))
	datatype=np.dtype('float32')
	dimensions=('y','x')
	FillValue=-99.0
	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
	var.setncattr('bilinear_offset_correction_details','%s applied to measured x offset (offset in pixels that was added to original match data - subtract x_pix_size*del_t*this_offset and divide by x_pix_size to get to original x offset in pixels)'%(offset_correction_type_applied))
	var.setncattr('grid_mapping',grid_mapping)
	var.setncattr('standard_name','x_offset_pixels')
	var.setncattr('long_name','x offset correction in pixels')
	var.setncattr('units','pixels')
	# var[:] = np.flipud(vx_nomask).astype('float32')
	var[:] = i_bilinear_offset_corr_arr.astype('float32')
	if not(args.nlf):
		log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))
	
	varname='applied_bilinear_y_offset_correction_in_pixels'
	print('writing %s'%(varname))
	datatype=np.dtype('float32')
	dimensions=('y','x')
	FillValue=-99.0
	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
	var.setncattr('bilinear_offset_correction_details','%s applied to measured y offset (offset in pixels that was added to original match data - subtract y_pix_size*del_t*this_offset and divide by y_pix_size to get to original y offset in pixels)'%(offset_correction_type_applied))
	var.setncattr('grid_mapping',grid_mapping)
	var.setncattr('standard_name','y_offset_pixels')
	var.setncattr('long_name','y offset correction in pixels')
	var.setncattr('units','pixels')
	# var[:] = np.flipud(vx_nomask).astype('float32')
	var[:] = j_bilinear_offset_corr_arr.astype('float32')
	if not(args.nlf):
		log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))


varname='vx'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','x_velocity')
var.setncattr('long_name','x component of velocity')
var.setncattr('units',del_t_speedunit_str)
# var[:] = np.flipud(vx_nomask).astype('float32')
var[:] = vx_nomask.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

varname='vy'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','y_velocity')
var.setncattr('long_name','y component of velocity')
var.setncattr('units',del_t_speedunit_str)
# var[:] = np.flipud(vy_nomask).astype('float32')
var[:] = vy_nomask.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

varname='vv'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','speed')
var.setncattr('long_name','magnitude of velocity')
var.setncattr('units',del_t_speedunit_str)
# var[:] = np.flipud(vv_nomask).astype('float32')
var[:] = vv_nomask.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

varname='vx_masked'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','x_velocity_masked')
var.setncattr('long_name','x component of velocity (masked)')
var.setncattr('units',del_t_speedunit_str)
var.setncattr('masking_info','masked_where(((del_corr_arr<%4.3f)&(corr_arr<%4.3f))|(corr_arr<%4.3f))'%(dcam,cam,cam1))
# var[:] = np.flipud(vx).astype('float32')
var[:] = vx.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

varname='vy_masked'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','y_velocity_masked')
var.setncattr('long_name','y component of velocity (masked)')
var.setncattr('units',del_t_speedunit_str)
var.setncattr('masking_info','masked_where(((del_corr_arr<%4.3f)&(corr_arr<%4.3f))|(corr_arr<%4.3f))'%(dcam,cam,cam1))
# var[:] = np.flipud(vy).astype('float32')
var[:] = vy.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

varname='vv_masked'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=vel_nodata_val
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','speed_masked')
var.setncattr('long_name','magnitude of velocity (masked)')
var.setncattr('units',del_t_speedunit_str)
var.setncattr('masking_info','masked_where(((del_corr_arr<%4.3f)&(corr_arr<%4.3f))|(corr_arr<%4.3f))'%(dcam,cam,cam1))
# var[:] = np.flipud(vv).astype('float32')
var[:] = vv.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))
	
varname='corr'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=corr_nodata_val
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','corr')
var.setncattr('long_name','peak correlation value')
var.setncattr('units',' ')
# var[:] = np.flipud(corr_arr).astype('float32')
var[:] = corr_arr.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

varname='del_corr'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=corr_nodata_val
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','del_corr')
var.setncattr('long_name','delta corr between primary and second peak')
var.setncattr('units',' ')
# var[:] = np.flipud(del_corr_arr).astype('float32')
var[:] = del_corr_arr.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

# lgo_mask_image_utm - added in v5p2
if args.offset_correction_lgo_mask:
	varname='lgo_mask'
	print('writing %s'%(varname))
	datatype=np.dtype('byte')
	dimensions=('y','x')
	if lgo_mask_nodata:
		FillValue=np.byte(lgo_mask_nodata)
	else:
		FillValue=np.byte(255)
	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
	var.setncattr('grid_mapping',grid_mapping)
	var.setncattr('standard_name','lgo_mask')
	var.setncattr('long_name','land(1) glacier(0) ocean(2) mask')
	var.setncattr('units',' ')
# 	var[:] = np.flipud(lgo_mask_image_utm).astype('byte')
	var[:] = lgo_mask_image_utm.astype('byte')
	if not(args.nlf):
		log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

varname='d2idx2'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=curvature_nodata_value
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','d2idx2')
var.setncattr('long_name','corr peak curvature in x direction')
var.setncattr('units',' ')
# var[:] = np.flipud(d2idx2).astype('float32')
var[:] = d2idx2.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

varname='d2jdx2'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=curvature_nodata_value
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','d2jdx2')
var.setncattr('long_name','corr peak curvature in y direction')
var.setncattr('units',' ')
# var[:] = np.flipud(d2jdx2).astype('float32')
var[:] = d2jdx2.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

#######################################################
# commented this out in version 5p8 because it can be calculated directly from d2idx2 and d2jdx2
#######################################################
#
# varname='d2dx2_mean'
# print 'writing %s'%(varname)
# datatype=np.dtype('float32')
# dimensions=('y','x')
# FillValue=curvature_nodata_value
# var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
# var.setncattr('grid_mapping',grid_mapping)
# var.setncattr('standard_name','d2dx2_mean')
# var.setncattr('long_name','mean corr peak curvature')
# var.setncattr('units',' ')
# # var[:] = np.flipud(d2dx2_mean).astype('float32')
# var[:] = d2dx2_mean.astype('float32')
# if not(args.nlf):
# 	log_output_lines.append('added  %s to netCDF file %d, %d \n'%(varname, out_pixels, out_lines))

varname='del_i'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=None
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','i_pixel_offset')
var.setncattr('long_name','i pixel offset (positive in image right direction, original image pixel size, no offset correction applied)')
var.setncattr('units','pixels')
# var[:] = np.flipud(vx_nomask).astype('float32')
var[:] = del_i.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))

varname='del_j'
print('writing %s'%(varname))
datatype=np.dtype('float32')
dimensions=('y','x')
FillValue=None
var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, zlib=True, complevel=2, shuffle=True)
var.setncattr('grid_mapping',grid_mapping)
var.setncattr('standard_name','j_pixel_offset')
var.setncattr('long_name','j pixel offset (positive in image down direction (if pix_y_m is negative), original image pixel size, no offset correction applied)')
var.setncattr('units','pixels')
# var[:] = np.flipud(vx_nomask).astype('float32')
var[:] = del_j.astype('float32')
if not(args.nlf):
	log_output_lines.append('added  %s to netCDF file %d, %d '%(varname, out_pixels, out_lines))


# enter end information into log so it will be in nc file - repeated to screen again at the real end of the program
if psutil_available and not(args.nlf):
	t_log('At end - psutil reports process %s using '%(args.out_name_base) + str(memory_usage_psutil()),outlogdisabled=args.nlf)
if resource_available and not(args.nlf):
	t_log('At end - resource module reports process %s using '%(args.out_name_base) + str(memory_usage_resource()),outlogdisabled=args.nlf)
	
if not(args.nlf):   # add log file output to nc file
	inlog=''.join(log_output_lines)
	nc_outfile.createDimension('chars',len(inlog))

	varname='chars'
	datatype=np.dtype('int32')
	dimensions=('chars')
	FillValue=False
	lsd = None
	# var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue, least_significant_digit=lsd,zlib=zlib,complevel=complevel,shuffle=shuffle,fletcher32=fletcher32)
	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
	var.setncattr('standard_name','character_index')
	var.setncattr('long_name','character index for log file')
	var.setncattr('units',' ')
	var[:] = range(len(inlog))


	varname='processing_log'
	datatype=np.dtype('S1')
	dimensions=('chars')
	FillValue=None
	lsd = None
	var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
	var.setncattr('standard_name','sc_pycorr_process_log')
	var.setncattr('long_name','log from feature tracking processing')
	for indx,character in enumerate(inlog):
		var[indx]=character


nc_outfile.sync() # flush data to disk
nc_outfile.close()



if psutil_available:
	print('At end - psutil reports process %s using '%(args.out_name_base),memory_usage_psutil())
if resource_available:
	print('At end - resource module reports process %s using '%(args.out_name_base),memory_usage_resource())

t_log("End Program",outlogdisabled=args.nlf)
