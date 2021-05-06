import numpy as np
import gdal
import os
# import subprocess as sp
import string
import sys
import time
import datetime as dt
import argparse
import osr
from scipy.ndimage.filters import gaussian_filter
# import multiprocessing
import glob

try:
	import resource
	resource_available=True
	print 'resource found'
	def memory_usage_resource():
		# return the memory usage in MB
		usage = resource.getrusage(resource.RUSAGE_SELF)
		return('%7.3f local %7.3f stack %7.3f max memory'%(float(usage.ru_idrss)/float(2 ** 20),float(usage.ru_isrss)/float(2 ** 20),float(usage.ru_maxrss)/float(2 ** 20)))
except:
	resource_available=False
	print 'resource package not found, proceeding without resource memory usage reports'


class GeoImg_noload:  # modified 9/7/2015 for LO8 fname -> date    modified to noload for sc application - don't read image on setup - will read image data in main code to hp filter, delete...and keep footprint small
	"""geocoded image input and info
		a=GeoImg(in_file_name,indir='.')
			a.img will contain image
			a.parameter etc..."""
	def __init__(self, in_filename,in_dir='.',datestr=None,datefmt='%m/%d/%y'):
		self.filename = in_filename
		self.in_dir_path = in_dir  #in_dir can be relative...
		self.in_dir_abs_path=os.path.abspath(in_dir)  # get absolute path for later ref if needed
		self.gd=gdal.Open(self.in_dir_path + os.path.sep + self.filename)
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
		if (datestr is not None):
			self.imagedatetime=dt.datetime.strptime(datestr,datefmt)
		elif ((self.filename.find('LC8') == 0) | (self.filename.find('LO8') == 0) | \
				(self.filename.find('LE7') == 0) | (self.filename.find('LT5') == 0) | \
				(self.filename.find('LT4') == 0)):	# looks landsat like - try parsing the date from filename (contains day of year)
			self.sensor=self.filename[0:3]
			self.path=int(self.filename[3:6])
			self.row=int(self.filename[6:9])
			self.year=int(self.filename[9:13])
			self.doy=int(self.filename[13:16])
			self.imagedatetime=dt.date.fromordinal(dt.date(self.year-1,12,31).toordinal()+self.doy)
		else:
			self.imagedatetime=None  # need to throw error in this case...or get it from metadata
# 		self.img=self.gd.ReadAsArray().astype(np.float32)   # works for L8 and earlier - and openCV correlation routine needs float or byte so just use float...
		self.srs=osr.SpatialReference(wkt=self.gd.GetProjection())
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
# 		self.img=self.gd.ReadAsArray().astype(np.uint8)		# L7 and earlier - doesn't work with plt.imshow...
# 		self.img_ov2=self.img[0::2,0::2]
# 		self.img_ov10=self.img[0::10,0::10]


# these need to be global - based on numpy limits for int16 datatype
int_limit_maxm1 = np.iinfo('int16').max - 1  # get max (and min) signed int16 values for rescaling float hp before converting back to signed int16
int_limit_minp2 = np.iinfo('int16').min + 2
int_nodata_val = np.iinfo('int16').min
format = "GTiff"
driver = gdal.GetDriverByName( format )


def make_hp_func(instuff):
	process_time=time.time()
	indir,b8_img_name,hpdir,hp1filename,hpargs=instuff
	if resource_available:
		print 'memory use at start of hp function: %s'%(memory_usage_resource())
# 	print 'working on image: %s'%(b8_img_name)
	img1=GeoImg_noload(b8_img_name,in_dir=indir)
	if resource_available:
		print 'memory use after noload: %s'%(memory_usage_resource())
	# 	t_log("open image")
	img1_img=img1.gd.ReadAsArray().astype(np.float32)
	print 'pixels,lines ',img1_img.shape
	if resource_available:
		print 'memory use after load as float32: %s'%(memory_usage_resource())
	hp_arr_1=np.zeros_like(img1_img,dtype=np.float32)
	if resource_available:
		print 'memory use after alloc for hp_arr_1 as float32: %s'%(memory_usage_resource())
	gaussian_filter(img1_img,hpargs.gfilt_sigma,output=hp_arr_1)
	if resource_available:
		print 'memory use after gaussian_filter: %s'%(memory_usage_resource())
	hp_arr_1 -= img1_img
	hp_arr_1 *= -1.0
# 	for j in range(img1_img.shape[0]):
# 		for i in range(img1_img.shape[1]):
# 			tmp=hp_arr_1[j,i]
# 			hp_arr_1[j,i]=img1_img[j,i]-tmp
# 	hp_arr_1=img1_img-gaussian_filter(img1_img,hpargs.gfilt_sigma,output=None)
	if resource_available:
		print 'memory use after hp: %s'%(memory_usage_resource())
	if img1.nodata_value:
		nodata_index = img1_img == img1.nodata_value
	else:
		nodata_index = img1_img == 0  # if there is no no-data value set, use 0 in the orginal image as no-data
	img1_img=None
	if resource_available:
		print 'memory use after del of img1_img: %s'%(memory_usage_resource())
	stddev=np.std(hp_arr_1)
	hp_min=np.min(hp_arr_1)
	hp_max=np.max(hp_arr_1)
	if resource_available:
		print 'memory use after minmaxstd: %s'%(memory_usage_resource())
# 	img_min=np.min(img1_img[img1_img!=0.0])
# 	img_max=np.max(img1_img)
	hp_larger=np.max([np.abs(hp_min), hp_max])
	if (hp_larger<=float(int_limit_maxm1)):  # scale out to +- hp_larger
# 		hp_arr_1=int_limit_maxm1 * (hp_arr_1/hp_larger)
		hp_arr_1 *= (float(int_limit_maxm1)/hp_larger)
		scaling=float(int_limit_maxm1)/hp_larger
# 		print 'scaling is +- %f mapped to +- %f'%(hp_larger, int_limit_maxm1)
	elif ((10.0*stddev)<=float(int_limit_maxm1)): # or scale out to +- 10*std
# 		hp_arr_1=int_limit_maxm1 * (hp_arr_1/(10.0*stddev))
		hp_arr_1 *= (float(int_limit_maxm1)/(10.0*stddev))
		scaling=float(int_limit_maxm1)/(10.0*stddev)
# 		print 'scaling is +- %f mapped to +- %f'%(10.0*stddev, int_limit_maxm1)
	else:
# 		hp_arr_1=hp_arr_1  # leave it alone and clip it  (min radiometric resolution in hp is 1 DN from the original image...)
		scaling=1.0
	if resource_available:
		print 'memory use after rescaling: %s'%(memory_usage_resource())
# 		print 'scaling is one to one DN image to DN hp, with hp clipped at +- %f'%(int_limit_maxm1)
# 	hp_arr_1=int_limit_maxm1 * (hp_arr_1/(2.0*stddev))
	print 'image %s hp std dev: %f  hp min %f hp max %f scaling %f hp DN : 1 img DN'%(b8_img_name,stddev,hp_min,hp_max,scaling)
	hp_arr_1[hp_arr_1<=float(int_limit_minp2)]=float(int_limit_minp2)  # clip to min value plus 2 - min value is used for no_data
	hp_arr_1[hp_arr_1>float(int_limit_maxm1)]=float(int_limit_maxm1)   # p2 and m1 are to help with float to int16 conversion issues (rounding won't produce wrapping)
	hp_arr_1=np.int16(hp_arr_1)
# 	if img1.nodata_value:
# 		hp_arr_1[img1_img == img1.nodata_value] = int_nodata_val
# 	else:
# 		hp_arr_1[img1_img == 0] = int_nodata_val  # if there is no no-data value set, use 0 in the orginal image as no-data
	hp_arr_1[nodata_index] = int_nodata_val
	if resource_available:
		print 'memory use just before tiff output: %s'%(memory_usage_resource())
	dict={} # set up tiff tag
	dict['TIFFTAG_IMAGEDESCRIPTION'] = "original image: %s"%(b8_img_name)+" python_script_parameters:" + " ".join(sys.argv) + " hp_scaling: %f "%(scaling) + ' Arg ' + str(hpargs)
	dst_filename = hpdir + '/' + hp1filename
	(out_lines,out_pixels)=hp_arr_1.shape
	out_bands=1
	compression_options=['COMPRESS=DEFLATE','PREDICTOR=1']
	dst_ds = driver.Create( dst_filename, out_pixels, out_lines, out_bands, gdal.GDT_Int16, options = compression_options)
	dst_ds.SetMetadata(dict)  # write tiff tag
	dst_ds.SetGeoTransform( img1.gt ) 
	dst_ds.SetProjection( img1.proj )
	dst_ds.GetRasterBand(1).SetNoDataValue( int_nodata_val )
	dst_ds.GetRasterBand(1).WriteArray( (hp_arr_1).astype('int16') )
	dst_ds = None # done, close the dataset
	# 		if psutil_available:
	# 			print 'img1_int16_hp written - using ',memory_usage_psutil(),memory_usage_resource()
# 	img1 = None
# 	img1_img = None
# 	hp_arr_1 = None
	if resource_available:
		print 'memory use just after tiff output: %s'%(memory_usage_resource())
	runtime=time.time()-process_time
# 	print 'done with image %s in %5.1f seconds\n'%(hp1filename, runtime)
	return runtime

# set up command line arguments
parser = argparse.ArgumentParser( \
    description="""high-pass filters Landsat band 8 image, writes out as signed int16,
    scaled to range of hp image, or if that gives output single DN > input single DN,
    to +-10 standard deviations (same restriction), or, if needed, clipped to int16 
    range (so one output DN = one input DN).
    
    output format: as a lossless (simple) compressed geotiff""",
    epilog='>>  <<',
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-B8_dir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='source dir for image [.]')
parser.add_argument('B8_file', 
                    action='store', 
                    type=str, 
                    default=None,
                    help='source B8 filename')
parser.add_argument('-hp_dir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='hp image output dir for both new, and possibly prior, hp images[.]')
parser.add_argument('-gfilt_sigma', 
                    action='store', 
                    type=float, 
                    default=3.0, 
                    help='gaussian filter sigma (standard deviation for Gaussian kernel) [%(default)f]')
# the rest of the parameters are flags - if raised, set to true, otherwise false
parser.add_argument('-v', 
                    action='store_true',  
                    default=False, 
                    help='verbose - extra diagnostic and image info put in log file [False if not raised]')
args = parser.parse_args()

outdir=args.hp_dir

if(args.v):
	print '#', args
	
indir=args.B8_dir
b8_img_name=args.B8_file
hpdir=args.hp_dir

if not(os.path.isdir(hpdir)):
	if not(os.path.isdir(os.path.dirname(hpdir))):
		print '>>>>>>>>>>>>>>>parent directory %s does not exist - halting'%(os.path.dirname(hpdir))
		sys.exit(1)
	else:
		os.makedirs(hpdir)
		print '>>>>>>>>>>>>>>>created directory %s '%(hpdir)
	
hp1filename=args.B8_file.replace('_B8.TIF','_B8_hp.tif')
deltime=make_hp_func([indir,b8_img_name,hpdir,hp1filename,args])

if resource_available:
	print 'file %s/%s processed in %f seconds using %s'%(hpdir,hp1filename,deltime,memory_usage_resource())
else:
	print 'file %s/%s processed in %f seconds'%(hpdir,hp1filename,deltime)



