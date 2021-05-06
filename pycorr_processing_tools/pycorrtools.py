import numpy as np
from scipy.interpolate import RectBivariateSpline as RBS
import cv2
import gdal
import gdalconst as gdc  # constants for gdal - e.g. GA_ReadOnly, GA_Update ( http://www.gdal.org )
import os
# import subprocess as sp
# import string
# import random
# import sys
# import time
import datetime as dt
import argparse
import osr
from scipy.ndimage.filters import gaussian_filter
# from matplotlib.pyplot import get_cmap
# import netCDF4
# import re

##########################################################################################
#
# pycorrtools
# 
# library of functions used in image-to-image offset tracking, with an ice focus
#
# built initially from code written by Mark Fahnestock 2013-2020
# modeled on functionality of IMCORR (M. Fahnestock 1992-93) (a C wrapper for greycorr fortran fft-based correlation routine for tiepointing)
#
##########################################################################################









##########################################################################################
# chip_corr - returns subpixel location of peak match for source chip within target chip
###############
# 
# corr_setup_params = {
#                         'peak_block_halfsize':peak_block_halfsize,
#                         'rbs_halfsize':rbs_halfsize,
#                         'rbs_order':rbs_order,
#                         'offset_one_tenth_j':offset_one_tenth_j,
#                         'offset_one_tenth_i':offset_one_tenth_i,
#                         'offset_one_hundredth_j':offset_one_hundredth_j,
#                         'offset_one_hundredth_i':offset_one_hundredth_i,
#                         'd2xdx2_spacing':d2xdx2_spacing,
#                         'corr_nodata_val':corr_nodata_val,
#                         'curvature_nodata_value':curvature_nodata_value
#                     }
#
#       corr_return_values = {
#                               corr,       - correlation value at sub-pixel peak
#                               del_corr,   - corr minus next highest peak (peak_block_halfsize blocked out around peak)
#                               del_i,      - offset in i direction (decimal pixels)
#                               del_j,      - offset in j direction (decimal pixels)
#                               d2idx2,     - curvature of peak in i direction
#                               d2jdx2      - curvature of peak in j direction
#                            }
#                           
#           
#
##########################################################################################
def chip_corr(chip_src, chip_tar, **corr_setup_params):
        """ docstring not here yet...see comment above
        """
        # new_cent_loc is the integer number of pixels into the correlation surface in both i and j 
        #       to the "0 offset" location (the middle of the correlation surface array)
        # note peak1_loc array indicies are i,j instead of j,i, so new_cent_loc values switched here 
        #      - also resultCv.shape[::-1] reversed to get proper order below 
        new_cent_loc = [(chip_tar.shape[1] - chip_src.shape[1])/2.0, (chip_tar.shape[0] - chip_src.shape[0])/2.0]         
        resultCv=cv2.matchTemplate(chip_src, chip_tar, cv2.TM_CCOEFF_NORMED)
        mml_1=cv2.minMaxLoc(resultCv)
        peak1_corr=mml_1[1]
        peak1_loc=np.array(mml_1[3])

        tempres_peak1_zeroed=resultCv.copy()    # del_corr is difference in peak height and next highest peak with area around first peak zeroed out - need copy to zero out peak
        tempres_peak1_zeroed[(peak1_loc[1]-peak_block_halfsize):(peak1_loc[1]+peak_block_halfsize+1),(peak1_loc[0]-peak_block_halfsize):(peak1_loc[0]+peak_block_halfsize+1)]=0
                
        corr = peak1_corr
        del_corr = peak1_corr - cv2.minMaxLoc(tempres_peak1_zeroed)[1]
        # if peak is far enough from the edges of the correlation surface, fit locale of peak with a spline and find sub-pixel offset at 1/100th of a pixel level
        if((np.array([n-rbs_halfsize for n in peak1_loc])>=np.array([0,0])).all() & 
           (np.array([(n+rbs_halfsize) for n in peak1_loc ])<np.array(list(resultCv.shape[::-1]))).all()):# 
            # at this point we have the i,j loc of integer pixel offset, need to chase peak at sub-pixel level
            # so set up RBS spline function to interpolate correlation surface locally
            rbs_p=RBS(list(range(-rbs_halfsize,rbs_halfsize+1)),
                      list(range(-rbs_halfsize,rbs_halfsize+1)),
                      resultCv[(peak1_loc[1]-rbs_halfsize):(peak1_loc[1]+rbs_halfsize+1),
                               (peak1_loc[0]-rbs_halfsize):(peak1_loc[0]+rbs_halfsize+1)],
                      kx=rbs_order,
                      ky=rbs_order)
            spoffi=0.0
            spoffj=0.0
            
            ###################################################
            # first walk surface to find peak to tenths of a pixel
            ###################################################
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
                else:
                    # peak not in center, move
                    mmi=mml[3][0]
                    mmj=mml[3][1]
                    peak1_corr=mml[1]
                    spoffi=spoffi+offset_one_tenth_i[mmj,mmi]
                    spoffj=spoffj+offset_one_tenth_j[mmj,mmi]
            ###################################################
            # now walk surface to find peak to hundredths of a pixel (inside the tenth just found)
            ###################################################
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
                else:
                    # peak not in center, move
                    mmi=mml[3][0]
                    mmj=mml[3][1]
                    peak1_corr=mml[1]
                    spoffi=spoffi+offset_one_hundredth_i[mmj,mmi]
                    spoffj=spoffj+offset_one_hundredth_j[mmj,mmi]
            ####################################################
            # end walk surface
            ####################################################

            corr = peak1_corr
            del_corr = peak1_corr - cv2.minMaxLoc(tempres_peak1_zeroed)[1]
            deli=(peak1_loc[0]+spoffi)-new_cent_loc[0] # see new_cent_loc note above
            delj=(peak1_loc[1]+spoffj)-new_cent_loc[1] # see new_cent_loc note above
            d2i_pts=rbs_p.ev([spoffj,spoffj,spoffj],[spoffi-d2xdx2_spacing,spoffi,spoffi+d2xdx2_spacing])
            d2j_pts=rbs_p.ev([spoffj-d2xdx2_spacing,spoffj,spoffj+d2xdx2_spacing],[spoffi,spoffi,spoffi])
            d2idx2_t=(2.0*d2i_pts[1] - (d2i_pts[0]+d2i_pts[2]))/(np.power(d2xdx2_spacing,2.0))
            d2jdx2_t=(2.0*d2j_pts[1] - (d2j_pts[0]+d2j_pts[2]))/(np.power(d2xdx2_spacing,2.0))
        else: # peak is too close to an edge of correlation surface, so didn't do sub-pixel fit.
            deli=peak1_loc[0]-new_cent_loc[0] # see new_cent_loc note above
            delj=peak1_loc[1]-new_cent_loc[1] # see new_cent_loc note above
            d2idx2_t=curvature_nodata_value
            d2jdx2_t=curvature_nodata_value

        corr_return_values = {
                          'corr':corr,
                          'del_corr':del_corr,
                          'del_i':deli,
                          'del_j':delj,
                          'd2idx2':d2idx2_t,
                          'd2jdx2':d2jdx2_t
                       }
        
        return corr_return_values

# end chip_corr
