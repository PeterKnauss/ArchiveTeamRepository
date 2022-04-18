# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 20:37:03 2022

@author: peter
"""

import pandas
import math
import astroquery
from astroquery.simbad import Simbad
import astropy.coordinates as coord
import splat.database
from astropy import units as u
import glob
import os
import sys
from astropy.io import fits
from openpyxl import load_workbook
col_top = ['PROG_ID','OBSERVER','DATE_OBS']
cols_new = {
    'Source Name': 'TCS_OBJ',
    'RA': 'TCS_RA',
    'Dec': 'TCS_DEC',
    'UT Time': 'TIME_OBS',
    'MJD': 'MJD_OBS',
    'HA': 'TCS_HA',
    'PA': 'POSANGLE',
    'Parallactic': 'TCS_PA',
    'Airmass': 'TCS_AM',
    'Integration': 'ITIME',
    'Coadds': 'CO_ADDS',
    'Type': 'DATATYPE',
    'Slit': 'SLIT',
    'Mode': 'GRAT',
    'Program': 'PROG_ID',
    'Observer': 'OBSERVER',
}

cols_old = {
    'Source Name': 'OBJECT',
    'RA': 'RA',
    'Dec': 'DEC',
    'UT Time': 'TIME_OBS',
    'HA': 'HA',
    'PA': 'POSANGLE',
    'Airmass': 'AIRMASS',
    'Integration': 'ITIME',
    'Coadds': 'CO_ADDS',
    'Slit': 'SLIT',
    'Mode': 'GRAT',
    'Observer': 'OBSERVER',
}

# basefolder is where you downloaded the fits files
# For example:
basefolder='C:\\Users\\peter\\Desktop\\Raw Data'
def makelog(date):
    folder = basefolder+'/{}/'.format(date)
    if os.path.isdir(folder) == False: 
        raise ValueError('Cannot find folder {}'.format(folder))
    #files = glob.glob(folder+'/data/*.fits')
    files = glob.glob(folder+'*.fits')
    if len(files) == 0: 
        raise ValueError('Cannot find any .fits data files in {}'.format(folder+'/data/'))
    dpc=pandas.DataFrame()
    dpk=pandas.DataFrame()
    dp = pandas.DataFrame()
    magnitudes=['B','V','J','H','K']
    dp['File'] = [f.split('/')[-1] for f in files]
#dp['File number'] = [int(f.split('.')[-3]) for f in files]

# select which columns to use
    hdu = fits.open(files[0])
    h = hdu[0].header
    hdu.close()
    if 'TCS_OBJ' in list(h.keys()): cols = cols_new
    else: cols = cols_old

# make log	
    for c in list(cols.keys()): 
        dp[c] = ['']*len(files)
    old_name = None
    for i,f in enumerate(files):
        hdu = fits.open(f)
        hdu.verify('silentfix')
        h = hdu[0].header
        hdu.close()
        
        for c in list(cols.keys()):
            dp.loc[i,c] = h[cols[c]]
        if 'arc' in f: 
            dp.loc[i,'Source Name'] = 'arclamp'
        if 'flat' in f: 
            dp.loc[i,'Source Name'] = 'flat field'

        name = str(dp.loc[i,'Source Name'])
        coordinate=str(dp.loc[i,'RA']+' '+ dp.loc[i,'Dec'])
        
  # find target magnitude or flux with proper coordinates. using 2MASS catalog for J,H,K flux and querySimbad for B,V.      

        if old_name!= name: # if old_name is not name, run this loop to find magnitude/flux.

            for mag in magnitudes:
                proper_coord=splat.properCoordinates(coordinate)
                query=splat.database.querySimbad(proper_coord)
                query2MASS=splat.database.queryVizier(proper_coord, catalog='2MASS', radius=30*u.arcsec, nearest=True)
                if len(query.columns) == 0:
                    dp.loc[i,'Object Type'] = 'N/A'
                    dp.loc[i,'Spectral Type'] = 'N/A'
                else:
                    if 'arc' in f:
                        dp.loc[i,'Object Type'] = 'None'
                        dp.loc[i,'Spectral Type'] = 'None'
                    elif 'flat' in f:
                        dp.loc[i,'Object Type'] = 'None'
                        dp.loc[i,'Spectral Type'] = 'None'
                    else:
                        dp.loc[i,'Object Type'] = query.loc[0,'OTYPE']
                        dp.loc[i,'Spectral Type'] = query.loc[0,'SP_TYPE']
                if 'arc' in f:
                    dp.loc[i,'%s Flux' %mag]='None'
                elif 'flat' in f:
                    dp.loc[i,'%s Flux' %mag]='None'
                else:
                    if mag in ['B','V']:
                        # using 'try' statement in case for queried objects not in database (i.e that returns an empty table, or raises an error)
                        try: 
                            flux=float(query['FLUX_%s' %mag])
                            dp.loc[i,'%s Flux' %mag]=flux
                        except:
                            if mag in ['J','H','K']:
                                flux=float(query2MASS['%smag' %mag])
                                dp.loc[i,'%s Flux' %mag]=flux
                            pass
                        
                     
        else: # For old_name the same as name run the loop below to repeat the known object type, spectral type, and magnitudes.
            for mag in magnitudes:
                dp.loc[i,'Object Type'] = dp.loc[int(i-1),'Object Type']
                dp.loc[i,'Spectral Type'] = dp.loc[int(i-1),'Spectral Type']
                if 'arc' in f:
                    dp.loc[i,'%s Flux' %mag]='None'
                elif 'flat' in f:
                    dp.loc[i,'%s Flux' %mag]='None'
                else:
                    if dp.loc[int(i-1),'%s Flux' %mag] == 'N/A':
                        dp.loc[i,'%s Flux' %mag] = 'N/A'
                    else:
                        flux = dp.loc[int(i-1),'%s Flux' %mag]
                        dp.loc[i,'%s Flux' %mag] = flux
            
        old_name = name
        
  # 2MASS catalog & Simbad query reference for target sources.
    for i,l in enumerate(dp['Source Name']):
        
        dpcopy=dp.copy(deep=True)
        dpcopy.sort_values('UT Time',inplace=True)
        dpcopy.reset_index(inplace=True,drop=True)
        coordinate=str(dpcopy.loc[i,'RA']+' '+ dpcopy.loc[i,'Dec'])
        proper=splat.properCoordinates(coordinate)
        query2MASS=splat.database.queryVizier(proper_coord, catalog='2MASS', radius=30*u.arcsec, nearest=True)
        dpk.loc[i,'Source Name']= dpcopy.loc[i,'Source Name']
        dpc.loc[i,'RA']=float(proper.ra.deg)
        dpc.loc[i,'DEC']=float(proper.dec.deg)
        
    dpc=splat.database.prepDB(dpc)
    dpq=splat.database.queryXMatch(dpc, catalog='2MASS', radius=30.*u.arcsec)
    dpj=splat.database.queryXMatch(dpq, catalog='Simbad', radius=30.*u.arcsec)

    for i, f in enumerate(dpj['DESIGNATION']):
        
        for cols in list(dpj.columns):
            
            dpk.loc[i, cols]=dpj.loc[i, cols]
            
            if 'flat' in dpk.loc[i, 'Source Name']:
                dpk.loc[i,cols]=''
            if 'arc' in dpk.loc[i, 'Source Name']:
                dpk.loc[i, cols]=''
# The for..if...if loop is a rough code trying to replace JHK mags by checking the error jhk-mags. 
# basically if the cell is empty then don't use the JHKmag there, if it's not empty then replace the JHK mag on the "Files" sheet with the 2MASS catalog
# This definitely needs some corrections for files with moving sources...
        for mag in ['J','H','K']:
            
            if ((dpk.loc[i,'2MASS_e_Jmag'])!= False) and (dpk.loc[i,'2MASS_e_Hmag']!= '') and (dpk.loc[i,'2MASS_e_Kmag']!= ''):
                
                if dpk.loc[i,'SIMBAD_main_id']==False:
                    dp.loc[i,'%s Flux' %mag]=dpk.loc[i,'2MASS_%smag' %mag]
            else:
                if dpk.loc[i,'SIMBAD_main_id']==False:
                    if dpk.loc[i,'Source Name'] == dpk.loc[i+1, 'Source Name']:                                   
                        dp.loc[i, '%s Flux' %mag] = dpk.loc[i+1, '%s Flux' %mag]
                    
    dp['Notes'] = ['']*len(files)


    # Write dataframe in an excel sheet. 

    with pandas.ExcelWriter(folder+'logs_{}.xlsx'.format(date)) as writer:
        dp.sort_values('UT Time',inplace=True)
        dp.reset_index(inplace=True,drop=True)
        dp.to_excel(writer,sheet_name='Files',index=False)
        dpk.reset_index(inplace=True,drop=True)
        dpk.to_excel(writer, sheet_name='Target 2MASS & Simbad Query' )
          
    print('log written to {}'.format(folder+'logs_{}.xlsx'.format(date)))
    return
#have code make score for source-callibrator pair
#turns RA into an angle in radians
def RA_angle(ra):
    x = float(ra[0:2]) * 15
    y = float(ra[3:5]) * 0.25
    z = float(ra[6:10]) * (1/240)
    angle = math.radians(x+y+z)
    return(angle)

#turns Dec into an angle in radians
def Dec_angle(dec):
    if dec[0] == "-":
        x = float(dec[0:3])
        y = float(dec[4:6]) * (1/60)
        z = float(dec[7:])* (1/3600)
        angle = math.radians(x-y-z)
        return(angle)
    else:
        x = float(dec[0:2])
        y = float(dec[3:5]) * (1/60)
        z = float(dec[6:]) * (1/3600)
        angle = math.radians(x+y+z)
        return(angle)

#Gets the part of the score from airmass and time
def partial_score(std_vals, obj_vals):
    score = 0
    scalers = [1,2]
    for i in range(len(scalers)):
        if std_vals[i] >= obj_vals[i]:
            aux = abs(std_vals[i]-obj_vals[i])
        else:
            aux = abs(obj_vals[i]-std_vals[i])
        aux = aux * scalers[i]
        score += aux
    return(score)

#turns Dec and Ra into distance
def distance(Dec_1, Dec_2, RA_1, RA_2):
    a = math.sin((Dec_2-Dec_1)/2)**2
    b = math.cos(Dec_1)*math.cos(Dec_2)*math.sin((RA_2-RA_1)/2)**2
    c = math.sqrt(a+b)
    dist = abs(2*math.asin(c))
    return(dist)

#turns time into hours
def hours(time):
    x = float(time[0:2])
    y = float(time[3:5]) * (1/60)
    z = float(time[6:14]) * (1/3600)
    tm = x+y+z
    return(tm)

#gets the final score of the pair
#and this is the function that will be used for the final calculation of score
#takes a list of 4 strings: ['RA', 'Dec','UT Time', 'airmass']
#So an example of a list it could take is ['10:16:13.70','29:18:41.9','05:52:16.323347','1.024']

def Score(std,obj):
    RA_std = RA_angle(std[0])
    Dec_std = Dec_angle(std[1])
    RA_obj = RA_angle(obj[0])
    Dec_obj = Dec_angle(obj[1])
    dist = distance(Dec_std, Dec_obj, RA_std, RA_obj)
    time_std = hours(std[2])
    time_obj = hours(obj[2])
    std_list = [time_std, float(std[3])]
    obj_list = [time_obj, float(obj[3])]
    score = partial_score(std_list, obj_list)
    score = score + dist
    return(score)

makelog(10/15/2001)

# external function call
#if __name__ == '__main__':
 #   if len(sys.argv) > 1: 
  #      makelog(sys.argv[1])
