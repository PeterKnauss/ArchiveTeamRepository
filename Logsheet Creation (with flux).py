# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 20:37:03 2022

@author: peter
"""

import pandas
import astroquery
from astroquery.simbad import Simbad
import astropy.coordinates as coord
import splat.database
from astropy import units as u
import glob
import os
import sys
from astropy.io import fits

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
basefolder='C:\\Users\\peter\\Desktop\\Raw Data'

def makelog(date):
    #folder = basefolder+'/{}/'.format(date)
    #if os.path.isdir(folder) == False: 
        #raise ValueError('Cannot find folder {}'.format(folder))
    #files = glob.glob(folder+'/data/*.fits')
    files = glob.glob(basefolder+'\\*.fits')
    if len(files) == 0: 
        raise ValueError('Cannot find any .fits data files in {}'.format(basefolder+'/data/'))

    dp = pandas.DataFrame()
    dpc=pandas.DataFrame()
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

        if old_name != name:
          coordinate=str(dp.loc[i,'RA']+' '+ dp.loc[i,'Dec'])
          for mag in magnitudes:
              proper_coord=splat.properCoordinates(coordinate)
              query=splat.database.querySimbad(proper_coord)
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
                  try:
                      flux=float(query['FLUX_%s' %mag])
                      dp.loc[i,'%s Flux' %mag]=flux
                  except:
                      dp.loc[i, '%s Flux' %mag]='N/A'
                      pass
        else:
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
    for i,l in enumerate(dp['Source Name']):
        coordinate=str(dp.loc[i,'RA']+' '+ dp.loc[i,'Dec'])
        proper=splat.properCoordinates(coordinate)

        if 'flat' in l:
            pass
        elif 'arc' in l:
            pass
        else:
            proper=splat.properCoordinates(coordinate)
            dpc.loc[i,'RA']=float(proper.ra.deg)
            dpc.loc[i,'DEC']=float(proper.dec.deg)
    dpc=splat.database.prepDB(dpc)
    dpq=splat.database.queryXMatch(dpc, catalog='2MASS', radius=30.*u.arcsec)
    dpk=splat.database.queryXMatch(dpq, catalog='Simbad', radius=30.*u.arcsec)    
    dp['Notes'] = ['']*len(files)    

    with pandas.ExcelWriter(folder+'logs_{}.xlsx'.format(date)) as writer:
        dp.sort_values('UT Time',inplace=True)
        dp.reset_index(inplace=True,drop=True)
        dp.to_excel(writer,sheet_name='Files',index=False)
        dpk.reset_index(inplace=True,drop=True)
        dpk.to_excel(writer, sheet_name='Target 2MASS & Simbad Query' )	
    print('log written to {}'.format(basefolder+'logs_{}.xlsx'.format(date)))

    return

makelog(10/15/2001)

# external function call
#if __name__ == '__main__':
 #   if len(sys.argv) > 1: 
  #      makelog(sys.argv[1])
