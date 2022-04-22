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
from openpyxl import load_workbook
import numpy as np

#-----------------------------------------------------------------------------#
# Values for new and old columns

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
###############################################################################

#Log Creation Subroutines#

###############################################################################
# Original Settings Required

# basefolder is where you downloaded the fits files
# For example:
basefolder='C:\\Users\\peter\\Desktop\\Raw Data'
magnitudes=['B','V','J','H','K']

#-----------------------------------------------------------------------------#
# Opens and extracts data from hdu files

def openfiles(date): #Opens and extracts data from hdu files
    folder = basefolder
    #####USED IN IDL#####folder = basefolder+'/{}/'.format(date)
    if os.path.isdir(folder) == False: 
        raise ValueError('Cannot find folder {}'.format(folder))
    ####USED IN IDL#####files = glob.glob(folder+'/data/*.fits')
    files = glob.glob(basefolder+'\\*.fits')
    if len(files) == 0: 
        raise ValueError('Cannot find any .fits data files in {}'
                         .format(folder+'/data/'))
    dpc = pandas.DataFrame()
    dp = pandas.DataFrame()
    dp['File'] = [f.split('/')[-1] for f in files]

# select which set of columns to use
    hdu = fits.open(files[0])
    h = hdu[0].header
    hdu.close()
    if 'TCS_OBJ' in list(h.keys()): cols = cols_new
    else: cols = cols_old
    
    return cols, files, dp, dpc

#-----------------------------------------------------------------------------#
# Finds target magnitude with proper coordinates using 2Mass catalog
# J,H,K flux and Simbad for B,V

def magnitude_get(i, dp, old_name, f):
        name = str(dp.loc[i,'Source Name'])
        coordinate=str(dp.loc[i,'RA']+' '+ dp.loc[i,'Dec'])
        
    # if old_name is not name, run this loop to find magnitude/flux.
        if old_name!= name: 

            for mag in magnitudes:
                proper_coord=splat.properCoordinates(coordinate)
                query=splat.database.querySimbad(proper_coord)
                query2MASS=splat.database.queryVizier(proper_coord, 
                            catalog='2MASS', radius=30*u.arcsec, nearest=True)
                if len(query.columns) == 0:
                    dp.loc[i,'Spectral Type'] = 'N/A'
                else:
                    if 'arc' in f:
                        dp.loc[i,'Spectral Type'] = 'None'
                    elif 'flat' in f:
                        dp.loc[i,'Spectral Type'] = 'None'
                    else:
                        dp.loc[i,'Spectral Type'] = query.loc[0,'SP_TYPE']
                if 'arc' in f:
                    dp.loc[i,'%s Flux' %mag]='None'
                elif 'flat' in f:
                    dp.loc[i,'%s Flux' %mag]='None'
                else:
                    if mag in ['B','V']:
  
   # Use 'try' statement in case for queried objects not in database 
   # (i.e that returns an empty table, or raises an error)
                        try: 
                            flux=float(query['FLUX_%s' %mag])
                            dp.loc[i,'%s Flux' %mag]=flux
                        except:
                            dp.loc[i, '%s Flux' %mag]='N/A'
                            pass
                    if mag in ['J','H','K']:
                        try:
                            flux=float(query2MASS['%smag' %mag])
                            dp.loc[i,'%s Flux' %mag]=flux
                        except:
                            dp.loc[i, '%s Flux' %mag]='N/A'
                            pass
                        
    # For old_name the same as name run the loop below to repeat the known 
    # object type, spectral type, and magnitudes.
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
        
        return dp

#-----------------------------------------------------------------------------#
# 2Mass catalog and Simbad query reference for target sources

def query_reference(dp, dpc):
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

    return dp, dpc, dpk

#-----------------------------------------------------------------------------#
# Write dataframes to an excel sheet

def writer(date, dp, dpk, dpsl):
    with pandas.ExcelWriter(basefolder+'logs_{}.xlsx'.format(date)) as writer:
        dp.sort_values('UT Time',inplace=True)
        dp.reset_index(inplace=True,drop=True)
        dp.to_excel(writer,sheet_name='Files',index=False)
        dpk.reset_index(inplace=True,drop=True)
        dpk.to_excel(writer, sheet_name='Target 2MASS & Simbad Query',index = False)
        #dpsl.to_excel(writer, sheet_name='Source List', index = False)
        dpsl2.reset_index(inplace=True,drop=True)
        dpsl2.to_excel(writer, sheet_name='Source List2', index=False )
          
    print('log written to {}'.format(basefolder+'logs_{}.xlsx'.format(date)))
    return

#-----------------------------------------------------------------------------#
# Create new DataFrame for source list

#def source_list(final, dp):
#    dpsl = pandas.DataFrame()
#    for i, sources in enumerate(final):
#        if 'flat' in sources:
#            pass
#        else:
#            dpsl.loc[i,'Source Name'] = sources
#            dpsl.loc[i,'RA'] = final[sources]['ra'][0]
#            dpsl.loc[i,'Dec'] = final[sources]['dec'][-1]

#    return dpsl
        
#-----------------------------------------------------------------------------#
# Moving vs Fixed

def moving_fixed(final, source):
    ra_first = final[source]['ra_first']
    ra_last = final[source]['ra_last']
    dec_first = final[source]['dec_first']
    dec_last = final[source]['dec_last']
    ra_first_prop = float(splat.properCoordinates(str(ra_first) + ' ' + str(dec_first)).ra.deg)
    ra_last_prop = float(splat.properCoordinates(str(ra_last) + ' ' + str(dec_last)).ra.deg)
    dec_first_prop = float(splat.properCoordinates(str(ra_first) + ' ' + str(dec_first)).dec.deg)
    dec_last_prop = float(splat.properCoordinates(str(ra_last) + ' ' + str(dec_last)).dec.deg)
        
    ra_diff = abs(ra_last_prop - ra_first_prop)
    dec_diff = abs(dec_last_prop - dec_first_prop)
    
    if ra_diff > 0.001 or dec_diff > 0.001:
        object_type = 'moving'
    else:
        object_type = 'fixed'
        
    #if source in ['2001 be10','1620','20790','2000 xl44','110','29']:
        #print(source)
        #print('ra_diff')
        #print(ra_diff)
        #print('dec_diff')
        #print(dec_diff)
        #print('------------------')
        
    return object_type
    
###############################################################################

# Batch Creation Subroutines

###############################################################################
# Create Batches

def add_batch(source, batch, final, prefix, airmass, ra_first, ra_last, dec_first, dec_last, uttime):

    # If we haven't see this Source before, create a place for it
    if not final.get(source):
        final[source] = {'types': {'calibrator': [], 'target': []}}

    # Add the data to the final
    for type in ['calibrator', 'target']:
        if batch[type]:
            final[source]['types'][type].append({'start': batch[type][0], 'end': batch[type][-1], 'airmass': airmass})
    final[source]['prefix'] = prefix
    final[source]['ra_last'] = ra_last
    final[source]['ra_first'] = ra_first
    final[source]['dec_last'] = dec_last
    final[source]['dec_first'] = dec_first
    final[source]['UT Time'] = uttime

#-----------------------------------------------------------------------------#
# Create Dictionaries for each batch

def create_dictionaries(mode, dp):

    # Set up some variables
    final = {}
    batch = {'calibration': [], 'calibrator': [], 'target': []}
    prefix_old = None
    source_old = None
    airmass_old = None
    ra_old = None
    dec_old = None
    ra_first = None
    dec_first = None
    uttime = None

    # Filter data by Mode
    data = []
    for index in np.arange(0, len(dp.index)):
        if mode == dp.iloc[index]['Mode']:
            data.append(dp.iloc[index])

    # If the filter left us with nothing, give up
    if not data:
        print('There is no data with mode %s' %mode)
        return None

    # Run through each row of the filtererd data
    for row in data:

        # Get the individual values, using part of File for Source if Source is generic
        prefix = row['File'][0:-11]
        number = row['File'][-11:-7]
        source = row['Source Name']
        if source == 'Object_Observed':
            source = row['File'][0:-11]
        if source in ['flat field', 'arclamp']:
            continue
            #source = 'flatlamp'
            #prefix = 'flat/arc'
        ra = row['RA']
        dec = row['Dec']
        integration = row['Integration']
        airmass = row['Airmass']

        # If this iteration of the loop has a new Source, process the values we have been saving
        if source_old is not None and source.lower() != source_old.lower():

            # Add this batch to the final result
            add_batch(source_old.lower(), batch, final, prefix_old, airmass_old, ra_first, ra_old, dec_first, dec_old, uttime)

            # Start a new batch
            batch = {'calibrator': [], 'target': []}
            ra_first = ra
            dec_first = dec
            uttime = row['UT Time']
       
        if ra_first == None:
            ra_first = ra
            dec_first = dec
            uttime = row['UT Time']

        # Remember the last Source we saw, so we can tell if it changes with the next line
        prefix_old = prefix
        source_old = source
        ra_old = ra
        dec_old = dec
        airmass_old = airmass

        # The actual smarts -- figure out what type the star is based on the Integration
        if integration <= 70.0:
            type = 'calibrator'
        else:
            type = 'target'

        # Add the Number to this batch as the particular Type
        batch[type].append(number)

    # Add the last batch to the final
    add_batch(source.lower(), batch, final, prefix, airmass, ra_first, ra, dec_first, dec, uttime)

    return final
##############################################################################################
def get_source_name_list(dp):
    source_name_list=[]
    for i,l in enumerate(dp['Source Name']):
        dpcopy=dp.copy(deep=True)
        dpcopy.sort_values('UT Time',inplace=True)
        dpcopy.reset_index(inplace=True,drop=True)
        coordinate=str(dpcopy.loc[i,'RA']+' '+ dpcopy.loc[i,'Dec'])
        proper=splat.properCoordinates(coordinate)
        dpc.loc[i.'RA']=float(proper.ra.deg)
        dpc.loc[i,'DEC']=float(proper.dec.deg)
        if 'flat' in dpcopy.loc[i, 'Source Name']:
            pass
        elif 'arc' in dpcopy.loc[i,'Source Name']:
            pass
        else:
            if  dpcopy.loc[i,'Source Name' ] in source_name_list:
                pass
           
            else:
                source_name_list.append(dpcopy.loc[i,'Source Name'])
    return source_name_list  

def get_avg(source_name_list):
    dpsl2=pandas.DataFrame()  
    avg_list=[]
    dpcc=pandas.DataFrame()
    num=0
    for name in source_name_list:
        ra_list=[]
        dec_list=[]
        dpsl2.loc[num,'Source Name']=name
        num=num+1
        for i,l in enumerate(dpcopy['Source Name']):
        
            if dpcopy.loc[i, 'Source Name'] in name:
                ra_list.append(dpc.loc[i,'RA'])
                dec_list.append(dpc.loc[i,'DEC'])
        
        ra_sum=sum(ra_list)
        dec_sum=sum(dec_list)
        ra_avg=ra_sum/(len(ra_list))
        dec_avg=dec_sum/(len(dec_list))
        dpcc=dpcc.append({'DEC':dec_avg, 'RA':ra_avg }, ignore_index=True)
    dpc=splat.database.preDB(dpcc)
    dpq=splat.database.queryXMatch(dpc, catalog='2MASS', radius=30.*u.arcsec)
    dpj=splat.database.queryXMatch(dpq, catalog='Simbad', radius=30.*u.arcsec)
    for i, f in enumerate(dpj['DESIGNATION']):
        for cols in list(dpj.columns):
            dpsl2.loc[i,cols]=dpj.loc[i,cols]
    return dpsl2
###############################################################################
# Actually make the log

def makelog(date):
    
    cols, files, dp, dpc = openfiles(date)
    for c in list(cols.keys()): 
        dp[c] = ['']*len(files)
    old_name = None
    object_type = None
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
        
    final = create_dictionaries('LowRes15',dp)
    dp['Object Type'] = ['']*len(files)
    
    print(final)
        
    for i,f in enumerate(files):
        source = dp.loc[i,'Source Name'].lower()
        if source == 'arclamp':
            dp.loc[i,'Object Type'] = 'Calibration'
        elif source == 'flat field':
            dp.loc[i,'Object Type'] = 'Calibration'
        else:
            object_type = moving_fixed(final, source)
            dp.loc[i,'Object Type'] = object_type
        if object_type == 'fixed':
            dp = magnitude_get(i, dp, old_name, f)
    snl=get_source_name_list(dp)
    get_avg(snl)
    #dp, dpc, dpk = query_reference(dp, dpc)

    dp['Notes'] = ['']*len(files)

    #dpsl = source_list(final, dp)
    
    
    writer(date, dp, dpsl2)
    
    return dp

#-----------------------------------------------------------------------------#

makelog(10/15/2001)
