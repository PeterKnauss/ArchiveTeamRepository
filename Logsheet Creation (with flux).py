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
import math
from astroquery.jplsbdb import SBDB

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
    ####USED IN IDL#####files = glob.glob(folder+ '.fits')
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
#Code and function used to get the scores

#Code used to find RA
def RA_angle(ra):
    x = float(ra[0:2]) * 15
    y = float(ra[3:5]) * 0.25
    z = float(ra[6:10]) * (1/240)
    angle = math.radians(x+y+z)
    return(angle)

#Code for Dec
def Dec_angle(dec):
    if dec[0] == "-":
        x = float(dec[0:3])
        y = float(dec[4:6]) * (1/60)
        z = float(dec[7:10])* (1/3600)
        angle = math.radians(x-y-z)
        return(angle)
    else:
        x = float(dec[0:2])
        y = float(dec[3:5]) * (1/60)
        z = float(dec[7:9]) * (1/3600)
        angle = math.radians(x+y+z)
        return(angle)

#Code to find distance between two objects
def distance(Dec_1, Dec_2, RA_1, RA_2):
    a = math.sin((Dec_2-Dec_1)/2)**2
    b = math.cos(Dec_1)*math.cos(Dec_2)*math.sin((RA_2-RA_1)/2)**2
    c = math.sqrt(a+b)
    dist = 0.8*abs(math.asin(c))
    return(dist)

#code used to turn UT time into hours
def hours(time):
    x = float(time[0:2])
    y = float(time[3:5]) * (1/60)
    z = float(time[6:14]) * (1/3600)
    tm = x+y+z
    return(tm)

#generates the score form airmass and time
def partial_score(std_vals, obj_vals):
    Peter = 0
    scalers = [3,0.3]
    for i in range(len(scalers)):
        if std_vals[i] >= obj_vals[i]:
            aux = abs(std_vals[i]-obj_vals[i])
        else:
            aux = abs(obj_vals[i]-std_vals[i])
        aux = aux * scalers[i]
        Peter += aux
    return Peter

#generates final score
def Score(std,obj):
    RA_std = RA_angle(std[0])
    Dec_std = Dec_angle(std[1])
    RA_obj = RA_angle(obj[0])
    Dec_obj = Dec_angle(obj[1])
    dist = distance(Dec_std, Dec_obj, RA_std, RA_obj)
    time_std = hours(std[2])
    time_obj = hours(obj[2])
    std_list = [time_std, float(std[4])]
    obj_list = [time_obj, float(obj[4])]
    score = partial_score(std_list, obj_list)
    score = score + dist
    return score

#gives 2D list of scores given some dictionary
def Get_Scores(dictionary):
    target = []
    calibrator = []
    #This for loop pulls RA, Dec, UT time, and airmass for each object
    #and puts it in a list, that then gets added to either the target
    #or calibrator list defined above
    for i in dictionary.keys():
        temp = []
        RA = dictionary[i]['RA'][0]
        Dec = dictionary[i]['Dec'][0]
        Time = dictionary[i]['UT Time'][0]
        temp.append(str(RA))
        temp.append(str(Dec))
        temp.append(str(Time))
        temp.append(str(i))

    #logic used to identify if something is a target or a calibrator
        try: 
            calibrator_diff = int(dictionary[i]['types']['calibrator'][0]['end']) - int(dictionary[i]['types']['calibrator'][0]['start'])
        except IndexError:
            calibrator_diff = None
        try: 
            target_diff = int(dictionary[i]['types']['target'][0]['end']) - int(dictionary[i]['types']['target'][0]['start'])
        except IndexError:
            target_diff = None
    
        #print(calibrator_diff)
        #print(target_diff)

        if calibrator_diff == None:
            airmass = dictionary[i]['types']['target'][0]['airmass']
            temp.append(str(airmass))
            target.append(temp)
        
        elif target_diff == None:
            airmass = dictionary[i]['types']['calibrator'][0]['airmass']
            temp.append(str(airmass))
            calibrator.append(temp)

        elif calibrator_diff > target_diff:
            airmass = dictionary[i]['types']['calibrator'][0]['airmass']
            temp.append(str(airmass))
            calibrator.append(temp)
            
        elif target_diff > calibrator_diff:
            airmass = dictionary[i]['types']['target'][0]['airmass']
            temp.append(str(airmass))
            target.append(temp)

    #Creates a 2D list of scores where the columns are calibrators and
    #the rows are targets
    Scores = []
    for i in target:
        row = []
        for j in calibrator:
            a = Score(i,j)
            row.append(a)
        Scores.append(row)
    
    #picks the best calibrator for each target based on its score
    Best = []
    for i in Scores:
        pair = []
        a = min(i)
        b = i.index(a)
        pair.append(a)
        pair.append(b)
        Best.append(pair)

    return Best, calibrator, target
#-----------------------------------------------------------------------------#
# Write dataframes to an excel sheet

def writer(date, dp, dpsl):
    with pandas.ExcelWriter(basefolder+'logs_{}.xlsx'.format(date)) as writer:
        dp.sort_values('UT Time',inplace=True)
        dp.reset_index(inplace=True,drop=True)
        dp.to_excel(writer,sheet_name='Files',index=False)
        dpsl.reset_index(inplace=True,drop=True)
        dpsl.to_excel(writer, sheet_name='Source List', index=False )
          
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
    ra_first = final[source]['RA'][0]
    ra_last = final[source]['RA'][-1]
    dec_first = final[source]['Dec'][0]
    dec_last = final[source]['Dec'][-1]
    ra_first_prop = float(splat.properCoordinates(str(ra_first) + ' ' + str(dec_first)).ra.deg)
    ra_last_prop = float(splat.properCoordinates(str(ra_last) + ' ' + str(dec_last)).ra.deg)
    dec_first_prop = float(splat.properCoordinates(str(ra_first) + ' ' + str(dec_first)).dec.deg)
    dec_last_prop = float(splat.properCoordinates(str(ra_last) + ' ' + str(dec_last)).dec.deg)
        
    ra_diff = abs(ra_last_prop - ra_first_prop)
    dec_diff = abs(dec_last_prop - dec_first_prop)
    
    first_time = hours(final[source]['UT Time'][0])
    last_time = hours(final[source]['UT Time'][-1])
    time_diff = last_time - first_time
    
    #print(time_diff)
    
    ra_time_diff = ra_diff / time_diff
    dec_time_diff = dec_diff / time_diff
    
    if ra_time_diff > 0.003 or dec_time_diff > 0.003:
        object_type = 'moving'
    else:
        object_type = 'fixed'
        
    #if source in ['2001 be10','1620','20790','2000 xl44','110','29']:
    #print(source)
    #print('ra_diff')
    #print(ra_time_diff)
    #print('dec_diff')
    #print(dec_time_diff)
    #print('------------------')
        
    return object_type
    
###############################################################################

# Batch Creation Subroutines

###############################################################################
# Create Batches

def add_batch(source, batch, final, prefix, airmass, ra_first, ra_last, dec_first, dec_last, ut_first, ut_last):

    # If we haven't see this Source before, create a place for it
    if not final.get(source):
        final[source] = {'types': {'calibrator': [], 'target': []}}

    # Add the data to the final
    for type in ['calibrator', 'target']:
        if batch[type]:
            final[source]['types'][type].append({'start': batch[type][0], 'end': batch[type][-1], 'airmass': airmass})
    final[source]['prefix'] = prefix
    final[source]['RA'] = [ra_first, ra_last]
    final[source]['Dec'] = [dec_first, dec_last]
    final[source]['UT Time'] = [ut_first, ut_last]

#-----------------------------------------------------------------------------#
# Create Dictionaries for each batch

def create_dictionaries(dp):

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
    ut_first = None
    ut_old = None

    # Filter data by Mode
    data = []
    for index in np.arange(0, len(dp.index)):
         data.append(dp.iloc[index])

    # If the filter left us with nothing, give up
    #if not data:
        #print('There is no data with mode %s' %mode)
        #return None # comment out this line if the data has empty cells in MODE or has no mode
                    # and add data.append(dp.iloc[index]) to take account into data with no mode

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
        ut = row['UT Time']

        # If this iteration of the loop has a new Source, process the values we have been saving
        if source_old is not None and source.lower() != source_old.lower():

            # Add this batch to the final result
            add_batch(source_old.lower(), batch, final, prefix_old, airmass_old, ra_first, ra_old, dec_first, dec_old, ut_first, ut_old)

            # Start a new batch
            batch = {'calibrator': [], 'target': []}
            ra_first = ra
            dec_first = dec
            ut_first = row['UT Time']
       
        if ra_first == None:
            ra_first = ra
            dec_first = dec
            ut_first = ut

        # Remember the last Source we saw, so we can tell if it changes with the next line
        prefix_old = prefix
        source_old = source
        ra_old = ra
        dec_old = dec
        airmass_old = airmass
        ut_old = ut

        # The actual smarts -- figure out what type the star is based on the Integration
        if integration < 70:
            type = 'calibrator'
        else:
            type = 'target'

        # Add the Number to this batch as the particular Type
        batch[type].append(number)

    # Add the last batch to the final
    add_batch(source.lower(), batch, final, prefix, airmass, ra_first, ra, dec_first, dec, ut_first, ut_old)

    return final

##############################################################################################

def get_source_list(dp):
    source_name_list=[]
    dpc=pandas.DataFrame()
    for i,l in enumerate(dp['Source Name']):
        dpcopy=dp.copy(deep=True)
        dpcopy.sort_values('UT Time',inplace=True)
        dpcopy.reset_index(inplace=True,drop=True)
        coordinate=str(dpcopy.loc[i,'RA']+' '+ dpcopy.loc[i,'Dec'])
        proper=splat.properCoordinates(coordinate)
        dpc.loc[i,'RA']=float(proper.ra.deg)
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

    dpsl=pandas.DataFrame()  
    #avg_list=[]
    dpcc=pandas.DataFrame()
    num=0
    for name in source_name_list:
        ra_list=[]
        dec_list=[]
        dpsl.loc[num,'Source Name']=name
        try:
            #Horizons works too. For example,
            #dpsl.loc[num , 'Moving object- Abs Mag ']= float(Horizons(id=dpsl.loc[num,'Source Name']).ephemerides()['H'])
            dpsl.loc[num,'Moving obj- Abs Mag ']= float(SBDB.query(dpsl.loc[num,'Source Name'], phys=True)['phys_par']['H'])
            dpsl.loc[num,'Moving obj- Diameter ']= SBDB.query(dpsl.loc[num,'Source Name'], phys=True)['phys_par']['diameter']
            dpsl.loc[num,'Moving obj- Diameter_sig']= SBDB.query(dpsl.loc[num,'Source Name'], phys=True)['phys_par']['diameter_sig']
            dpsl.loc[num,'Moving obj- Rotation period ']= SBDB.query(dpsl.loc[num,'Source Name'], phys=True)['phys_par']['rot_per']
        except:
            dpsl.loc[num,'Moving obj- Abs Mag ']='N/A'                    
            dpsl.loc[num,'Moving obj- Diameter ']= 'N/A'
            dpsl.loc[num,'Moving obj- Diameter_sig']= 'N/A'
            dpsl.loc[num,'Moving obj- Rotation period ']= 'N/A'
            
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
    dpc=splat.database.prepDB(dpcc)
    dpq=splat.database.queryXMatch(dpc, catalog='2MASS', radius=30.*u.arcsec)
    dpj=splat.database.queryXMatch(dpq, catalog='Simbad', radius=30.*u.arcsec)
    for i, f in enumerate(dpj['DESIGNATION']):
        for cols in list(dpj.columns):
            dpsl.loc[i,cols]=dpj.loc[i,cols]
            
    return dpsl

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
        
    final = create_dictionaries(dp)
    dp['Object Type'] = ['']*len(files)
    
    #print(final)
        
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
            
    dpsl=get_source_list(dp)
    
    best, calibrator, target = Get_Scores(final)
    #print('Best')
    #print(best)
    #print('calibrator')
    #print(calibrator)
    #print('target')
    #print(target)
    
    for number, lists in enumerate(target):
        name = lists[3]
        rows = dp.loc[(dp['Source Name'].str.lower() == name) & (dp['Object Type'] == 'moving')]
        if rows.empty:
            #INSERT STUFF FOR MOVING TARGET IDENTIFICATION HERE
            pass
        else:
            prefix = final[name]['prefix']
            start_of_target = final[name]['types']['target'][0]['start']
            end_of_target = final[name]['types']['target'][0]['end']
            calibrator_index = best[number][1]
            calibrator_name = calibrator[calibrator_index][3]
            start_of_calibrator = final[calibrator_name]['types']['calibrator'][0]['start']
            end_of_calibrator = final[calibrator_name]['types']['calibrator'][0]['end']
            calibrator_rows = dp.loc[(dp['Source Name'].str.lower() == calibrator_name)]
            B_mag = calibrator_rows['B Flux'][0]
            V_mag = calibrator_rows['V Flux'][0]
            
            print('Prefix'
                  '%s'
                  
                  'Object Range'
                  '%s - %s' 
                  
                  'Standard Range'
                  '%s - %s' 
                  
                  'Standard B Mag'
                  '%s'
                  
                  'Standard V Mag'
                  '%s'
                  
                  'File Name'
                  'spex_prism_%s_%s'
                  
                  (prefix, start_of_target, end_of_target, start_of_calibrator, end_of_calibrator, B_mag, V_mag, name, date))
                  
                  
    #dp, dpc, dpk = query_reference(dp, dpc)

    dp['Notes'] = ['']*len(files)

    #dpsl = source_list(final, dp)
    
    
    writer(date, dp, dpsl)
    
    return dp

#-----------------------------------------------------------------------------#

makelog(10/15/2001)
