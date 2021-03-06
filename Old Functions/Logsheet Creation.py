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
from astropy.io import ascii
import numpy as np
import math
import re
from astroquery.jplsbdb import SBDB
from sbident import SBIdent # pip install git+https://github.com/bengebre/sbident
from astropy.time import Time
from astropy.coordinates import SkyCoord
from prettytable import PrettyTable, NONE, HEADER #pip install prettytable
from tkinter import *
from tkinter import filedialog

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

#-----------------------------------------------------------------------------#
# Opens and extracts data from hdu files

def openfiles(date): #Opens and extracts data from hdu files
    folder = basefolder.replace('/','//')
    #####USED IN IDL#####folder = basefolder+'/{}/'.format(date)
    if os.path.isdir(folder) == False: 
        raise ValueError('Cannot find folder {}'.format(folder))
    ####USED IN IDL#####files = glob.glob(folder+ '.fits')
    files = glob.glob(basefolder+'\\*.fits')
    files = sorted(files, key=lambda _: re.sub(r'[^\d]', '', _)[-4:])
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
                        dp.loc[i,'Spectral Type'] = ''
                    elif 'flat' in f:
                        dp.loc[i,'Spectral Type'] = ''
                    else:
                        dp.loc[i,'Spectral Type'] = query.loc[0,'SP_TYPE']
                if 'arc' in f:
                    dp.loc[i,'%s Flux' %mag]=''
                elif 'flat' in f:
                    dp.loc[i,'%s Flux' %mag]=''
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
        z = float(dec[7:])* (1/3600)
        angle = math.radians(x-y-z)
        return(angle)
    else:
        x = float(dec[0:2])
        y = float(dec[3:5]) * (1/60)
        z = float(dec[6:]) * (1/3600)
        angle = math.radians(x+y+z)
        return(angle)

#Code to find distance between two objects
def distance(Dec_1, Dec_2, RA_1, RA_2):
    a = math.sin((Dec_2-Dec_1)/2)**2
    b = math.cos(Dec_1)*math.cos(Dec_2)*math.sin((RA_2-RA_1)/2)**2
    c = math.sqrt(a+b)
    dist = (10/math.pi)*abs(math.asin(c))
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
    scalers = [2/11,20]
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
def Get_Scores(dictionary, dp):
    target = []
    calibrator = []
    cals = []
    #This for loop pulls RA, Dec, UT time, and airmass for each object
    #and puts it in a list, that then gets added to either the target
    #or calibrator list defined above. Only does it for fixed sources
    for i in dictionary.keys():
        rows = dp.loc[(dp['Source Name'].str.lower() == i) & (dp['Object Type'] == 'fixed')]
        if rows.empty:
            if 'flatlamp' in i:
                temp = []
                cals_time = dictionary[i]['UT Time'][0]
                cals_hr = hours(cals_time)
                temp.append(cals_hr)
                temp.append(i)
                cals.append(temp)
            else:
                temp = []
                RA = dictionary[i]['RA'][0]
                Dec = dictionary[i]['Dec'][0]
                Time = dictionary[i]['UT Time'][0]
                temp.append(str(RA))
                temp.append(str(Dec))
                temp.append(str(Time))
                temp.append(str(i))
                try:
                    airmass = dictionary[i]['types']['target'][0]['airmass']
                except IndexError:
                    airmass = dictionary[i]['types']['calibrator'][0]['airmass']
                temp.append(str(airmass))
                target.append(temp)
                pass
        else:
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
    
    #picks the best calibrator for each target based on its score and adds a warning is the score is too high
    Best = []
    for i in Scores:
        pair = []
        a = min(i)
        b = i.index(a)
        pair.append(a)
        pair.append(b)
        if a >= 10:
            pair.append(True)
        else:
            pair.append(False)
        Best.append(pair)

    select_cals = []
    for i in target:
        diffs = []
        temp = []
        for j in cals:
            aux = abs(hours(i[2]) - j[0])
            diffs.append(aux)
        a = min(diffs)
        b = diffs.index(a)
        cals_name = cals[b][1]
        temp.append(b)
        temp.append(cals_name)
        temp.append(i[3])
        select_cals.append(temp)

    return Best, calibrator, target, select_cals
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

def add_batch(source, batch, final, prefix, airmass, calibration_number, ra_first, ra_last, dec_first, dec_last, ut_first, ut_last):

    # If we haven't see this Source before, create a place for it
    if not final.get(source):
        final[source] = {'types': {'calibration': [], 'calibrator': [], 'target': []}}

    # Add the data to the final
    for type in ['calibration', 'calibrator', 'target']:
        if batch[type]:
            final[source]['types'][type].append({'start': batch[type][0], 'end': batch[type][-1], 'airmass': airmass})
    final[source]['prefix'] = prefix
    final[source]['RA'] = [ra_first, ra_last]
    final[source]['Dec'] = [dec_first, dec_last]
    final[source]['UT Time'] = [ut_first, ut_last]
    if len(final[source]['types']['calibration']) != 0:
        calibration_number = calibration_number + 1

    return calibration_number
    
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
    calibration_number = 1

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
            source = 'flatlamp %s' % calibration_number
            prefix = 'flat/arc %s' % calibration_number
        ra = row['RA']
        dec = row['Dec']
        integration = row['Integration']
        airmass = row['Airmass']
        ut = row['UT Time']

        # If this iteration of the loop has a new Source, process the values we have been saving
        if source_old is not None and source.lower() != source_old.lower():

            # Add this batch to the final result
            calibration_number = add_batch(source_old.lower(), batch, final, prefix_old, airmass_old, calibration_number, ra_first, ra_old, dec_first, dec_old, ut_first, ut_old)
            
            # Start a new batch
            batch = {'calibration': [], 'calibrator': [], 'target': []}
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
        if integration < 1:
            type = 'calibration'
        elif integration < 70:
            type = 'calibrator'
        else:
            type = 'target'

        # Add the Number to this batch as the particular Type
        batch[type].append(number)

    # Add the last batch to the final
    add_batch(source.lower(), batch, final, prefix, airmass, calibration_number, ra_first, ra, dec_first, dec, ut_first, ut_old)

    return final

##############################################################################################

def get_source_list(dp,date):
    source_name_list=[]
    dpc=pandas.DataFrame()
    time_list=[]
    coord_list=[]
    index_list_source=[]
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
                time_list.append(str(date[0:4])+'-'+str(date[4:6])+'-'+str(date[6:])+' '+str(dpcopy.loc[i,'UT Time'])) 
                coord_list.append(str(dpcopy.loc[i,'RA'])+' '+str(dpcopy.loc[i,'Dec']))
                index_list_source.append(i)

    dpsl=pandas.DataFrame()  
    #avg_list=[]
    dpcc=pandas.DataFrame()
    mpc_obs='568'
    num=0
    for name in source_name_list:
        ra_list=[]
        dec_list=[]
        dpsl.loc[num,'Source Name']=name
        
        if dpcopy.loc[index_list_source[num], 'Object Type']== 'moving':
            
            try:
                SBDB.query(dpsl.loc[num, 'Source Name'])['message']
                time=Time(time_list[num])
                proper_coord=splat.properCoordinates(coord_list[num])
                center=SkyCoord(int(proper_coord.ra.deg), int(proper_coord.dec.deg), unit='deg')
                sbid=SBIdent(mpc_obs, time, center, maglim=20, hwidth=1)
                obj_name=sbid.results[0]['Object name']
                id_name=str(re.findall(r'\d+',obj_name )[0])
                try:
                    dpsl.loc[num,'Moving obj- Abs Mag ']= SBDB.query(id_name, phys=True)['phys_par']['H']
                except:
                    dpsl.loc[num,'Moving obj- Abs Mag ']='N/A'
                try:
                    dpsl.loc[num,'Moving obj- Diameter ']= SBDB.query(id_name, phys=True)['phys_par']['diameter']
                except:
                    dpsl.loc[num,'Moving obj- Diameter ']='N/A'
                try:
                    dpsl.loc[num,'Moving obj- Diameter_sig']= SBDB.query(id_name, phys=True)['phys_par']['diameter_sig']
                except:
                    dpsl.loc[num,'Moving obj- Diameter_sig']='N/A'
                try:
                    dpsl.loc[num,'Moving obj- Rotation period ']= SBDB.query(id_name, phys=True)['phys_par']['rot_per']
                except:
                    dpsl.loc[num,'Moving obj- Rotation period ']='N/A'
                try:
                    dpsl.loc[num,'Moving obj- SPK-ID']=SBDB.query(id_name)['object']['spkid']
                except:
                    dpsl.loc[num,'Moving obj- SPK-ID']='N/A'
            except:
                try:
                    dpsl.loc[num,'Moving obj- Abs Mag ']= SBDB.query(dpsl.loc[num,'Source Name'], phys=True)['phys_par']['H']
                except:
                    dpsl.loc[num,'Moving obj- Abs Mag ']='N/A'
                try:
                    dpsl.loc[num,'Moving obj- Diameter ']= SBDB.query(dpsl.loc[num,'Source Name'], phys=True)['phys_par']['diameter']
                except:
                    dpsl.loc[num,'Moving obj- Diameter ']= 'N/A'
                try:
                    dpsl.loc[num,'Moving obj- Diameter_sig']= SBDB.query(dpsl.loc[num,'Source Name'], phys=True)['phys_par']['diameter_sig']
                except:
                    dpsl.loc[num,'Moving obj- Diameter_sig']= 'N/A'
                try:
                    dpsl.loc[num,'Moving obj- Rotation period ']= SBDB.query(dpsl.loc[num,'Source Name'], phys=True)['phys_par']['rot_per']
                except:
                    dpsl.loc[num,'Moving obj- Rotation period ']= 'N/A'
                try:
                    dpsl.loc[num,'Moving obj- SPK-ID']=SBDB.query(dpsl.loc[num,'Source Name'])['object']['spkid']
                except:
                    dpsl.loc[num,'Moving obj- SPK-ID']='N/A'

            
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
            
    dpsl=get_source_list(dp, str(date))
    
    best, calibrators, targets, cals = Get_Scores(final, dp)
    #print('Best:', best)
    #print('calibrator:', calibrators)
    #print('target:', targets)
    #print('Selected Cals:', cals)
    
    for number, lists in enumerate(targets):
        name = lists[3]
        #rows = dp.loc[(dp['Source Name'].str.lower() == name) & (dp['Object Type'] == 'fixed')]
        #if rows.empty:
            #INSERT STUFF FOR MOVING TARGET IDENTIFICATION HERE
            #pass
        #else:
        prefix = final[name]['prefix']
        try:
            start_of_target = final[name]['types']['target'][0]['start']
        except IndexError:
            start_of_target = final[name]['types']['calibrator'][0]['start']
        try:
            end_of_target = final[name]['types']['target'][0]['end']
        except IndexError:
            end_of_target = final[name]['types']['calibrator'][0]['end']
        calibrator_index = best[number][1]
        calibrator_name = calibrators[calibrator_index][3]
        start_of_calibrator = final[calibrator_name]['types']['calibrator'][0]['start']
        end_of_calibrator = final[calibrator_name]['types']['calibrator'][0]['end']
        calibrator_rows = dp.loc[(dp['Source Name'].str.lower() == calibrator_name)]
        B_mag = round(calibrator_rows['B Flux'].iloc[0],2)
        V_mag = round(calibrator_rows['V Flux'].iloc[0],2)
        calibration_name = cals[number][1]
        start_of_calibration = final[calibration_name]['types']['calibration'][0]['start']
        end_of_calibration = final[calibration_name]['types']['calibration'][0]['end']
        calibration_range = '{0}-{1}'.format(start_of_calibration, end_of_calibration)
        calibrator_range = '{0}-{1}'.format(start_of_calibrator, end_of_calibrator)
        target_range = '{0}-{1}'.format(start_of_target, end_of_target)
        
        #print('---------------------------------')
        #print('Source Name:', name)
        #print('Calibrator Name:', calibrator_name)
        #print('Calibration Name:', calibration_name)
        #print('Prefix:', prefix)
        #print('Object Range:', start_of_target, '-', end_of_target)
        #print('Standard Range:',start_of_calibrator, '-', end_of_calibrator)
        #print('Calibration Range:', start_of_calibration, '-', end_of_calibration)
        #print('Standard B Mag:', B_mag)
        #print('Standard V Mag:', V_mag)
        #print('File Name:', 'spex_prism_%s_%s' % (name, date))
        
        if number == 0:   
            final_table = PrettyTable()
            final_table.field_names = ['cals', 'prefix1', 'obj', 'prefix2', 'std', 'ext. params',
                       'obj scl range', 'obj shape flag', 'std scl range', 
                       'std shape flag', 'std b mag', 'std v mag', 'shift flag', 
                       'shift range', 'filename', 'force flag']
        
        final_table.add_row([calibration_range, prefix, target_range, prefix, 
                             calibrator_range, '2.5,2,2.2,2,0', '1.4-1.8','0',
                             '1-2','1', B_mag, V_mag, '1', '1.75-2.05', 
                             'spex_prism_{0}_{1}'.format(name,date), '0'])
    
    #make final table look like how we need
    final_table.header = False
    final_table.hrules = HEADER
    
    input_file = 'C:\\Users\\peter\\Desktop\\input.txt'
    with open(input_file, 'w') as file:
        if cols == cols_old: 
            instrument = 'SpeX'
        if cols == cols_new: 
            instrument = 'uSpeX'
        lines = ['#','# instrument = '+instrument,'# rawpath = ','# calpath = ','# procpath = ','#']
        file.write('# \n')
        file.write('# instrument = '+instrument+'\n')
        file.write('# rawpath = \n')
        file.write('# calpath = \n')
        file.write('# procpath = \n')
        file.write('# \n')
        file.write(str(final_table))
        
    print('input file written to {}'.format(input_file))
    

    #dp, dpc, dpk = query_reference(dp, dpc)

    dp['Notes'] = ['']*len(files)

    #dpsl = source_list(final, dp)
    
    writer(date, dp, dpsl)
    
    return dp

#-----------------------------------------------------------------------------#

magnitudes=['B','V','J','H','K']

dirs = []
root = Tk()
root.attributes('-topmost',True)
root.geometry('200x100')
def select():
    root.directory = filedialog.askdirectory()
    dirs.append(root.directory)
    return root.directory
select_btn = Button(root, text='Select',command=select)
close_btn = Button(root, text = 'Exit', command=root.destroy)
select_btn.place(relx = 0.5, rely = 0.35, anchor=CENTER)
close_btn.place(relx = 0.5, rely = 0.65, anchor=CENTER)
root.mainloop()
print(dirs)

for directory in dirs:
    basefolder = str(directory)
    print('')
    print(basefolder)
    makelog(20011015)    
    
