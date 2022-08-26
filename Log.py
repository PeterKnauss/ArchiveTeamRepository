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
import errno
import shutil
from pathlib import Path
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import math
import re
from astroquery.jplsbdb import SBDB
from sbident import SBIdent # pip install git+https://github.com/bengebre/sbident
from astropy.time import Time
from astropy.coordinates import SkyCoord
from prettytable import PrettyTable, NONE, HEADER #pip install prettytable
#from tkinter import *
#from tkfilebrowser import askopendirnames # pip install tkfilebrowser
import traceback

#-----------------------------------------------------------------------------#
# Values for new and old columns

col_top = ['PROG_ID','OBSERVER','DATE_OBS']
cols_new = {
    'Source Name': 'TCS_OBJ',
    'RA': 'TCS_RA',
    'Dec': 'TCS_DEC',
    'UT Date': 'DATE_OBS',
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
    'UT Date': 'DATE_OBS',
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
    folder = basefolder
    #####USED IN REMOTE MACHINE#####
    #basefolder='/data/SpeX'
    #folder = basefolder+'/{}/'.format(date)
    if os.path.isdir(folder) == False: 
        raise ValueError('Cannot find folder {}'.format(folder))
    ####USED IN REMOTE MACHINE#####files = glob.glob(folder+ '*.fits')
    files = glob.glob(basefolder+'/*.fits')
    files = sorted(files, key=lambda _: re.sub(r'[^\d]', '', _)[-4:])
    if len(files) == 0: 
        raise ValueError('Cannot find any .fits data files in {}'
                         .format(folder+'/data/'))
    dpc = pandas.DataFrame()
    dp = pandas.DataFrame()
    dp['File'] = [f.split('/')[-1] for f in files]

# select which set of columns to use
    hdu = fits.open(files[0], ignore_missing_end=True)
    h = hdu[0].header
    hdu.close()
    if 'TCS_OBJ' in list(h.keys()):
        cols = cols_new
        global Cols_New_Label
        Cols_New_Label = True
    else: cols = cols_old
    
    return cols, files, dp, dpc

#-----------------------------------------------------------------------------#
# Finds target magnitude with proper coordinates using 2Mass catalog
# J,H,K flux and Simbad for B,V

def magnitude_get(i, dp, old_name, f):
        name = str(dp.loc[i,'Source Name'])
        coordinate=str(dp.loc[i,'RA']+' '+ dp.loc[i,'Dec'])
    # if old_name is not name, run this loop to find magnitude/flux.
        if old_name == None or old_name!= name: 
            for mag in magnitudes:
                proper_coord=splat.properCoordinates(coordinate)
                query=splat.database.querySimbad(proper_coord, nearest=True)
                query2MASS=splat.database.queryVizier(proper_coord, 
                            catalog='2MASS', radius=30*u.arcsec, nearest=True)
                if len(query.columns) == 0:
                    dp.loc[i,'Spectral Type'] = 'N/A'
                else:
                    dp.loc[i,'Spectral Type'] = query.loc[0,'SP_TYPE']
                    if dp.loc[i,'Spectral Type'] == '':
                        query = splat.database.querySimbad(name, nearest=True, isname=True)
                        if len(query.columns) == 0:
                            dp.loc[i,'Spectral Type'] = 'N/A'
                        else:
                            dp.loc[i,'Spectral Type'] = query.loc[0,'SP_TYPE']
  
   # Use 'try' statement in case for queried objects not in database 
   # (i.e that returns an empty table, or raises an error)
                if mag in ['B','V']:
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
                #dp.loc[i,'Object Type'] = dp.loc[int(i-1),'Object Type']
                dp.loc[i,'Spectral Type'] = dp.loc[int(i-1),'Spectral Type']
                if dp.loc[int(i-1),'%s Flux' %mag] == 'N/A':
                    dp.loc[i,'%s Flux' %mag] = 'N/A'
                else:
                    flux = dp.loc[int(i-1),'%s Flux' %mag]
                    dp.loc[i,'%s Flux' %mag] = flux
            
        old_name = name
        
        return dp, old_name

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
    score_partial = 0
    scalers = [2,20]
    for i in range(len(scalers)):
        if std_vals[i] >= obj_vals[i]:
            aux = abs(std_vals[i]-obj_vals[i])
        else:
            aux = abs(obj_vals[i]-std_vals[i])
        aux = aux * scalers[i]
        score_partial += aux
    return score_partial

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
def Get_Scores(dictionary, dp, date, proc_path, format_input, reduction):
    target = []
    calibrator = []
    cals = []
    #This for loop pulls RA, Dec, UT time, and airmass for each object
    #and puts it in a list, that then gets added to either the target
    #or calibrator list defined above. Only does it for fixed sources
    for i in dictionary.keys():
        if i == 'empty':
            i = ''
            rows = dp.loc[(dp['Source Name'].str.lower() == i) & (dp['Object Type'] == 'fixed')]
            i = 'empty'
        else: rows = dp.loc[(dp['Source Name'].str.lower() == i) & (dp['Object Type'] == 'fixed')]
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
                    calibrator_diff = int(dictionary[i]['types']['calibrator'][0]['end']) - int(dictionary[i]['types']['calibrator'][0]['start'])
                except IndexError:
                    calibrator_diff = None
                try: 
                    target_diff = int(dictionary[i]['types']['target'][0]['end']) - int(dictionary[i]['types']['target'][0]['start'])
                except IndexError:
                    target_diff = None

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

    #check to see if there are no targets or no calibrators in the dataset
    if len(target) == 0:
        dpsl_notarget = get_source_list(dp, str(date))
        writer(proc_path, date, dp, dpsl_notarget, format_input, reduction)
        raise Exception('There are no targets in dataset {}'.format(date))
    if len(calibrator) == 0:
        dpsl_nocalibs = get_source_list(dp, str(date))
        writer(proc_path, date, dp, dpsl_nocalibs, format_input, reduction)
        raise Exception('There are no calibrators in dataset {}'.format(date))

    #Creates a 2D list of scores where the columns are calibrators and
    #the rows are targets
    scores = []
    for i in target:
        row = []    
        for j in calibrator:
            a = Score(i,j)
            row.append(a)
        scores.append(row)

    #picks the best calibrator for each target based on its score and adds a warning is the score is too high
    best = []
    for i in scores:
        pair = []
        mode_same_check_best = False
        mode_best_index = 0
        while mode_same_check_best == False:
            scores_sorted = sorted(i)
            a = scores_sorted[mode_best_index]
            b = i.index(a)
            target_scores_index = scores.index(i)
            target_scores_name = target[target_scores_index][3]
            calibrator_scores_name = calibrator[b][3]
            mode_best_target = dictionary[target_scores_name]['Mode']
            mode_best_calibrator = dictionary[calibrator_scores_name]['Mode']
            print(target_scores_name, calibrator_scores_name, mode_best_target, mode_best_calibrator)
            if mode_best_target == mode_best_calibrator:
                mode_same_check_best = True
            else:
                mode_same_check_best = False
                mode_best_index += 1
        pair.append(a)
        pair.append(b)
        if a >= 10:
            pair.append(True)
        else:
            pair.append(False)
        best.append(pair)

    #picks best calibration set for each target
    select_cals = []

    for i in target:
        mode_same_check_cals = False
        mode_cals_index = 0
        temp = []
        while mode_same_check_cals == False:
            diffs = []
            for j in cals:
                aux = abs(hours(i[2]) - j[0])
                diffs.append(aux)
            diffs_sorted = sorted(diffs)
            a = diffs_sorted[mode_cals_index]
            b = diffs.index(a)
            cals_name = cals[b][1]
            mode_cals_target = dictionary[i[3]]['Mode']
            mode_cals_cal = dictionary[cals_name]['Mode']
            if mode_cals_target == mode_cals_cal:
                mode_same_check_cals = True
            else:
                mode_same_check_cals = False
                mode_cals_index += 1
        temp.append(b)
        temp.append(cals_name)
        temp.append(i[3])
        select_cals.append(temp)

    return best, calibrator, target, select_cals
#-----------------------------------------------------------------------------#
# Write dataframes to an excel sheet

def writer(proc_path, date, dp, dpsl, format_input, reduction_path):
    if format_input.lower() == 'excel':
        with pandas.ExcelWriter(proc_path+'/logs_{}.xlsx'.format(date)) as writer:
            dp.reset_index(inplace=True,drop=True)
            dp.to_excel(writer,sheet_name='Files',index=False)
            dpsl=dpsl.drop(['DEC','RA'], axis=1)
            dpsl.reset_index(inplace=True,drop=True)
            dpsl.to_excel(writer, sheet_name='Source List', startrow=1, header=False, index=False )
            print('log written to {}'.format(proc_path+'/logs_{}.xlsx'.format(date)))

        if os.path.isdir('/data/SpeX/LOGS'):
            try:
                shutil.copy(proc_path+'/logs_{}.xlsx'.format(date), '/data/SpeX/LOGS/logs_{}.xlsx'.format(date))
                print('log written to /data/SpeX/LOGS as xlsx')
            except PermissionError:
                try:
                    shutil.copy(proc_path+'/logs_{}.xlsx'.format(date),reduction_path+'/LOGS/logs_{}.xlsx'.format(date))
                    print('log written to {} as xlsx'.format(reduction_path+'/LOGS'))
                except:
                    raise Exception('You do not have permission to add files to /data/SpeX/LOGS and do not have a LOGS folder at {}'.format(reduction_path))
    
    else:
        dp.reset_index(inplace=True,drop=True)
        dp.to_csv(proc_path+'/logs_{}.csv'.format(date),index=False)
        dpsl=dpsl.drop(['DEC','RA'], axis=1)
        dpsl.reset_index(inplace=True,drop=True)
        dpsl.to_csv(proc_path+'/source_list_{}.csv'.format(date), index=False )
        print('log and source list written to {} as csv'.format(proc_path))

        if os.path.isdir('/data/SpeX/LOGS'):
            try:
                shutil.copy(proc_path+'/logs_{}.csv'.format(date),'/data/SpeX/LOGS/logs_{}.csv'.format(date))
                shutil.copy(proc_path+'/source_list_{}.csv'.format(date),'/data/SpeX/LOGS/source_list_{}.csv'.format(date))
                print('log and source list written to /data/SpeX/LOGS as csv')
            except PermissionError:
                try:
                    shutil.copy(proc_path+'/logs_{}.csv'.format(date),reduction_path+'/LOGS/logs_{}.csv'.format(date))
                    shutil.copy(proc_path+'/source_list_{}.csv'.format(date),reduction_path+'/LOGS/source_list_{}.csv'.format(date))
                    print('log and source list written to {} as csv'.format(reduction_path+'/LOGS'))
                except:
                    raise Exception('You do not have permission to add files to /data/SpeX/LOGS and do not have a LOGS folder at {}'.format(reduction_path))
        
    #Dropped columns : 'DEC' and 'RC' - since 'Coordinates' column already covered this info.
          
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
    
    try:
        ra_time_diff = ra_diff / time_diff
        dec_time_diff = dec_diff / time_diff
    except ZeroDivisionError:
        return False
    
    total_time_diff = math.sqrt(float(ra_time_diff)**2 + float(dec_time_diff)**2)

    if total_time_diff > 0.015:
        object_type = 'moving'
    else:
        object_type = 'fixed'
        
    return object_type
    
###############################################################################

# Batch Creation Subroutines

###############################################################################
# Create Batches

def add_batch(source, batch, final, prefix, airmass, calibration_number, ra_first, ra_last, dec_first, dec_last, ut_first, ut_last, mode):

    # Check if source already exists, and if it does, check that its mode matches the new source's mode (edge case)
    # Protects against sources of different modes with the same name
    if 'flatlamp' in source:
        source = 'flatlamp %s' % calibration_number
        prefix = 'flat/arc %s' % calibration_number
    elif final.get(source) and final[source]['Mode'] != mode:
            source = source+'_'+mode

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
    final[source]['Mode'] = mode
    if 'flatlamp' in source:
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
    row_old = None
    calibration_number = 1
    spectral_type_B_check = False

    data = []
    for index in np.arange(0, len(dp.index)):
         data.append(dp.iloc[index])

    # Run through each row of the filtererd data
    for row in data:

        # Get the individual values, using part of File for Source if Source is generic
        if Cols_New_Label == True:
            prefix = row['File'][0:-12]
        else:
            prefix = row['File'][0:-11]
        if Cols_New_Label == True:
            number = row['File'][-12:-7].lstrip('0')
        else:
            number = row['File'][-11:-7].lstrip('0')
        source = row['Source Name']
        if source == 'Object_Observed':
            source = row['File'][0:-11]
        if source == '':
            source = 'Empty'
        if source in ['flat field', 'arclamp']:
            source = 'flatlamp %s' % calibration_number
            prefix = 'flat/arc %s' % calibration_number
        ra = row['RA']
        dec = row['Dec']
        integration = row['Integration']
        airmass = row['Airmass']
        ut = row['UT Time']
        mode = row['Mode']

        # If this iteration of the loop has a new Source, process the values we have been saving
        if source_old is not None and source.lower() != source_old.lower():          
            # Add this batch to the final result
            calibration_number = add_batch(source_old.lower(), batch, final, prefix_old, airmass_old, calibration_number, ra_first, ra_old, dec_first, dec_old, ut_first, ut_old, mode_old)
            
            # Start a new batch
            batch = {'calibration': [], 'calibrator': [], 'target': []}
            ra_first = ra
            dec_first = dec
            ut_first = row['UT Time']
        
        elif source_old is not None and any(x in source for x in ['flat','arc']) and any(x in source_old for x in ['flat','arc']):
            if mode != mode_old:
                # Add this batch to the final result
                calibration_number = add_batch(source_old.lower(), batch, final, prefix_old, airmass_old, calibration_number, ra_first, ra_old, dec_first, dec_old, ut_first, ut_old, mode_old)               

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
        if 'flatlamp' in source:
            source_old = 'flatlamp %s' % calibration_number
            prefix_old = 'flat/arc %s' % calibration_number
        else:
            prefix_old = prefix
            source_old = source
        if '.b.' in row['File']:
            try:
                ra_old = row_old['RA']
                dec_old = row_old['Dec']   
            except TypeError:   
                ra_old = ra
                dec_old = dec
        else:
            ra_old = ra
            dec_old = dec
        airmass_old = airmass
        ut_old = ut
        row_old = row
        mode_old = mode

        # The actual smarts -- figure out what type the star is
        spectral_type = 'N/A'
        spectral_type_flag = False
        if 'flatlamp' in source:
            type = 'calibration'
        else:
            if ra == '' or dec == '' or airmass == '':
                raise Exception("There is a source with blank coordinates or airmass. Please remove or manually reduce.")
            try:
                spectral_type = splat.database.querySimbad(source, nearest=True, isname=True)['SP_TYPE'][0]
                spectral_type_flag = True
            except KeyError:
                try:
                    proper_coord=splat.database.properCoordinates(str(ra)+' '+str(dec))
                    spectral_type = splat.database.querySimbad(proper_coord, nearest=True)['SP_TYPE'][0]
                    spectral_type_flag = True
                except KeyError:
                    pass
            if spectral_type != 'N/A' and any(_ in spectral_type for _ in ['A', 'F', 'G','B']) and spectral_type_flag == False:
                if 'B' in spectral_type and spectral_type_B_check == False:
                    print('Spectral Type B used as a calibrator')
                    spectral_type_B_check = True #Makes this print only happen once
                type = 'calibrator'        
            else:
                type = 'target'

        # Add the Number to this batch as the particular Type
        batch[type].append(number)

    # Add the last batch to the final
    add_batch(source_old.lower(), batch, final, prefix_old, airmass_old, calibration_number, ra_first, ra_old, dec_first, dec_old, ut_first, ut_old, mode_old)

    return final
#---------------------------------------------------------------------------------------------
#Create cals and proc folders

# This function is designed to run on the remote machine. os.getcwd will get the current directory.
# Assuming we ran this code outside any folder (ie home/user/) we want to direct ourselves to reductions folder and create
# a folder with the input date. Then create proc and cals folder in it. 



# READ, if you want to know the detail instruction
#--------------------------------------------------------------------------------#

# That is : login to remote machine --> open terminal --> run this code to create folders

# The method I used was to add a New  'Untitled.py' file from one of the python editor files (can find one through 'Assessories'
# from 'Application' icon on the top left corner --> click 'MATE Search Tool' and just search 'python', or 'editor') 
# then I copied and pasted my code from jupyter notebook onto the .py file I just created (don't have to email the codes to yourself every time now).

# Once you have the code --> save it as some 'filename.py' in the user directory or elsewhere if you prefer, then in the terminal type 'cd' to get to 
# the start directory (or you cd to the path if you save the python code elsewhere like 'Desktop')

# I ran python on the remote machine terminal by typing 'python'.
# The format to run the code : you have to import 'filename', then type  'filename'.'function name'(arguments)
# For example : i saved my codes as createfolder.py with function create_folder(date) 
# Then in terminal i'd type : python --> import createfolder --> createfolder.create_folder(20001010) 

#--------------------------------------------------------------------------------#

def create_folder(path, date, path_input, overwrite_input):
    if path_input == '/data/SpeX/':
        proc_directory='proc'
        cals_directory='cals'
        reduction_directory=os.path.join(str(Path.home()), 'reductions')
        LOGS_directory=os.path.join(reduction_directory, 'LOGS')
        parental_directory= os.path.join(reduction_directory, str(date) )
    
        #should be something like : /home/user/reductions/(the date)/
    
        #Check if 'reductions' is created. if not print a comment and ask for one
        #and then checks if the 'date' folder already exist in reductions folder
        try:
            os.mkdir(LOGS_directory)
        except OSError as err:
            if err.errno == errno.ENOENT:
                print('FileNotFound, Create a reductions folder')
            if err.errno == errno.EEXIST:
                pass
            else:
                raise
        try:
            os.mkdir(parental_directory)
        except OSError as err:
            if err.errno == errno.ENOENT:
                print('FileNotFound, Create a reduction folder')
            if err.errno == errno.EEXIST:
                print('FileExistError: %s file already exist in reductions folder' %date)
            else:
                raise
        
        proc_path=os.path.join(parental_directory, proc_directory)
        cals_path=os.path.join(parental_directory, cals_directory)
        raw_path=os.path.join('/data/SpeX/', str(date))
        
        # Check if data exist in raw_path
        try: 
            os.mkdir(raw_path)
        except OSError as err:
            if err.errno == errno.ENOENT:
                print('FileNotFound: %s does not exist in raw folder' %date)
            if err.errno == errno.EEXIST:
                print('FileExist: Good ! File exist in raw folder.')
            else:
               raise

        # Making the folders
        try:
            os.mkdir(proc_path)
        except OSError as err:
            if err.errno == errno.ENOENT:
                print('FileNotFound : Missing some folders, please check your path')
            if err.errno == errno.EEXIST:
                if overwrite_input.lower() == 'true' or overwrite_input.lower() == 'yes':
                    proc_input = 'yes'
                elif overwrite_input.lower() == 'false' or overwrite_input.lower() == 'no':
                    proc_input == 'no'
                else:
                    print('FileExistError : File exist, Do you want to OVERWRITE proc folder? Please enter yes or no: ')
                    proc_input=input()
                if 'yes' in proc_input.lower():
                    shutil.rmtree(proc_path)
                    os.makedirs(proc_path)
                    print('File Renewed')
                else:
                    print('Nothing Changed')
            else:
               raise

        try:
            os.mkdir(cals_path)
        except OSError as err:
            if err.errno == errno.ENOENT:
                print('FileNotFound : Missing some folders, please check your path')
            if err.errno == errno.EEXIST:
                if overwrite_input.lower() == 'true' or overwrite_input.lower() == 'yes':
                    cals_input = 'yes'
                elif overwrite_input.lower() == 'false' or overwrite_input.lower() == 'no':
                    cals_input = 'no'
                else:
                    print('FileExistError : File exist, Do you want to OVERWRITE cals folder? Please enter yes or no: ')
                    cals_input=input()
                if 'yes' in cals_input.lower():
                    shutil.rmtree(proc_path)
                    os.makedirs(proc_path)
                    print('File Renewed')
                else:
                    print('Nothing Changed')
            else:
               raise
    
    else:
        reduction_directory = ''
        raw_path = path
        cals_path = os.path.dirname(path) + '/{} cals'.format(date)
        proc_path = os.path.dirname(path) + '/{} proc'.format(date)
        os.mkdir(cals_path)
        os.mkdir(proc_path)
    

    return raw_path, cals_path, proc_path, reduction_directory
#---------------------------------------------------------------------------------------------
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
        try:
            proper=splat.properCoordinates(coordinate)
            dpc.loc[i,'RA']=float(proper.ra.deg)
            dpc.loc[i,'DEC']=float(proper.dec.deg)
        except ValueError:
            dpc.loc[i,'RA'] = ''
            dpc.loc[i,'DEC'] = ''
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
    dpcc=pandas.DataFrame()
    mpc_obs='568'
    num=0
    for name in source_name_list:
        ra_list=[]
        dec_list=[]
        dpsl.loc[num,'Source Name']=name
        dpsl.loc[num,'Object Type']=''
        dpsl.loc[num,'Date']=str(date[0:4])+'-'+str(date[4:6])+'-'+str(date[6:])
        for i,k in enumerate(dpcopy['Source Name']):
            if name == dpcopy.loc[i,'Source Name']:
                if bool(dpsl.loc[num,'Object Type'])== False:
                    dpsl.loc[num, 'Object Type']=dpcopy.loc[i,'Object Type']

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

def makelog(raw_path, cals_path, proc_path, date, format_input, reduction):
    
    cols, files, dp, dpc = openfiles(date)
    for c in list(cols.keys()): 
        dp[c] = ['']*len(files)
    old_name = None
    object_type = None
    for i,f in enumerate(files):
        hdu = fits.open(f,ignore_missing_end=True)
        hdu.verify('fix')
        h = hdu[0].header
        hdu.close()
        for c in list(cols.keys()):
            try:
                dp.loc[i,c] = h[cols[c]]
            except KeyError:
                dp.loc[i,c] = ''
        if 'arc' in f: 
            dp.loc[i,'Source Name'] = 'arclamp'
        if 'flat' in f: 
            dp.loc[i,'Source Name'] = 'flat field'  

    #Create a new dataframe that includes only prefixes and indices of all rows in our dataframe.
    #Then check this new dataframe for duplicates to check for duplicate indices, which spextoolette cant handle
    dp_duplicate_check = pandas.DataFrame(dp['File'].str[0:-7])
    dp_duplicate_check.drop_duplicates(subset='File',inplace=True)
    #If this new duplicated dp and the original dp have the same number of rows, then nothing was dropped and all is good
    #If not, raise an exception
    if len(dp.index) != len(dp_duplicate_check.index):
        raise Exception('Prefix and Index duplicate detected, please fix or reduce manually') 
        
    #Sort dataframe by UT Data and Time before creating the final dictionary
    dp.sort_values(by=['UT Date','UT Time'],inplace=True)
    dp.reset_index(inplace=True,drop=True)

    final = create_dictionaries(dp)
    dp['Object Type'] = ['']*len(files)
    
    #print(final)
    
    #Iterate through all the files and set object type and magnitude
    for i,f in enumerate(files):
        source = dp.loc[i,'Source Name'].lower()
        if 'arc' in source or 'flat' in source:
            dp.loc[i,'Object Type'] = 'Calibration'
            object_type = 'calibration'
        else:
            if source == '':
                source = 'empty'
            if source == 'object_observed':
                source = dp.loc[i,'File'][0:-11].lower()
            object_type = moving_fixed(final, source)
            #If there is a single picture error, print all but the last seven columns (Object/Spectral Type, Fluxes)
            if object_type == False:
                dp_limited = dp.drop(columns=dp.columns[-7:])
                dpsl_limited = get_source_list(dp, str(date))
                dpsl_limited.drop(columns='Object Type',inplace=True)
                writer(proc_path, date, dp_limited, dpsl_limited, format_input, reduction)
                raise Exception('Single picture detected. Limited spreadsheet printed, please manually reduce.')
            dp.loc[i,'Object Type'] = object_type

        if object_type == 'fixed':
            dp, old_name = magnitude_get(i, dp, old_name, f)
        else:
            old_name = str(dp.loc[i,'Source Name'])
            
    #Makes source list based on main dataframe
    dpsl=get_source_list(dp, str(date))
    
    #Matches up targets with correct calibrators and calibration sets
    best, calibrators, targets, cals = Get_Scores(final, dp, date, proc_path, format_input, reduction)

    #Check for sources that are logged as "fixed" but don't show up on Simbad. If there is a source 
    #that is logged as "moving" on the same night, change the fixed targets to moving.
    if 'moving' in dp['Object Type'].unique():
        fixed_target_indices = dp.loc[(dp['Object Type'] == 'fixed') & (dp['Spectral Type'] == 'N/A')].index.values
        for index in fixed_target_indices:
            dp['Object Type'].loc[index] = 'moving'
            dp['Spectral Type'].loc[index] = ''
            dp['B Flux'].loc[index] = ''
            dp['V Flux'].loc[index] = ''
            dp['J Flux'].loc[index] = ''
            dp['H Flux'].loc[index] = ''
            dp['K Flux'].loc[index] = ''
    
    #Initializing the dataframe with the columns that we need 
    data_table = pandas.DataFrame(columns = ['cals', 'prefix1', 'obj', 'prefix2', 'std', 'ext. params',
                            'obj scl range', 'obj shape flag', 'std scl range', 
                            'std shape flag', 'std b mag', 'std v mag', 'scalelines flag', 'shift flag', 
                            'shift range', 'filename', 'force flag']) 

    for target_index, lists in enumerate(targets,0):
        name = lists[3]
        for index_number in range(len(final[name]['types']['target'])):   
            if index_number == 0:
                display_name = name.replace(' ','')
            else:
                name_number = int(index_number+1)
                display_name = name.replace(' ','')+'-{}'.format(name_number)
            prefix = final[name]['prefix']
            mode = final[name]['Mode']
            start_of_target = final[name]['types']['target'][index_number]['start']
            end_of_target = final[name]['types']['target'][index_number]['end']
            calibrator_index = best[target_index][1]
            calibrator_name = calibrators[calibrator_index][3]
            calibrator_prefix = final[calibrator_name]['prefix']
            start_of_calibrator = final[calibrator_name]['types']['calibrator'][0]['start']
            end_of_calibrator = final[calibrator_name]['types']['calibrator'][0]['end']
            calibrator_mode = str(final[calibrator_name]['Mode'])
            if calibrator_mode in calibrator_name:
                calibrator_name = calibrator_name.rsplit('_')[0]
            calibrator_rows = dp.loc[(dp['Source Name'].str.lower() == calibrator_name) & (dp['Mode'] == calibrator_mode)]
            try:
                V_mag = round(float(calibrator_rows['V Flux'].iloc[0]),2)
            except ValueError:
                V_mag = 'N/A'
            try:
                B_mag = round(float(calibrator_rows['B Flux'].iloc[0]),2)
            except ValueError:
                B_mag = 'N/A'
            if pandas.isnull(V_mag) or V_mag == 'N/A' or V_mag == '':
                V_mag = B_mag
            calibration_name = cals[target_index][1]
            start_of_calibration = final[calibration_name]['types']['calibration'][0]['start']
            end_of_calibration = final[calibration_name]['types']['calibration'][0]['end']
            calibration_range = '{0}-{1}'.format(start_of_calibration, end_of_calibration)
            calibrator_range = '{0}-{1}'.format(start_of_calibrator, end_of_calibrator)
            target_range = '{0}-{1}'.format(start_of_target, end_of_target)        
            
            if 'long' in mode.lower() or 'short' in mode.lower() or 'lxd' in mode.lower() or 'sxd' in mode.lower():
                data_table.loc[len(data_table.index)] = ['# '+ calibration_range, prefix, target_range, calibrator_prefix, 
                                 calibrator_range, '2.5,2,2.2,2,0', '1.0-1.7','0',
                                 '1-2','1', B_mag, V_mag, '0', '1', '1.75-2.05', 
                                 'spex-prism_{0}_{1}'.format(display_name,date), '1']
            else:
                data_table.loc[len(data_table.index)] = [calibration_range, prefix, target_range, calibrator_prefix, 
                                     calibrator_range, '2.5,2,2.2,2,0', '1.0-1.7','0',
                                     '1-2','1', B_mag, V_mag, '0', '1', '1.75-2.05', 
                                     'spex-prism_{0}_{1}'.format(display_name,date), '1']
    
    input_file = proc_path+'/input_{}.txt'.format(date)
    with open(input_file, 'w') as file:
        if cols == cols_old: 
            instrument = 'SpeX'
        if cols == cols_new: 
            instrument = 'uSpeX'
        file.write('# \n')
        file.write('# instrument = '+instrument+'\n')
        file.write('# rawpath = {} \n'.format(raw_path))
        file.write('# calpath = {} \n'.format(cals_path))
        file.write('# procpath = {} \n'.format(proc_path))
        file.write('# \n')
    
    data_table = data_table.applymap(lambda x:str(x).center(40))
    data_table.to_csv(input_file,'|', mode= 'a', index = False, header = False)
        
    print('input file written to {}'.format(input_file))

    dp['Notes'] = ['']*len(files)

    #dp, dpc, dpk = query_reference(dp, dpc)
    #dpsl = source_list(final, dp)
    
    writer(proc_path, date, dp, dpsl, format_input, reduction)
    
    return dp

#-----------------------------------------------------------------------------#

magnitudes=['B','V','J','H','K']

if __name__ == '__main__':

    path_input = '/data/SpeX/'
    dates_input = ''
    format_input = ''
    overwrite_input = ''
    
    for argument in sys.argv:
        if argument[0:4] == 'date':
            dates_input = argument[5:].split(',')
        if argument[0:4].lower() == 'path':
            path_input = argument[5:]
        if argument[0:6].lower() == 'format':
            format_input = argument[7:]
        if argument[0:9].lower() == 'overwrite':
            overwrite_input = argument[10:]

    if dates_input == '':
        print('Please input a date or set of dates (using format 200101*) you would like to run:')
        dates_input = input()
        # Code used for a "Selectable"input
        #root = Tk()
        #root.title('Select Files')
        #root.attributes('-topmost',True)
        #root.geometry('250x100')
        #def select():
        #    root.directories = askopendirnames(initialdir=path_input ,title='Select Files')
        #    return root.directories
        #select_btn = Button(root, text='Select',command=select)
        #close_btn = Button(root, text = 'Exit', command=root.destroy)
        #select_btn.place(relx = 0.5, rely = 0.35, anchor=CENTER)
        #close_btn.place(relx = 0.5, rely = 0.65, anchor=CENTER)
        #root.mainloop()
        #input_directories = root.directories

    #else:
    input_directories = []
    for dates in dates_input:
        if dates == '**':
            input_directory = glob.glob(os.path.join(path_input, dates))
            length = len(input_directory)
            print('You are about to run {} files, are you sure you want to do this? yes/no:'.format(length))
            assuredness = input()
            if 'yes' in assuredness:
                input_directories.extend(input_directory)
                pass
            else:
                raise Exception('Process Canceled')
        else:
            input_directory = glob.glob(os.path.join(path_input, dates))
            if len(input_directory) == 0:
                pass
            else:
                input_directories.extend(input_directory)

    if len(input_directories) == 0:
        raise Exception('There are no folders with those dates in {}'.format(path_input))
    
    print(input_directories)
    for directory in input_directories:
            Cols_New_Label = False
            basefolder = str(directory)
            date = os.path.basename(basefolder)
            print(date)
            try:
                raw, cals, proc, reduction = create_folder(basefolder, date, path_input, overwrite_input)
                makelog(raw, cals, proc, date, format_input, reduction)
            except Exception as err:
                if reduction == '':
                    print(traceback.format_exc())
                else:
                    ErrorLog = open(os.path.join(reduction,'Error Log.txt'),'a+')
                    ErrorLog.write(date+'\n')
                    ErrorLog.write(traceback.format_exc())
                    ErrorLog.write('\n')
                    ErrorLog.close()
                    print(traceback.format_exc())
        
        
    

