# -*- coding: utf-8 -*-

"""
Created on Fri Feb 18 16:18:26 2022
@author: peter

"""

import pandas as pd
import numpy as np
import sys
try:
    import tkinter as tk
    from tkinter import filedialog as fd
except Exception:
    pass


# #####################################################################################################################

# Select and open up the file

try:
    tk.Tk().withdraw()
    file = fd.askopenfilename()
except Exception:
    file = sys.argv[1] if len(sys.argv) == 2 else 'default.csv'

df = None
try:
    # XLSX using the default Excel library
    df = pd.read_excel(file)
except Exception:
    # CSV
    df = pd.read_csv(file)

# If we couldn't load anything, fail
if df is None:
    print('*** Could not read "%s"; exiting' % file)
    sys.exit(-1)

# #####################################################################################################################

# Subroutines

# Add a batch of data to the final results
def add_batch(source, batch, final, prefix, airmass, ra, dec):

    # If we haven't see this Source before, create a place for it
    if not final.get(source):
        final[source] = {'types': {'calibration': [], 'calibrator': [], 'target': []}}

    # Add the data to the final
    for type in ['calibration', 'calibrator', 'target']:
        if batch[type]:
            final[source]['types'][type].append({'start': batch[type][0], 'end': batch[type][-1], 'airmass': airmass})
    final[source]['prefix'] = prefix
    final[source]['ra'] = ra
    final[source]['dec'] = dec


# Create a set of final results
def create_final(mode):

    # Set up some variables
    final = {}
    batch = {'calibration': [], 'calibrator': [], 'target': []}
    prefix_old = None
    source_old = None
    airmass_old = None
    ra_old = None
    dec_old = None

    # Filter data by Mode
    data = []
    for index in np.arange(0, len(df.index)):
        if mode == df.iloc[index]['Mode']:
            data.append(df.iloc[index])

    # If the filter left us with nothing, give up
    if not data:
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
            source = 'flatlamp'
            prefix = 'flat/arc'
        ra = row['RA'][0:-6]
        dec = row['Dec'][0:-5]
        integration = row['Integration']
        airmass = row['Airmass']

        # If this iteration of the loop has a new Source, process the values we have been saving
        if source_old is not None and source.lower() != source_old.lower():

            # Add this batch to the final result
            add_batch(source_old.lower(), batch, final, prefix_old, airmass_old, ra_old, dec_old)

            # Start a new batch
            batch = {'calibration': [], 'calibrator': [], 'target': []}

        # Remember the last Source we saw, so we can tell if it changes with the next line
        prefix_old = prefix
        source_old = source
        ra_old = ra
        dec_old = dec
        airmass_old = airmass

        # The actual smarts -- figure out what type the record is based on the Integration
        if integration <= 1.8:
            type = 'calibration'
        elif integration <= 70.0:
            type = 'calibrator'
        else:
            type = 'target'

        # Add the Number to this batch as the particular Type
        batch[type].append(number)

    # Add the last batch to the final
    add_batch(source.lower(), batch, final, prefix, airmass, ra, dec)

    return final


# Print the final results
def print_final(final, mode):

    if final is None:
        print('  There is no data for Mode "%s"' % (mode))
        return

    # Loop through all the sources
    for source in final.keys():

        # Display the current Source and Coordinates
        print('  %s (prefix: %s):' % (source, final[source]['prefix']))
        print('%15s: %s, %s' % ('coordinates', final[source]['ra'], final[source]['dec']))

        # Loop through all the types for this source
        for type in final[source]['types'].keys():

            # Loop through all the ranges for this type
            for index, range in enumerate(final[source]['types'][type]):

                # Print type type and range
                if range['start'] == range['end']:
                    print('%15s: %s' % (type, range['start']), end='')
                else:
                    print('%15s: %s - %s' % (type, range['start'], range['end']), end='')

                print(' (%s: %s)' % ('airmass', range['airmass']))
        print()

# #####################################################################################################################

# Main code

# Run through each mode
for mode in ['LowRes15', 'ShortXD']:
    title = '%s Data Sets' % (mode)
    print('%s\n%s' % (title, '=' * len(title)))
   
    data = create_final(mode)

    print_final(data, mode)
