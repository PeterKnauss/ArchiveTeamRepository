# Log Creation Repository

Repository for the Log Creation team on the Spex Archive team at UCSD Cool Star Lab.

# Requirements

To install all required libraries, run "pip install -r requirements.txt" from CMD while in the project directory.

`[User@splat ~]$ pip install -r requirements.txt`

# Inputs

There are multiple inputs that can change how the code runs or make it easier to run. The order of inputs does not matter.

All inputs follow the format of a keyword, followed by an equal sign, followed by the input.

## Date Input

To run a specific date or a specific set of dates, you must add "date=" then your input.

To run just one date, insert the date you would like to run in the 8 number, year/month/day format .

e.g. --> `[User@splat ArchiveTeamRepository]$ python Log.py date=20010101`

To run multiple dates, the input should follow the \* notation, preceeded by the month or year that you would like to run.

e.g. for all dates in 2001; `date=2001*`

e.g. for all dates in January 2001; `date= 200101*`

To run multiple sets of dates, follow all above rules but separate inputs with commas (eg. 2001*,20010101,2002*...)

## Path Input

The default path that dates are grabbed from is `/data/SpeX`, used on the guacamole virtual machine.

To set a different default path, add an input starting with "path=" then the full path to the folder where the date folders that you would like to run are.

## Format Input

The default format of the output files is as a .csv file.

To change the output to .xlsx file, add "format=excel" as an input.

## Overwrite Input

The default function is that if a cals/proc folder already exists for a date that is being run, a yes/no input will be required to overwrite the folder or not.

To automatically choose to either not overwrite or overwrite these folders, add "overwrite=no" or "overwrite=yes" respectively.

`[User@splat ArchiveTeamRepository]$ python Log.py date=20010101 overwrite=yes`
