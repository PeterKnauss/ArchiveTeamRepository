import os
import errno

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

def create_folder(date):    
    proc_directory='proc'
    cals_directory='cals'
    reduction_directory=os.path.join(os.getcwd(), 'reductions')
    parental_directory= os.path.join(reduction_directory, str(date) )
    
    #should be something like : /home/user/reductions/(the date)/
    
    #Check if 'reductions' is created. if not print a comment and ask for one
    #and then checks if the 'date' folder already exist in reductions folder
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
            print('FileExistError: File exist in raw folder')
        else:
            raise
    # Making the folders
    try:
        os.mkdir(proc_path)
    except OSError as err:
        if err.errno == errno.ENOENT:
            print('FileNotFound : Missing some folders, please check your path')
        if err.errno == errno.EEXIST:
            print('FileExistError : File exist, might overwrite !!')
        else:
            raise

    try:
        os.mkdir(cals_path)
    except OSError as err:
        if err.errno == errno.ENOENT:
            print('FileNotFound : Missing some folders, please check your path')
        if err.errno == errno.EEXIST:
            print('FileExistError : File exist, might overwrite !!')
        else:
            raise

