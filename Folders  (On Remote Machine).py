import os

# import shutil
# import shutil if you want to create a 'data' folder inside our specific observation 'date' folder

# This function is designed to run on the remote machine. os.getcwd will get the current directory.
# Assuming we ran this code outside any folder (ie home/user/) we want to direct ourselves to reductions folder and create
# a folder with the input date. Then create proc and cals folder in it. 
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
    except FileNotFoundError:
        print('Missing folder : "reductions". Create a "reductions" folder !!')
    except FileExistError:
        print("File/Folder already exist. Be careful, you might overwrite !!" )
        
    proc_path=os.path.join(parental_directory, proc_directory)
    cals_path=os.path.join(parental_directory, cals_directory)
    raw_path=os.path.join('/data/SpeX/', str(date))
    #Check if data exist in raw_path and if not prints a comment
    try: 
        os.mkdir(raw_path)
    except FileExistError:
        print('File exists in raw folder, Good')
    except FileNotFoundError:
        print('No such file, Check your input date or file path again ')
    #----------------------------------------------------------------------
    #This part of the code is for creating/copying raw'data' files in the 'date' folder if needed on remote machine.
    #Most likely don't need this part
    
    #data_path=os.path.join(parental_directory, 'data')
    #data=os.makedirs(data_path)
    #shutil.copyfile(raw_path, data)
    
    #----------------------------------------------------------------------
    # using mkdir() so it doesn't create folders by accident 
    try:
        os.mkdir(proc_path)
        os.mkdir(cals_path)
    except FileNotFoundError:
        print('File not Found, Something went wrong, Please check file path or previous comments for missing folders !')
    except FileExistError:
        print('File already exist, you might overwrite !!')
