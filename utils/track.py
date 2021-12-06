################################################################################
# Tracking in .txt file
################################################################################

import os

def track(message, track_name):
    '''
    Objective:
        - Write results and tracking of the script in .txt file.
    Input:
        - message: message to write in .txt file.
        - track_name: name of .txt file to write results.
    '''
    
    ## Define path to .txt file
    base_path = '/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/output'
    file = '{}/{}'.format(base_path, track_name)
    
    ## If it does not already exist, create directory that will contain the results and tracking of the script
    try:
        os.makedirs(base_path)
    except OSError:
        pass

    ## Open the file in append & read mode ('a+')
    with open(file, "a+") as file_object:
        ## Move read cursor to the start of the file
        file_object.seek(0)
        
        ## If the file is not empty, then append new line '\n'
        data = file_object.read()
        if len(data) > 0 :
            file_object.write("\n")
        
        ## Append text at the end of file
        file_object.write(message)