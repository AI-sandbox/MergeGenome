################################################################################
# Tracks the results of a script in a .txt file
# Similar to print() but saving the results in a .txt file
################################################################################

import os

def track(message, track_path):
    '''
    Objective:
        - Write tracking message in .txt file stored at track_path. If the .txt file does not exist, it is created.
    Input:
        - message: message to write at the end of .txt file.
        - track_path: path to output .txt file that will store the tracking of the script.
    '''
       
    ## If it does not already exist, create directory that will contain the .txt with the tracking of the script
    try:
        ## First, remove the name of the .txt file to check if the directory exists
        last = [pos for pos, char in enumerate(track_path) if char == '/']
        os.makedirs(track_path[:last[-1]+1])
    except OSError:
        pass

    ## Open the file in append & read mode ('a+')
    with open(track_path, "a+") as file_object:
        ## Move read cursor to the start of the file
        file_object.seek(0)
        data = file_object.read()
        if len(data) > 0 :
            ## If the file is not empty...
            # Append new line '\n'
            file_object.write("\n")
        
        ## Append text message at the end of the file
        file_object.write(message)
