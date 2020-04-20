import multiprocessing as mp

def multiprocess_guard(): #return False if code is being executed from LPM's child processes
    return mp.current_process().name!='LPM_child' #built to be implemented as a more proper if __name__=='__main__' guard for Windows multiprocessing