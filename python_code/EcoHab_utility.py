from EcoHab import *
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io
import itertools

def combine_tunnel(labels_wt):
    """to combine the tunnel state to the next box"""
    labels_wot = np.array(labels_wt)
    if np.ndim(labels_wt) == 1:
        for i in np.flip(range(np.size(labels_wot))):
            if labels_wot[i] == 0:
                labels_wot[i] = labels_wot[i+1]
    else:
        for imice in range(np.shape(labels_wot)[1]):
            for i in np.flip(range(np.shape(labels_wot)[0])):
                if labels_wot[i, imice] == 0:
                    labels_wot[i, imice] = labels_wot[i+1, imice]
            
    return labels_wot


def change_address(mice_address, mouseId, boxId, in_out_flag):
    boxId = int(boxId)
    in_out_flag = int(in_out_flag)
    mouseId = int(mouseId)
    mice_address_new = np.array(mice_address, dtype=int)
    if in_out_flag > 0:
        mice_address_new[mouseId] = boxId;
    else:
        mice_address_new[mouseId] = 5;
        
    return mice_address_new


def change_config(pk, boxId, in_out_flag):
    '''
    Update P(k) given time_stamp info (we can upgrade this function to update other observables)
    '''
    boxId_ = int(boxId)
    in_out_flag = int(in_out_flag)
    pk_new = np.array(pk, dtype=int)
    pk_new[boxId] += in_out_flag
    return pk_new


# ---------
# Generate P(k) from microstates
# ---------
def config_to_pk(mice_address):
    nbox = 4
    pk = np.empty((nbox,),dtype=int)
    for boxId in np.arange(1,nbox+1):
        pk[boxId-1] = sum(mice_address==boxId)
        
    return pk