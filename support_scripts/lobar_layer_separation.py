import pandas as pd
import nibabel as nib
import numpy as np
import os
import glob
import sys
import argparse
import getopt

dict_lobes_mapping={
    'FL':[1],
    'FR':[2],
    'PL':[3],
    'PR':[4],
    'OL':[5],
    'OR':[6],
    'TL':[7],
    'TR':[8],
    'BG':[9],
    'IT':[10],
    'F':[1,2],
    'P':[3,4],
    'O':[5,6],
    'T':[7,8],
    'BGIT':[9,10],
    'L':[1,3,5,7],
    'R':[2,4,6,8]
}

dict_layers_mapping={
    '1':[1],
    '2':[2],
    '3':[3],
    '4':[4,5],
    'Periv':[1,2],
    'Mid':[2,3],
    'Juxta':[3,4,5]
}

def process_layers(layers_data, les_data):
    new_dict = {}
    list_keys = dict_layers_mapping.keys()
    for k in list_keys:
        overlap = 0
        reg = 0
        for val in dict_layers_mapping[k]:
            overlap += np.sum(les_data[layers_data==val])
            reg += np.asarray(np.where(layers_data==val)).shape[1] 
        freq = overlap/reg
        dist = overlap/np.sum(les_data)
        new_dict['LesLayers%s'%k] = overlap
        new_dict['RegLayers%s'%k] = reg
        new_dict['FreqLayers%s'%k] = freq
        new_dict['DistLayers%s'%k] = dist
    return new_dict

def process_lobes_layers(lobes_data,layers_data,les_data):
    new_dict = {}
    list_lobes = dict_lobes_mapping.keys()
    list_layers = dict_layers_mapping.keys()
    for ko in list_lobes:
        lobes_tmp = np.zeros_like(lobes_data)
        for val in dict_lobes_mapping[ko]:
            lobes_tmp = np.where(lobes_data==val,np.ones_like(lobes_data),lobes_tmp)
        for ka in list_layers:
            layers_tmp = np.zeros_like(layers_data)
            for val in dict_layers_mapping[ka]:
                layers_tmp = np.where(layers_data==val,np.ones_like(layers_data),layers_tmp)
            mix_tmp = lobes_tmp * layers_tmp
            overlap = np.sum(les_data[mix_tmp==1])
            reg = np.sum(mix_tmp)
            freq = overlap / reg
            dist = overlap / np.sum(les_data)
            new_dict['LesLobes%s_Layers%s'%(ko,ka)] = overlap
            new_dict['RegLobes%s_Layers%s'%(ko,ka)] = reg
            new_dict['FreqLobes%s_Layers%s'%(ko,ka)] = freq
            new_dict['DistLobes%s_Layers%s'%(ko,ka)] = dist
    return new_dict


def process_lobes(lobes_data, les_data):
    new_dict = {}
    list_keys = dict_lobes_mapping.keys()
    for k in list_keys:
        overlap = 0
        reg = 0
        for val in dict_lobes_mapping[k]:
            overlap += np.sum(les_data[lobes_data==val])
            reg += np.asarray(np.where(lobes_data==val)).shape[1] 
        freq = overlap/reg
        dist = overlap/np.sum(les_data)
        new_dict['LesLobes%s'%k] = overlap
        new_dict['RegLobes%s'%k] = reg
        new_dict['FreqLobes%s'%k] = freq
        new_dict['DistLobes%s'%k] = dist
    return new_dict



def get_id_tp(lobes_name):
    basename = os.path.split(lobes_name)[1]
    id = basename.split('_')[0]
    tp = basename.split('_')[1]
    return id, tp

def main(argv):

    parser = argparse.ArgumentParser(description='Lobar and Layers separation')
    parser.add_argument('-lesion', dest='lesion', metavar='file with the database of results',
                        type=str, required=True,
                        help='File to read the data from')
    parser.add_argument('-lobes', dest='lobes', metavar='file with the database of results',
                        type=str, required=True,
                        help='File to read the data from')
    parser.add_argument('-layers', dest='layers', metavar='file with the database of results',
                        type=str, required=True,
                        help='File to read the data from')       
    parser.add_argument('-id', dest='id', default='',
                        type=str)
    parser.add_argument('-name_out', dest='name_out', type=str, help='name of the output file', default='Completed')

    

    try:
        args = parser.parse_args(argv)

    except getopt.GetoptError:
        print('lobar_layer_separation.py -lesion <ids with ref> -layers <file to complete>  -lobes <file for lobes> -path_out <path where to save completed output> -name_out ')
        sys.exit(2)

    lesion = nib.load(args.lesion).get_fdata()
    layers = nib.load(args.layers).get_fdata()
    lobes = nib.load(args.lobes).get_fdata()
    dict_lobes = process_lobes(lobes,lesion)
    print('Lobes processed')
    dict_layers = process_layers(layers, lesion)
    print('Layers processed')
    dict_layers_lobes = process_lobes_layers(lobes, layers,lesion)
    print('Layers and Lobes processed')
    dict_fin = {}
    dict_fin.update(dict_lobes)
    dict_fin.update(dict_layers)
    dict_fin.update(dict_layers_lobes)
    df_fin = pd.DataFrame.from_dict([dict_fin])
    df_fin['id'] = args.id
    df_fin['les_tot'] = np.sum(lesion)
    df_fin.to_csv(args.name_out)

if __name__ == "__main__":
    main(sys.argv[1:])

