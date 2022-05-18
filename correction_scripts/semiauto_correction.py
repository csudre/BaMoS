








import pandas as pd
import nibabel as nib

import argparse
import os
import sys
import numpy as np
from scipy.ndimage.morphology import binary_dilation
from scipy.ndimage import label as cc

LIST_ARTEFACTS=[24,31,60,59,48,49,4,90, 32,33, 62,63,76,77]
LIST_EXTEND=[125,126,145,147]


def find_associated_label(point, connect):
    value = -1
    if connect[point[0],point[1],point[2]] >0:
        value =  connect[point[0],point[1],point[2]]
    else:
        print('Not direct, creating search zone')
        seg_search = create_search_zone(point,connect.shape)
        values = np.unique(seg_search*connect)
        if np.max(values>0):
            value = np.max(values)
    return value

def create_allowed_zone(point, shape, range_x=-1, range_z=-1, range_y=-1):
    allowed_x = np.ones(shape)
    allowed_y = np.ones(shape)
    allowed_z = np.ones(shape)
    if range_x >= 0:
        allowed_x = np.zeros(shape)
        allowed_x[point[0]-range_x:point[0]+range_x,:,:] = 1
    if range_y >= 0:
        allowed_y = np.zeros(shape)
        allowed_y[:,point[1]-range_y:point[1]+range_y,:] = 1
    if range_z >= 0:
        allowed_y = np.zeros(shape)
        allowed_y[:,:,point[2]-range_z:point[2]+range_z] = 1
    return allowed_x * allowed_y * allowed_z

def create_search_zone(point, shape, zone=2):
    seg_zone = np.zeros(shape)
    for x in range(-zone,zone+1):
        for y in range(-zone,zone+1):
            for z in range(-zone,zone+1):
                seg_zone[point[0]+x,point[1]+y,point[2]+z] = 1
    return seg_zone

def associated_seg(lesion, connect, value):
    seg_connect = np.where(connect==value, np.ones_like(connect),
                           np.zeros_like(connect))
    lesion_new = lesion * seg_connect
    return lesion_new

def correct_lesion(lesion_tocorr, correction, type_corr):
    print('Correction volume is ', np.sum(correction))
    print('Lesion to corr volume is', np.sum(lesion_tocorr))
    if type_corr == 'add' or type_corr=='de novo':
        #lesion_corr = lesion_tocorr + correction
        lesion_corr = np.where(lesion_tocorr > 0.4, lesion_tocorr, correction)
	return lesion_corr
    if type_corr == 'remove':
        lesion_corr = lesion_tocorr - correction
        return lesion_corr
    if type_corr == 'extend':
        lesion_corr = np.where(lesion_tocorr>0.4,lesion_tocorr,correction)
        return lesion_corr
    
    print('Unknown correction type ', type_corr)
    return lesion_tocorr

def process_full_case(lesion_tocorr, lesion_init, connect, df_point
                      ):
    lesion_corr = lesion_tocorr
    for i in range(df_point.shape[0]):
        case = df_point.iloc[i]
        point = [case['x'], case['y'], case['z']]
        type_corr = case['type']
        label = find_associated_label(point, connect)
        print('Associated label is %d' % int(label))
        if label > 0:
            correction = associated_seg(lesion_init, connect, label)

            lesion_corr = correct_lesion(lesion_corr, correction, type_corr)
	    print(np.sum(lesion_corr))
    return lesion_corr

def associated_mahalseg(mahal, connect, value, thresh=3.0):
    seg_connect = np.where(connect == value, np.ones_like(connect),
                           np.zeros_like(connect))
    mahal_new = mahal[...,-1]
    mahal_thresh = mahal_new / thresh
    lesion_new = mahal_thresh * seg_connect
    lesion_new = np.where(lesion_new<1, lesion_new, np.ones_like(lesion_new))
    lesion_new = np.where(lesion_new>=0.5, lesion_new, np.zeros_like(
        lesion_new))
    return lesion_new

def associated_mahalnovo(mahal, point, thresh):
    seg_connect = create_search_zone(point, mahal.shape[0:3], 2)
    mahal_new = mahal[..., -1]
    mahal_thresh = mahal_new / thresh
    lesion_new = mahal_thresh * seg_connect
    lesion_new = np.where(lesion_new < 1, lesion_new, np.ones_like(lesion_new))
    lesion_new = np.where(lesion_new >= 0.5, lesion_new, np.zeros_like(
        lesion_new))
    return lesion_new

def process_full_mahal(lesion_tocorr, mahal, connect, parcellation, df_point,
                       region=10,
                       thresh=3 ):
    lesion_corr = lesion_tocorr
    artefact_zone = np.where(parcellation>98, np.ones_like(
        parcellation), np.zeros_like(parcellation))
    artefact_extend = np.zeros_like(parcellation)

    for f in LIST_ARTEFACTS:
        artefact_zone = np.where(parcellation==f, np.ones_like(parcellation),
                                 artefact_zone)
    for f in LIST_EXTEND:
        artefact_extend = np.where(parcellation==f, np.ones_like(
            parcellation), artefact_extend)
    artefact_extend = binary_dilation(artefact_extend,iterations=3)

    artefact_zone += artefact_extend
    for i in range(df_point.shape[0]):
	print(df_point.iloc[i])
        range_x = -1
        range_y = -1
        range_z = -1
        case = df_point.iloc[i]
        point = [case['x'], case['y'], case['z']]
        if 'range_x' in df_point.columns:
            range_x = case['range_x']
        if 'range_y' in df_point.columns:
            range_y = case['range_y']
        if 'range_z' in df_point.columns:
            range_z = case['range_z']
        allowed_zone = create_allowed_zone(point, shape=parcellation.shape,
                                           range_x=range_x, range_y=range_y,
                                           range_z=range_z)
        type_corr = case['type']
        label = find_associated_label(point, connect)
        print('Associated label is %d' % int(label))
        if label > 0:
            correction = associated_mahalseg(mahal, connect, label,
                                             thresh)
            correction *= allowed_zone
            if case['artefact']:
                correction = np.where(artefact_zone>=1, np.zeros_like(
                    correction), correction)
                cc_corr, numb_cc = cc(correction)
                if numb_cc > 1:
                    cc_value = find_associated_label(point, cc_corr)
                    correction = np.where(cc_corr==cc_value,correction,
                                          np.zeros_like(parcellation))

            lesion_corr = correct_lesion(lesion_corr, correction, type_corr)
        elif label==-1 and type_corr=='de novo':
            print('Creating de novo')
            correction = associated_mahalnovo(mahal, point, thresh)
            correction *= allowed_zone
            if case['artefact']:
                correction = np.where(artefact_zone >= 1, np.zeros_like(
                    correction), correction)
                cc_corr, numb_cc = cc(correction)
                if numb_cc > 1:
                    cc_value = find_associated_label(point, cc_corr)
                    correction = np.where(cc_corr == cc_value, correction,
                                          np.zeros_like(parcellation))
            lesion_corr = np.where(lesion_corr > 0.4, lesion_corr,correction)
	    
        else:
	    print('impossible to find anything to correct')
	print(np.sum(lesion_corr))
	print(np.sum(lesion_corr)-np.sum(lesion_tocorr))
    lesion_corr = np.where(lesion_corr < 0, np.zeros_like(lesion_corr), lesion_corr)
    
    return lesion_corr


def main(argv):

    parser = argparse.ArgumentParser(description='Correct specific lesion '
                                                 'choices')
    parser.add_argument('-csv', dest='csv', metavar='csv with location of '
                                                    'specific points 4 '
                                                    'columns x y z and')
    parser.add_argument('-connect', dest='connect', metavar='seg pattern',
                        type=str, required=False,
                        help='RegExp pattern for the segmentation files')
    parser.add_argument('-init', dest='init', action='store', type=str,
                        help='RegExp pattern for the reference files')
    parser.add_argument('-tocorr', dest='tocorr', action='store',
                        default='TrainingLongData_3004.csv')
    parser.add_argument('-save_name', dest='save_name', action='store',
                        default='', help='name to save results')
    parser.add_argument('-mahal', dest='mahal', action='store',
                        default='TrainingLongData_3004.csv',type=str)
    parser.add_argument('-parc', dest='parc', action='store',
                        default='TrainingLongData_3004.csv', type=str)
    parser.add_argument('-init_connect', dest='init_connect', action='store',
                        type=str)
    parser.add_argument('-region', dest='region', action='store',type=int,
                        default=10)
    parser.add_argument('-thresh', dest='thresh', action='store', type=float,
                        default=3)
    parser.add_argument('-range_z', dest='range_z',action='store',type=float,
                        default=-1)
    parser.add_argument('-range_x', dest='range_x', action='store', type=float,
                        default=-1)
    parser.add_argument('-range_y', dest='range_y', action='store', type=float,
                        default=-1)


    try:
        args = parser.parse_args()
        # print(args.accumulate(args.integers))
    except argparse.ArgumentTypeError:
        print('semiauto_correction.py -csv <csv_file> -init '
              '<initial lesion segmentation> -connect <initial connected '
              'component> -tocorr '
              '<file to correct> '
              '-save_name <name for saving>   ')
    df_corr = pd.read_csv(args.csv)
    lesion_tocorr = nib.load(args.tocorr).get_data()
    if args.connect is not None:
        init = nib.load(args.init).get_data()
        connect = nib.load(args.connect).get_data()
        lesion_corr = process_full_case(lesion_tocorr,init,connect,df_corr)
    else:
        mahal = nib.load(args.mahal).get_data()
        parcellation = nib.load(args.parc).get_data()
        connect = nib.load(args.init_connect).get_data()
        lesion_corr = process_full_mahal(lesion_tocorr, mahal, connect,
                                         parcellation,
                                         df_corr, args.thresh)
    nii_base = nib.load(args.tocorr)
    new_name = args.tocorr.lstrip('.nii.gz') +'_'+args.save_name +'.nii.gz'
    new_nii = nib.Nifti1Image(lesion_corr, nii_base.affine)
    nib.save(new_nii, new_name)

if __name__ == "__main__":
     main(sys.argv[1:])
