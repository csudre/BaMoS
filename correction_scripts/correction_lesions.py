import nibabel as nib
import pandas as pd
from Parsing_LesionFile_args import parse
import os
import argparse
import numpy as np
import sys
from scipy.ndimage.morphology import binary_dilation
from scipy.ndimage.morphology import binary_erosion
from scipy.ndimage import label


def create_ventricles(parcellation):
    ventr = np.where(np.logical_or(parcellation==52, parcellation==53),
                     np.ones_like(parcellation), np.zeros_like(parcellation))
    return ventr

def create_border_thalamusventr(parcellation):
    ventr = np.where(np.logical_or(parcellation == 52, parcellation == 53),
                     np.ones_like(parcellation), np.zeros_like(parcellation))
    thal = np.where(np.logical_or(parcellation == 60, parcellation == 61),
                     np.ones_like(parcellation), np.zeros_like(parcellation))
    border_thal = binary_dilation(ventr, iterations=1) * thal
    return border_thal


def slicewise_choroid(lesion, ventricles, border_thal):
    intersect = lesion * (ventricles+border_thal)

    indices = np.asarray(np.where(intersect>0)).T
    slices = np.unique(indices[:,2])
    choroid = np.zeros_like(lesion)
    choroid_temp = np.zeros_like(lesion)
    for sl in slices:
        full_slice = intersect[..., sl]
        bin_slice = np.asarray(full_slice, dtype=bool)
        labelled, numb_labels = label(bin_slice)
        for lab in range(1, numb_labels+1):
            lab_temp = np.where(labelled == lab, np.ones_like(labelled),
                                np.zeros_like(labelled))
            dil_bord = binary_dilation(lab_temp, iterations=1) - lab_temp
            ventrinter = dil_bord * ventricles[..., sl]
            ratio = np.sum(ventrinter)*1.0/np.sum(dil_bord)
            # print(ratio, sl, lab, np.sum(lab_temp))
            if ratio > 0.75:
                choroid[..., sl] = np.where(lab_temp, lesion[...,sl],
                                            choroid[...,sl])
            if ratio >= 0.25:
                choroid_temp[..., sl] = np.where(lab_temp, lesion[..., sl],
                                            choroid_temp[..., sl])
    aggreg_choroid = np.squeeze(np.sum(choroid, axis=2))
    ind_aggreg = np.where(aggreg_choroid>0)
    bin_aggreg = np.where(aggreg_choroid>0, np.ones_like(aggreg_choroid),
                          np.zeros_like(aggreg_choroid))
    dilated_aggreg = binary_dilation(bin_aggreg, iterations=1)
    # print(np.sum(choroid), np.sum(choroid_temp))
    aggreg_3d = np.zeros_like(choroid_temp)
    for sl in slices:
        aggreg_3d[..., sl] = aggreg_choroid
        if np.sum(choroid_temp[..., sl]) > 0:
            bin_slice = np.where(choroid_temp[..., sl]>0, np.ones_like(
                aggreg_choroid), np.zeros_like(aggreg_choroid))
            labelled, numb_labels = label(bin_slice)
            for lab in range(1, numb_labels+1):
                temp_lab = np.where(labelled == lab, np.ones_like(bin_slice),
                                    np.zeros_like(bin_slice))

                if np.sum(dilated_aggreg * temp_lab) > 0:
                    choroid[...,sl] = np.where(temp_lab>0, choroid_temp[...,
                                                                        sl],
                                               choroid[...,sl])
    print('Choroid correction is' ,np.sum(choroid))
    return choroid, aggreg_3d

def ventr_corr(lesion, parcellation, dil=1):
    ventr = create_ventricles(parcellation)
    #struct = np.zeros((3,3,3))

    dil_ventr = binary_dilation(ventr,iterations=dil)
    corr_ventr  = lesion*dil_ventr
    print("Ventr correction is ", np.sum(corr_ventr))
    return corr_ventr

def cortical_corr(file, lesion, connected):
    df_parse = pd.DataFrame.from_dict(parse(file, 2))
    to_remove = np.zeros_like(lesion)
    df_toremove = df_parse[df_parse['ot']=='16']
    for lab in df_toremove['label']:
        label = int(lab)
        to_remove = np.where(connected == label, lesion, to_remove)
    print("Cortical 16 volume is ", np.sum(to_remove))
    return to_remove


def hypot1_corr(file, lesion, connected):
    df_parse = pd.DataFrame.from_dict(parse(file,2))
    to_remove = np.zeros_like(lesion)
    df_toremove = df_parse[df_parse['mahalmean1'].astype(float)<-5.5]
    for lab in df_toremove['label']:
        label = int(lab)
        to_remove = np.where(connected==label, lesion, to_remove)
    print("To remove because of T1 hypo is ", np.sum(to_remove))
    return to_remove


def cortical_sheet(lesion, parcellation, thresh_vol=10, ratio_thresh=0.55):
    shape = lesion.shape
    gm = (parcellation>98)
    correction_axial = np.zeros_like(lesion)
    structure = np.ones([3,3])
    for s in range(0, shape[2]):

        lesion_temp = lesion[:,:,s]
        correct_temp = np.zeros_like(lesion_temp)
        gm_temp = gm[:,:,s]
        bin_temp = (lesion_temp > 0.5)
        if np.sum(bin_temp) > thresh_vol:
            conn_temp, num_lab = label(bin_temp, structure=structure)
            for val in range(1, num_lab+1):
                seg_temp = (conn_temp==val).astype(float)
                value_temp = np.where(seg_temp, lesion_temp, np.zeros_like(
                    lesion_temp))
                if np.sum(seg_temp) > thresh_vol/2:
                    border_ext = binary_dilation(seg_temp) - seg_temp
                    gm_border = border_ext * gm_temp
                    ratio_gmborder = np.sum(gm_border) / np.sum(seg_temp)
                    if ratio_gmborder > ratio_thresh:
                        correct_temp += value_temp
            correction_axial[:,:,s] = correct_temp
    correction_coronal = np.zeros_like(lesion)
    for s in range(0, shape[1]):
        lesion_temp = lesion[:, s, :]
        correct_temp = np.zeros_like(lesion_temp)
        gm_temp = gm[:, s, :]
        bin_temp = (lesion_temp > 0.5)
        if np.sum(bin_temp) > thresh_vol:
            conn_temp, num_lab = label(bin_temp, structure=structure)
            for val in range(1, num_lab + 1):
                seg_temp = (conn_temp == val).astype(float)
                value_temp = np.where(seg_temp, lesion_temp, np.zeros_like(
                    lesion_temp))
                if np.sum(seg_temp) > thresh_vol/2:
                    border_ext = binary_dilation(seg_temp) - seg_temp
                    gm_border = border_ext * gm_temp
                    ratio_gmborder = np.sum(gm_border) / np.sum(seg_temp)
                    if ratio_gmborder > ratio_thresh:
                        correct_temp += value_temp
            correction_coronal[:, s, :] = correct_temp

    correction_sagittal = np.zeros_like(lesion)
    for s in range(0, shape[0]):
        lesion_temp = lesion[s, :, :]
        correct_temp = np.zeros_like(lesion_temp)
        gm_temp = gm[s, :, :]
        bin_temp = (lesion_temp > 0.5)
        if np.sum(bin_temp) > thresh_vol:
            conn_temp, num_lab = label(bin_temp, structure=structure)
            for val in range(1, num_lab + 1):
                seg_temp = (conn_temp == val).astype(float)
                value_temp = np.where(seg_temp, lesion_temp, np.zeros_like(
                    lesion_temp))
                if np.sum(seg_temp) > thresh_vol/2:
                    border_ext = binary_dilation(seg_temp) - seg_temp
                    gm_border = border_ext * gm_temp
                    ratio_gmborder = np.sum(gm_border) / np.sum(seg_temp)
                    if ratio_gmborder > ratio_thresh:
                        correct_temp += value_temp
            correction_sagittal[s, :, :] = correct_temp
    correction_fin = correction_coronal + correction_axial + correction_sagittal
    print("Sagittal, Axial, Coronal correction:%f %f %f" % (np.sum(
        correction_sagittal), np.sum(correction_axial), np.sum(
        correction_coronal)))
    correction_fin = np.minimum(lesion, correction_fin)
    return correction_fin





def main(argv):

    parser = argparse.ArgumentParser(description='Correction of lesions')
    parser.add_argument('-les', dest='lesion', metavar='lesion file',
                        type=str, required=True,
                        help='Path to lesion file')
    parser.add_argument('-parc', dest='parc', action='store',
                        default='', type=str,
                        help='path to parcellation file')
    parser.add_argument('-connect', dest='connect', action='store',
                        default='',
                        type=str)
    parser.add_argument('-dil', dest='dil_ventr', action='store', type=int,
                       default=1, help='Number of dilation iterations for '
                                       'ventricle correction')
    parser.add_argument('-label', dest='label', action='store', type=str,
                        default='')
    parser.add_argument('-id', dest='id', action='store',
                        default='', help='name to save results')
    parser.add_argument('-corr',nargs='+', dest='corr', action='store',type=str,
                        default='all',choices=['choroid', 'cortex',
                                                'sheet', 'choroid-cortex',
                                                'choroid-sheet',
                                                'cortical-sheet',
                                                'all','hypot1','ventr'])

    try:
        args = parser.parse_args()
        # print(args.accumulate(args.integers))
    except argparse.ArgumentTypeError:
        print('compare_segmentation.py -s <segmentation_pattern> -r '
              '<reference_pattern> -t <threshold> -a <analysis_type> '
              '-save_name <name for saving> -save_maps  ')

        sys.exit(2)


    lesion_nii = nib.load(args.lesion)
    lesion_data  = lesion_nii.get_data()
    print('Lesion range before corr ',np.min(lesion_data), np.max(lesion_data))
    path_results = os.path.split(args.lesion)[0]
    if args.connect != '':
        connect_nii = nib.load(args.connect)
    choroid = None
    cortical = None
    sheet = None
    ventr_correction = None
    hypo_corr = None
    if args.corr == 'choroid' or args.corr == 'all' or args.corr \
            =='choroid-sheet' or args.corr == 'choroid-cortex':
        parc_nii = nib.load(args.parc)
        ventricles = create_ventricles(parc_nii.get_data())
        border_thal = create_border_thalamusventr(parc_nii.get_data())
        choroid, aggreg_3d = slicewise_choroid(lesion_data, ventricles,
                                               border_thal)
        nii_choroid = nib.Nifti1Image(choroid, lesion_nii.affine)
        nib.save(nii_choroid, os.path.join(path_results,
                                           'ChoroidCorr_'+args.id+'.nii.gz'))
        nii_agg = nib.Nifti1Image(aggreg_3d, lesion_nii.affine)
        nib.save(nii_agg, os.path.join(path_results,
                                           'AggregChor_' + args.id + '.nii.gz'))

    if 'all' in args.corr or 'cortex' in args.corr or 'choroid-cortex' in \
            args.corr or 'cortical-sheet' in args.corr:
        cortical = cortical_corr(args.label, lesion_data,
                       connect_nii.get_data())
        nii_cortical = nib.Nifti1Image(cortical, lesion_nii.affine)
        nib.save(nii_cortical, os.path.join(path_results,
                                           'CorticalCorr_' + args.id +
                                            '.nii.gz'))

    if 'all' in args.corr or 'sheet' in args.corr or 'choroid-sheet' in \
                args.corr or 'cortical-sheet' in args.corr:
        parc_nii = nib.load(args.parc)
        parc_data = parc_nii.get_data()
        sheet = cortical_sheet(lesion_data, parc_data)
        nii_cortical = nib.Nifti1Image(cortical, lesion_nii.affine)
        nib.save(nii_cortical, os.path.join(path_results,
                                            'SheetCorr_' + args.id +
                                            '.nii.gz'))
    if 'hypot1' in args.corr:
        hypo_corr = hypot1_corr(args.label, lesion_data, connect_nii.get_data())
        nii_hypot1 = nib.Nifti1Image(hypo_corr, lesion_nii.affine)
        nib.save(nii_hypot1, os.path.join(path_results,
                                        'HypoCorr' + args.id +
                                        '.nii.gz'))
    if 'ventr' in args.corr:
        parc_data = nib.load(args.parc).get_data()
        ventr_correction = ventr_corr(lesion_data, parc_data,
                                      dil=args.dil_ventr)
        nii_ventr = nib.Nifti1Image(ventr_correction, lesion_nii.affine)
        nib.save(nii_ventr, os.path.join(path_results,
                                          'VentrCorr' + args.id +
                                          '.nii.gz'))

    if choroid is not None:
        lesion_data -= choroid
        print('Lesion range after choroid is',np.min(lesion_data), np.max(lesion_data))
    if cortical is not None:
        lesion_data -= cortical
        print('Lesion range after cortical is', np.min(lesion_data), np.max(lesion_data))

    if sheet is not None:
        lesion_data -= sheet
        print('Lesion range after sheet is', np.min(lesion_data), np.max(lesion_data))

    if hypo_corr is not None:
        lesion_data -= hypo_corr
        print('Lesion range after hypo corr is', np.min(lesion_data), np.max(lesion_data))

    if ventr_correction is not None:
        lesion_data -= ventr_correction
        print('Lesion range after ventr is', np.min(lesion_data), np.max(lesion_data))

    nii_lesion = nib.Nifti1Image(lesion_data, lesion_nii.affine)
    nib.save(nii_lesion, os.path.join(path_results,
                                      'CorrectLesion_'+args.id+'.nii.gz'))
if __name__ == '__main__':
    main(sys.argv[1:])







