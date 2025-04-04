import nibabel as nib
import pandas as pd
from scipy.ndimage import label
import os
import glob
import sys
import argparse
import numpy as np

struct1=np.zeros([3,3])
struct1[1,1,1] = 1
struct1[0,1,1] = 1
struct1[2,1,1]=1
struct1[1,0,1]=1
struct1[1,2,1] =1
struct1[1,1,0]=1
struct1[1,1,2]=1

struct2=np.ones([3,3])
struct2[0,0,0] = 0
struct2[2,2,2] = 0
struct2[0,0,2] = 0
struct2[0,2,0] = 0
struct2[2,0,0] = 0
struct2[2,2,0] = 0
struct2[2,0,2] = 0
struct2[0,2,2] = 0

struct3 = np.ones([3,3])

def count_lesions(file):
    data = nib.load(file).get_fdata()
    nlab1, labels1 = label(data, struct1)
    nlab2, labels2 = label(data, struct2)
    nlab3, labels3 = label(data, struct3)

    return nlab1, nlab2, nlab3


def main(argv):

    parser = argparse.ArgumentParser(description='Assess number of lesions in file across files')
    parser.add_argument('-out', dest='output', metavar='Output file')
    parser.add_argument('-reg_exp', dest='reg_exp', metavar='seg pattern',
                        type=str, required=True,
                        help='RegExp pattern for the segmentation files')
 


    try:
        args = parser.parse_args()
        # print(args.accumulate(args.integers))
    except argparse.ArgumentTypeError:
        print('count_lesions.py -reg_exp XXX -out XXX ')

    list_lesions = glob.glob(args.reg_exp)
    print(len(list_lesions),' to process for count')
    list_dict = []
    for k in list_lesions:
        dict_new = {}
        dict_new['file'] = k
        n1, n2, n3 = count_lesions(k)
        dict_new['nles_1'] = n1
        dict_new['nles_2'] = n2
        dict_new['nles_3'] = n3
        list_dict.append(dict_new)
    
    df_count = pd.DataFrame.from_dict(list_dict)
    df_count.to_csv(args.output)
    

if __name__ == "__main__":
     main(sys.argv[1:])

