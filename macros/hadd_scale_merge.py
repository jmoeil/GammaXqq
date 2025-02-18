import ROOT
import os
import sys
import argparse
import numpy as np

from xs_mc import xs_gjets

def main():
    parser = argparse.ArgumentParser(
        description='''Scale all histos in a file to integrated lumi
        ''',
        usage='use "%(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    

    parser.add_argument("-d", "--directories", dest="directories", help="Input directories", nargs='+', type=str, default='')
    parser.add_argument("--norm", dest="norm", help="Name of the histo used to find the number of events", type=str, default='h_nvtx')
    parser.add_argument("-l", "--lumi", dest="lumi", help="Integrated lumi (in /pb)", type=str, default='1.0')
    parser.add_argument("-o", "--output", dest="output", help="Ouput file", type=str, default='')

    args = parser.parse_args()
    
    str_listdir = ''
    for i in args.directories:
        #Merge all root files in directory
        os.system('cd '+i+'; hadd all.root *.root; cd -')
        xs = 0
        for key in xs_gjets.keys():
            if i.find(key)>=0:
                xs =  xs_gjets[key]
        #Scale them
        print('python3 /user/jmoeil/GammaXqq/macros/scaletointegratedlumi.py -i all.root --norm '+ args.norm + ' -l ' + args.lumi + ' --xs ' + format(xs))
        os.system('cd '+i+'; python3 /user/jmoeil/GammaXqq/macros/scaletointegratedlumi.py -i all.root --norm '+ args.norm + ' -l ' + args.lumi + ' --xs ' + format(xs))
        str_listdir = str_listdir+i+"/all_rescaled.root "
    
    #Merge all files together
    print("hadd "+ args.output +" "+ str_listdir)
    os.system("hadd "+ args.output + " "+str_listdir )
    

if __name__ == '__main__':
    main()




    
