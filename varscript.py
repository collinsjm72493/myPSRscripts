"""
varscript.py
Justin Collins
July 8, 2014

A minor script for creating variance images from rectified files.
"""

import numpy
import os, sys, glob
import pyfits
import re

print 'Variance image creation routine.'
print '------------------------------------------'
cont=0
while cont==0:
    prodir=raw_input('Enter the product directory as a string: ')
    if os.path.exists(prodir):
        cont=1
        #Get the read noise and gain for each product file
        main,prod=prodir.split('product')
        rawdir=(main+'raw/')
        rawlist=glob.glob(rawdir+'*.fits')
        for k in range(len(rawlist)):
            if re.search('_var',rawlist[k])!=None:
                os.remove(rawlist[k])
        rawlist=glob.glob(rawdir+'*.fits')
        infile_list=glob.glob(prodir+'xmcxgbp*.fits')
        for m in range(len(infile_list)):
            if re.search('_var',infile_list[m])!=None:
                print 'Warning: File name %s already exists.' % infile_list[m]
                con=0
                while con!=1:
                    ow=raw_input('Would you like to overwrite the file? (y/n, case sensitive): ')
                    if ow=='y':
                        print 'Overwriting file %s' % infile_list[m]
                        os.remove(infile_list[m])
                        con=1
                    elif ow=='n':
                        print 'File will not be overwritten.'
                        con=1
                    else:
                        print 'Invalid input.  Try again.'
        for m in range(len(infile_list)):
            try:
                if re.search('_var',infile_list[m])!=None:
                    infile_list.remove(infile_list[m])
            except:
                pass
        if len(rawlist)==len(infile_list):
            for n in range(len(infile_list)):
                p=pyfits.open(infile_list[n],mode='update')
                r=pyfits.open(rawlist[n])
                framenoise=[]
                for w in range(len(r)):
                    if w!=0:
                        framenoise.append(r[w].header['RDNOISE'])
                rdnoise=numpy.mean(framenoise)
                gain=1
                p['SCI'].header.set('RDNOISE',rdnoise)
                p['SCI'].header.set('GAIN',gain)
                p.flush()
        else:
            sys.exit('Error: Not all raw files have product files.')
                
        #Create variance files
        print 'List of files to build from: ',infile_list
        print '-------------------------------'
        varlist=[]
        for file in infile_list:
            f=pyfits.open(file)
            for k in range(len(f)):
                if k!=0:
                    rdnoise=f['SCI'].header['RDNOISE']
                    gain=f['SCI'].header['GAIN']
                    if (f[k].data<=0).any():
                        j=numpy.where(f[k].data>0)
                        min_pos=f[k].data[j].min()
                        j=numpy.where(f[k].data<0)
                        f[k].data[j]=min_pos
                        j=numpy.where(f[k].data==0)
                        f[k].data[j]=1000000
                    f[k].data=(f[k].data+(rdnoise/gain)**2)
            root,ext=os.path.splitext(file)
            varfile=(root+'_var'+ext)
            print 'Attempting to write variance file: ',varfile
            if os.path.exists(varfile):
                print 'Warning: File %s was not overwritten.  New file will not be created.' % varfile
            else:
                f.writeto(varfile)
                print 'File %s successfully written' % varfile
                varlist.append(varfile)
            print '----------------------------------------'
        varfiles=','.join(['%s' % n for n in varlist])
        print 'The variance files which were written are: ',varfiles
    else: #prodir path does not exist
        print 'Error: Path %s does not exist.  Try again.' % prodir
        print '------------------------------------'
