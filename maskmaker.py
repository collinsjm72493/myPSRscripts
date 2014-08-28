"""
maskmaker.py
Justin Collins
July 8, 2014

A simple routine for creating masked arrays of rectified files.
"""

import os,sys,glob
import pyfits
import numpy
import numpy.ma
import matplotlib
import re
import pylab
from scipy import sqrt
import time

print 'Masked array creation and sky subtraction routine.'
print '------------------------------------------'

pylab.ion()
pylab.close('all')
cont=0
while cont!=1:
    workdir=raw_input('Enter the directory that you would like to work in: ')
    if os.path.exists(workdir):
        cont=1
    else:
        print 'Error: Directory %s does not exist.  Try again.' % workdir
        print '-------------------------------------------'
con=0
while con!=1:
    img=raw_input('Enter the file that you would like to create a masked array for: ')
    img=(workdir+img)
    root,ext=os.path.splitext(img)
    if os.path.exists(img):
        con=1
        f=pyfits.open(img)
        #Look for corresponding variance image
        varfile=(root+'_var'+ext)
        if os.path.exists(varfile):
            print 'Pulling variance information from file %s.' % varfile
            v=pyfits.open(varfile)
        else:
            print 'Error: File %s does not have a corresponding variance file.' % img
            sys.exit('No variance file')
        for k in range(len(f)):
            if k!=0:
                #Check whether input file is a variance file
                if re.search('_var',img)!=None or re.search('_sky',img)!=None or re.search('_skysub',img)!=None or re.search('xmcxgbp',img)==None:
                    sys.exit('Sky fit to only be applied to wavelength-calibrated files.')
                else:
                    print 'Working image is not a variance, sky, or sky subtracted image.'

                #Identify the size of the incoming file and mask the zero counts
                x=numpy.arange(f[k].data.shape[1])
                y=numpy.arange(f[k].data.shape[0])
                X,Y=numpy.meshgrid(x,y)
                maskzero=numpy.ma.array(f[k].data,mask=f[k].data==0)
                pylab.figure(); pylab.imshow(maskzero); pylab.suptitle('Zeros Masked')
                oned=maskzero.mean(axis=1)

                #Find the primary object
                primval=oned[475:565].max()
                primloc=numpy.where(oned==primval)[0]
                print 'Primary object is at row %s.' % int(primloc[0])

                #Find the secondary object
                ave1=oned.mean()
                secvalone=oned.data[335:380].max()
                secone=numpy.where(oned.data==secvalone)[0]
                if len(secone)!=1:
                    secvalone=0
                secvaltwo=oned.data[660:730].max()
                sectwo=numpy.where(oned.data==secvaltwo)[0]
                if len(sectwo)!=1:
                    secvaltwo=0
                try:
                    int(secvalone)
                except:
                    secvalone=0
                try:
                    int(secvaltwo)
                except:
                    secvaltwo=0
                secval=max(secvalone,secvaltwo)
                if secvalone>secvaltwo:
                    secloc=secone
                    trimlow=int(secloc-85)
                    trimhigh=int(primloc+100)
                else:
                    secloc=sectwo
                    trimlow=int(primloc-95)
                    trimhigh=int(secloc+75)
                print 'Secondary object is at row %s.' % int(secloc)
                pylab.figure(); pylab.plot(y,oned,primloc,primval/2,'g*',secloc,secval/2,'r*'); pylab.suptitle('One-Dimensional Compression of Masked Zeroes')

                #Mask the object lines
                mask=numpy.abs(Y-primloc[0])<=10
                mask=(mask) | (numpy.abs(Y-secloc)<=10)
                maskobj=numpy.ma.array(maskzero.data,mask=(mask) | (maskzero.mask))
                pylab.figure(); pylab.imshow(maskobj); pylab.suptitle('Objects Masked')
                onedtwo=maskobj.mean(axis=1)
                pylab.figure(); pylab.plot(onedtwo); pylab.suptitle('One-Dimensional Compression of Masked Objects')
                

                #Identify other object(s)
                nz=numpy.where(onedtwo!=0)[0]
                sdo=numpy.std(onedtwo[nz])
                lineloc=[]
                for i in range(len(onedtwo)):
                    m1=i-50
                    if m1<0:
                        m1=0
                    m2=i+25
                    if m2>len(onedtwo):
                        m2=len(onedtwo)
                    t1=i-150
                    if t1<0:
                        t1=0
                    t2=i+150
                    if t2>len(onedtwo):
                        t2=len(onedtwo)
                    sd=(numpy.std(onedtwo[m1:m2]))
                    if abs(onedtwo[i]-onedtwo[m1:m2].mean())>=(2.25*sd) or abs(onedtwo[i]-onedtwo[t1:t2].mean())>=(1.75*sdo):
                        lineloc.append(i)
                linepts=[]
                grp=5
                for n in range(len(lineloc)):
                    l1=lineloc[n]-grp
                    l2=lineloc[n]+grp
                    mn=onedtwo[l1:l2].mean()
                    em=onedtwo[l1:l2].max()
                    ap=onedtwo[l1:l2].min()
                    if abs(em-mn)>abs(ap-mn):
                        linepts.append(em)
                    else:
                        linepts.append(ap)
                linept=[]
                for i in linepts:
                    if i not in linept:
                        linept.append(i)
                linelocs=[]
                for b in range(len(linept)):
                    line=numpy.where(onedtwo.data==linept[b])
                    line=int(line[0])
                    linelocs.append(line)
                print 'Other object lines are at rows: %s' % linelocs
                pylab.figure(); pylab.plot(onedtwo); pylab.suptitle('Non-Major Object Locations')
                for i in range(len(linelocs)):
                    pylab.plot(linelocs[i],onedtwo.mean(),'r*')
                    pylab.vlines(linelocs[i],onedtwo.min(),onedtwo.max())

                #Manual Check
                mcon=0
                while mcon!=1:
                    yn=raw_input('Were any other object lines missed or wrongfully marked? (y/n, case sensitive) ')
                    if yn=='y':
                        mode=raw_input('Do you wish to add a line that was missed (+) or remove a line that was incorrectly marked (-)? ')
                        if mode=='+':
                            try:
                                line=int(input('Where is the line that was missed? '))
                                linelocs.append(line)
                                pylab.clf(); pylab.plot(onedtwo); pylab.suptitle('Non-Major Object Line Identification')
                                for i in range(len(linelocs)):
                                    pylab.plot(linelocs[i],onedtwo.mean(),'r*')
                                    pylab.vlines(linelocs[i],onedtwo.min(),onedtwo.max())
                            except:
                                print('Invalid input.  Try again.')
                        elif mode=='-':
                            try:
                                line=int(input('Where is the line that should not have been added? '))
                                linelocs.remove(line)
                                pylab.clf(); pylab.plot(onedtwo); pylab.suptitle('Non-Major Object Line Identification')
                                for i in range(len(linelocs)):
                                    pylab.plot(linelocs[i],onedtwo.mean(),'r*')
                                    pylab.vlines(linelocs[i],onedtwo.min(),onedtwo.max())
                            except:
                                print('Invalid input.  Try again.')
                        else:
                            print('Invalid input.  Try again.')
                    elif yn=='n':
                        print('Moving on to object masking.')
                        mcon=1
                    else:
                        print('Invalid input.  Try again.')
                linelocs=numpy.array(linelocs)
                linelocs=numpy.sort(linelocs,axis=None)
                linelocs=linelocs.tolist()
                print 'Non-major object lines were confirmed at row(s): %s' % linelocs

                #Mask other object lines
                maskwidth=5
                mask2=numpy.abs(Y-linelocs[0])<=maskwidth
                for i in range(len(linelocs)-1):
                    mask2=(mask2) | (numpy.abs(Y-linelocs[i+1])<=maskwidth)
                masklines=numpy.ma.array(maskobj.data,mask=(mask2) | (maskobj.mask))
                pylab.figure(); pylab.imshow(masklines); pylab.suptitle('Object Lines Masked')
                onedthree=masklines.mean(axis=1)
                pylab.figure(); pylab.plot(onedthree); pylab.suptitle('One-Dimensional Compression of Masked Object Lines')
                
                #Trim the data
                trimdata=masklines[trimlow:trimhigh,:]
                pylab.figure(); pylab.imshow(trimdata); pylab.suptitle('Trimmed Region of Masked Data')
                onedfour=trimdata.mean(axis=1)
                pylab.figure(); pylab.plot(onedfour); pylab.suptitle('One-Dimensional Compression of Trimmed Region')

                #Mask around the wavelength gaps
                zerocols=(numpy.where(f[k].data[int(secloc),:]==0))[0]
                zeromask=zerocols.tolist()
                for i in zerocols:
                    lowcol=i-7
                    modcol=i-5
                    highcol=i+7
                    if lowcol<0:
                        lowcol=0
                    if modcol<0:
                        modcol=0
                    if highcol>zerocols.max():
                        highcol=zerocols.max()
                    if lowcol not in zeromask:
                        zeromask.append(lowcol)
                    if modcol not in zeromask:
                        zeromask.append(modcol)
                    if highcol not in zeromask:
                        zeromask.append(highcol)
                zeromask=numpy.sort(zeromask,axis=None)
                gapmask=trimdata.mask
                gapmask[:,zeromask]=True
                maskgap=numpy.ma.array(trimdata,mask=gapmask)
                pylab.figure(); pylab.imshow(maskgap); pylab.suptitle('Wavelength Gaps Masked')
                maskgap1d=maskgap.mean(axis=1)
                pylab.figure(); pylab.plot(maskgap1d); pylab.suptitle('One-D Compression of Masked Gaps')
            
                #Fit the sky data
                sky=numpy.zeros_like(masklines.data)
                bin=int(input('What would you like the bin size to be? '))
                fitorder=int(input('What order should the fit be? '))
                for xx in x:
                    x1=xx-bin
                    x2=xx+bin
                    if x1<0:
                        x1=0
                    if x2>x.max():
                        x2=x.max()                
                    try:
                        if bin!=0:
                            fit=numpy.ma.polyfit(y[trimlow:trimhigh],masklines[trimlow:trimhigh,x1:x2].mean(axis=1),fitorder)
                            sky[:,xx]=numpy.polyval(fit,y)
                        else:
                            fit=numpy.ma.polyfit(y[trimlow:trimhigh],masklines[trimlow:trimhigh,xx],fitorder)
                            sky[:,xx]=numpy.polyval(fit,y)
                    except:
                        sky[:,xx]=0
                sky=sky[trimlow:trimhigh,:]
                sky[:,zeromask]=trimdata[:,zeromask]
                pylab.figure(); pylab.imshow(sky); pylab.suptitle('Sky Fit')
                pylab.figure(); pylab.plot(maskgap1d); pylab.plot(sky.mean(axis=1)); pylab.suptitle('One-D of Sky Fit')

                #Check the sky fit
                nz2=numpy.where(maskgap1d.data!=0)
                rem=maskgap-sky
                remoned=rem[nz2[0],:].mean(axis=1)
                pylab.figure(); pylab.plot(remoned); pylab.suptitle('Sky Fit Applied')
                pylab.show()
                
                #Get the variance data and finish checking the fit
                vardata=numpy.ma.array(v[k].data,mask=(mask2) | (maskobj.mask))
                vardata=vardata[trimlow:trimhigh,:]
                if vardata.shape==sky.shape:
                    chisqr=(((maskgap-sky)**2)/vardata).sum()
                    dof=(((1-maskgap.mask).sum())-(fitorder*len(x)))
                    print 'The Chi-Squared value is %s' % chisqr
                    print 'The DOF is %s' %dof
                    cod=chisqr/dof
                    print 'The Chi-Squared divided by DOF is then %s' % cod
                    acon=0
                    while acon!=1:
                        acc=raw_input('Is this Chi-Squared/DOF acceptable? (y/n, case sensitive) ')
                        if acc=='y':
                            print 'Moving on to writing sky file.'
                            acon=1
                        elif acc=='n':
                            acon=1
                            print 'Unfortunate. Terminating script'
                            sys.exit('Unacceptable Chi-Squared/DOF')
                        else:
                            print 'Invalid input.  Try again.'
                else:
                    sys.exit('Warning: Cannot statistically evaluate data.  Variance file and sky fit have different shape')

                #Write the sky file
                runtime=time.strftime("%c")
                skyfile=(root+'_sky'+ext)
                f[k].data=sky
                f[k].header['CRTIME']=(runtime,'Time at which this file was written')
                f[k].header['BIN']=(2*bin,'No. columns in moving boxcar fit')
                f[k].header['FITORDER']=(fitorder,'Order of the polynomial for sky fit')
                if os.path.exists(skyfile):
                    print 'Warning: File name %s already exists' % skyfile
                    ovw=0
                    while ovw!=1:
                        acon=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
                        if acon=='y':
                            print 'Overwriting file %s' % skyfile
                            os.remove(skyfile)
                            f.writeto(skyfile)
                            ovw=1
                        elif acon=='n':
                            ovw=1
                            print 'File %s will not be written' % skyfile
                        else:
                            print 'Invalid input.  Try again.'
                else:
                    print 'Writing file %s' % skyfile
                    f.writeto(skyfile)

                #Write the sky subtracted file
                skysubfile=(root+'_skysub'+ext)
                f[k].data=pyfits.open(img)[k].data[trimlow:trimhigh,:]-sky
                f[k].header['CRTIME']=(runtime,'Time at which this file was written')
                f[k].header['CHISQR']=(chisqr,'Chi-Squared value calculated')
                f[k].header['DOF']=(dof,'Degrees of freedom')
                f[k].header['NUMMASK']=(((1-trimdata.mask).sum()),'Number of unmasked data points in the fit')
                f[k].header['COD']=(cod,'Ratio of Chi-Squared to DOF')
                if os.path.exists(skysubfile):
                    print 'Warning: File name %s already exists' % skysubfile
                    ov=0
                    while ov!=1:
                        ocon=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
                        if ocon=='y':
                            print 'Overwriting file %s' % skysubfile
                            os.remove(skysubfile)
                            f.writeto(skysubfile)
                            ov=1
                        elif ocon=='n':
                            print 'File %s will not be written.' % skysubfile
                            ov=1
                        else:
                            print 'Invalid input.  Try again.'
                else:
                    print 'Writing file %s' % skysubfile
                    f.writeto(skysubfile)

    else:
        print 'Error: File %s does not exist.  Try again.' % img
        print '--------------------------------------'
