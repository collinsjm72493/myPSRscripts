"""
specextract.py
Written by: Justin Collins
Started on: July 18, 2014

This script will extract a one-dimensional spectrum from sky-subtracted files.
"""

import pyfits
import pylab
import scipy
import numpy
import numpy.ma
import os,sys,glob
import re
import time
import astropy
from astropy import wcs

pylab.ion()
pylab.close('all')

#Spectral lines; to be used only for comparison in this step
lines=[4226.728,4341.42,4481.2,4383.5447,4250.7866,4260.4741,4254.35,4271.76,4668,4144,4404.5,4863,6563,3968.5,3933.7,5891,5167.3,5172.7,5183.6,5889.95,5895.92]

#Operate on a directory containing sky-subtracted files
cont=0
while cont!=1:
    prodir=raw_input('Please indicate the directory (must contain sky-subtracted files) you would like to operate in: ')
    if os.path.exists(prodir):
        cont=1
        main,prod=prodir.split('product')
        rawdir=main+'raw/'
        skysublist=glob.glob(prodir+'*_skysub.fits')
        if len(skysublist)==0:
            sys.exit('Directory '+prodir+' does not contain sky-subtracted files.  Hint: Ensure ending forward slash is present.')
        else:
            print 'Sky-subtracted files to work with are: %s' % skysublist

        #Establish the location of the primary and secondary objects, as well as data parameters
        img1=skysublist[0]
        testfits=pyfits.open(img1)
        testdat=testfits[1].data
        x,y=(numpy.arange(testdat.shape[1]),numpy.arange(testdat.shape[0]))
        X,Y=numpy.meshgrid(x,y)
        test1d=testdat.mean(axis=1)
        primloc=int(numpy.where(test1d==test1d.max())[0])
        print 'The primary object is at (trimmed) row %s' % primloc
        primval=test1d[primloc]
        secloc1=int(numpy.where(test1d==test1d[:primloc-10].max())[0])
        secloc2=int(numpy.where(test1d==test1d[primloc+10:].max())[0])
        if test1d[secloc1]>test1d[secloc2]:
            secloc=secloc1
        else:
            secloc=secloc2
        print 'The secondary object is at (trimmed) row %s' % secloc
        secval=test1d[secloc]
        pylab.figure(); pylab.plot(test1d); pylab.plot(primloc,primval/2,'g*',secloc,secval/2,'r*')
        pylab.suptitle('Object Locations'); pylab.ylabel('Counts'); pylab.xlabel('Pixel Row')
        pylab.vlines(primloc,test1d.min(),test1d.max()); pylab.vlines(secloc,test1d.min(),test1d.max())

        #Identify the flux region for the objects
        #Defaults +/- 10
        primy1=primloc-10
        primy2=primloc+10
        secy1=secloc-10
        secy2=secloc+10
        ave1=test1d.mean()
        pylab.plot(primy1,ave1,'bx',primy2,ave1,'bx',secy1,ave1,'bx',secy2,ave1,'bx')
        pylab.vlines(primy1,test1d.min(),test1d.max()); pylab.vlines(primy2,test1d.min(),test1d.max())
        pylab.vlines(secy1,test1d.min(),test1d.max()); pylab.vlines(secy2,test1d.min(),test1d.max())
        pylab.text(primy1-5,test1d.max()+5,'1, '+str(primy1),rotation='vertical',fontsize=9)
        pylab.text(primy2-5,test1d.max()+5,'2, '+str(primy2),rotation='vertical',fontsize=9)
        pylab.text(secy1-5,test1d.max()+5,'3, '+str(secy1),rotation='vertical',fontsize=9)
        pylab.text(secy2-5,test1d.max()+5,'4, '+str(secy2),rotation='vertical',fontsize=9)
        icon=0
        while icon!=1:
            fix=raw_input('Would you like to modify the object flux regions? (y/n, case sensitive) ')
            if fix=='y':
                icon=1
                fcon=0
                while fcon!=1:
                    obj=raw_input('Which line would you like to move? (1 for primary left, 2 for primary right, 3 for secondary left, 4 for secondary right) ')
                    try:
                        if obj=='1':
                            primy1=int(input('Where should the primary left line be located? '))
                        elif obj=='2':
                            primy2=int(input('Where should the primary right line be located? '))
                        elif obj=='3':
                            secy1=int(input('Where should the secondary left line be located? '))
                        elif obj=='4':
                            secy2=int(input('Where should the secondary right line be located? '))
                        else:
                            print('Invalid input.  Try again.')
                        pylab.clf(); pylab.plot(test1d); pylab.plot(primloc,primval/2,'g*',secloc,secval/2,'r*')
                        pylab.suptitle('Object Locations'); pylab.ylabel('Counts'); pylab.xlabel('Pixel Row')
                        pylab.vlines(primloc,test1d.min(),test1d.max()); pylab.vlines(secloc,test1d.min(),test1d.max())
                        pylab.plot(primy1,ave1,'bx',primy2,ave1,'bx',secy1,ave1,'bx',secy2,ave1,'bx')
                        pylab.vlines(primy1,test1d.min(),test1d.max()); pylab.vlines(primy2,test1d.min(),test1d.max())
                        pylab.vlines(secy1,test1d.min(),test1d.max()); pylab.vlines(secy2,test1d.min(),test1d.max())
                        pylab.text(primy1-5,test1d.max()+5,'1, '+str(primy1),rotation='vertical',fontsize=9)
                        pylab.text(primy2-3,test1d.max()+5,'2, '+str(primy2),rotation='vertical',fontsize=9)
                        pylab.text(secy1-5,test1d.max()+5,'3, '+str(secy1),rotation='vertical',fontsize=9)
                        pylab.text(secy2-3,test1d.max()+5,'4, '+str(secy2),rotation='vertical',fontsize=9)
                        mcon=0
                        while mcon!=1:
                            imore=raw_input('Would you like to make another correction? (y/n, case sensitive) ')
                            if imore=='y':
                                mcon=1
                                fcon=0
                            elif imore=='n':
                                mcon=1
                                fcon=1
                            else:
                                print('Invalid input.  Try again.')
                                mcon=0
                                fcon=0
                    except:
                        print('Invalid input.  Returning to parent loop.')
                        fcon=0
            elif fix=='n':
                icon=1
                print('Moving on, then!')
            else:
                icon=0
                print('Invalid input. Try again.')

        print 'The primary flux will be taken between rows %s and %s.' % (primy1,primy2)
        print 'The secondary flux will be taken between rows %s and %s.' % (secy1,secy2)

        #Identify the flux region for the variance image (row numbers were not maintained when trim was applied)
        wcimg1,ext=img1.split('_skysub')
        varimg1=wcimg1+'_var'+ext
        testvarfits=pyfits.open(varimg1)
        testvardat=numpy.ma.array(testvarfits[1].data,mask=testvarfits[1].data>=10000)
        testvar1d=testvardat.mean(axis=1)
        varprimloc=int(numpy.where(testvar1d==testvar1d.max())[0])
        print 'The primary object for the variance file is at (untrimmed) row %s' % varprimloc
        varprimval=testvar1d[varprimloc]
        varsecloc1=int(numpy.where(testvar1d==testvar1d[:varprimloc-10].max())[0])
        varsecloc2=int(numpy.where(testvar1d==testvar1d[varprimloc+10:].max())[0])
        if testvar1d[varsecloc1]>testvar1d[varsecloc2]:
            varsecloc=varsecloc1
        else:
            varsecloc=varsecloc2
        print 'The secondary object for the variance file is at (untrimmed) row %s' % varsecloc
        varsecval=testvar1d[varsecloc]
        pylab.figure(); pylab.plot(testvar1d); pylab.plot(varprimloc,varprimval/2,'g*',varsecloc,varsecval/2,'r*'); pylab.suptitle('Object Locations (Variance File)'); pylab.ylabel('Counts'); pylab.xlabel('Pixel Row')
        pylab.vlines(varprimloc,testvar1d.min(),testvar1d.max()); pylab.vlines(varsecloc,testvar1d.min(),testvar1d.max())

        diffprim1=primloc-primy1
        diffprim2=primy2-primloc
        diffsec1=secloc-secy1
        diffsec2=secy2-secloc

        varprimy1=varprimloc-diffprim1
        varprimy2=varprimloc+diffprim2
        varsecy1=varsecloc-diffsec1
        varsecy2=varsecloc+diffsec2
        varave1=testvar1d.mean()
        pylab.plot(varprimy1,varave1,'bx',varprimy2,varave1,'bx',varsecy1,varave1,'bx',varsecy2,varave1,'bx')
        pylab.vlines(varprimy1,testvar1d.min(),testvar1d.max()-25); pylab.vlines(varprimy2,testvar1d.min(),testvar1d.max()-25)
        pylab.vlines(varsecy1,testvar1d.min(),testvar1d.max()-25); pylab.vlines(varsecy2,testvar1d.min(),testvar1d.max()-25)
        pylab.text(varprimy1-5,testvar1d.max()-5,'1, '+str(varprimy1),rotation='vertical',fontsize=9)
        pylab.text(varprimy2-5,testvar1d.max()-5,'2, '+str(varprimy2),rotation='vertical',fontsize=9)
        pylab.text(varsecy1-5,testvar1d.max()-5,'3, '+str(varsecy1),rotation='vertical',fontsize=9)
        pylab.text(varsecy2-5,testvar1d.max()-5,'4, '+str(varsecy2),rotation='vertical',fontsize=9)
        print 'The primary flux for the variance file will be taken between rows %s and %s.' % (varprimy1,varprimy2)
        print 'The secondary flux for the variance file will be taken between rows %s and %s.' % (varsecy1,varsecy2)

        #Extract the spectrum from each image individually
        for i in range(len(skysublist)):
            #Build a list of corresponding images to get
            path,imgno=skysublist[i].split('xmcxgbp')
            imgno,ext=os.path.splitext(imgno)
            imgno=imgno.split('_skysub')[0]
            img=path+'xmcxgbp'+imgno+'_skysub'+ext
            skyimg=path+'xmcxgbp'+imgno+'_sky'+ext
            rawimg=rawdir+imgno+ext
            wcimg=path+'xmcxgbp'+imgno+ext
            varimg=path+'xmcxgbp'+imgno+'_var'+ext
            #The raw image is required for important header data
            
            #Open all relevant images
            wc=pyfits.open(wcimg)
            f=pyfits.open(img)
            s=pyfits.open(skyimg)
            r=pyfits.open(rawimg)
            v=pyfits.open(varimg)

            #Pull relevant header info
            obj=str(r[0].header['OBJECT'])
            ra=str(r[0].header['RA'])
            pmra=str(r[0].header['PM-RA'])
            dec=str(r[0].header['DEC'])
            pmdec=str(r[0].header['PM-DEC'])
            date=str(r[0].header['DATE-OBS'])
            obstime=str(r[0].header['TIME-OBS'])
            jd=str(r[0].header['JD'])
            exptime=str(r[0].header['EXPTIME'])
            airmass=str(r[0].header['AIRMASS'])
            moonang=str(r[0].header['MOONANG'])
            pixscale=str(r[0].header['PIXSCALE'])
            ccdsum=str(r[0].header['CCDSUM'])
            grating=str(r[0].header['GRATING'])
            grtilt=str(r[0].header['GRTILT'])
            filt=str(r[0].header['FILTER'])
            gain=str(r[1].header['GAIN'])
            rdnoise=str(r[1].header['RDNOISE'])
            binned=str(f[1].header['BIN']+1)
            fito=str(f[1].header['FITORDER'])
            cod=str(f[1].header['COD'])
            nummask=str(f[1].header['NUMMASK'])
            wavescale=str(wc[1].header['CD1_1'])

            #Actual extraction portion of script
            if len(f)==len(s):
                for k in range(len(f)):
                    if k!=0:
                        print 'Current file: %s' % img
                        print 'Current extension: %s' % k
                        #Get wavelengths from wavelength-calibrated file
                        w=wcs.WCS(wc[k].header)
                        wave=[]
                        for m in range(len(x)):
                            wave.append(numpy.array([x[m],primloc]))
                        waves=w.wcs_pix2world(wave,1)[:,0]

                        #Simple sum unidimensional spectrum
                        primdat=numpy.ma.array(f[k].data[primy1:primy2,:],mask=v[k].data[varprimy1:varprimy2,:]>=100000)
                        skyprimdat=numpy.ma.array(s[k].data[primy1:primy2,:],mask=v[k].data[varprimy1:varprimy2,:]>=100000)
                        varprimdat=numpy.ma.array(v[k].data[varprimy1:varprimy2,:],mask=v[k].data[varprimy1:varprimy2,:]>=100000)
                        subprimct=primdat.sum(axis=0)
                        skyprimct=skyprimdat.sum(axis=0)
                        varprimct=varprimdat.sum(axis=0)
                        secdat=numpy.ma.array(f[k].data[secy1:secy2,:],mask=v[k].data[varsecy1:varsecy2,:]>=100000)
                        skysecdat=numpy.ma.array(s[k].data[secy1:secy2,:],mask=v[k].data[varsecy1:varsecy2,:]>=100000)
                        varsecdat=numpy.ma.array(v[k].data[varsecy1:varsecy2,:],mask=v[k].data[varsecy1:varsecy2,:]>=100000)
                        subsecct=secdat.sum(axis=0)
                        skysecct=skysecdat.sum(axis=0)
                        varsecct=varsecdat.sum(axis=0)

                        #Prompt to print a table of counts at each column
                        pcon=0
                        while pcon!=1:
                            pr=raw_input('Would you like to print out the tables? (y/n, case-sensitive) ')
                            if pr=='y':
                                pcon=1
                                print 'Primary Object:'
                                print '{0:^20} {1:^20} {2:^20} {3:^20} {4:^20}'.format('Pixel Column', 'Wavelength', 'Sky Counts', 'Sky-sub Counts', 'Variance Counts')
                                for n in range(len(subprimct)):
                                    try:
                                        print '{0:^20} {1:^20.3f} {2:^20.3f} {3:^20.3f} {4:^20.3f}'.format(x[n], waves[n], skyprimct[n], subprimct[n], varprimct[n])
                                        print '-------------------------------------------------------------------------------------------------------'
                                    except:
                                        print '{0:^20} {1:^20.3f} {2:^20} {3:^20} {4:^20}'.format(x[n], waves[n], skyprimct[n], subprimct[n], varprimct[n])
                                        print '-------------------------------------------------------------------------------------------------------'
                                print '-------------------------------------------------------------------------------------------------------'
                                print '-------------------------------------------------------------------------------------------------------'
                                print 'Secondary Object:'
                                print '{0:^20} {1:^20} {2:^20} {3:^20} {4:^20}'.format('Pixel Column', 'Wavelength', 'Sky Counts', 'Sky-sub Counts', 'Variance Counts')
                                for n in range(len(subsecct)):
                                    try:
                                        print '{0:^20} {1:^20.3f} {2:^20.3f} {3:^20.3f} {4:^20.3f}'.format(x[n], waves[n], skysecct[n], subsecct[n], varsecct[n])
                                        print '-------------------------------------------------------------------------------------------------------'
                                    except:
                                        print '{0:^20} {1:^20} {2:^20} {3:^20} {4:^20}'.format(x[n], waves[n], skysecct[n], subsecct[n], varsecct[n])
                                        print '-------------------------------------------------------------------------------------------------------'
                            elif pr=='n':
                                pcon=1
                                print 'Moving on, then!'
                            else:
                                print 'Invalid input.  Try again.'

                        #Plot the extraction and prompt to save the figures
                        primfig='../Desktop/Screenshots/SpectrumExtractionPics/'+main+'primextr_xmcxgbp'+imgno+'_skysubext'+str(k)+'.png'
                        pylab.figure(); pylab.plot(waves,subprimct,label='Sky-Subtracted'); pylab.suptitle('Primary Object Extraction for '+img+', Extension No. '+str(k))
                        pylab.ylabel('Counts'); pylab.xlabel('Wavelength (A)')
                        for w in lines:
                            if waves.min()<w<waves.max():
                                pylab.text(w-5,subprimct.max()-15,str(w),rotation='vertical',fontsize=9)
                                pylab.vlines(w,0,subprimct.max()-50,linestyles='dotted')
                        primscon=0
                        while primscon!=1:
                            primsave=raw_input('Would you like to save this figure? (y/n, case-sensitive) ')
                            if primsave=='y':
                                primscon=1
                                primadcon=0
                                while primadcon!=1:
                                    primad=raw_input('Would you like to adjust this figure prior to saving? (y/n, case-sensitive) ')
                                    if primad=='y':
                                        primwait=raw_input('Hit <return> when you are finished adjusting the figure.')
                                        primadcon=1
                                    elif primad=='n':
                                        primadcon=1
                                    else:
                                        print 'Invalid input.  Try again.'
                                if os.path.exists(primfig):
                                    print 'Warning: File name %s already exists.  This file will be overwritten.' % primfig
                                    os.remove(primfig)
                                    pylab.savefig(primfig)
                                else:
                                    pylab.savefig(primfig)
                            elif primsave=='n':
                                print 'This figure will not be saved.'
                                primscon=1
                            else:
                                print 'Invalid input.  Try again.'

                        secfig='../Desktop/Screenshots/SpectrumExtractionPics/'+main+'secextr_xmcxgbp'+imgno+'_skysubext'+str(k)+'.png'
                        pylab.figure(); pylab.plot(waves,subsecct,label='Sky-Subtracted'); pylab.suptitle('Secondary Object Extraction for '+img+', Extension No.'+str(k))
                        pylab.ylabel('Counts'); pylab.xlabel('Wavelength (A)')
                        for w in lines:
                            if waves.min()<w<waves.max():
                                pylab.text(w-5,subsecct.max()-15,str(w),rotation='vertical',fontsize=9)
                                pylab.vlines(w,0,subsecct.max()-50,linestyles='dotted')
                        secscon=0
                        while secscon!=1:
                            secsave=raw_input('Would you like to save this figure? (y/n, case-sensitive) ')
                            if secsave=='y':
                                secscon=1
                                secadcon=0
                                while secadcon!=1:
                                    secad=raw_input('Would you like to adjust this figure prior to saving? (y/n, case-sensitive) ')
                                    if secad=='y':
                                        secwait=raw_input('Hit <return> when you are finished adjusting the figure.')
                                        secadcon=1
                                    elif secad=='n':
                                        secadcon=1
                                    else:
                                        print 'Invalid input.  Try again.'
                                if os.path.exists(secfig):
                                    print 'Warning: File name %s already exists.  This file will be overwritten.' % primfig
                                    os.remove(secfig)
                                    pylab.savefig(secfig)
                                else:
                                    pylab.savefig(secfig)
                            elif secsave=='n':
                                print 'This figure will not be saved.'
                                secscon=1
                            else:
                                print 'Invalid input.  Try again.'
                        
                        #Build the table(s) to be written into a txt file
                        crtime=time.strftime("%c")
                        outtable=path+imgno+'_extrcol'+str(k)+'.txt'
                        if os.path.exists(outtable):
                            print 'Warning: File %s already exists.' % outtable
                            owcon=0
                            while owcon!=1:
                                ov=raw_input('Would you like to overwrite this file? (y/n, case-sensitive) ')
                                if ov=='y':
                                    owcon=1
                                    t=open(outtable,'w')
                                    t.write('# This table contains data for the extraction of the sky-subtracted file '+img+' \n# as well as relevant header data from corresponding images. \n')
                                    t.write('# The following header data were deemed relevant: \n')
                                    t.write('# OBJECT: '+obj+'\n')
                                    t.write('# Target RA: '+ra+'\n')
                                    t.write('# Proper motion in RA (mas/yr): '+pmra+'\n')
                                    t.write('# Target Declination: '+dec+'\n')
                                    t.write('# Proper motion in Dec (mas/yr): '+pmdec+'\n')
                                    t.write('# Date of observation: '+date+'\n')
                                    t.write('# Time of observation: '+obstime+'\n')
                                    t.write('# Julian Day: '+jd+'\n')
                                    t.write('# Duration of exposure (s): '+exptime+'\n')
                                    t.write('# Airmass: '+airmass+'\n')
                                    t.write('# Angle between the Moon and pointing (degrees): '+moonang+'\n')
                                    t.write('# Raw image pixel scale (arcseconds/unbinned pixel): '+pixscale+'\n')
                                    t.write('# On-chip CCD summation: '+ccdsum+'\n')
                                    t.write('# Grating used: '+grating+'\n')
                                    t.write('# Grating angle: '+grtilt+'\n')
                                    t.write('# Filter used: '+filt+'\n')
                                    t.write('# Nominal CCD gain (e/ADU): '+gain+'\n')
                                    t.write('# Nominal readout noise (e): '+rdnoise+'\n')
                                    t.write('# Wavelength Dispersion (A/pixel): '+wavescale+'\n')
                                    t.write('# Number of pixel columns used in "moving boxcar" for sky fit: '+binned+'\n')
                                    t.write('# Order of polynomial used to fit the sky: '+fito+'\n')
                                    t.write('# Value of Chi-Squared divided by DOF in application of sky fit: '+cod+'\n')
                                    t.write('# Number of unmasked pixels used in the creation of the sky fit: '+nummask+'\n')
                                    t.write('# \n')
                                    t.write('# Object Locations (with respect to the trimmed images produced in sky subtraction): \n')
                                    t.write('# Primary Object: '+str(primloc)+'\n')
                                    t.write('# Secondary Object: '+str(secloc)+'\n')
                                    t.write('# Flux regions: \n')
                                    t.write('# Primary, rows '+str(primy1)+' to '+str(primy2)+'\n')
                                    t.write('# Secondary, rows '+str(secy1)+' to '+str(primy2)+'\n')
                                    t.write('# \n')
                                    t.write('# \n')
                                    t.write('# \n')
                                    t.write('# The following data were added to this file at: '+crtime+'\n')
                                    t.write('# \n')
                                    t.write('# Object extraction from: '+img+'\n')
                                    t.write('# \n')
                                    t.write('# {0:^18} {1:^45} {2:^45} \n'.format(' ','Primary Object','Secondary Object'))
                                    t.write('# {0:^8} {1:^10} {2:^15} {3:^15} {4:^15} {5:^15} {6:^15} {7:^15} \n'.format('Pix Col','Wavelength','Prim Sky','Prim Sky-sub','Prim Var','Sec Sky','Sec Sky-Sub','Sec Var'))
                                    t.write('# {0:^8} {1:^10} {2:^15} {3:^15} {4:^15} {5:^15} {6:^15} {7:^15} \n'.format('(pix)','(A)','(cts)','(cts)','(cts**2)','(cts)','(cts)','(cts**2)'))
                                    for n in range(len(subprimct)):
                                        try:
                                            t.write('{0:^10} {1:^10.3f} {2:^15.3f} {3:^15.3f} {4:^15.3f} {5:^15.3f} {6:^15.3f} {7:^15.3f} \n'.format(x[n], waves[n], skyprimct[n], subprimct[n], varprimct[n], skysecct[n], subsecct[n], varsecct[n]))
                                        except:
                                            t.write('{0:^10} {1:^10} {2:^15} {3:^15} {4:^15} {5:^15} {6:^15} {7:^15} \n'.format(x[n], waves[n], skyprimct[n], subprimct[n], varprimct[n], skysecct[n], subsecct[n], varsecct[n]))
                                        t.write('#-------------------------------------------------------------------------------------------------------\n')
                                    print 'File %s has been written.' % outtable
                                elif ov=='n':
                                    owcon=1
                                    apcon=0
                                    while apcon!=1:
                                        ap=raw_input('Would you like to append to this file? (y/n, case-sensitive) ')
                                        if ap=='y':
                                            apcon=1
                                            t=open(outtable,'a')
                                            t.write('# \n')
                                            t.write('# \n')
                                            t.write('# The following data were added to this file at: '+crtime+'\n')
                                            t.write('# \n')
                                            t.write('# Object extraction from: '+img+'\n')
                                            t.write('# \n')
                                            t.write('# {0:^8} {1:^10} {2:^10} {3:^10} {4:^10} {5:^10} {6:^10} {7:^10} \n'.format('Pix Col','Wavelength','Prim Sky','Prim Sky-sub','Prim Var','Sec Sky','Sec Sky-sub','Sec Var'))
                                            t.write('{0:^10} {1:^10} {2:^10} {3:^10} {4:^10} {5:^10} {6:^10} {7:^10} \n'.format('(pix)','(A)','(cts)','(cts)','(cts**2)','(cts)','(cts)','(cts**2)'))
                                            for n in range(len(subprimct)):
                                                try:
                                                    t.write('{0:^10} {1:^10.3f} {2:^15.3f} {3:^15.3f} {4:^15.3f} {5:^15.3f} {6:^15.3f} {7:^15.3f} \n'.format(x[n], waves[n], skyprimct[n], subprimct[n], varprimct[n], skysecct[n], subsecct[n], varsecct[n]))
                                                except:
                                                    t.write('{0:^10} {1:^10} {2:^15} {3:^15} {4:^15} {5:^15} {6:^15} {7:^15} \n'.format(x[n], waves[n], skyprimct[n], subprimct[n], varprimct[n], skysecct[n], subsecct[n], varsecct[n]))
                                                t.write('#-------------------------------------------------------------------------------------------------------\n')
                                        elif ap=='n':
                                            apcon=1
                                            print 'File %s will not be updated.' % outtable
                                        else:
                                            print 'Invalid input.  Try again.'
                                else:
                                    print 'Invalid input.  Try again.'
                        else:
                            os.open(outtable,os.O_CREAT)
                            t=open(outtable,'w')
                            t.write('# This table contains data for the extraction of the sky-subtracted file '+img+' \n# as well as relevant header data from corresponding images. \n')
                            t.write('# The following header data were deemed relevant: \n')
                            t.write('# OBJECT: '+obj+'\n')
                            t.write('# Target RA: '+ra+'\n')
                            t.write('# Proper motion in RA (mas/yr): '+pmra+'\n')
                            t.write('# Target Declination: '+dec+'\n')
                            t.write('# Proper motion in Dec (mas/yr): '+pmdec+'\n')
                            t.write('# Date of observation: '+date+'\n')
                            t.write('# Time of observation: '+obstime+'\n')
                            t.write('# Julian Day: '+jd+'\n')
                            t.write('# Duration of exposure (s): '+exptime+'\n')
                            t.write('# Airmass: '+airmass+'\n')
                            t.write('# Angle between the Moon and pointing (degrees): '+moonang+'\n')
                            t.write('# Raw image pixel scale (arcseconds/unbinned pixel): '+pixscale+'\n')
                            t.write('# On-chip CCD summation: '+ccdsum+'\n')
                            t.write('# Grating used: '+grating+'\n')
                            t.write('# Grating angle: '+grtilt+'\n')
                            t.write('# Filter used: '+filt+'\n')
                            t.write('# Nominal CCD gain (e/ADU): '+gain+'\n')
                            t.write('# Nominal readout noise (e): '+rdnoise+'\n')
                            t.write('# Wavelength Dispersion (A/pixel): '+wavescale+'\n')
                            t.write('# Number of pixel columns used in "moving boxcar" for sky fit: '+binned+'\n')
                            t.write('# Order of polynomial used to fit the sky: '+fito+'\n')
                            t.write('# Value of Chi-Squared divided by DOF in application of sky fit: '+cod+'\n')
                            t.write('# Number of unmasked pixels used in the creation of the sky fit: '+nummask+'\n')
                            t.write('# \n')
                            t.write('# \n')
                            t.write('# Object Locations (with respect to the trimmed images produced in sky subtraction): \n')
                            t.write('# Primary Object: '+str(primloc)+'\n')
                            t.write('# Secondary Object: '+str(secloc)+'\n')
                            t.write('# Flux regions: \n')
                            t.write('# Primary, rows '+str(primy1)+' to '+str(primy2)+'\n')
                            t.write('# Secondary, rows '+str(secy1)+' to '+str(primy2)+'\n')
                            t.write('# \n')
                            t.write('# \n')
                            t.write('# The following data were added to this file at: '+crtime+'\n')
                            t.write('# Object extraction from: '+img+'\n')
                            t.write('# \n')
                            t.write('# {0:^18} {1:^25} {2:^45} \n'.format(' ','Primary Object','Secondary Object'))
                            t.write('# {0:^8} {1:^10} {2:^15} {3:^15} {4:^15} {5:^15} {6:^15} {7:^15} \n'.format('Pix Col', 'Wavelength', 'Prim Sky', 'Prim Sky-sub', 'Prim Var', 'Sec Sky', 'Sec Sky-Sub', 'Sec Var'))
                            t.write('# {0:^8} {1:^10} {2:^15} {3:^15} {4:^15} {5:^15} {6:^15} {7:^15} \n'.format('(pix)', '(A)', '(cts)', '(cts)', '(cts**2)', '(cts)', '(cts)', '(cts**2)'))
                            for n in range(len(subprimct)):
                                try:
                                    t.write('{0:^10} {1:^10.3f} {2:^15.3f} {3:^15.3f} {4:^15.3f} {5:^15.3f} {6:^15.3f} {7:^15.3f} \n'.format(x[n], waves[n], skyprimct[n], subprimct[n], varprimct[n], skysecct[n], subsecct[n], varsecct[n]))
                                except:
                                    t.write('{0:^10} {1:^10} {2:^15} {3:^15} {4:^15} {5:^15} {6:^15} {7:^15} \n'.format(x[n], waves[n], skyprimct[n], subprimct[n], varprimct[n], skysecct[n], subsecct[n], varsecct[n]))
                                t.write('#-------------------------------------------------------------------------------------------------------\n')
                            print 'File %s has been written.' % outtable

                      
    else:
        print 'Path %s does not exist.  Try again.' % prodir
