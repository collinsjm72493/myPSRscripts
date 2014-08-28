"""
fluxcal_gem.py
Aug. 18, 2014
Justin Collins

Use the calibration from Daniel's calibrated file to flux-calibrate the remainder of the data
"""

import os,sys,glob
import numpy as np
import pylab
import pyfits

pylab.ion()
pylab.close('all')

c=pyfits.open('cetgsgS20130613S0126.fits')
cdat=c['SCI'].data[0]
cwaves=[]
for i in range(len(cdat)):
    w=3874.37963867188+(0.23520332574844*i)
    cwaves.append(w)
pylab.figure(); pylab.plot(cwaves,cdat); pylab.suptitle('Spectrum from cetgsgS20130613S0126.fits')
pylab.xlabel('Wavelength (A)'); pylab.ylabel('Flux (erg/s/cm**2/A)')

e=pyfits.open('etgsgS20130613S0126.fits')
edat=e['SCI'].data
edat=edat[0]

con=0
while con!=1:
    fitorder=input('What would you like the fit order to be? ')
    if fitorder==abs(int(fitorder)):
        con=1
    else:
        print 'Invalid input.  Try again.'

polyfit=np.polyfit(cwaves,cdat/edat,fitorder)
polycal=np.polyval(polyfit,cwaves)
pylab.figure(); pylab.plot(cwaves,polycal,'r',label='Calibration Polynomial')
pylab.suptitle('Polynomial Fit of Calibrated / Non-Calibrated Quotient'); pylab.xlabel('Wavlength (A)')
pylab.plot(cwaves,cdat/edat,'b',label='Quotient for Polynomial Fit')
pylab.legend(loc='best')

infiles=glob.glob('etgsg*.fits')
infiles.remove('etgsgS20130613S0126.fits')
for n in range(len(infiles)):
    imgno=infiles[n].split('etgsgS')[1].split('.fits')[0]
    ef=pyfits.open(infiles[n])
    efdat=ef['SCI'].data
    if efdat.shape[0]!=1:
        efdat1=np.zeros((1,efdat.shape[1]))
        for c in range(efdat.shape[1]):
            efdats=[]
            for r in range(efdat.shape[0]):
                efdats.append(efdat[r,c])
            efdat1[0,c]=np.mean(efdats)
        efcal=efdat1*polycal
    else:
        efcal=efdat*polycal
    wave1=ef['SCI'].header['CRVAL1']
    dwave=ef['SCI'].header['CD1_1']
    waves=[]
    for m in range(len(efdat[0])):
        w=wave1+(dwave*m)
        waves.append(w)
    pylab.figure(); pylab.subplot(2,1,1); pylab.plot(waves,efdat[0]); pylab.xlabel('Wavelength (A)'); pylab.ylabel('Counts')
    pylab.subplot(2,1,2); pylab.plot(waves,efcal[0]); pylab.xlabel('Wavelength (A)'); pylab.ylabel('Flux (erg/s/cm**2/A)')
    pylab.suptitle('Flux Calibration for File: '+infiles[n])
    
    #Prompt to save figures
    scon=0
    while scon!=1:
        fsave=raw_input('Would you like to save this figure? (y/n, case sensitive) ')
        if fsave=='y':
            scon=1
            wait=raw_input('Press <return> when ready to save the figure.')
            figfile=imgno+'_fluxcal.png'
            if os.path.exists(figfile):
                print 'Warning: File name %s already exists.' % figfile
                ocon=0
                while ocon!=1:
                    ow=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
                    if ow=='y':
                        os.remove(figfile)
                        ocon=1
                        pylab.savefig(figfile)
                    elif ow=='n':
                        print 'File %s will not be overwritten.' % figfile
                        ocon=1
                    else:
                        print 'Invalid input.  Try again.'
            else:
                print 'Saving figure to file %s' % figfile
                pylab.savefig(figfile)
        elif fsave=='n':
            scon=1
            print 'Figure will not be saved.'
        else:
            print 'Invalid input.  Try again.'
    
    #Save calibrated file to FITS
    outfits='cetgsgS'+imgno+'.fits'
    ef['SCI'].data=efcal
    if os.path.exists(outfits):
        print 'Warning: The output file %s already exists.' % outfits
        wcon=0
        while wcon!=1:
            owr=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
            if owr=='y':
                wcon=1
                ef.writeto(outfits,clobber=True)
            elif owr=='n':
                wcon=1
                print 'File %s will not be overwritten.'
            else:
                print 'Invalid input.  Try again.'
    else:
        print 'Writing new FITS file %s' % outfits
        ef.writeto(outfits)
    
