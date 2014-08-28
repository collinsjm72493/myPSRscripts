"""
calcmjd.py
Justin Collins
Aug 13, 2014

A simple routine to search through subdirectories and insert the MJD
to every sky-subtracted file
"""

import astropy
from astropy.time import Time
import os,sys,glob
import re
import pyfits
from astropy.table import Table

dirlist=os.listdir('./') # Builds a list of directories and files within the current directory
## Note: this script was written with the basic assumption of starting in my salt_reduced/ directory
ufiles=[]
mjdvals=[]
for dirx in dirlist:
    if os.path.isdir(dirx)==True:
        prodir=dirx+'/product/' # Again, this is specific to my path-naming nomenclature
        ssfiles=glob.glob(prodir+'*_skysub.fits')
        for x in ssfiles:
            print 'Current file being updated: '+x
            ufiles.append(x)
            f=pyfits.open(x)
            dt=f[0].header['DATE-OBS']+' '+f[0].header['TIME-OBS'] # Get the observation date/time in ISO format
            t=Time(dt,format='iso',scale='utc')
            mjd=t.mjd # The Time object contains the MJD as an attribute
            exptime=f[0].header['EXPTIME']
            mjd=mjd+((exptime/2)/3600/24) # Add half of the exposure time, converted from s to days
            mjdvals.append(mjd)
            print '{0:^20} {1:^20}'.format('Input Date/Time','Output MJD')
            print '{0:^20} {1:^20}'.format(dt,mjd)
            # Update the FITS header with the MJD
            f[0].header['GMJD']=(mjd,'The geocentric Modified Julian Day at the midpoint of the exposure') 
            # To make the change permanent, the file needs to be written
            f.writeto(x,clobber=True)
            f.close()

# Add these mjd values to the hdf5 table of the radial velocities for each observation
if len(ufiles)==len(mjdvals):
    rv=Table.read('rvtable2.hdf5',path='OVA')
    for i in range(len(ufiles)):
        obsno=ufiles[i].split('xmcxgbp')[1].split('_skysub.fits')[0]
        for k in range(len(rv)):
            if rv['Observation'][k]==obsno:
                rv['GMJD'][k]=mjdvals[i]
    rv.write('rvtable2.hdf5',path='OVA',overwrite=True)
    print 'The updated radial velocity table:'
    print rv
else:
    print 'Number of updated files does not match number of MJD values.'
