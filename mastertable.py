"""
mastertable.py
Aug. 29, 2014
Justin Collins

Pull all SALT and Gemini data and build a master reference
table for the collective results.
"""

import os,sys,glob
import pyfits
import pylab
import astropy
from astropy.table import Table

# The data I am collecting are the Observation No. (file name),
# telescope/instrument used, observation date and time, 
# geocentric MJD, barycentric MJD, grating used, grating angle,
# seeing (in arcsec), and pulsar RV (+/- uncertainty); also 
# putting in orbital phase (for completeness)

# UPDATE: Also adding wavelength resolution and slit width

# Initialize lists
mobsno=[]
mtelins=[]
mdate=[]
mtime=[]
mgmjd=[]
mbmjd=[]
mgrating=[]
mgrangle=[]
mseeing=[]
mrv=[]
munc=[]
mphase=[]
mwres=[]
msw=[]

# Starting with SALT
salttable='salt_reduced/phaservtable.hdf5'
st=Table.read(salttable,path='OVA')
for i in range(len(st)):
    # Pull what data I can from the 'salttable'
    sobsno=st['Observation'][i]
    print 'The current SALT observation number is: '+sobsno
    telins='SALT RSS'
    sgmjd=st['GMJD'][i]
    sbmjd=st['BMJD'][i]
    sphase=st['Phase'][i]
    srv=st['Prim RV'][i]
    sunc=st['Prim Unc'][i]

    # Pull some other data from the rectified image's FITS header
    if sobsno=='P201305120029' or sobsno=='P201305120030':
        sdir='512a/'
    elif sobsno=='P201305120054' or sobsno=='P201305120055':
        sdir='512b/'
    elif sobsno=='P201306160046' or sobsno=='P201306160047':
        sdir='616v2a/'
    elif sobsno=='P201306160112' or sobsno=='P201306160113':
        sdir='616v2b/'
    else:
        sdir=sobsno.split('P20130')[1][:3]
        if sdir=='607':
            sdir='607v2a/'
        elif sdir=='518' or sdir=='520' or sdir=='617' or sdir=='618':
            sdir=sdir+'v2/'
        else:
            sdir=sdir+'/'
    sprodir='salt_reduced/'+sdir+'product/'
    srectfile=sprodir+'xmcxgbp'+sobsno+'.fits'
    print 'The corresponding SALT rectified image is: '+srectfile
    s=pyfits.open(srectfile)
    sdate=s[0].header['DATE-OBS']
    stime=s[0].header['TIME-OBS']
    sgrating=s[0].header['GRATING']
    sgrangle=s[0].header['GRTILT']
    swres=s[1].header['CD1_1']
    maskid=s[0].header['MASKID']
    if maskid=='PL0100N002':
        ssw=1.0
    elif maskid=='PL0150N001':
        ssw=1.5
    else:
        sys.exit('Unsupported MASKID: '+maskid)
    
    # Pull the seeing from the '_model.hdf5' file
    smod=sprodir+sobsno+'_model.hdf5'
    sm=Table.read(smod,path='Written082614')
    sseeing=sm['Prim Seeing (A,as)'][1]
    adat=[sobsno,telins,sdate,stime,sgmjd,sbmjd,sgrating,sgrangle,sseeing,srv,sunc,sphase]
    print 'The data being added to the master table are: '
    print '{0:<20} {1:<20}'.format('OBSNO:',sobsno)
    print '{0:<20} {1:<20}'.format('TEL/INS:',telins)
    print '{0:<20} {1:<20}'.format('OBSDATE:',sdate)
    print '{0:<20} {1:<20}'.format('OBSTIME:',stime)
    print '{0:<20} {1:<20.4f}'.format('GMJD:',sgmjd)
    print '{0:<20} {1:<20.4f}'.format('BMJD:',sbmjd)
    print '{0:<20} {1:<20}'.format('GRATING:',sgrating)
    print '{0:<20} {1:<20}'.format('GR ANGLE:',sgrangle)
    print '{0:<20} {1:<20.4f}'.format('SEEING (as):',sseeing)
    print '{0:<20} {1:<20}'.format('RV:',srv)
    print '{0:<20} {1:<20}'.format('RV UNC:',sunc)
    print '{0:<20} {1:<20.4f}'.format('ORB PHASE:',sphase)

    mobsno.append(sobsno)
    mtelins.append(telins)
    mdate.append(sdate)
    mtime.append(stime)
    mgmjd.append(sgmjd)
    mbmjd.append(sbmjd)
    mgrating.append(sgrating)
    mgrangle.append(sgrangle)
    mseeing.append(sseeing)
    mrv.append(srv)
    munc.append(sunc)
    mphase.append(sphase)
    mwres.append(swres)
    msw.append(ssw)

    print '--------------------------------------------------'

# Now get the Gemini data
gemdir='gemini_psrj1723/Gemini/'
gtable=gemdir+'phaservtable_gem.hdf5'
gt=Table.read(gtable,path='OVA')
for k in range(len(gt)):
    gtelins='Gemini GMOS-S'
    # Pull the data available in gtable
    gobsno=gt['Observation'][k]
    print 'The current Gemini observation number is: '+gobsno
    ggmjd=gt['GMJD'][k]
    gbmjd=gt['BMJD'][k]
    grv=gt['RV'][k]
    gunc=gt['ERR'][k]
    gphase=gt['Phase'][k]

    # Pull other data from a extracted image's FITS header
    gimg=gemdir+'cetgsgS20130'+gobsno+'.fits'
    print 'The corresponding Gemini extracted image is: '+gimg
    g=pyfits.open(gimg)
    gdate=g[0].header['DATE-OBS']
    gtime=g[0].header['TIME-OBS']
    ggrating=g[0].header['GRATING']
    ggrangle=g[0].header['GRTILT']
    gwres=g['SCI'].header['CD1_1']
    gsw=g[0].header['MASKNAME'].split('arcsec')[0]
    
    # Pull the seeing from the '_model.hdf5' table
    gmodtable=gemdir+gobsno+'_model.hdf5'
    gm=Table.read(gmodtable,path='Written082714')
    gseeing=gm['Seeing (A,as)'][1]

    idat=[gobsno,gtelins,gdate,gtime,ggmjd,gbmjd,ggrating,ggrangle,gseeing,grv,gunc,gphase]
    print 'The data being added to the master table are: '
    print '{0:<20} {1:<20}'.format('OBSNO:',gobsno)
    print '{0:<20} {1:<20}'.format('TEL/INS:',gtelins)
    print '{0:<20} {1:<20}'.format('OBSDATE:',gdate)
    print '{0:<20} {1:<20}'.format('OBSTIME:',gtime)
    print '{0:<20} {1:<20.4f}'.format('GMJD:',ggmjd)
    print '{0:<20} {1:<20.4f}'.format('BMJD:',gbmjd)
    print '{0:<20} {1:<20}'.format('GRATING:',ggrating)
    print '{0:<20} {1:<20}'.format('GR ANGLE:',ggrangle)
    print '{0:<20} {1:<20.4f}'.format('SEEING (as):',gseeing)
    print '{0:<20} {1:<20}'.format('RV:',grv)
    print '{0:<20} {1:<20}'.format('RV UNC:',gunc)
    print '{0:<20} {1:<20.4f}'.format('ORB PHASE:',gphase)

    mobsno.append(gobsno)
    mtelins.append(gtelins)
    mdate.append(gdate)
    mtime.append(gtime)
    mgmjd.append(ggmjd)
    mbmjd.append(gbmjd)
    mgrating.append(ggrating)
    mgrangle.append(ggrangle)
    mseeing.append(gseeing)
    mrv.append(grv)
    munc.append(gunc)
    mphase.append(gphase)
    mwres.append(gwres)
    msw.append(gsw)

    print '--------------------------------------------------'

# Now write all the results to the master table
mdat=[mobsno,mtelins,mdate,mtime,mgmjd,mbmjd,mgrating,mgrangle,mseeing,mwres,msw,mrv,munc,mphase]
mt=Table(mdat,names=('Observation','Telescope/Instrument','Observation Date',
                           'Observation Time','GMJD','BMJD','Grating','Grating Angle',
                           'Seeing','Wavelength Res','Slit Width (as)','RV','RV Uncertainty','Orbital Phase'))
print 'The summary table is now:'
print mt
mastertable='summarytable.hdf5'
if os.path.exists(mastertable):
    print 'Warning: File name '+mastertable+' already exists.'
    con=0
    while con!=1:
        o=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
        if o=='y':
            con=1
            print 'Overwriting the summary table.'
            mt.write(mastertable,path='OVA',overwrite=True)
        elif o=='n':
            con=1
            print 'The summary table will not be updated.'
        else:
            print 'Invalid input.  Try again.'
else:
    print 'Writing the summary table.'
    mt.write(mastertable,path='OVA')
