import os,sys,glob
import numpy as np
import scipy.ndimage as nd
from astropy.constants import c
import astropy.units as u
from astropy.units import Quantity
from astropy.table import Table, Column
import warnings
from numpy import genfromtxt
from scipy import interpolate
from astropy.io import ascii
import pylab
import pyfits
import time
import jdcal
from scipy.optimize import curve_fit
import re

pylab.ion()
pylab.close('all')

# TODO: use unit stuff more now it is better developped
def vbyc(v):
    return (Quantity(v, unit=u.km/u.s)/c).to(1).value


def doppler(v):
    return 1.+vbyc(v)


def doppler_undo(v):
    return 1./doppler(v)

doppler_logw = vbyc


class Fit1(object):
    def __init__(self, w, f, e, mw, mf, npol, doppler):
        self.w = w
        self.f = f
        self.e = e
        self.mw = mw
        self.mf = mf
        self.npol = npol
        self.doppler = doppler
        # pre-calculate f/e and polynomial bases
        self.fbye = f/self.e
        self.wmean = w.mean()
        # Vp[:,0]=1, Vp[:,1]=w, .., Vp[:,npol]=w**npol
        self.Vp = np.polynomial.polynomial.polyvander(w/self.wmean-1., npol)
        self.rcond = len(f)*np.finfo(f.dtype).eps
        self.sol = None

    def __call__(self, v):
        mfi = np.interp(self.w, self.mw*self.doppler(v), self.mf)
        # V[:,0]=mfi/e, Vp[:,1]=mfi/e*w, .., Vp[:,npol]=mfi/e*w**npol
        V = self.Vp*(mfi/self.e)[:,np.newaxis]
        # normalizes different powers
        scl = np.sqrt((V*V).sum(0))
        sol, resids, rank, s = np.linalg.lstsq(V/scl, self.fbye, self.rcond)
        self.sol = (sol.T/scl).T
        if rank != self.npol:
            msg = "The fit may be poorly conditioned"
            warnings.warn(msg)
        fit = np.dot(V, self.sol)*self.e
        chi2 = np.sum(((self.f-fit)/self.e)**2)
        return chi2, fit, mfi


def fittable(obs, model, *args, **kwargs):
    return fit(obs['w'], obs['flux'], obs['err'], model['Wavelength'], model['Scaled Model'],
               *args, **kwargs)


def fit(w, f, e, mw, mf, vgrid, npol,
        sigrange=None, vrange=None, doppler=doppler, plot=False):

    vgrid = Quantity(vgrid, u.km/u.s)
    chi2 = Table([vgrid.value, np.zeros_like(vgrid.value)], names=['v','chi2'])
    chi2['v'].units = vgrid.unit

    fit1 = Fit1(w, f, e, mw, mf, npol, doppler)

    chi2['chi2'] = np.array([fit1(v)[0] for v in vgrid])

    chi2.meta['ndata'] = len(f)
    chi2.meta['npar'] = npol+1+1
    chi2.meta['ndof'] = chi2.meta['ndata']-chi2.meta['npar']

    if plot:
        import matplotlib.pylab as plt
        plt.scatter(chi2['v'], chi2['chi2'])

    if vrange is None and sigrange is None or len(vgrid) < 3:
        ibest = chi2['chi2'].argmin()
        vbest, bestchi2 = chi2[ibest]
        chi2.meta['vbest'] = vbest
        chi2.meta['verr'] = 0.
        chi2.meta['bestchi2'] = bestchi2
    else:
        vbest, verr, bestchi2 = minchi2(chi2, vrange, sigrange, plot=plot)

    _, fit, mfi = fit1(vbest)
    chi2.meta['wmean'] = fit1.wmean
    chi2.meta['continuum'] = fit1.sol
    return chi2, fit, mfi


def minchi2(chi2, vrange=None, sigrange=None,
            fitcol='chi2fit', plot=False):
    assert vrange is not None or sigrange is not None
    if sigrange is None:
        sigrange = 1e10
    if vrange is None:
        vrange = 1e10

    iminchi2 = chi2['chi2'].argmin()
    ndof = float(chi2.meta['ndof'])
    iok = np.where((chi2['chi2'] <
                    chi2['chi2'][iminchi2]*(1.+sigrange**2/ndof)) &
                   (abs(chi2['v']-chi2['v'][iminchi2]) <= vrange))

    p = np.polynomial.Polynomial.fit(chi2['v'][iok], chi2['chi2'][iok],
                                     2, domain=[])

    if plot:
        import matplotlib.pylab as plt
        plt.scatter(chi2['v'][iok], chi2['chi2'][iok], c='g')
        plt.plot(chi2['v'], p(chi2['v']))

    vbest = -p.coef[1]/2./p.coef[2]
    # normally, get sigma from where chi2 = chi2min+1, but best to scale
    # errors, so look for chi2 = chi2min*(1+1/ndof) ->
    # a verr**2 = chi2min/ndof -> verr = sqrt(chi2min/ndof/a)
    verr = np.sqrt(p(vbest)/p.coef[2]/ndof)
    chi2.meta['vbest'] = vbest
    chi2.meta['verr'] = verr
    chi2.meta['bestchi2'] = p(vbest)
    if fitcol is not None:
        if fitcol in chi2.colnames:
            chi2[fitcol] = p(chi2['v'])
        else:
            chi2.add_column(Column(data=p(chi2['v']), name=fitcol))

    return chi2.meta['vbest'], chi2.meta['verr'], chi2.meta['bestchi2']


def observe(model, wgrid, slit, seeing, overresolve, offset=0.):
    """Convolve a model with a seeing profile, truncated by a slit, & pixelate

    Parameters
    ----------
    model: Table (or dict-like)
       Holding wavelengths and fluxes in columns 'w', 'flux'
    wgrid: array
       Wavelength grid to interpolate model on
    slit: float
       Size of the slit in wavelength units
    seeing: float
       FWHM of the seeing disk in wavelength units
    overresolve: int
       Factor by which detector pixels are overresolved in the wavelength grid
    offset: float, optional
       Offset of the star in the slit in wavelength units (default 0.)

    Returns
    -------
    Convolved model: Table
       Holding wavelength grid and interpolated, convolved fluxes
       in columns 'w', 'flux'
    """
    # make filter
    wgridres = np.min(np.abs(np.diff(wgrid)))
    filthalfsize = np.round(slit/2./wgridres)
    filtgrid = np.arange(-filthalfsize,filthalfsize+1)*wgridres
    # sigma ~ seeing-fwhm/sqrt(8*ln(2.))
    filtsig = seeing/np.sqrt(8.*np.log(2.))
    filt = np.exp(-0.5*((filtgrid-offset)/filtsig)**2)
    filt /= filt.sum()
    # convolve with pixel width
    filtextra = int((overresolve-1)/2+0.5)
    filt = np.hstack((np.zeros(filtextra), filt, np.zeros(filtextra)))
    filt = nd.convolve1d(filt, np.ones(overresolve)/overresolve)
    mint = np.interp(wgrid, model['Wavelength'], model['Model Flux'])
    mconv = nd.convolve1d(mint, filt)
    return Table([wgrid, mconv], names=('w','flux'), meta={'filt': filt}),mconv,filt

########################################################################
# All code below written by Justin Collins

#Gaussian function definition, for calling by curve_fit
def gauss(x,a,b,x0,sigma):
    return a+(b*np.exp(-(x-x0)**2/(2*(sigma**2))))

pylab.ion()
pylab.close('all')

#Prompt the user to input the fit order
con2=0
while con2!=1:
    fitorder=input('Please indicate the order of the polyfit to be used in scaling: ')
    if fitorder==int(fitorder):
        con2=1
    else:
        print 'Input must be an integer.  Try again.'

#Build list of files to work on from a directory
infiles=glob.glob('cetgsg*.fits')
if len(infiles)==0:
    sys.exit('The indicated path '+prodir+' contains no extracted data files.')
for k in range(len(infiles)):
    #Get corresponding raw and wavelength-calibrated files
    dat=infiles[k]
    datno=dat.split('cetgsgS20130')[1].split('.fits')[0]
    if re.search('613',datno)!=None:
        setno='set2/'
    elif re.search('712',datno)!=None:
        setno='set3/'
    elif re.search('713',datno)!=None or re.search('714',datno)!=None or re.search('715',datno)!=None:
        setno='set4/'
    elif re.search('810',datno)!=None:
        setno='set5/'
    else:
        sys.exit('The image '+dat+' does not have a corresponding rectified image.')
    rectimg=setno+'tgsgS20130'+datno+'.fits'
    if not os.path.exists(rectimg):
        sys.exit('The file '+dat+' does not have a corresponding rectified image '+rectimg)
    r=pyfits.open(rectimg)
    
    #Get the extracted data
    f=pyfits.open(dat)
    fdat=f['SCI'].data
    fdat1=fdat[0]
    wave1=f['SCI'].header['CRVAL1']
    dwave=f['SCI'].header['CD1_1']
    waves=[]
    for c in range(len(fdat[0])):
        waves.append(wave1+(dwave*c))
    waves=np.array(waves)
    if f['VAR'].data.shape[0]!=1:
        try:
            vdat=[]
            for q in range(f['VAR'].data.shape[1]):
                vs=[]
                for b in range(f['VAR'].data.shape[0]):
                    vs.append(f['VAR'].data[b,q])
                vdat.append(np.mean(vs))
        except:
            vdat=f['VAR'].data
    else:
        vdat=f['VAR'].data
        vdat=vdat[0]
    vdat=np.array(vdat)
    sdat=np.sqrt(vdat)
                        
    #Establish values for the parameters to be used in "observe"
    model=ascii.ui.read('G5model.dat')
    wgres=dwave #(A/unbinned pixel)
    picscale=float(f[0].header['PIXSCALE']) #(as/unbinned pixel)
    resolution=wgres/picscale #(A/as)

    #Calculate the seeing
    r1d=r['SCI'].data.mean(axis=1)
    primloc=int(np.where(r1d==r1d[200:800].max())[0])
    pleft=primloc-10
    pright=primloc+10
    prows=np.arange(pleft,pright+1)
    base=r1d[pleft:pright].min()
    pmax=r1d[primloc]

    #Fit a Gaussian to the extraction of the primary object
    primopt,primcov=curve_fit(gauss,prows,r1d[prows],p0=[base,pmax-base,primloc,1])
    psig=np.abs(primopt[-1])
    pfwhm=2.*np.sqrt(2.*np.log(2.))*psig
    halfmax=base+(0.5*(pmax-base))
    pylab.figure(); pylab.plot(r1d); pylab.suptitle('One-D Compression of Rectified Image '+rectimg)
    pylab.xlabel('Binned Pixel Row'); pylab.ylabel('Mean Counts')
    pylab.plot(primloc,halfmax,'r*')
    pylab.plot(prows,gauss(prows,*primopt),'ro:')
    pylab.text(primloc,pmax+2,'Gaussian Fit',rotation='vertical',fontsize=9)
    pylab.vlines(primloc-(0.5*pfwhm),base-5,pmax+5); pylab.vlines(primloc+(0.5*pfwhm),base-5,pmax+5)
    pylab.text(primloc,base-2,'Prim FWHM',rotation='vertical',fontsize=9)
    seeing=pfwhm #unbinned pix
    asseeing=seeing*picscale #as
    seeing=asseeing*resolution #(A)

    # Prepare a table column to save the seeing
    pseeing=np.empty_like(waves)
    pseeing[:]=np.nan
    pseeing[0]=seeing
    pseeing[1]=asseeing
    pseecol=Column(data=pseeing,name='Seeing (A,as)')

    #Establish the slit width
    sw=float(f[0].header['MASKNAME'].split('arcsec')[0]) # (arcsec)
    sw=sw*resolution #(A)

    #Establish the value for overresolve
    grating=str(f[0].header['GRATING'])
    if grating=='B1200+_G5321':
        R=3744 # Value obtained from: http://www.gemini.edu/?q=node/10375
    else:
        sys.exit('Unrecognized grating '+grating) 
    overresolve=np.around((waves.mean())/(R*wgres))
            
    #Execute the "observe"
    #For primary:
    print 'Parameters going in to "observe" for the primary object:'
    print(' {0:^20} {1:^20} {2:^20} {3:^20}'.format('Slit Width (A)','Seeing (as)','Seeing (A)','Overresolve'))
    print(' {0:^20.5f} {1:^20.5f} {2:^20.5f} {3:^20}'.format(sw,asseeing,seeing,overresolve))
    outtable1,mconv,filt=observe(model,waves,sw,seeing,overresolve)
    pylab.figure(); pylab.plot(filt); pylab.suptitle('Convolution Kernel for Primary, Seeing = '+str(seeing))
    
    #Fit a polynomial to the data divided by the model
    #Primary:
    mdiv=fdat1/mconv
            
    #Fit the data divided by the model
    # Primary:
    mfit=np.ma.polyfit(waves[-np.isnan(mdiv)],mdiv[-np.isnan(mdiv)],fitorder)
    normmod=np.polyval(mfit,waves)
    # Make a column of the fitting parameters to be saved to .hdf5 table
    mfitcol=np.empty_like(waves)
    mfitcol[:]=np.nan
    for q in range(len(mfit)):
        mfitcol[q]=mfit[q]
    polycol=Column(data=mfitcol,name='Poly Coeff')

    #Multiply the fit by the model
    # Primary:
    mod=normmod*mconv

    #View the now-scaled primary model
    pylab.figure(); pylab.subplot(2,1,1); pylab.plot(waves,mconv,label='Convolved & Interpolated Model')
    pylab.suptitle('Checking Model (Primary Object)')
    pylab.plot(waves,fdat1,label='Extracted Data for Primary Object'); pylab.xlabel('Wavelength (A)')
    pylab.legend(loc='best')
    pylab.subplot(2,1,2); pylab.plot(waves,fdat1,label='Extracted Data for Primary Object')
    pylab.plot(waves,mod,label='Normalized & Fitted Model'); pylab.suptitle('Checking Scaled Model (Primary Object)')
    pylab.xlabel('Wavelength (A)'); pylab.ylabel('Flux (erg/s/cm**2/A)'); pylab.legend(loc='best')

    #Prompt to save the image
    figfile=datno+'_primscaledmodel.png'
    con3=0
    while con3!=1:
        imsave=raw_input('Would you like to save the image of the scaled model? (y/n, case sensitive) ')
        if imsave=='y':
            con3=1
            wait=raw_input('Press <return> when ready to save.')
            if os.path.exists(figfile):
                print 'Warning: File %s already exists.' % figfile
                con4=0
                while con4!=1:
                    ow=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
                    if ow=='y':
                        con4=1
                        os.remove(figfile)
                        pylab.savefig(figfile)
                    elif ow=='n':
                        con4=1
                        print 'File %s will not be replaced.' % figfile
                    else:
                        print 'Invalid input.  Try again.'
            else:
                print 'Writing file %s' % figfile
                pylab.savefig(figfile)
        elif imsave=='n':
            con3=1
            print 'The image will not be saved.'
        else:
            print 'Invalid input.  Try again.'

    #Create a table of the model to be written to an hdf5 file
    outtable2=Table([waves,mod],names=('Wavelength','Prim Scaled Model'),meta={'name':'Model Data'})
    outtable2.add_column(pseecol)
    outtable2.add_column(polycol)
    outfile=datno+'_model.hdf5'
    month,day,year=time.strftime('%x').split('/')
    if os.path.exists(outfile):
        print 'Warning: File %s already exists.' % outfile
        con5=0
        while con5!=1:
            ovw=raw_input('Would you like to append to this file? (y/n, case sensitive) ')
            if ovw=='y':
                con5=1
                try:
                    outtable2.write(outfile,path='Written'+month+day+year,append=True)
                except:
                    outtable2.write(outfile,path='Written'+month+day+year,overwrite=True)
            elif ovw=='n':
                con5=1
                print 'File %s will not be updated.' % outfile
            else:
                print 'Invalid input.  Try again.'
    else:
        outtable2.write(outfile,path='Written'+month+day+year)

    #Velocity fit, following Daniel's example
    #Paper shows roughly -150 < v < 200, will try from -300 to +300
    v=range(-300,300)
    for jj in range(len(v)):
        v[jj]=float(v[jj])
    if grating=='B1200+_G5321':
        ind1=np.argmin(np.abs(waves-4000))
        ind2=np.argmin(np.abs(waves-5200))
    else:
        sys.exit('Unsupported Grating type: '+grating)
    c=ind1
    d=ind2
    waves2=waves[c:d]
    fprim=fdat1[c:d]
    ch=np.zeros_like(v)
    err=sdat[c:d]
    for j in range(len(v)):
        #Calculate chi**2 for an array of velocities
        q=interpolate.InterpolatedUnivariateSpline(waves+(waves*v[j]/(2.9979E5)),mod)
        zz=np.polyfit(waves2[-np.isnan(fprim)],(fprim[-np.isnan(fprim)]/q(waves2[-np.isnan(fprim)])),fitorder)
        if np.isnan(zz).any():
            sys.exit('Polyfit returned NaN values.  Polyfit coefficients are: '+str(zz))
        pp=np.poly1d(zz)
        res=np.zeros_like(waves)
        for m in range(c,d):           
            if fdat1[m]>=0 and pp(waves[m])>=0 and err[m-c]>0:
                res[m-c]=float((((fdat1[m]-q(waves[m])*pp(waves[m]))/(err[m-c]))**2))
        ch[j]=float(res.sum())

    #Ding!
    print('\a')

    # Plot the distribution of chi**2 results
    pylab.figure(); pylab.plot(v,ch); pylab.suptitle('Chi**2 Results for Primary, File '+dat)
    pylab.xlabel('Velocity (km/s)'); pylab.ylabel('Chi**2')
    pchifile='chi2'+datno+'_'+month+day+year+'.png'
    if os.path.exists(pchifile):
        print 'Warning: File name %s already exists.' % pchifile
        chcon=0
        while chcon!=1:
            chwr=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
            if chwr=='y':
                chcon=1
                os.remove(pchifile)
                pylab.savefig(pchifile)
            elif chwr=='n':
                chcon=1
                print 'Primary Chi**2 figure will not be saved.'
            else:
                print 'Invalid input.  Try again.'
    else:
        pylab.savefig(pchifile)
            
    #Fit a parabola around the minimum of the chi**2 array
    ymin1=v[np.argmin(ch)]
    para=np.polyfit(v[np.argmin(ch)-6:np.argmin(ch)+7],ch[np.argmin(ch)-6:np.argmin(ch)+7],2)
    print 'The velocity corresponding to the minimum chi**2 value for the primary object from file '+dat+' is %s' % ymin1
    xx=np.arange(v[np.argmin(ch)-6],v[np.argmin(ch)+6],1)
    cc=para[0]*xx**2+para[1]*xx+para[2]
    ymin=xx[np.argmin(cc)]

    #Plot results and calculate uncertainty
    pylab.figure(); pylab.plot(xx,cc); pylab.suptitle('Parabolic Fit of Chi**2 Function for Primary Object from File '+dat)
    pylab.xlabel('Velocity (km/s)'); pylab.ylabel('Chi**2')
    pylab.plot(xx,[cc[np.argmin(cc)]+1]*len(xx))
    unc=np.abs((xx[np.argmin(np.abs(cc-(cc[np.argmin(cc)]+1)))])-(xx[np.argmin(cc)]))
    print 'The radial velocity for the primary object in the '+datno+' observation is '+str(ymin)+' +/- %s km/s' % unc
            
    #Analyze the shift at this velocity
    shiftwaves=waves+(waves*(ymin/2.9979E5))
    q2=interpolate.InterpolatedUnivariateSpline(shiftwaves,fdat1)
    if np.isnan(q2(waves[-np.isnan(fdat1)])).any():
        sys.exit('Interpolation of shifted data returned NaN values.')
    zz2=np.polyfit(waves,fdat1/q2(waves),fitorder)
    if np.isnan(zz2).any():
        sys.exit('Polyfit of shifted array returned NaN values.  Polyfit coefficients are: '+str(zz2))
    pp2=np.poly1d(zz2)
    res2=np.zeros_like(waves)
    for n in range(len(waves)):
        if sdat[n]>0:
            res2[n]=((fdat1[n]-q2(waves[n])*pp2(waves[n]))/(sdat[n]))
    ch2=np.sum(res2**2)
    pylab.figure(); pylab.subplot(2,1,1); pylab.plot(waves,fdat1,label='Extracted Data')
    pylab.plot(shiftwaves,mod,label='Shifted, Fitted & Scaled Model'); pylab.xlabel('Wavelength (A)'); pylab.ylabel('Flux (erg/s/cm**2/A)')
    pylab.suptitle('Residuals for Shifted Primary Model, Obs. Date '+datno+', Velocity = '+str(ymin)); pylab.legend(loc='best')
    pylab.subplot(2,1,2); pylab.plot(waves,res2,'.'); pylab.xlabel('Wavelength (A)'); pylab.ylabel('Residuals')
    print 'The value of Chi**2 for the velocity of the primary object is %s ' % ch2
    
    #Prompt to save the residuals figure for the primary object
    presfile='vshift'+datno+'.png'
    con5=0
    while con5!=1:
        pressave=raw_input('Would you like to save the image of the primary residuals? (y/n, case sensitive) ')
        if pressave=='y':
            con5=1
            preswait=raw_input('Press <return> when ready to save.')
            if os.path.exists(presfile):
                print 'Warning: File %s already exists.' % presfile
                con6=0
                while con6!=1:
                    prow=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
                    if prow=='y':
                        con6=1
                        os.remove(presfile)
                        pylab.savefig(presfile)
                    elif prow=='n':
                        con6=1
                        print 'File %s will not be replaced.' % presfile
                    else:
                        print 'Invalid input.  Try again.'
            else:
                print 'Writing file %s' % presfile
                pylab.savefig(presfile)
        elif pressave=='n':
            con5=1
            print 'The image will not be saved.'
        else:
            print 'Invalid input.  Try again.'

    #Write results to hdf5 table for all observations
    print 'Writing results to rvtable.hdf5:'
    et=float(f[0].header['EXPTIME'])
    gmjd=float(f[0].header['MJD-OBS'])+(0.5*et)
    rvdata=[datno,gmjd,ymin,unc]
    rvtable='rvtable_gem.hdf5'
    if os.path.exists(rvtable):
        rv=Table.read(rvtable,path='OVA')
        for d in range(len(rv['Observation'])):
            try:
                if rv['Observation'][d]==datno:
                    print 'This observation already exists at Table Row index: '+str(d)
                    print 'Overwriting row entry for '+datno
                    rv.remove_row(d)
            except:
                pass
        rv.add_row(rvdata)
    else:
        rvdata=[[datno],[gmjd],[ymin],[unc]]
        rv=Table(data=rvdata,names=('Observation','GMJD','RV','ERR'),meta={'name':'Gemini RV Table'})
    print rv
    rv.write(rvtable,path='OVA',overwrite=True)
