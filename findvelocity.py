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

##################################################################################
# All code below written by Justin Collins

#Gaussian function definition, for calling by curve_fit
def gauss(x,a,b,x0,sigma):
    return a+(b*np.exp(-(x-x0)**2/(2*(sigma**2))))


#Build list of files to work on from a directory
con1=0
while con1!=1:
    prodir=raw_input('Please indicate the directory (must contain extracted data) you would like to work in: ')
    main,junk=prodir.split('product')
    rawdir=main+'raw/'
    if os.path.exists(prodir):
        con1=1
        infiles=glob.glob(prodir+'*_extrcol1.txt')
        if len(infiles)==0:
            sys.exit('The indicated path '+prodir+' contains no extracted data files.')
        for k in range(len(infiles)):
            #Get corresponding raw and wavelength-calibrated files
            dat=infiles[k]
            datno=dat.split(prodir)[1].split('_extrcol1.txt')[0]
            rawimg=main+'raw/'+datno+'.fits'
            if not os.path.exists(rawimg):
                sys.exit('The file '+dat+' does not have a corresponding raw image '+rawimg)
            r=pyfits.open(rawimg)
            wcimg=prodir+'xmcxgbp'+datno+'.fits'
            ssimg=prodir+'xmcxgbp'+datno+'_skysub.fits'
            if not os.path.exists(wcimg):
                sys.exit('The file '+dat+' does not have a corresponding wavelength-calibrated image '+wcimg)
            w=pyfits.open(wcimg)

            #Get the extracted data
            skysubdata=genfromtxt(dat)
            skysubwaves=np.array(skysubdata[:,1])
            skysubprim=np.ma.array(skysubdata[:,3],mask=(1-(skysubdata[:,3]>0)))
            skysubprimvar=np.ma.array(skysubdata[:,4],mask=(1-(skysubdata[:,4]>0)))
            skysubprimerr=np.ma.sqrt(skysubprimvar)
            skysubsec=np.ma.array(skysubdata[:,6],mask=(1-(skysubdata[:,6]>0)))
            skysubsecvar=np.ma.array(skysubdata[:,7],mask=(1-(skysubdata[:,7]>0)))
            skysubsecerr=np.ma.sqrt(skysubsecvar)
                        
            #Establish values for the parameters to be used in "observe"
            spacebin=4
            specbin=2
            model=ascii.ui.read('modeltable.dat')
            wgres=float(w[1].header['CD1_1']) #(A/binned pixel)
            wgres2=wgres/specbin #(A/unbinned pixel)
            picscale=float(w[0].header['PIXSCALE']) #(as/unbinned pixel)
            resolution=wgres2/picscale #(A/as)

            #Calculate the seeing by fitting a Gaussian
            r1d=r[1].data.mean(axis=1)
            oblocs=Table.read('objectlocations.hdf5',path='Written 08-26')
            ind=int(np.where(oblocs['Data Set']==main.split('/')[0])[0])
            primloc=int(oblocs['Primary'][ind])
            primloc=int(np.where(r1d==r1d[primloc-3:primloc+3].max())[0])
            pmax=r1d[primloc]
            pleft=primloc-5
            pright=primloc+5
            prows=np.arange(pleft,pright+1)
            base=r1d[pleft:pright].min()
            halfmax=base+(0.5*(pmax-base))
            pylab.figure(); pylab.plot(r1d); pylab.suptitle('One-D Compression of Raw Image '+rawimg)
            pylab.xlabel('Binned Pixel Row'); pylab.ylabel('Mean Counts')
            pylab.plot(primloc,halfmax,'r*')
            pylab.vlines(pleft,base,pmax); pylab.vlines(pright,base,pmax)
            pylab.text(pleft,pmax+5,'P1',fontsize=9)
            pylab.text(pright,pmax+5,'P2',fontsize=9)
            primopt,primcov=curve_fit(gauss,prows,r1d[prows],p0=[base,pmax-base,primloc,1])
            psig=np.abs(primopt[-1])
            pfwhm=2.*np.sqrt(2.*np.log(2.))*psig
            pylab.plot(prows,gauss(prows,*primopt),'ro:')
            pylab.text(prows.mean(),pmax+5,'Gaussian Fit',rotation='vertical',fontsize=9)
            pylab.vlines(prows.mean()-(0.5*pfwhm),base-5,pmax+5); pylab.vlines(prows.mean()+(0.5*pfwhm),base-5,pmax+5)
            pylab.text(prows.mean(),base-2,'Prim FWHM',rotation='vertical',fontsize=9)

            # Interactively identify the Gaussian range
            conq=0
            while conq!=1:
                fixg=raw_input('Would you like to move one of the bounds for the Gaussian fit? (y/n, case sensitive) ')
                if fixg=='y':
                    conq=1
                    fcon=0
                    while fcon!=1:
                        rfix=raw_input('Which line would you like to move? (P1 or P2) ')
                        if rfix=='P1':
                            rcon=0
                            while rcon!=1:
                                rnew=input('What row should the line be moved to? ')
                                if rnew==int(rnew) and np.abs(primloc-rnew)<=20:
                                    rcon=1
                                    pleft=rnew
                                else:
                                    print 'Requires input of an integer within 20 rows of the primary location.'
                        elif rfix=='P2':
                            r2con=0
                            while r2con!=1:
                                r2new=input('What row should the line be moved to? ')
                                if r2new==int(r2new) and np.abs(primloc-r2new)<=20:
                                    r2con=1
                                    pright=r2new
                                else:
                                    print 'Requires input of an integer within 20 rows of the primary object.'
                        else:
                            print 'Invalid input. Codeloc2.  Try again.'
                        prows=np.arange(pleft,pright+1)
                        base=r1d[pleft:pright].min()
                        pylab.clf(); pylab.plot(r1d); pylab.suptitle('One-D Compression of Raw Image '+rawimg)
                        pylab.xlabel('Binned Pixel Row'); pylab.ylabel('Mean Counts')
                        pylab.plot(primloc,halfmax,'r*')
                        pylab.vlines(pleft,base,pmax); pylab.vlines(pright,base,pmax)
                        pylab.text(pleft,pmax+5,'P1',fontsize=9)
                        pylab.text(pright,pmax+5,'P2',fontsize=9)
                        primopt,primcov=curve_fit(gauss,prows,r1d[prows],p0=[base,pmax-base,primloc,1])
                        psig=np.abs(primopt[-1])
                        pfwhm=2.*np.sqrt(2.*np.log(2.))*psig
                        pylab.plot(prows,gauss(prows,*primopt),'ro:')
                        pylab.text(prows.mean(),pmax+5,'Gaussian Fit',rotation='vertical',fontsize=9)
                        pylab.vlines(prows.mean()-(0.5*pfwhm),base-5,pmax+5); pylab.vlines(prows.mean()+(0.5*pfwhm),base-5,pmax+5)
                        pylab.text(prows.mean(),base-2,'Prim FWHM',rotation='vertical',fontsize=9)
                        dcon=0
                        while dcon!=1:
                            domore=raw_input('Would you like to modify another primary line? (y/n, case sensitive) ')
                            if domore=='y':
                                dcon=1
                                fcon=0
                                rcon=1
                            elif domore=='n':
                                dcon=1
                                fcon=1
                                rcon=1
                                print 'Moving on to the secondary object'
                            else:
                                print 'Invalid input. Codeloc1.  Try again.'

                        
                elif fixg=='n':
                    conq=1
                    print 'The points will not be moved.'
                else:
                    print 'Invalid input. Codeloc3.  Try again.'
            
            
            # Calculate the seeing resultant from the Gaussian fit
            seeing=pfwhm #binned pix
            asseeing=seeing*spacebin*picscale #as
            seeing=asseeing*resolution #(A)
            # Prepare a table column to save the seeing
            pseeing=np.empty_like(skysubwaves)
            pseeing[:]=np.nan
            pseeing[0]=seeing
            pseeing[1]=asseeing
            pseecol=Column(data=pseeing,name='Prim Seeing (A,as)')
            
            #Again for secondary
            secloc=int(oblocs['Secondary'][ind])
            secloc=int(np.where(r1d==r1d[secloc-3:secloc+3].max())[0])
            sbase=r1d[secloc-10:secloc+10].min()
            smax=r1d[secloc]
            shalfmax=sbase+(0.5*(smax-sbase))
            sleft=secloc-5
            sright=secloc+5
            srows=np.arange(sleft,sright+1)
            secopt,seccov=curve_fit(gauss,srows,r1d[srows],p0=[sbase,smax-sbase,secloc,1])
            ssig=np.abs(secopt[-1])
            sfwhm=2.*np.sqrt(2.*np.log(2.))*ssig
            pylab.plot(secloc,shalfmax,'g*')
            pylab.plot(srows,gauss(srows,*secopt),'ro:')
            pylab.text(srows.mean(),smax+7,'Gaussian Fit, Secondary',rotation='vertical',fontsize=9)
            pylab.vlines(srows.mean()-(0.5*sfwhm),sbase-5,smax+5); pylab.vlines(srows.mean()+(0.5*sfwhm),sbase-5,smax+5)
            pylab.text(srows.mean(),sbase-2,'Sec FWHM',rotation='vertical',fontsize=9)
            pylab.vlines(sleft,sbase,smax); pylab.vlines(sright,sbase,smax)
            pylab.text(sleft,smax+5,'S1',fontsize=9)
            pylab.text(sright,smax+5,'S2',fontsize=9)

            # Interactive part of the Gaussian fit
            sconq=0
            while sconq!=1:
                sfixg=raw_input('Would you like to move one of the bounds for the secondary Gaussian fit? (y/n, case sensitive) ')
                if sfixg=='y':
                    sconq=1
                    sfcon=0
                    while sfcon!=1:
                        srfix=raw_input('Which line would you like to move? (S1 or S2) ')
                        if srfix=='S1':
                            srcon=0
                            while srcon!=1:
                                srnew=input('What row should the line be moved to? ')
                                if srnew==int(srnew) and np.abs(secloc-srnew)<=20:
                                    srcon=1
                                    sleft=srnew
                                else:
                                    print 'Requires input of an integer within 20 rows of the primary location.'
                        elif srfix=='S2':
                            sr2con=0
                            while sr2con!=1:
                                sr2new=input('What row should the line be moved to? ')
                                if sr2new==int(sr2new) and np.abs(secloc-sr2new)<=20:
                                    sr2con=1
                                    sright=sr2new
                                else:
                                    print 'Requires input of an integer within 20 rows of the primary object.'
                        else:
                            print 'Invalid input.  Codeloc 4.  Try again.'
                        srows=np.arange(sleft,sright+1)
                        secopt,seccov=curve_fit(gauss,srows,r1d[srows],p0=[sbase,smax-sbase,secloc,1])
                        ssig=np.abs(secopt[-1])
                        sfwhm=2.*np.sqrt(2.*np.log(2.))*ssig

                        # Plot the primary stuff again; don't want to lose that
                        pylab.clf(); pylab.plot(r1d); pylab.suptitle('One-D Compression of Raw Image '+rawimg)
                        pylab.xlabel('Binned Pixel Row'); pylab.ylabel('Mean Counts')
                        pylab.plot(primloc,halfmax,'r*')
                        pylab.vlines(pleft,base,pmax); pylab.vlines(pright,base,pmax)
                        pylab.text(pleft,pmax+5,'P1',fontsize=9)
                        pylab.text(pright,pmax+5,'P2',fontsize=9)
                        primopt,primcov=curve_fit(gauss,prows,r1d[prows],p0=[base,pmax-base,primloc,1])
                        psig=np.abs(primopt[-1])
                        pfwhm=2.*np.sqrt(2.*np.log(2.))*psig
                        pylab.plot(prows,gauss(prows,*primopt),'ro:')
                        pylab.text(prows.mean(),pmax+5,'Gaussian Fit',rotation='vertical',fontsize=9)
                        pylab.vlines(prows.mean()-(0.5*pfwhm),base-5,pmax+5)
                        pylab.vlines(prows.mean()+(0.5*pfwhm),base-5,pmax+5)
                        pylab.text(prows.mean(),base-2,'Prim FWHM',rotation='vertical',fontsize=9)

                        # Plot the new secondary things
                        pylab.plot(secloc,shalfmax,'g*')
                        pylab.plot(srows,gauss(srows,*secopt),'ro:')
                        pylab.text(srows.mean(),smax+7,'Gaussian Fit, Secondary',rotation='vertical',fontsize=9)
                        pylab.vlines(srows.mean()-(0.5*sfwhm),sbase-5,smax+5)
                        pylab.vlines(srows.mean()+(0.5*sfwhm),sbase-5,smax+5)
                        pylab.text(srows.mean(),sbase-2,'Sec FWHM',rotation='vertical',fontsize=9)
                        pylab.vlines(sleft,sbase,smax); pylab.vlines(sright,sbase,smax)
                        pylab.text(sleft,smax+5,'S1',fontsize=9)
                        pylab.text(sright,smax+5,'S2',fontsize=9)
                        sdcon=0
                        while sdcon!=1:
                            sdomore=raw_input('Would you like to modify another secondary line? (y/n, case sensitive) ')
                            if sdomore=='y':
                                sdcon=1
                                sfcon=0
                            elif sdomore=='n':
                                sdcon=1
                                sfcon=1
                                print 'Moving on, then.'
                            else:
                                print 'Invalid input.  Codeloc 5.  Try again.'

                elif sfixg=='n':
                    sconq=1
                    print 'The points will not be moved.'
                else:
                    print 'Invalid input.  Codeloc 6.  Try again.'
            
            # Calculate the seeing resultant from these parameters
            secseeing=sfwhm #binned pix
            assecseeing=secseeing*spacebin*picscale #as
            secseeing=assecseeing*resolution #A
            # Prepare a table column to save the seeing
            secsee=np.empty_like(skysubwaves)
            secsee[:]=np.nan
            secsee[0]=secseeing
            secsee[1]=assecseeing
            secseecol=Column(data=pseeing,name='Sec Seeing (A,as)')

            #Establish the slit width
            maskid=str(w[0].header['MASKID'])
            if maskid=='PL0100N002':
                sw=1.0 #(as)
            elif maskid=='PL0150N001':
                sw=1.5 #(as)
            else:
                sys.exit('Unrecognized slit MASKID '+maskid)
            sw=sw*resolution #(A)

            #Establish the value for overresolve
            grating=str(w[0].header['GRATING'])
            if grating=='PG2300':
                R=2850
            elif grating=='PG0900':
                R=1000
            else:
                sys.exit('Unrecognized grating '+grating) 
            overresolve=np.around((skysubwaves.mean())/(R*wgres))
            
            #Execute the "observe"
            #For primary:
            print 'Parameters going in to "observe" for the primary object:'
            print(' {0:^20} {1:^20} {2:^20} {3:^20}'.format('Slit Width (A)','Seeing (as)','Seeing (A)','Overresolve'))
            print(' {0:^20.5f} {1:^20.5f} {2:^20.5f} {3:^20}'.format(sw,asseeing,seeing,overresolve))
            outtable1,mconv,filt=observe(model,skysubwaves,sw,seeing,overresolve)
            #For secondary:
            print 'Parameters going in to "observe" for the secondary object:'
            print(' {0:^20} {1:^20} {2:^20} {3:^20}'.format('Slit Width (A)','Seeing (as)','Seeing (A)','Overresolve'))
            print(' {0:^20.5f} {1:^20.5f} {2:^20.5f} {3:^20}'.format(sw,assecseeing,secseeing,overresolve))
            souttable1,smconv,sfilt=observe(model,skysubwaves,sw,secseeing,overresolve)

            pylab.figure(); pylab.plot(filt); pylab.suptitle('Convolution Kernel for Primary, Seeing = '+str(seeing))
            pylab.figure(); pylab.plot(sfilt); pylab.suptitle('Convolution Kernel for Secondary, Seeing = '+str(secseeing))

            #Fit a polynomial to the data divided by the model
            #Primary:
            mdiv=skysubprim/mconv
            #Secondary:
            smdiv=skysubsec/smconv
            
            #Prompt the user to input the fit order
            con2=0
            while con2!=1:
                fitorder=input('Please indicate the order of the polyfit to be used in scaling: ')
                if fitorder==int(fitorder):
                    con2=1
                else:
                    print 'Input must be an integer.  Try again.'
            
            #Fit the data divided by the model
            # Primary:
            mfit=np.ma.polyfit(skysubwaves[-np.isnan(mdiv)],mdiv[-np.isnan(mdiv)],fitorder)
            normmod=np.polyval(mfit,skysubwaves)
            # Secondary:
            smfit=np.ma.polyfit(skysubwaves[-np.isnan(smdiv)],smdiv[-np.isnan(smdiv)],fitorder)
            snormmod=np.polyval(smfit,skysubwaves)

            #Multiply the fit by the model
            # Primary:
            mod=normmod*mconv
            # Secondary:
            smod=snormmod*smconv

            #View the now-scaled primary model
            pylab.figure(); pylab.subplot(2,1,1)
            pylab.plot(skysubwaves,mconv,label='Convolved & Interpolated Model')
            pylab.suptitle('Checking Model (Primary Object)')
            pylab.plot(skysubwaves,skysubprim,label='Extracted Data for Primary Object')
            pylab.xlabel('Wavelength (A)')
            pylab.legend(loc='best')
            pylab.subplot(2,1,2); pylab.plot(skysubwaves,skysubprim,label='Extracted Data for Primary Object')
            pylab.plot(skysubwaves,mod,label='Normalized & Fitted Model')
            pylab.suptitle('Checking Scaled Model (Primary Object)')
            pylab.xlabel('Wavelength (A)'); pylab.legend(loc='best')

            #Prompt to save the image
            figfile=prodir+datno+'_primscaledmodel.png'
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
            
            #View the now-scaled secondary model
            pylab.figure(); pylab.subplot(2,1,1); pylab.plot(skysubwaves,smconv,label='Convolved & Interpolated Model')
            pylab.suptitle('Checking Model (Secondary Object)')
            pylab.plot(skysubwaves,skysubsec,label='Extracted Data for Secondary Object'); pylab.xlabel('Wavelength (A)')
            pylab.legend(loc='best')
            pylab.subplot(2,1,2); pylab.plot(skysubwaves,skysubsec,label='Extracted Data for Secondary Object')
            pylab.plot(skysubwaves,smod,label='Normalized & Fitted Model')
            pylab.suptitle('Checking Scaled Model (Secondary Object)')
            pylab.xlabel('Wavelength (A)'); pylab.legend(loc='best')

            #Prompt to save this figure as well
            sfigfile=prodir+datno+'_secscaledmodel.png'
            scon3=0
            while scon3!=1:
                simsave=raw_input('Would you like to save the image of the secondary scaled model? (y/n, case sensitive) ')
                if simsave=='y':
                    scon3=1
                    swait=raw_input('Press <return> when ready to save the secondary model.')
                    if os.path.exists(sfigfile):
                        print 'Warning: File %s already exists.' % sfigfile
                        scon4=0
                        while scon4!=1:
                            sow=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
                            if sow=='y':
                                scon4=1
                                os.remove(sfigfile)
                                pylab.savefig(sfigfile)
                            elif sow=='n':
                                scon4=1
                                print 'File %s will not be replaced.' % sfigfile
                            else:
                                print 'Invalid input.  Try again.'
                    else:
                        print 'Writing file %s' % sfigfile
                        pylab.savefig(sfigfile)
                elif simsave=='n':
                    scon3=1
                    print 'The image will not be saved.'
                else:
                    print 'Invalid input.  Try again.'

            #Create a table of the model to be written to an hdf5 file
            outtable2=Table([skysubwaves,mod,smod],names=('Wavelength','Prim Scaled Model',
                                                          'Sec Scaled Model'),meta={'name':'Model Data'})
            outtable2.add_column(pseecol,index=2)
            outtable2.add_column(secseecol)
            outfile=prodir+datno+'_model.hdf5'
            # Add polyfit coefficients to the table (for normalization of composite spectra)
            mfitcol=np.empty_like(skysubwaves)
            mfitcol[:]=np.nan
            for q in range(len(mfit)):
                mfitcol[q]=mfit[q]
            smfitcol=np.empty_like(skysubwaves)
            smfitcol[:]=np.nan
            for l in range(len(smfit)):
                smfitcol[l]=smfit[l]
            polycol=Column(data=mfitcol,name='Prim Poly Coeff')
            spolycol=Column(data=smfitcol,name='Sec Poly Coeff')
            outtable2.add_column(polycol,index=3)
            outtable2.add_column(spolycol)
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
            if grating=='PG2300':
                ind1=np.argmin(np.abs(skysubwaves-4500))
                ind2=np.argmin(np.abs(skysubwaves-4900))
            elif grating=='PG0900':
                ind1=np.argmin(np.abs(skysubwaves-4800))
                ind2=np.argmin(np.abs(skysubwaves-5700))
            c=ind1
            d=ind2
            waves=skysubwaves[c:d]
            fprim=skysubprim[c:d]
            fsec=skysubsec[c:d]
            ch=np.zeros_like(v)
            sch=np.zeros_like(v)
            for j in range(len(v)):
                #Calculate chi**2 for an array of velocities
                q=interpolate.InterpolatedUnivariateSpline(skysubwaves+(skysubwaves*v[j]/(2.9979E5)),mod)
                sq=interpolate.InterpolatedUnivariateSpline(skysubwaves+(skysubwaves*v[j]/(2.9979E5)),smod)
                zz=np.polyfit(waves[-np.isnan(fprim)],(fprim[-np.isnan(fprim)]/q(waves[-np.isnan(fprim)])),fitorder)
                szz=np.polyfit(waves[-np.isnan(fsec)],(fsec[-np.isnan(fsec)]/sq(waves[-np.isnan(fsec)])),fitorder)
                if np.isnan(zz).any() or np.isnan(szz).any():
                    sys.exit('Polyfit returned NaN values.  Polyfit coefficients are: '+str(zz)+
                             ' for primary and: '+str(szz)+' for secondary.')
                pp=np.poly1d(zz)
                spp=np.poly1d(szz)
                res=np.zeros_like(waves)
                sres=np.zeros_like(waves)
                for m in range(c,d):
                    if skysubprim[m]>=0 and pp(skysubwaves[m])>=0 and skysubprimerr[m]>0:
                        res[m-c]=float((((skysubprim[m]-q(skysubwaves[m])*pp(skysubwaves[m]))/(skysubprimerr[m]))**2))
                    if skysubsec[m]>=0 and spp(skysubwaves[m])>=0 and skysubsecerr[m]>0:
                        sres[m-c]=float((((skysubsec[m]-sq(skysubwaves[m])*spp(skysubwaves[m]))/(skysubsecerr[m]))**2))
                ch[j]=float(res.sum())
                sch[j]=float(sres.sum())
            pylab.figure(); pylab.plot(v,ch); pylab.suptitle('Chi**2 Results for Primary, File '+dat)
            pylab.xlabel('Velocity (km/s)'); pylab.ylabel('Chi**2')
            pchifile=prodir+'chi2'+datno+'_'+month+day+year+'.png'
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
            pylab.figure(); pylab.plot(v,sch); pylab.suptitle('Chi**2 Results for Secondary, File '+dat)
            pylab.xlabel('Velocity (km/s)'); pylab.ylabel('Chi**2')
            schifile=prodir+'chi2'+datno+'secondary_'+month+day+year+'.png'
            if os.path.exists(schifile):
                print 'Warning: File name %s already exists.' % schifile
                schcon=0
                while schcon!=1:
                    schwr=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
                    if schwr=='y':
                        schcon=1
                        os.remove(schifile)
                        pylab.savefig(schifile)
                    elif schwr=='n':
                        schcon=1
                        print 'Primary Chi**2 figure will not be saved.'
                    else:
                        print 'Invalid input.  Try again.'
            else:
                pylab.savefig(schifile)
            
            #Fit a parabola around the minimum of the chi**2 array
            ymin1=v[np.argmin(ch)]
            symin1=v[np.argmin(sch)]
            para=np.polyfit(v[np.argmin(ch)-6:np.argmin(ch)+7],ch[np.argmin(ch)-6:np.argmin(ch)+7],2)
            spara=np.polyfit(v[np.argmin(sch)-6:np.argmin(sch)+7],sch[np.argmin(sch)-6:np.argmin(sch)+7],2)
            print('The velocity corresponding to the minimum chi**2 value for the primary object from file '
                  +dat+' is %s' % ymin1)
            print('The velocity corresponding to the minimum chi**2 value for the secondary object from file '+dat+
                  ' is %s' % symin1)
            xx=np.arange(v[np.argmin(ch)-6],v[np.argmin(ch)+6],1)
            sxx=np.arange(v[np.argmin(sch)-6],v[np.argmin(sch)+6],1)
            cc=para[0]*xx**2+para[1]*xx+para[2]
            scc=spara[0]*sxx**2+spara[1]*sxx+spara[2]
            ymin=xx[np.argmin(cc)]
            symin=sxx[np.argmin(scc)]

            #Plot results and calculate uncertainty
            pylab.figure(); pylab.plot(xx,cc)
            pylab.suptitle('Parabolic Fit of Chi**2 Function for Primary Object from File '+dat)
            pylab.xlabel('Velocity (km/s)'); pylab.ylabel('Chi**2')
            pylab.plot(xx,[cc[np.argmin(cc)]+1]*len(xx))
            unc=np.abs((xx[np.argmin(np.abs(cc-(cc[np.argmin(cc)]+1)))])-(xx[np.argmin(cc)]))
            print 'The radial velocity for the primary object in the '+datno+' observation is '+str(ymin)+' +/- %s km/s' % unc

            #Again for secondary
            pylab.figure(); pylab.plot(sxx,scc)
            pylab.suptitle('Parabolic Fit of Chi**2 Funciton for Secondary Object from File '+dat)
            pylab.xlabel('Velocity (km/s)'); pylab.ylabel('Chi**2')
            pylab.plot(sxx,[scc[np.argmin(scc)]+1]*len(sxx))
            sunc=np.abs((sxx[np.argmin(np.abs(scc-(scc[np.argmin(scc)]+1)))])-(sxx[np.argmin(scc)]))
            print('The radial velocity for the secondary object in the '+datno+' observation is '
                  +str(symin)+' +/- %s km/s' % sunc)
            
            #Analyze the shift at this velocity
            shiftwaves=skysubwaves+(skysubwaves*(ymin/2.9979E5))
            shiftwaves2=skysubwaves+(skysubwaves*(symin/2.9979E5))
            primnotmask=skysubprim[-skysubprim.mask]
            secnotmask=skysubsec[-skysubsec.mask]
            corrwaves=skysubwaves[-skysubprim.mask]
            scorrwaves=skysubwaves[-skysubsec.mask]
            correrr=skysubprimerr[-skysubprim.mask]
            scorrerr=skysubsecerr[-skysubsec.mask]
            q2=interpolate.InterpolatedUnivariateSpline(shiftwaves[-skysubprim.mask],primnotmask)
            sq2=interpolate.InterpolatedUnivariateSpline(shiftwaves2[-skysubsec.mask],secnotmask)
            if np.isnan(q2(skysubwaves[-np.isnan(skysubprim)])).any():
                sys.exit('Interpolation of shifted data returned NaN values.')
            if np.isnan(sq2(skysubwaves[-np.isnan(skysubsec)])).any():
                sys.exit('Interpolation for secondary object of shifted data returned NaN values.')
            zz2=np.polyfit(corrwaves,primnotmask/q2(corrwaves),fitorder)
            szz2=np.polyfit(scorrwaves,secnotmask/sq2(scorrwaves),fitorder)
            if np.isnan(zz2).any():
                sys.exit('Polyfit of shifted array returned NaN values.  Polyfit coefficients are: '+str(zz2))
            if np.isnan(szz2).any():
                sys.exit('Polyfit for secondary object of shifted array returned NaN values.  Polyfit coefficients are: '+str(szz2))
            pp2=np.poly1d(zz2)
            spp2=np.poly1d(szz2)
            res2=np.zeros_like(corrwaves)
            sres2=np.zeros_like(scorrwaves)
            for n in range(len(corrwaves)):
                res2[n]=((primnotmask[n]-q2(corrwaves[n])*pp2(corrwaves[n]))/(correrr[n]))
            for sn in range(len(scorrwaves)):
                sres2[sn]=((secnotmask[sn]-sq2(scorrwaves[sn])*spp2(scorrwaves[sn]))/(scorrerr[sn]))
            ch2=np.sum(res2**2)
            sch2=np.sum(sres2**2)
            pylab.figure(); pylab.subplot(2,1,1); pylab.plot(skysubwaves,skysubprim,label='Extracted Data')
            pylab.plot(shiftwaves,mod,label='Shifted, Fitted & Scaled Model')
            pylab.xlabel('Wavelength (A)'); pylab.ylabel('Counts')
            pylab.suptitle('Residuals for Shifted Primary Model, Obs. Date '+datno+', Velocity = '+str(ymin))
            pylab.legend(loc='best')
            pylab.subplot(2,1,2); pylab.plot(corrwaves,res2,'.'); pylab.xlabel('Wavelength (A)'); pylab.ylabel('Residuals')
            print 'The value of Chi**2 for the velocity of the primary object is %s ' % ch2
            
            #Ding!
            print('\a')

            #Prompt to save the residuals figure for the primary object
            presfile=prodir+'vshift'+datno+'.png'
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

            #Plot the residuals for the secondary object
            pylab.figure(); pylab.subplot(2,1,1); pylab.plot(skysubwaves,skysubsec,label='Shifted Secondary Data')
            pylab.plot(shiftwaves2,smod,label='Shifted, Fitted & Scaled Secondary Model')
            pylab.xlabel('Wavelength (A)'); pylab.ylabel('Counts')
            pylab.suptitle('Residuals for Shifted Secondary Model, Obs. Date '+datno+', Velocity = '+str(symin))
            pylab.legend(loc='best')
            pylab.subplot(2,1,2); pylab.plot(scorrwaves,sres2,'.'); pylab.xlabel('Wavelength (A)'); pylab.ylabel('Residuals')
            print 'The value of Chi**2 for the velocity of the secondary object is %s ' % sch2

            #Prompt to save the residuals figure for the primary object
            sresfile=prodir+'vshiftsec'+datno+'.png'
            con7=0
            while con7!=1:
                sressave=raw_input('Would you like to save the image of the secondary residuals? (y/n, case sensitive) ')
                if sressave=='y':
                    con7=1
                    sreswait=raw_input('Press <return> when ready to save.')
                    if os.path.exists(sresfile):
                        print 'Warning: File %s already exists.' % sresfile
                        con8=0
                        while con8!=1:
                            srow=raw_input('Would you like to overwrite this file? (y/n, case sensitive) ')
                            if srow=='y':
                                con8=1
                                os.remove(sresfile)
                                pylab.savefig(sresfile)
                            elif srow=='n':
                                con8=1
                                print 'File %s will not be replaced.' % sresfile
                            else:
                                print 'Invalid input.  Try again.'
                    else:
                        print 'Writing file %s' % sresfile
                        pylab.savefig(sresfile)
                elif sressave=='n':
                    con7=1
                    print 'The image will not be saved.'
                else:
                    print 'Invalid input.  Try again.'
            
            #Write results to hdf5 table for all observations
            print 'Writing results to rvtable.hdf5:'
            ss=pyfits.open(ssimg)
            hjd=ss[0].header['JD']
            rvdata=[datno,hjd,ymin,unc,symin,sunc]
            rvtable='rvtable.hdf5'
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
            print rv
            rv.write(rvtable,path='OVA',overwrite=True)


    else:
        print 'Path %s does not exist.  Try again.' % prodir
