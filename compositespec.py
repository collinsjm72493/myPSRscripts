"""
compositespec.py
Started on Aug. 21, 2014
Justin Collins

This script takes in all of the extracted spectra from the the SALT and Gemini observations
of the pulsar PSR J1723-2837 and creates a weighted composite spectrum for all the
observations from each of the two telescopes (i.e. two composite spectra).
"""

import numpy as np
from numpy import genfromtxt
from astropy.table import Table,Column
import os,sys,glob
import pyfits
import re
import pylab
import scipy
from scipy.interpolate import interp1d

pylab.ion()
pylab.close('all')

# Bring in the RV Data from the two hdf5 tables
saltrvtable='salt_reduced/rvtable2.hdf5'
saltdat=Table.read(saltrvtable,path='OVA')
gemrvtable='gemini_psrj1723/Gemini/rvtable_gem.hdf5'
gemdat=Table.read(gemrvtable,path='OVA')
rvgem=gemdat['RV']
rvsalt=saltdat['Prim RV']


# Gather all the extracted spectra, in their various data forms
saltfiles=[]
saltdirs=os.listdir('./salt_reduced/')
for i in saltdirs:
    if os.path.isdir('./salt_reduced/'+i):
        if re.search('819',i)==None and re.search('819mod',i)==None and \
           re.search('819test',i)==None and re.search('909',i)==None and re.search('909test',i)==None:
            try:
                saltextr=glob.glob('./salt_reduced/'+i+'/product/*_extrcol1.txt')
                for n in range(len(saltextr)):
                    saltfiles.append(saltextr[n])
            except:
                pass
#print 'The SALT files to be composited are:'
#print saltfiles
#print '\n'

if len(saltfiles)!=len(saltdat):
    sys.exit('Not all of the extracted SALT spectra are in the corresponding hdf5 RV table: '+saltrvtable)

gemfiles=glob.glob('gemini_psrj1723/Gemini/cetgsg*.fits')
#print 'The Gemini files to be composited are: '
#print gemfiles
#print '\n'

if len(gemfiles)!=len(gemdat):
    sys.exit('Not all of the extracted Gemini spectra are in the corresponding hdf5 RV table: '+gemrvtable)


# Build a "master" wavelength array for the SALT data
saltmwaves=np.arange(4200,5200.01,0.01)

# Reference spectral lines, to check alignment
refspeclines=[4226.728,4341.42,4383.5447,4250.7866,4260,4741,4254.35,4271.760,4668,4404.5,4861.35,5167.3,5172.7,5183.6]

# Prompt for plotting
con1=0
while con1!=1:
    plotsalt=raw_input('Would you like to plot the stacked SALT spectra? (y/n, case sensitive) ')
    if plotsalt=='y':
        con1=1
        print 'Stacked SALT spectra will be plotted.'
    elif plotsalt=='n':
        con1=1
        print 'Stacked SALT spectra will NOT be plotted.'
    else:
        print 'Invalid input.  Try again.'

# Shift the SALT spectra by the calculated velocity
if plotsalt=='y':
    pylab.figure(num=1)
    pylab.suptitle('Shifted SALT Spectra')
    pylab.xlabel('Wavelength (A)')

    pylab.figure(num=2)
    pylab.suptitle('Shifted SALT Spectra')
    pylab.xlabel('Wavelength (A)')

    pylab.figure(num=3)
    pylab.suptitle('Shifted SALT Spectra')
    pylab.xlabel('Wavelength (A)')
    for k in range(len(saltfiles)):
        obsno=saltfiles[k].split('/')[-1].split('_')[0]
        print 'Current observation number: '+obsno
        for n in range(len(saltdat)):
            if saltdat['Observation'][n]==obsno:
                ind=n
        saltprv=saltdat['Prim RV'][ind]
        saltpdrv=saltdat['Prim Unc'][ind]
        saltobs=saltdat['Observation'][ind]
        saltgmjd=saltdat['GMJD'][ind]
        print 'The corresponding RV data are:'
        print '{0:^15} {1:^15} {2:^10} {3:^10}'.format('Observation','GMJD','Prim RV','Prim Unc')
        print '{0:^15} {1:^15.5f} {2:^10.1f} {3:^10.1f}'.format(saltobs,saltgmjd,saltprv,saltpdrv)

        saltexdat=genfromtxt(saltfiles[k])
        saltwaves=np.array(saltexdat[:,1])
        saltpcts=np.ma.array(saltexdat[:,3],mask=(saltexdat[:,3]<=0))
        saltshiftwaves=saltwaves-(saltwaves*(saltprv/(2.9979E5)))
    
        if re.search('512',obsno)!=None or re.search('518',obsno)!=None or re.search('520',obsno)!=None:
            pylab.figure(num=1)
            pylab.plot(saltshiftwaves,saltpcts+(100*k),label=obsno)
            pylab.text(saltwaves.max()+10,np.mean(saltpcts[-np.isnan(saltpcts)])+(100*k),obsno,fontsize=9)
            pylab.text(saltwaves.min()-250,np.mean(saltpcts[-np.isnan(saltpcts)])+(100*k)-150,
                       'Velocity: '+str(saltprv)+' km/s',fontsize=9)
        elif re.search('712',obsno)!=None:
            pylab.figure(num=2)
            pylab.plot(saltshiftwaves,saltpcts+(100*k),label=obsno)
            pylab.text(saltwaves.max()+10,np.mean(saltpcts[-np.isnan(saltpcts)])+(100*k),obsno,fontsize=9)
            pylab.text(saltwaves.min()-250,np.mean(saltpcts[-np.isnan(saltpcts)])+(100*k)-50,
                       'Velocity: '+str(saltprv)+' km/s',fontsize=9)
        else:
            pylab.figure(num=3)
            pylab.plot(saltshiftwaves,saltpcts+(100*k),label=obsno)
            pylab.text(saltwaves.max()+10,np.mean(saltpcts[-np.isnan(saltpcts)])+(100*k),obsno,fontsize=9)
            pylab.text(saltwaves.min()-250,np.mean(saltpcts[-np.isnan(saltpcts)])+(100*k)-450,
                       'Velocity: '+str(saltprv)+' km/s',fontsize=9)

    # Plot the G5 star model here as well (one of the scaled models)
    for n in range(1,4):
        if n==1:
            modelfile='salt_reduced/512a/product/P201305120029_model.hdf5'
            model=Table.read(modelfile,path='Written081114')
        elif n==2:
            modelfile='salt_reduced/712/product/P201307120021_model.hdf5'
            model=Table.read(modelfile,path='Written081214')
        elif n==3:
            modelfile='salt_reduced/918/product/P201309180006_model.hdf5'
            model=Table.read(modelfile,path='Written081314')
        
        modwaves=model['Wavelength']
        modflux=model['Prim Scaled Model']

        wave1=saltwaves.min()
        wave2=saltwaves.max()
        ind1=np.argmin(np.abs(modwaves-wave1))
        ind2=np.argmin(np.abs(modwaves-wave2))

        pylab.figure(num=n)
        if n==1:
            pylab.plot(modwaves[ind1:ind2],modflux[ind1:ind2]+(1000))
            pylab.text(modwaves[ind1:ind2].max()+10,np.mean(modflux[ind1:ind2])+(1000),'G5 Star Model',fontsize=9)
        elif n==2:
            pylab.plot(modwaves[ind1:ind2],modflux[ind1:ind2]+2800)
            pylab.text(modwaves[ind1:ind2].max()+10,np.mean(modflux[ind1:ind2])+(2800),'G5 Star Model',fontsize=9)
        elif n==3:
            pylab.plot(modwaves[ind1:ind2],modflux[ind1:ind2]+3500)
            pylab.text(modwaves[ind1:ind2].max()+10,np.mean(modflux[ind1:ind2])+(3500),'G5 Star Model',fontsize=9)

    # Plot some reference spectral lines
    for x in refspeclines:
        pylab.figure(num=1)
        pylab.vlines(x,0,1350)
        pylab.text(x-5,1450,str(x),rotation='vertical',fontsize=9)
        
        pylab.figure(num=2)
        pylab.vlines(x,2300,2900)
        pylab.text(x-5,2950,str(x),rotation='vertical',fontsize=9)

        pylab.figure(num=3)
        pylab.vlines(x,1000,4250)
        pylab.text(x-5,4600,str(x),rotation='vertical',fontsize=9)

# Now again for the Gemini data
# Prompt first
con2=0
while con2!=1:
    plotgem=raw_input('Would you like to plot the stacked Gemini spectra? (y/n, case sensitive) ')
    if plotgem=='y':
        con2=1
        print 'The stacked Gemini spectra will be plotted.'
    elif plotgem=='n':
        con2=1
        print 'The stacked Gemini spectra will NOT be plotted.'
    else:
        print 'Invalid input.  Try again.'

# If prompt was yes, plot the spectra
if plotgem=='y':
    pylab.figure(num=4); pylab.suptitle('Shifted Gemini Spectra'); pylab.xlabel('Wavelength (A)')
    for c in range(len(gemfiles)):
        gobsno=gemfiles[c].split('/')[-1].split('20130')[1].split('.fits')[0]
        for g in range(len(gemdat)):
            if gemdat['Observation'][g]==gobsno:
                gind=g
    
        gemprv=gemdat['RV'][gind]
        gempdrv=gemdat['ERR'][gind]
        gemobs=gemdat['Observation'][gind]
        gemgmjd=gemdat['GMJD'][gind]
    
        gemfits=pyfits.open(gemfiles[c])
        gemexdat=gemfits['SCI'].data
        gemexdat=np.ma.array(gemexdat[0],mask=(gemexdat[0]==0))
        gemwave1=gemfits['SCI'].header['CRVAL1']
        gemdwave=gemfits['SCI'].header['CD1_1']
        gemwaves=[]
        for l in range(len(gemexdat)):
            gemwaves.append(gemwave1+(gemdwave*l))
        if len(gemwaves)!=len(gemexdat):
            sys.exit('Length of Gemini wavelength array does not match length of array of extracted data.')
        gemwaves=np.array(gemwaves)
        gemshiftwaves=gemwaves-(gemwaves*(gemprv/(2.9979E5)))

        if c==0:
            low=gemexdat.min()
        if c==len(gemfiles)-1:
            high=((1E-14)*c)+gemexdat.max()

        pylab.plot(gemshiftwaves,gemexdat+((5E-15)*c))
        pylab.text(gemshiftwaves.max()+10,np.mean(gemexdat)+((5E-15)*c),gobsno,fontsize=9)
        pylab.text(gemshiftwaves.min()-10,np.mean(gemexdat)+((5E-15)*c),'Velocity: '+str(gemprv)+' km/s',fontsize=9)

    for x in refspeclines:
        pylab.figure(num=4)
        pylab.vlines(x,low,high-(7E-14))
        pylab.text(x-5,high-(5.5E-14),str(x),rotation='vertical',fontsize=9)

    # Plot the G5 star model here as well (one of the scaled models)
    gemmodelfile='gemini_psrj1723/Gemini/810S0056_model.hdf5'
    gemmodel=Table.read(gemmodelfile,path='Written082014')

    gemmodwaves=gemmodel['Wavelength']
    gemmodflux=gemmodel['Prim Scaled Model']

    gwave1=gemwaves.min()
    gwave2=gemwaves.max()
    gind1=np.argmin(np.abs(gemmodwaves-gwave1))
    gind2=np.argmin(np.abs(gemmodwaves-gwave2))

    pylab.figure(num=4)
    pylab.plot(gemmodwaves[gind1:gind2],gemmodflux[gind1:gind2]+(9E-14))
    pylab.text(gemmodwaves[gind1:gind2].max()+10,np.mean(gemmodflux[gind1:gind2])+(9E-14),'G5 Star Model',fontsize=9)

##############################################

# Prompt to plot all spectra together
con3=0
while con3!=1:
    plotall=raw_input('Would you like to plot all of the stacked spectra? (y/n, case sensitive) ')
    if plotall=='y':
        con3=1
        print 'All spectra will be plotted.'

        # Plot all spectra together
        step=1500
        pylab.figure(num=5); pylab.suptitle('All Gemini & SALT Spectra'); pylab.xlabel('Wavelength (A)')
        for q in range(len(saltfiles)):
            obsno=saltfiles[q].split('/')[-1].split('P20130')[1].split('_')[0]
            for n in range(len(saltdat)):
                if saltdat['Observation'][n]=='P20130'+obsno:
                    sind=n
            saltprv=saltdat['Prim RV'][sind]
            saltpdrv=saltdat['Prim Unc'][sind]

            saltexdat=genfromtxt(saltfiles[q])
            saltwaves=np.array(saltexdat[:,1])
            saltpcts=np.ma.array(saltexdat[:,3],mask=(saltexdat[:,3]<=0))
            saltshiftwaves=saltwaves-(saltwaves*(saltprv/(2.9979E5)))

            if q==0:
                low2=np.min(saltprv)-50
    
            if re.search('512',obsno)!=None or re.search('518',obsno)!=None or re.search('520',obsno)!=None:
                scale=0.125
            elif re.search('712',obsno)!=None:
                scale=0.0625
            else:
                scale=0.25
            pylab.plot(saltshiftwaves,(saltpcts/scale)+(step*q),label=obsno)
            pylab.text(saltwaves.max()+10,(np.mean(saltpcts[-np.isnan(saltpcts)])/scale)+(step*q),'SALT'+obsno,fontsize=9)
            pylab.text(saltwaves.min()-150,(np.mean(saltpcts[-np.isnan(saltpcts)])/scale)+(step*q)-150,
                       'Velocity: '+str(saltprv)+' km/s',fontsize=9)

        # Now plot the model
        q=q+2
        modelfile='salt_reduced/512a/product/P201305120029_model.hdf5'
        model=Table.read(modelfile,path='Written081114')
        modwaves=model['Wavelength']
        modflux=model['Prim Scaled Model']

        modscale=0.125
        wave1=saltwaves.min()
        wave2=saltwaves.max()
        ind1=np.argmin(np.abs(modwaves-wave1))
        ind2=np.argmin(np.abs(modwaves-wave2))
        pylab.plot(modwaves[ind1:ind2],(modflux[ind1:ind2]/modscale)+(step*q),color='k')
        pylab.text(modwaves[ind1:ind2].max()+10,np.mean((modflux[ind1:ind2]/modscale))+(step*q),
                   'G5 Star Model',fontsize=14)

        # Now plot the Gemini data
        q=q+1
        for c in range(len(gemfiles)):
            q=q+1
            gobsno=gemfiles[c].split('/')[-1].split('20130')[1].split('.fits')[0]
            for g in range(len(gemdat)):
                if gemdat['Observation'][g]==gobsno:
                    gind=g
    
            gemprv=gemdat['RV'][gind]
            gempdrv=gemdat['ERR'][gind]
            gemobs=gemdat['Observation'][gind]
            gemgmjd=gemdat['GMJD'][gind]
    
            gemfits=pyfits.open(gemfiles[c])
            gemexdat=gemfits['SCI'].data
            gemexdat=np.ma.array(gemexdat[0],mask=(gemexdat[0]==0))
            gemwave1=gemfits['SCI'].header['CRVAL1']
            gemdwave=gemfits['SCI'].header['CD1_1']
            gemwaves=[]
            for l in range(len(gemexdat)):
                gemwaves.append(gemwave1+(gemdwave*l))
            if len(gemwaves)!=len(gemexdat):
                sys.exit('Length of Gemini wavelength array does not match length of array of extracted data.')
            gemwaves=np.array(gemwaves)
            gemshiftwaves=gemwaves-(gemwaves*(gemprv/(2.9979E5)))

            if c==len(gemfiles)-1:
                high2=np.max(gemexdat/(1.5E-17))+(step*q)

            if re.search('712',gobsno)!=None:
                if re.search('712S0130',gobsno)!=None:
                    gscale=(2E-19)
                else:
                    gscale=(5E-19)
            elif re.search('613',gobsno)!=None:
                gscale=(1E-18)
            else:
                gscale=(2E-18)
            pylab.plot(gemshiftwaves,(gemexdat/gscale)+(step*q))
            pylab.text(gemshiftwaves.max()+10,(np.mean(gemexdat)/gscale)+(step*q),'GEM'+gobsno,fontsize=9)
            pylab.text(gemshiftwaves.min()-10,(np.mean(gemexdat)/gscale)+(step*q),
                           'Velocity: '+str(gemprv)+' km/s',fontsize=9)

        # Now plot the reference spectral lines
        for x in refspeclines:
            pylab.vlines(x,low2,high2,linestyles='dashed')
            pylab.text(x-5,high2+(2.5*step),str(x),rotation='vertical',fontsize=9)

    elif plotall=='n':
        con3=1
        print 'The overall stacked spectra will not be plotted.'
    else:
        print 'Invalid input.  Try again.'

###################################################################

# Make the composite spectra

# Start with SALT

print '------------------------------------------------------------------------------'

# Shift to "best" spectrum; for SALT, I prefer the 5120029 observation

# Pull data from the 5120029 observation
tx1=genfromtxt(saltfiles[0]) # Data contained in the noted .txt file
waves1=tx1[:,1] # wavelengths
# Chop off the ends of the wavelength array to avoid problems with interpolation
w1=np.argmin(np.abs(waves1-4200))
w2=np.argmin(np.abs(waves1-5200))
waves1=waves1[w1:w2]
pcts1=tx1[:,3] # extracted counts
pcts1=pcts1[w1:w2]
pcts1=np.ma.array(pcts1,mask=(pcts1<=20)) # Mask the zeros
pvar1=tx1[:,4] # extracted variance
pvar1=pvar1[w1:w2]
pvar1=np.ma.array(pvar1,mask=pcts1.mask)

modtable1=saltfiles[0].split('extrcol1.txt')[0]+'model.hdf5' # The hdf5 table containing normalization parameters
t1=Table.read(modtable1,path='Written082614') # Open said table
# Pull out the normalization parameters (polynomial coefficients)
vals1=t1['Prim Poly Coeff'][-np.isnan(t1['Prim Poly Coeff'])]
# Pull the radial velocity of the 5120029 observation; for velocity correction of other obs:
rv1=saltdat['Prim RV'][0]

fig6=pylab.figure(num=6) # Figure 6 for comparing individual spectra to the "best" spectrum
fig7=pylab.figure(num=7) # Figure 7 for the continually-building composite spectrum

# Initialize list for incrementing composite spectrum
complist=[None]*len(saltfiles)
complist[0]=pcts1
compvarlist=[None]*len(saltfiles)
compvarlist[0]=pvar1
nspec=[None]*len(saltfiles)
nspec[0]=np.ones_like(waves1)

for x in range(len(saltfiles)):
    if not saltfiles[x]=='./salt_reduced/512a/product/P201305120029_extrcol1.txt':
        # Pull similar data from the other SALT observations
        obsno=saltdat['Observation'][x].split('P20130')[1]
        print 'Current observation: '+obsno
        tx2=genfromtxt(saltfiles[x])
        waves2=tx2[:,1]
        pcts2=tx2[:,3]
        pcts2=np.array(pcts2)
        pvar2=tx2[:,4]

        modtable2=saltfiles[x].split('extrcol1.txt')[0]+'model.hdf5'
        t2=Table.read(modtable2,path='Written082614')
        vals2=t2['Prim Poly Coeff'][-np.isnan(t2['Prim Poly Coeff'])]
        rv2=saltdat['Prim RV'][x]

        # Correction factors for normalization and velocity shift
        rvcor=rv2-rv1 # Velocity of the pulsar at this observation relative to the 5120029 observation
        shiftwaves2=waves2-(waves2*(rvcor/(2.9979E5))) # Shift the wavelengths using that relative velocity
        polycor=np.polyval(vals1,waves2)/np.polyval(vals2,waves2) # Normalization is quotient of two polynomials
        pvar2prop=(polycor**2)*pvar2 # Error propagation; Var(aX)=(a**2)*Var(X)
        pcts2norm=pcts2*polycor # Incoming flux is normalized to the 5120029 observation
        pcts2norm=np.ma.array(pcts2norm,mask=(pcts2norm<=0)) # Mask the zeros (or less)
        pvar2prop=np.ma.array(pvar2prop,mask=(pcts2norm.mask)) # Create a matching mask
        
        # Interpolate the normalized and shifted flux and variance
        method='linear' # The interpolation method to be used (see interp1d documentation)
        fcts=interp1d(shiftwaves2[-pcts2norm.mask],
                      pcts2norm[-pcts2norm.mask],kind=method) # Interpolation function for the flux
        fvar=interp1d(shiftwaves2[-pvar2prop.mask],
                      pvar2prop[-pvar2prop.mask],kind=method) # Interpolation function for the variance
        pcts2interp=fcts(waves1) # Interpolate onto same wavelength grid as 5120029
        mapcts2interp=np.ma.array(pcts2interp,mask=(pcts1.mask))
        pvar2interp=fvar(waves1)
        mapvar2interp=np.ma.array(pvar2interp,mask=(pcts1.mask))

        # Add these data to the list
        complist[x]=mapcts2interp
        compvarlist[x]=mapvar2interp
        nspec[x]=(x+1)*np.ones_like(waves1)
        
        pylab.figure(num=6); pylab.clf()
        pylab.plot(waves1,pcts1,'k',label='5120029 Spectrum, Unmodified')
        pylab.plot(waves1,pcts2interp,'r',label=obsno+' Spectrum, Normalized, Shifted, & Interpolated')
        pylab.suptitle('Checking Spectrum Shift')
        pylab.xlabel('Wavelength (A)'); pylab.ylabel('Flux (Arbitrary)')
        pylab.legend(loc='best')
        
        acon=0
        while acon!=1:
            acc=raw_input('Do the spectra sufficiently match (check Figure 6, will be buried) ? (y/n, case sensitive) ')
            if acc=='y':
                print 'Great! Spectrum will be added to the composite.'
                acon=1
            elif acc=='n':
                acon=1
                print 'That is unfortunate.  Quitting script.'
                sys.exit('Unacceptable spectrum normalization and/or fitting.')
            else:
                print 'Invalid input (line 426).  Try again.'

        # Save the figure (I will avoid prompting here for the sake of streamlining)
        comparefig='./salt_reduced/checkspec'+obsno+'.png'
        if os.path.exists(comparefig):
            os.remove(comparefig)
        wait=raw_input('Press <return> when you are ready to save the comparison figure.')
        pylab.savefig(comparefig)
        
        print '-------------------------------------------------------------'

print 'Moving on to creating the composite spectrum.'
# The composite spectrum is a weighted average of all the individual spectra.
# The spectra are weighted by their variance at every wavelength.
# As follows:
# NewFlux = ((F_o / var_o) + sum[(F_i / var_i)]) / ((1 / var_o) + sum[1 / var_i])

complist=np.ma.array(complist)
compvarlist=np.ma.array(compvarlist)
nspec=np.array(nspec)

compflux=np.ma.cumsum((complist/compvarlist),axis=0)/np.ma.cumsum((1/compvarlist),axis=0)
compvar=np.ma.cumsum((compvarlist/nspec),axis=0)

for n in range(len(compflux)):
    pylab.figure(num=7); pylab.clf()
    pylab.plot(waves1,compflux[n])
    pylab.suptitle('Composite Spectrum with '+str(n)+' Spectra Composited')
    pylab.xlabel('Wavelength (A)')
    pylab.ylabel('Flux (Arbitrary)')

    # Save the figure
    compositefig='./salt_reduced/composite'+str(n+1)+'spectraSALT.png'
    if os.path.exists(compositefig):
        os.remove(compositefig)
    wait2=raw_input('Press <return> when ready to save the composite figure.')
    pylab.savefig(compositefig)
