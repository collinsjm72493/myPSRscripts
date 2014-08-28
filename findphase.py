from astropy.table import Table,Column
import ephemeris
import de421
import numpy as np

# replace the .par file here with the .par file for your pulsar
orbit=ephemeris.ELL1Ephemeris('1723-2837.par')
jpleph = ephemeris.JPLEphemeris(de421)

# this is just an example
# you want some columns to use to calculate time/velocity
data=Table.read('rvtable2.hdf5',path='OVA')

# this will allow barycentering
# the important thing here is to use the right column from the data
# and get the coordinates from the .par file
# you also need to get the site right.  MK = Mauna Kea (Gemini)
# you can also use MP for Palomar
# or SAAO for SALT
obs=ephemeris.Observation(data['GMJD'],ra=orbit['RAJ'],dec=orbit['DECJ'],
                          site='SAAO')

# this actually calculates the barycentric delay and velocity
delay,rv=obs.calculate()

# as an example, this makes a new column for barycentric MJD
BMJD=data['GMJD']+delay/86400
data.add_column(Column(data=BMJD,name='BMJD',units='days',format='%.7f'),index=3)

# the rv should be applied to the radial velocities
# be careful about the sign

# now do the phase
t=orbit.mean_anomaly(BMJD)/2/np.pi
t-=np.floor(t)
phase=t
data.add_column(Column(data=phase,name='Phase',units='cycles',format='%.6f'),index=4)

data.write('phaservtable2.hdf5',path='OVA')
print 'The phase table: '
print data
