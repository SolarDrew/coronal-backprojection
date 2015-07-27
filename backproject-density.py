import sunpy
from sunpy.map import Map
import numpy as np
from matplotlib import cm, _cm, rc
rc('savefig', bbox='tight', pad_inches=0.5)
import matplotlib.pyplot as plt
import sys
from os.path import join, expanduser
sys.path.append(expanduser('~/CoronaTemps/'))
from temperature import TemperatureMap as TMap
from astropy import units as u
from scipy.interpolate import RegularGridInterpolator as rgi
from itertools import product

# Load a temperature map and get the EM data out of it
tmap = TMap('2011-02-15', maps_dir=expanduser('~/coronal-backprojection/'),
            n_params=3)
EMmap = Map(10.0**tmap.emission_measure, tmap.meta.copy())
EMlog = Map(tmap.emission_measure, tmap.meta.copy())
print np.nanmin(EMlog.data), np.nanmax(EMlog.data)

# Define ranges of reconstruction space
rsun = EMmap.rsun_arcseconds
xrange = (EMmap.xrange[0] / rsun, EMmap.xrange[1] / rsun)
yrange = (EMmap.yrange[0] / rsun, EMmap.yrange[1] / rsun)
zrange = [min([xrange[0], yrange[0]]), max([xrange[1], yrange[1]])]
print xrange, yrange, zrange

side = 1024

try:
  # Attempt to load distribution and coordinates files
  model = np.memmap(filename='distmodel', mode='r', dtype=np.float64,
                    shape=(side, side, side))
  x, y, z = np.memmap(filename='coords', mode='r', dtype=np.float64,
                      shape=(3, side, side, side))
  print 'Distribution model loaded'
except:
  print 'Calculating 3D distribution'
  # Creste a distribution file
  model = np.memmap(filename='distmodel', mode='w+', dtype=np.float64,
                    shape=(side, side, side))

  # Calculate array of distance from solar centre
  coords = np.memmap(filename='coords', mode='w+', dtype=np.float64,
                     shape=(3, side, side, side))
  coords[:] = np.mgrid[xrange[0]:xrange[1]:(xrange[1]-xrange[0])/side,
                       yrange[0]:yrange[1]:(yrange[1]-yrange[0])/side,
                       zrange[0]:zrange[1]:(zrange[1]-zrange[0])/side]
  coords.flush()
  print coords.shape
  x, y, z = coords
  projected_r = np.sqrt((x ** 2.0) + (y ** 2.0))
  r = np.sqrt((x ** 2.0) + (y ** 2.0) + (z ** 2.0))

  # Discount solar interior
  model[:] = r.copy()
  model[r < 1] = np.nan
  # Discount behind the sun
  model[(projected_r < 1) * (z > 0)] = np.nan
  # Discount excessive distance
  model[r > 1.25] = np.nan
  # Impose radial distribution based on emission measure scale height
  em_sclht = (23 * u.Mm).to(u.solRad).value
  model[:] = np.exp(-(model-1)/em_sclht)
  # Normalise each LOS
  model /= np.nansum(model, axis=2)[:, :, None]

  model.flush()
  print 'Distribution model calculated'

# Distribute EM over LOS according to model
EMmap = EMmap.superpixel((4096./side, 4096./side), method='average')
EMlog = EMlog.superpixel((4096./side, 4096./side), method='average')
print np.nanmin(EMlog.data), np.nanmax(EMlog.data)
dz = ((z[0, 0, 1] - z[0, 0, 0]) * u.solRad).to(u.cm).value
model = (model * EMmap.data[:, :, None]) / dz
model = np.log10(np.sqrt(model))
np.save('density', model[:].copy())
print 'Density cube calculated'# and saved'

emmapcubehelix = _cm.cubehelix(s=2.8, r=-0.7, h=1.4, gamma=1.0)
cm.register_cmap(name='emhelix', data=emmapcubehelix)
emcm = cm.get_cmap('emhelix')

#for j in np.linspace(0, model.shape[1]-1, 7):
fig = plt.figure(figsize=(32, 24))
print side/2, x[side/2]
slice = model[:, side/2, :]
#slice = 10.0**model[:, side/2, :]
plt.imshow(slice, cmap=emcm, interpolation='none',
           #vmin=10.0**np.nanmin(model), vmax=10.0**np.nanmax(model))
           vmin=np.nanmin(slice), vmax=np.nanmax(slice))
           #vmin=0, vmax=1)
plt.colorbar()
plt.savefig('slicemap')#_{}'.format(int(j)))
plt.close()

fig = plt.figure(figsize=(48, 24))
fig.add_subplot(1, 2, 1)
EMlog.plot(cmap=emcm, vmin=20, vmax=30)
plt.colorbar()

fig.add_subplot(1, 2, 2)
totmap = Map(np.log10(np.nansum((10.0**model)**2.0, axis=2) * dz), EMmap.meta.copy())
print np.nanmin(totmap.data), np.nanmax(totmap.data)
totmap.plot(cmap=emcm, vmin=20, vmax=30)
plt.colorbar()

plt.savefig('comparison')
plt.close()

# Convert cartesian density data into spherical so it can be plotted
x1 = np.linspace(xrange[0], xrange[1], side)
y1 = np.linspace(yrange[0], yrange[1], side)
z1 = np.linspace(zrange[0], zrange[1], side)
intf = rgi((x1, y1, z1), 10.0**model[:].copy())
rad = np.linspace(1.0, 1.25, 21)
lat = np.linspace(-np.pi/2, np.pi/2, 180)
lon = np.linspace(-np.pi, np.pi, 360)
print lon.shape, lat.shape, rad.shape
radi, lati, loni = np.array([i for i in product(rad, lat, lon)]).T
print loni.shape, lati.shape, radi.shape
yi = radi * np.cos(lati) * np.sin(loni)
xi = np.sin(lati)
zi = -radi * np.cos(lati) * np.cos(loni)
print xi.shape, yi.shape, zi.shape
print np.array([xi, yi, zi]).shape
modelsph = intf(np.array([xi, yi, zi]).T)
print modelsph.shape, np.nanmin(modelsph), np.nanmax(modelsph)
modelsph = modelsph.reshape(21, 180, 360)#[:, ::-1, :]

# Plot average density with height
avgdens = np.log10(np.nanmean(modelsph, axis=(1, 2)))[1:]
#em_sclht = (23 * u.Mm).to(u.solRad).value
#modeld = np.exp(-(rad-1)/em_sclht)
fig = plt.figure(figsize=(16, 12))
#plt.plot(rad-1, modeld, color='black', linestyle='--')
plt.plot(rad[1:]-1, avgdens, color='black')
plt.title('Average estimated density vs height')
plt.xlabel('Height above solar surface ($R_{\odot}$)')
plt.ylabel('Average density')
plt.savefig('avgd-vs-ht')
plt.close()

# Plot maps of density at constant height
modelsph = np.log10(modelsph)
for i, radslice in enumerate(modelsph):
    #print rad[i]
    fig = plt.figure(figsize=(32, 24))
    plt.imshow(radslice, cmap=emcm, interpolation='none',
               #vmin=np.nanmin(modelsph), vmax=np.nanmax(modelsph),
               vmin=np.nanmin(radslice), vmax=np.nanmax(radslice),
               extent=[-180, 180, -90, 90], origin='lower')
    plt.grid(linewidth=2)
    #plt.axhline(0, color='black')
    #plt.axvline(0, color='black')
    plt.colorbar(orientation='horizontal')
    plt.title('Estimated coronal density at h={}'.format(rad[i]-1))
    plt.xlabel('Longitude ($^{\circ}$)')
    plt.ylabel('Latitude ($^{\circ}$)')
    plt.savefig('radslice_h={}.png'.format(rad[i]-1))
    plt.close()

print 'Plots complete'
