import sunpy
from sunpy.map import Map
import numpy as np
from matplotlib import cm, _cm, rc, patches
rc('savefig', bbox='tight', pad_inches=0.5)
import matplotlib.pyplot as plt
import sys
from os.path import join, expanduser
sys.path.append(expanduser('~/CoronaTemps/'))
from temperature import TemperatureMap as TMap
from astropy import units as u
from scipy.interpolate import RegularGridInterpolator as rgi
from itertools import product

emmapcubehelix = _cm.cubehelix(s=2.8, r=-0.7, h=1.4, gamma=1.0)
cm.register_cmap(name='emhelix', data=emmapcubehelix)
emcm = cm.get_cmap('emhelix')

# Load a temperature map and get the EM data out of it
tmap = TMap('2011-02-15', maps_dir=expanduser('~/coronal-backprojection/'),
            n_params=3)
EMmap = Map(10.0**tmap.emission_measure, tmap.meta.copy())
EMlog = Map(tmap.emission_measure, tmap.meta.copy())
print np.nanmin(EMlog.data), np.nanmax(EMlog.data)
fig = plt.figure(figsize=(32, 24))
EMlog.plot(cmap=emcm)
plt.colorbar()
plt.savefig('input-em')
plt.close()

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
  fig = plt.figure(figsize=(32, 24))
  ax = fig.add_subplot(111, axisbg='gray')
  plt.imshow(np.log10(model[:, side/2, :]), cmap='afmhot', interpolation='none',
             vmin=0.0, vmax=1.0,
             extent=[zrange[0], zrange[1], xrange[1], xrange[0]])
  ax.add_artist(patches.Circle((0, 0), 1.0, color='black', fill=False))
  plt.colorbar()
  plt.savefig('rad-dist-model')
  plt.close()
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

#for j in np.linspace(0, model.shape[1]-1, 7):
fig = plt.figure(figsize=(32, 24))
ax = fig.add_subplot(111, axisbg='gray')
print side/2, x[side/2]
slice = model[:, side/2, :]
#slice = 10.0**model[:, side/2, :]
plt.imshow(slice, cmap=emcm, interpolation='none',
           #vmin=10.0**np.nanmin(model), vmax=10.0**np.nanmax(model))
           vmin=np.nanmin(slice), vmax=np.nanmax(slice),
           extent=[zrange[0], zrange[1], xrange[1], xrange[0]])
           #vmin=0, vmax=1)
ax.add_artist(patches.Circle((0, 0), 1.0, color='black', fill=False))
ax.add_artist(patches.Circle((0, 0), 1.1, color='black', fill=False, ls='dashed'))
plt.colorbar()
plt.savefig('slicemap')#_{}'.format(int(j)))
plt.close()

"""fig = plt.figure(figsize=(48, 24))
fig.add_subplot(1, 2, 1)
EMlog.plot(cmap=emcm, vmin=20, vmax=30)
plt.colorbar()

fig.add_subplot(1, 2, 2)
totmap = Map(np.log10(np.nansum((10.0**model)**2.0, axis=2) * dz), EMmap.meta.copy())
print np.nanmin(totmap.data), np.nanmax(totmap.data)
totmap.plot(cmap=emcm, vmin=20, vmax=30)
plt.colorbar()

plt.savefig('comparison')
plt.close()"""

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
modelsph = modelsph.reshape(21, 180, 360)

# Define lon and lat bottom corners for AR, QS and CH regions, in that order
regioncoords = np.array([[0, -25],
                         [15, 15],
                         [30, -45]])
crads = regioncoords.T * (np.pi / 180)
crads2 = (regioncoords+20).T * (np.pi / 180)
#print crads.shape
print crads
# Convert lon and lat to array indices
regioninds = np.array([[np.where(abs(lon - c) == abs(lon - c).min()) for c in crads[0]],
                       [np.where(abs(lat - c) == abs(lat - c).min()) for c in crads[1]]])
regioninds = regioninds.reshape((2, 3)).T
regioninds2 = np.array([[np.where(abs(lon - c) == abs(lon - c).min()) for c in crads2[0]],
                        [np.where(abs(lat - c) == abs(lat - c).min()) for c in crads2[1]]])
regioninds2 = regioninds2.reshape((2, 3)).T
#print regioninds
#print regioninds.shape
print lon[regioninds[:, 0]], lat[regioninds[:, 1]]
regions = []
labels = ['Active region', 'Quiet sun', 'Coronal hole']
for i, j in zip(regioninds, regioninds2):
    block = modelsph[1:, i[1]:j[1], i[0]:j[0]]
    print block.shape
    radavg = np.nanmean(block, axis=(1, 2))#modelsph[:, i[1]:j[1], i[0]:j[0]], axis=(1, 2))
    regions.append(np.log10(radavg))
#print regions

# Plot average density with height
avgdens = np.nanmean(modelsph, axis=(1, 2))[1:]
em_sclht = (46 * u.Mm).to(u.solRad).value
modeld = np.exp(-(rad-1)/em_sclht) * (10**(regions[0].max()))
fig = plt.figure(figsize=(16, 12))
plt.plot(rad-1, np.log10(modeld), color='black', linestyle='--', label='Model')
plt.plot(rad[1:]-1, np.log10(avgdens), color='black', label='Full')
for reg, label in zip(regions, labels):
    plt.plot(rad[1:]-1, reg, label=label)
plt.legend()
plt.ylim(7.0, 9.5)
plt.title('Average estimated density vs height')
plt.xlabel('Height above solar surface ($R_{\odot}$)')
plt.ylabel('Average density')
plt.savefig('avgd-vs-ht')
plt.close()

# Plot map of density at constant height
modelsph = np.log10(modelsph)
for i, radslice in enumerate(modelsph):
    if rad[i] != 1.1: continue
    print rad[i],
    fig = plt.figure(figsize=(32, 24))
    ax = fig.add_subplot(111, axisbg='gray')
    plt.imshow(radslice, cmap=emcm, interpolation='none',
               #vmin=np.nanmin(modelsph), vmax=np.nanmax(modelsph),
               vmin=np.nanmin(radslice), vmax=np.nanmax(radslice),
               extent=[-180, 180, -90, 90], origin='lower')
    plt.grid(linewidth=2)
    for coords in regioncoords:
        print coords
        rect = patches.Rectangle(coords, 20, 20, color='white', fill=False)
        ax.add_artist(rect)
    #plt.axhline(0, color='black')
    #plt.axvline(0, color='black')
    plt.colorbar(orientation='horizontal')
    plt.title('Estimated coronal density at h={}'.format(rad[i]-1))
    plt.xlabel('Longitude ($^{\circ}$)')
    plt.ylabel('Latitude ($^{\circ}$)')
    plt.savefig('radslice_h={}.png'.format(rad[i]-1))
    plt.close()

print 'Plots complete'
