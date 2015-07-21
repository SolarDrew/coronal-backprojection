import sunpy
from sunpy.map import Map
import numpy as np
from matplotlib import cm, _cm
import matplotlib.pyplot as plt
import sys
from os.path import join, expanduser
sys.path.append(expanduser('~/CoronaTemps/'))
from temperature import TemperatureMap as TMap
from astropy import units as u

# Load a temperature map and get the EM data out of it
tmap = TMap('2011-02-15', maps_dir=expanduser('~/coronal-backprojection/'),
            n_params=3)
EMmap = Map(10.0**tmap.emission_measure, tmap.meta.copy())

try:
  # Attempt to load distribution and coordinates files
  model = np.memmap(filename='distmodel', mode='r', dtype=np.float64,
                    shape=(1024, 1024, 1024))
  x, y, z = np.memmap(filename='coords', mode='r', dtype=np.float64,
                      shape=(3, 1024, 1024, 1024))
except:
  print 'Calculating 3D distribution'
  # Creste a distribution file
  model = np.memmap(filename='distmodel', mode='w+', dtype=np.float64,
                    shape=(1024, 1024, 1024))

  # Define ranges of reconstruction space
  rsun = EMmap.rsun_arcseconds
  xrange = (EMmap.xrange[0] / rsun, EMmap.xrange[1] / rsun)
  yrange = (EMmap.yrange[0] / rsun, EMmap.yrange[1] / rsun)
  zrange = [min([xrange[0], yrange[0]]), max([xrange[1], yrange[1]])]
  print xrange, yrange, zrange

  # Calculate array of distance from solar centre
  coords = np.memmap(filename='coords', mode='w+', dtype=np.float64,
                     shape=(3, 1024, 1024, 1024))
  coords[:] = np.mgrid[xrange[0]:xrange[1]:(xrange[1]-xrange[0])/1024,
                       yrange[0]:yrange[1]:(yrange[1]-yrange[0])/1024,
                       zrange[0]:zrange[1]:(zrange[1]-zrange[0])/1024]
  coords.flush()
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

# Distribute EM over LOS according to model
EMmap = EMmap.superpixel((4, 4), method='average')
dz = ((z[0, 0, 1] - z[0, 0, 0]) * u.solRad).to(u.cm).value
model = (model * EMmap.data[:, :, None]) / dz
model = np.log10(np.sqrt(model))

emmapcubehelix = _cm.cubehelix(s=2.8, r=-0.7, h=1.4, gamma=1.0)
cm.register_cmap(name='emhelix', data=emmapcubehelix)
emcm = cm.get_cmap('emhelix')

print 'Calculations successful'
#for j in np.linspace(0, model.shape[1]-1, 7):
fig = plt.figure(figsize=(32, 24))
#slice = model[:, int(j), :]
slice = model[:, int(model.shape[2]/2), :]
plt.imshow(slice, cmap=emcm, interpolation='none',
           vmin=np.nanmin(model), vmax=np.nanmax(model))
           #vmin=np.nanmin(slice), vmax=np.nanmax(slice))
           #vmin=0, vmax=1)
plt.colorbar()
plt.savefig('slicemap')#_{}'.format(int(j)))
plt.close()

print 'Plots complete'
