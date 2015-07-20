import sunpy
from sunpy.map import Map
import numpy as np
import matplotlib.pyplot as plt
import sys
from os.path import join, expanduser
sys.path.append(expanduser('~/CoronaTemps/'))
from temperature import TemperatureMap as TMap
from astropy import units as u

# Load a temperature map and get the EM data out of it
tmap = TMap('2011-02-15', maps_dir=expanduser('~/coronal-backprojection/'),
            n_params=3)
EMmap = Map(10**tmap.emission_measure, tmap.meta.copy())

# Define ranges of reconstruction space
xrange = EMmap.xrange.value / EMmap.rsun_obs.value
yrange = EMmap.yrange.value / EMmap.rsun_obs.value
zrange = [min([xrange[0], yrange[0]]), max([xrange[1], yrange[1]])]
print xrange, yrange, zrange

# Calculate array of distance from solar centre
x, y, z = np.mgrid[xrange[0]:xrange[1]:0.01,#(xrange[1]-xrange[0])/512,
                   yrange[0]:yrange[1]:0.01,#(yrange[1]-yrange[0])/512,
                   zrange[0]:zrange[1]:0.01]#(zrange[1]-zrange[0])/512]
projected_r = np.sqrt((x ** 2.0) + (y ** 2.0))
r = np.sqrt((x ** 2.0) + (y ** 2.0) + (z ** 2.0))

# Discount solar interior
model = r.copy()
model[r < 1] = np.nan
# Discount behind the sun
model[(projected_r < 1) * (z > 0)] = np.nan
# Discount excessive distance
model[r > 1.25] = np.nan
# Impose radial distribution based on emission measure scale height
em_sclht = (23 * u.Mm).to(u.solRad).value
model = np.exp(-(model-1)/em_sclht)
# Normalise each LOS
model /= np.nansum(model, axis=2)[:, :, None]

# Distribute EM over LOS according to model
#EMmap = EMmap.superpixel((8, 8))
#projection = model * EMmap.data
#projection = np.sqrt(projection)

print np.nanmin(model), np.nanmax(model)
for j in np.linspace(0, model.shape[1]-1, 5):
  fig = plt.figure(figsize=(32, 24))
  slice = model[:, int(j), :]
  plt.imshow(slice, cmap='hot', interpolation='none',
             #vmin=np.nanmin(model), vmax=np.nanmax(model))
             vmin=np.nanmin(slice), vmax=np.nanmax(slice))
             #vmin=0, vmax=1)
  plt.colorbar()
  plt.savefig('slicemap_{}'.format(int(j)))
  plt.close()
