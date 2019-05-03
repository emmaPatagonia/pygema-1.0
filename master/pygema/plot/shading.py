# fuente: http://rnovitsky.blogspot.com/2010/04/using-hillshade-image-as-intensity.html
#The difference in the shading colors derived from the method used to produce it. While the matplotlib method uses "hard light" method I use a "soft light" method. the matplotlib is converting the RGB colors to HSV and then calculate the new saturation and value according to the intensity. I use a formula based on the description of ImageMagick's pegtop_light.which is much faster as it is a single formula. Another advantage is the option to use a separate layer as the intensity and another as the data used for colors.

#!/bin/env python
from pylab import *
def set_shade(a,intensity=None,cmap=cm.jet,scale=10.0,azdeg=290.0,altdeg=45.0):

  if intensity is None:
# hilshading the data
    intensity = hillshade(a,scale=scale,azdeg=azdeg,altdeg=altdeg)
  else:
# or normalize the intensity
    intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())
# get rgb of normalized data based on cmap
  rgb = cmap((a-a.min())/float(a.max()-a.min()))[:,:,:3]
# form an rgb eqvivalent of intensity
  d = intensity.repeat(3).reshape(rgb.shape)
# simulate illumination based on pegtop algorithm.
  rgb = 2*d*rgb+(rgb**2)*(1-2*d)
  return rgb

def hillshade(data,scale=10.0,azdeg=165.0,altdeg=45.0):

  # convert alt, az to radians
  az = azdeg*pi/180.0
  alt = altdeg*pi/180.0
  # gradient in x and y directions
  dx, dy = gradient(data/float(scale))
  slope = 0.5*pi - arctan(hypot(dx, dy))
  aspect = arctan2(dx, dy)
  intensity = sin(alt)*sin(slope) + cos(alt)*cos(slope)*cos(-az - aspect - 0.5*pi)
  intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())
  return intensity

