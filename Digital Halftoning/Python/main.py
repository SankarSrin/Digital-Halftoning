#Digital Halftoning
#Mr.Sankarasrinivasan S
#Ref:  Digital Halftone Database (DHD): A Comprehensive Analysis on Halftone Types


import cv2
import numpy as np
import math
import halftoning
img = cv2.imread('peppers.tif')
img = cv2.resize(img, dsize=(512, 512), interpolation=cv2.INTER_CUBIC)
im=img[:,:,1]

## Bayer's Halftoning
halftoning.bayers(img)

##Ulichney Halftoning
halftoning.Ulichney(img)

##Error Diffusion
out=halftoning.errordiff(img)




