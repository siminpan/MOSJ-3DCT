# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import skimage.measure
import pydicom
import vtk
from vtk.util import numpy_support
import os
import matplotlib.pyplot as plt
import cv2

os.chdir('C:/Users/span/Documents/3DSlicerTutorial/CNN.test/')

# load data
# https://pyscience.wordpress.com/2014/09/08/dicom-in-python-importing-medical-image-data-into-numpy-with-pydicom-and-vtk/

PathDicom = "./23_8/"
reader = vtk.vtkDICOMImageReader()
reader.SetDirectoryName(PathDicom)
reader.Update()

# Load dimensions using `GetDataExtent`
_extent = reader.GetDataExtent()
ConstPixelDims = [_extent[1]-_extent[0]+1, _extent[3]-_extent[2]+1, _extent[5]-_extent[4]+1]

# Load spacing values
ConstPixelSpacing = reader.GetPixelSpacing()

# Get the 'vtkImageData' object from the reader
imageData = reader.GetOutput()
# Get the 'vtkPointData' object from the 'vtkImageData' object
pointData = imageData.GetPointData()
# Ensure that only one array exists within the 'vtkPointData' object
assert (pointData.GetNumberOfArrays()==1)
# Get the `vtkArray` (or whatever derived type) which is needed for the `numpy_support.vtk_to_numpy` function
arrayData = pointData.GetArray(0)

# Convert the `vtkArray` to a NumPy array
ArrayDicom = numpy_support.vtk_to_numpy(arrayData)
# Reshape the NumPy array to 3D using 'ConstPixelDims' as a 'shape'
ArrayDicom = ArrayDicom.reshape(ConstPixelDims, order='F')

# contour lines
a = ArrayDicom.copy()
plt.gray()

plt.imshow(a[:, :, 233])

d = a[:, :, 233].copy()
contours = skimage.measure.find_contours(d,200)
fig, ax = plt.subplots()
ax.imshow(d, cmap=plt.cm.gray)

for n, contour in enumerate(contours):
    ax.plot(contour[:, 1], contour[:, 0], linewidth=2)

ax.axis('image')
ax.set_xticks([])
ax.set_yticks([])
plt.show()


plt.figure(dpi=300)
plt.axes().set_aspect('equal', 'datalim')
plt.set_cmap(plt.gray())
plt.pcolormesh(np.flipud(a[:, :, 80]))

# Adaptive Thresholding
from skimage.filters import threshold_local

a = ArrayDicom.copy()
image = a[:, :, 233].copy()

block_size = 25
func = lambda arr: arr.mean()
binary_adaptive = threshold_local(image, block_size, 'generic',
                                        param=func, offset=0)
plt.imshow(binary_adaptive)


block_size = 23
binary_image1 = image > (threshold_local(image, block_size,offset=0, method='mean')+20)
plt.imshow(binary_image1)

block_size = 25
func = lambda arr: arr.mean()
binary_image2 = image > threshold_local(image, block_size, 'generic',
                                        param=func, offset=0)


plt.imshow(binary_image2)


from skimage.filters import sobel
elevation_map = sobel(binary_image1)
plt.imshow(elevation_map)


# local threshold 23_5 or 23_8
block_size = 23
binary_image1 = image > (threshold_local(image, block_size,offset=0, method='mean')+20)
plt.imshow(binary_image1)
# fill holes
from scipy import ndimage as ndi
fill_coins = ndi.binary_fill_holes(binary_image1)
plt.imshow(fill_coins)
# edge
from skimage.feature import canny
binary_image3 = canny(fill_coins, sigma=1)
plt.imshow(binary_image3)
# thicker the edge
from skimage import morphology
binary_image4 = morphology.dilation(binary_image3, morphology.disk(radius=15))
plt.imshow(binary_image4)
# remove the edge
binary_image5 = binary_image1 > binary_image4
plt.imshow(binary_image5)

# Append to 3d aray
# binary_image6 =  np.dstack((binary_image5, binary_image1))
# binary_image6.shape

result = np.empty([596, 596])
for i in range(a.shape[2]):
    image = a[:, :, i].copy()
    block_size = 23
    binary_image1 = image > (threshold_local(image, block_size,offset=0, method='mean')+20)
    fill_coins = ndi.binary_fill_holes(binary_image1)
    binary_image3 = canny(fill_coins, sigma=1)
    binary_image4 = morphology.dilation(binary_image3, morphology.disk(radius=15))
    binary_image5 = binary_image1 > binary_image4
    result =  np.dstack((result, binary_image5))

result.shape

# NumPy to VTK
# http://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/

# To remove small objects due to the segmented foreground noise, you may also consider trying skimage.morphology.remove_objects().
# Comparing edge-based and region-based segmentation
# https://scikit-image.org/docs/0.13.x/auto_examples/xx_applications/plot_coins_segmentation.html