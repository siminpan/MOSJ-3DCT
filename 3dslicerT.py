a = slicer.util.array('spgr')

print(a)

print(a.min(),a.max())

# Defines a python function ----
# Telling Slicer that the image data for node 'n' has been
# modified causes the slice view windows to refresh
def toggle():
 n = slicer.util.getNode('spgr')
 a = slicer.util.array('spgr')
 a[:] = a.max() - a
 n.Modified()
 print('Toggled')
toggle()

# show a Toggle button
c = qt.QPushButton('Toggle')
c.connect('clicked()',toggle)
b = qt.QPushButton('Toggle2')
b.connect('clicked()',toggle)
b.show()
c.show()


# max pooling ----
import numpy as np
import skimage.measure

a = slicer.util.array('Volume_1 Copy')
b = skimage.measure.block_reduce(a, (1,2,2), np.max)

volumeNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode')
volumeNode.CreateDefaultDisplayNodes()
updateVolumeFromArray(volumeNode, b)
setSliceViewerLayers(background=volumeNode)

# kernel for edge detection ----
import numpy as np
from scipy import ndimage

a = slicer.util.array('Volume_1 Copy')

a = slicer.util.array('Output Volume 0.05')

# trim based on kernel size
d = np.delete(a, [0,a.shape[2]-1], 2)
d = np.delete(d, [0,d.shape[1]-1], 1)

# kernel for Edge detection
kernel_laplace = np.array([np.array([1, 0, -1]), np.array([0, 0, 0]), np.array([-1, 0, 1])])

kernel_laplace = np.array([np.array([0, -1, 0]), np.array([-1, 4, -1]), np.array([0, -1, 0])])

kernel_laplace = np.array([np.array([1, 1, 1]), np.array([1, -8, 1]), np.array([1, 1, 1])])

# fit kernel for volme dim
kernel_laplace = np.expand_dims(kernel_laplace, axis=0)

b = ndimage.convolve(d, kernel_laplace, mode='reflect')

volumeNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode')
volumeNode.CreateDefaultDisplayNodes()
updateVolumeFromArray(volumeNode, b)
setSliceViewerLayers(background=volumeNode)

n = slicer.util.getNode('Volume_2')
n.Modified()


# try rewrite a existing volume ----
n = slicer.util.getNode('Volume_1 Copy')
a = slicer.util.array('Volume_1 Copy')
a = skimage.measure.block_reduce(a, (1,2,2), np.max)
n.Modified()

# load pythong output ----
b = np.loadtxt('C:/Users/span/Documents/3DSlicerTutorial/CNN.test/23_o_res_test.txt')
b = b.reshape((596, 596, 563))
b = np.moveaxis(b, -1, 0)
b = np.delete(b, 0, axis=0)
b = np.moveaxis(b, 1, 2)
b.shape

volumeNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode')
volumeNode.CreateDefaultDisplayNodes()
updateVolumeFromArray(volumeNode, b)
setSliceViewerLayers(background=volumeNode)

# IJK to RAS
n = slicer.util.getNode('Output Volume 0.05')
mat = vtk.vtkMatrix3x3()
n.GetIJKToRASDirections(mat)

b = slicer.util.getNode('Volume_1')
b.SetIJKToRASDirections(mat)
b.Modified()
