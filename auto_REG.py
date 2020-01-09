# reset all variables 
# too clean evern slicer.** can't work
import sys

sys.modules[__name__].__dict__.clear()

this = sys.modules[__name__]
for n in dir():
   if n[0]!='_': delattr(this, n) 
 
# delete individual names
del x
# remove them from the globals() object:
for name in dir():
    if not name.startswith('_'):
        del globals()[name]
        
for name in dir():
    if not name.startswith('_'):
        del locals()[name]
        
# list all Data item ----
import numpy
import pandas as pd
list1 = numpy.array(slicer.util.getNodes())

numpy.info(list1)
list1.shape


list2 = slicer.util.getNodes()

for key, value in list2.items():
    print(key, value)

key = 'Segmentation'
if key in list2:
    print(key, value)

print(list(list2.keys())[list(list2.values()).index(LinearTransform_3)])

keylist2 = list(list2.keys())

print(keylist2[2])


list(list2.items()).index(??)

list(list2.items())[2]
list(list2.items())[2][0]
keyword = 'Interaction'
for keyword in list(list2.items()):
    list(list2.items()).index(output)

output = [item for item in list(list2.items())
 if item[0] == 'Interaction']

[i for i, v in enumerate(list(list2.items())) if v[0] == 'ModelStorage_\d']

list(list2.keys())[2]
list(list2.values())[2][0]

# https://discourse.slicer.org/t/export-segmentation-to-labelmap-using-python/6801
# https://www.slicer.org/wiki/Documentation/Nightly/ScriptRepository#Export_labelmap_node_from_segmentation_node
#  ('Segmentation', (MRMLCorePython.vtkMRMLSegmentationNode)0000022A48FA3408)
seg = getNode('Segmentation')
exportedModelsNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelHierarchyNode')
slicer.modules.segmentations.logic().ExportAllSegmentsToModelHierarchy(seg, exportedModelsNode)

# ('ModelHierarchy_2', (MRMLCorePython.vtkMRMLModelHierarchyNode)0000022A48FA36A8)


# Femur segmentation using masked region growing in 3D Slicer
# https://www.youtube.com/watch?v=8Nbi1Co2rhY&t=31s

# find the left-t Right-t ----


# get the Model if not already ----


# mirror the Right-t -> Right-t-Mx ----


# Model Registration (Linear transformation_3) ----


# get Linear transformation_4 ----


# clone right volume ----


# apply Linear transformation_4 under Linear transformation_3 ----


# apply right volume copy under transformation_4 ----


# harden the transformation ----


# General Registration (Slicer Linear Transformation) ----


# apply Linear transformation_3 under Slicer Linear Transformation ----


# clone Right-t-Mx (Right-t-Mx Copy) ----


# apply Right-t-Mx Copy under transformation_3 ----


# harden the transformation ----


# get segmentation of Left-t and Right-t-Mx Copy (left volume)----


# manual cut ----


# Model to Model Distance (singed distance)----
