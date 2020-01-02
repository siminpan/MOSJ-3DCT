# list all Data item ----
import numpy
import pandas as pd
list1 = numpy.array(slicer.util.getNodes())

numpy.info(list1)
list1.shape


list2 = slicer.util.getNodes()

for key, value in list2.items():
    print(key, value)

key = "LinearTransform_3"
if key in list2:
    print(key, value)

keylist2 = list(list2.keys())

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
