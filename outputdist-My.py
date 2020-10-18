############## Working directory ###################
# https://stackoverflow.com/questions/5137497/find-current-directory-and-files-directory
import os
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
print(cwd)
# change the current working directory
os.chdir(path)


############ 3D Slicer directory ###########################
# https://discourse.slicer.org/t/working-directory-of-slicer/6248
# all nodes are saved relative to this path
slicer.mrmlScene.GetRootDirectory()
# write-able folder, you can use this to store any temporary data
slicer.app.temporaryPath
# Slicer core binary folder
slicer.app.slicerHome
# Slicer extensions folder
slicer.app.extensionsInstallPath
# path of a scripted module (in this example: Sample Data module)
slicer.modules.sampledata.path
############################################################

################# working process ############################################
n = slicer.util.getNode('VTK Output File_3')
a = slicer.util.array('VTK Output File_2')

dir(n)
dir(a)

name = 'VTK Output File_2'

modelNode = getNode(name)

############### Two model linked polydata ##########
newModel1 = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', name + "_invert")
newModel1.SetAndObservePolyData(modelNode.GetPolyData())
newModel1.CreateDefaultDisplayNodes()
newModel1.GetDisplayNode().Copy(modelNode.GetDisplayNode())
############### Linked model #################
newModel2 = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', name + "_linkedmodel")
newModel2.Copy(modelNode)
############### Two model sprate polydata ##########
# working progess ----
import copy
import numpy
newModel2 =  slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', name + "_Sinvert")
newModel2.SetAndObservePolyData(copy.deepcopy(modelNode.GetPolyData()))
newModel2.CreateDefaultDisplayNodes()
newModel2.GetDisplayNode().Copy(copy.deepcopy(modelNode.GetDisplayNode()))


newModel2 = copy.deepcopy(modelNode)
newModel2 = numpy.copy(modelNode)

newModel2.Modified()
#######################################################

############### 3D points ##########################
# slicer IO classes
name = 'VTK Output File_2'
modelNode = getNode(name)

outpath = 'C:/Users/span/Documents/MOSJ-TumorInjection/05.5_distances/'
# outpath = 'C:/Users/span/Documents/MOSJ-TumorInjection/05.6_3Dpoints/'
outfolder = 'NT-t/'
ainmalid = '08_slicer_IO'

# fileIn  = path + name + '.vtk'
fileOut  = outpath + outfolder + ainmalid + '.csv'

space3d = slicer.util.arrayFromModelPoints(modelNode)
np.savetxt(fileOut, space3d, delimiter=",")

################ 3D points and Signed all included ############################
# VTK IO classes
# https://discourse.vtk.org/t/save-polydata-to-csv-using-vtkdelimitedtextwriter/986/4
# working progess ----
import vtk
import csv


inpath = 'C:/Users/span/Documents/MOSJ-TumorInjection/05_3DSlicer/'
outpath = 'C:/Users/span/Documents/MOSJ-TumorInjection/05.6_3Dpoints/'
infolder = '39_20190401_VDRB2_RightHole_Left_DV/'
name = 'VTK Output File_2'
outfolder = 'DkkMoDRB-t/'
outcsv = '39_VTK_IO'

fileIn  = inpath + infolder + name + '.vtk'
fileOut  = outpath + outfolder + outcsv + '.csv'

reader = vtk.vtkGenericDataObjectReader()
reader.SetFileName(fileIn)
reader.Update()

point_obj = reader.GetOutput()
points = point_obj.GetPoints()
# vtkPoints
# Detailed Description
# represent and manipulate 3D points
# vtkPoints represents 3D points. The data model for vtkPoints is an array of vx-vy-vz triplets accessible by (point or cell) id.

table = vtk.vtkDataObjectToTable()
table.SetInputData(point_obj)
table.Update()
table.GetOutput().AddColumn(points.GetData())
table.Update()

writer = vtk.vtkDelimitedTextWriter()
writer.SetInputConnection(table.GetOutputPort())
writer.SetFileName(fileOut)
writer.Update()
writer.Write()

##################invert Color#####################################


newModel1 = getNode(name + "_Sinvert")
distanceArrayName = "Signed"
b = newModel1.GetPolyData().GetPointData().GetArray(distanceArrayName)
# https://stackoverflow.com/questions/7666981/how-to-set-data-values-on-a-vtkstructuredgrid
for i in range(b.GetNumberOfValues()):
    b.SetValue(i, -b.GetValue(i))

newModel1.Modified()

distances2 = vtk.util.numpy_support.vtk_to_numpy(newModel1.GetPolyData().GetPointData().GetArray(distanceArrayName))
distances2 = -distances2
c = b.GetData
dir(b)


b.GetValue(1)

##############################################################################

######## capture viewers / screenshot ######################
# https://www.slicer.org/wiki/Documentation/Nightly/ScriptRepository
# Capture a single view:
import ScreenCapture
l=ScreenCapture.ScreenCaptureLogic()
l.captureImageFromView(l.viewFromNode(slicer.util.getNode('vtkMRMLSliceNodeRed')), 'c:/tmp/red.png')

# Capture all the views save it into a file:
import ScreenCapture
cap = ScreenCapture.ScreenCaptureLogic()
cap.showViewControllers(False)
cap.captureImageFromView(None, slicer.mrmlScene.GetRootDirectory() + '/VTK Output File_2_BR.png')
cap.showViewControllers(True)

############################################################

##### get model to model distance by typing 5 lines of code: ########
# https://discourse.slicer.org/t/model-to-model-distance-module-questions/3481/6
# https://vtk.org/doc/nightly/html/classvtkDistancePolyDataFilter.html
# Get two model nodes that we want to compute distances of
m1=getNode('Segment_1')
m2=getNode('Segment_2')

# Compute distance
distanceFilter = vtk.vtkDistancePolyDataFilter()
distanceFilter.SetInputData(0, m1.GetPolyData())
distanceFilter.SetInputData(1, m2.GetPolyData())
distanceFilter.Update()
m1.SetAndObservePolyData(distanceFilter.GetOutput())

# Use the computed distance to color the node
m1.GetDisplayNode().SetActiveScalarName('Distance')
m1.GetDisplayNode().SetScalarVisibility(True)

##############################################################




# Start

name = 'VTK Output File_2'
modelNode = getNode(name)
distanceArrayName = "Signed"

# Get distances from point data
import numpy as np
distances = vtk.util.numpy_support.vtk_to_numpy(modelNode.GetPolyData().GetPointData().GetArray(distanceArrayName))

outpath = 'C:/Users/span/Documents/MOSJ-TumorInjection/05.5_distances/'
outfolder = 'NT-t/'
ainmalid = '08_distances'

fileOut  = outpath + outfolder + ainmalid + '.csv'
# write csv
# import numpy as np
np.savetxt(fileOut, distances, delimiter=",")

# Print basic stats
print(distances)
print("Minimum distance: %f" % min(distances))
print("25%% Percentile distance: %f" % np.percentile(distances, 25))
print("Mean distance: %f" % np.mean(distances))
print("Median distance: %f" % np.median(distances))
print("75%% Percentile distance: %f" % np.percentile(distances, 75))
print("Maximum distance: %f" % max(distances))

# For Slicer 4.10 I would probably suggest opening a terminal window to the location of PythonSlicer and doing
######################
# PythonSlicer -m pip install scipy
######################
#The error you are getting is because the “pip” package was updated and “main” was moved out. You could try some work around as was suggested earlier in this discussion such as

# from pip._internal import main as pipmain
# pipmain(['install','pandas'])

# pipmain(['uninstall','scipy'])

# max pooling
# scikit-image

# import os
# os.system('PythonSlicer -m pip install pandas')

# import pandas as pd
# df1 = pd.DataFrame(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
#                    columns=['a', 'b', 'c'])

# data = np.array([['','Min', '25th', 'Mean', 'Median', '75th', 'Max'],
#                 ['Row1', min(distances), np.percentile(distances, 25), np.mean(distances), np.median(distances), np.percentile(distances, 75), max(distances)]
#                 ])
data = np.array([[min(distances), np.percentile(distances, 25), np.mean(distances), np.median(distances), np.percentile(distances, 75), max(distances)]])
np.savetxt("C:/Users/span/Documents/MOSJ-TumorInjection/05_3DSlicer/distances/\
VDRB1_2Left-t-Ds.csv", data, delimiter=",")

# Compute histogram
histogram = np.histogram(distances, bins=100)

# Save results to a new table node
tableNode=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", modelNode.GetName() + " histogram-t")
updateTableFromArray(tableNode, histogram)
tableNode.GetTable().GetColumn(0).SetName("Count")
tableNode.GetTable().GetColumn(1).SetName("Intensity")


# slicer.util.saveNode(tableNode(), "C:/Users/span/Documents/MOSJ-TumorInjection/05_3DSlicer/distances/\
# vdrb2_righthole-DH.csv")

###########################################################
###########################################################
 #Projection plot serie
# projectionPlotSeries = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", "PCA projection")
# projectionPlotSeries.SetAndObserveTableNodeID(projectionTableNode.GetID())
# projectionPlotSeries.SetXColumnName("Count")
# projectionPlotSeries.SetYColumnName("Intensity")
# projectionPlotSeries.SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatter)
# projectionPlotSeries.SetLineStyle(slicer.vtkMRMLPlotSeriesNode.LineStyleNone)
# projectionPlotSeries.SetUniqueColor()
###########################################################
###########################################################

# show plot https://discourse.slicer.org/t/use-full-power-of-python-in-slicer/7162/5

# Create plot
plotSeriesNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", modelNode.GetName() + ' histogram-t')
plotSeriesNode.SetAndObserveTableNodeID(tableNode.GetID())
plotSeriesNode.SetXColumnName("Intensity")
plotSeriesNode.SetYColumnName("Count")
plotSeriesNode.SetPlotType(plotSeriesNode.PlotTypeScatterBar)
plotSeriesNode.SetColor(0, 0.6, 1.0)

# Create chart and add plot
plotChartNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode")
plotChartNode.AddAndObservePlotSeriesNodeID(plotSeriesNode.GetID())
plotChartNode.YAxisRangeAutoOff()
plotChartNode.SetYAxisRange(0, 500000)

# Show plot in layout
slicer.modules.plots.logic().ShowChartInLayout(plotChartNode)

# Show plot in layout

layoutManager = slicer.app.layoutManager()
layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpPlotView)
plotWidget = layoutManager.plotWidget(0)

plotViewNode = plotWidget.mrmlPlotViewNode()
plotViewNode.SetPlotChartNodeID(plotChartNode.GetID())
