modelNode = getNode('VTK Output File')
distanceArrayName = "Signed"

# Get distances from point data
import numpy as np
distances = vtk.util.numpy_support.vtk_to_numpy(modelNode.GetPolyData().GetPointData().GetArray(distanceArrayName))

# Print basic stats
print("Minimum distance: %f" % min(distances))
print("Maximum distance: %f" % max(distances))
print("Mean distance: %f" % np.mean(distances))

# Compute histogram
histogram = np.histogram(distances, bins=100)

# Save results to a new table node
tableNode=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", modelNode.GetName() + " histogram")
updateTableFromArray(tableNode, histogram)
tableNode.GetTable().GetColumn(0).SetName("Count")
tableNode.GetTable().GetColumn(1).SetName("Intensity")

# Create plot

plotDataNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotDataNode")
plotDataNode.SetAndObserveTableNodeID(tableNode.GetID())
plotDataNode.SetXColumnName("Intensity")
plotDataNode.SetYColumnName("Count")

plotChartNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode")
plotChartNode.AddAndObservePlotDataNodeID(plotDataNode.GetID())
plotChartNode.SetAttribute("Type", "Bar") # delete this line for line plot

# Show plot in layout

layoutManager = slicer.app.layoutManager()
layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpPlotView)
plotWidget = layoutManager.plotWidget(0)

plotViewNode = plotWidget.mrmlPlotViewNode()
plotViewNode.SetPlotChartNodeID(plotChartNode.GetID())
