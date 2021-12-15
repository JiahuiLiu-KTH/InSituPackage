from paraview.simple import *
import numpy as np
import vtk
import time
import math
from vtkmodules.vtkCommonCore import vtkMath
from mpi4py import MPI

paraview.simple._DisableFirstRenderCameraReset()

# Get input dataset
input = FindSource('input')

#Get uncleand grid
grid = input.GetClientSideObject().GetOutputDataObject(0)

# remove duplicate points
cleanGrid = CleantoGrid(Input=input)
cleanGrid.UpdatePipeline()

# Get the actual unstructured grid
uGrid = cleanGrid.GetClientSideObject().GetOutputDataObject(0)

points = grid.GetPoints()

poly = vtk.vtkPolyData()
poly.SetPoints(points)

vFilter = vtk.vtkVertexGlyphFilter()
cFilter = vtk.vtkCleanPolyData()

vFilter.SetInputData(poly)
vFilter.Update()
cFilter.SetInputData(vFilter.GetOutput())
cFilter.Update()
newPoly = cFilter.GetOutput()
print(newPoly.GetNumberOfPoints())

p0 = uGrid.GetPoint(0)
p1 = uGrid.GetPoint(300)
distSquared = vtkMath.Distance2BetweenPoints(p0, p1)
dist = math.sqrt(distSquared)

result = vtk.vtkIdList()
pl = vtk.vtkPointLocator()
pl.SetDataSet(newPoly)
pl.AutomaticOn()
pl.BuildLocator()
pl.FindPointsWithinRadius(dist, p0, result)

for i in range(result.GetNumberOfIds()):
    print(result.GetId(i))