from paraview.simple import *
import numpy as np
import vtk
import time
from mpi4py import MPI

ELE = 1

def FindNeighborPoints(id, uGrid):
    cellList = vtk.vtkIdList()
    neighborPoints = set()
    uGrid.GetPointCells(id, cellList)
    for i in range(cellList.GetNumberOfIds()):
        cellId = cellList.GetId(i)
        cell = uGrid.GetCell(cellId)
        for edgeId in range(cell.GetNumberOfEdges()):
            edge = cell.GetEdge(edgeId)
            if edge.GetPointId(0) == id:
                neighborPoints.add(edge.GetPointId(1))
                continue
            elif edge.GetPointId(1) == id:
                neighborPoints.add(edge.GetPointId(0))
                continue

    return neighborPoints


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
# input_14_2vtu = XMLUnstructuredGridReader(FileName=['/home/jiahui/InSituPackage/Test/SmallWing/data/input_14/input_14_2.vtu'])
# input_14_2vtu.PointArrayStatus = ['pressure', 'velocity', 'temperature']
# RenameSource('input', input_14_2vtu)

""" ------------------------- Real code starts here ------------------------- """
print("------------------")
tic = time.perf_counter()

# Get input dataset
input = FindSource('input')

# remove duplicate points
cleanGrid = CleantoGrid(Input=input)
cleanGrid.UpdatePipeline()

# Get the actual unstructured grid
uGrid = cleanGrid.GetClientSideObject().GetOutputDataObject(0)

numPoints = uGrid.GetNumberOfPoints()
print(f"Number of Points: {numPoints}")
# Find neighbours
neighborsBook = []

for pointId in range(numPoints):
    neighbors = FindNeighborPoints(pointId, uGrid)
    neighborsBook.append(neighbors)

toc = time.perf_counter()
print(f"Execution Time1 is {toc - tic} seconds")

# Test readings
# print(uGrid)
# print(sGrid)
# print(uGrid.GetPointData().GetArray(1).GetTuple(0))
# print(uGrid.GetPointData().GetArray(1).GetTuple(0)[1])

tic = time.perf_counter()

# Get input dataset
input = FindSource('input')

# # remove duplicate points
# cleanGrid = CleantoGrid(Input=input)
# cleanGrid.UpdatePipeline()

# Get the actual unstructured grid
uGrid = cleanGrid.GetClientSideObject().GetOutputDataObject(0)

numPoints = uGrid.GetNumberOfPoints()

maxima = []
minima = []

for pointId in range(numPoints):
    isMaxima = True
    isMinima = True

    currentValue = uGrid.GetPointData().GetArray(1).GetTuple(pointId)[ELE]
    neighbors = neighborsBook[pointId]

    for ids in neighbors:
        tempVal = uGrid.GetPointData().GetArray(1).GetTuple(ids)[ELE]
        if tempVal >= currentValue:
            isMaxima = False
            break

    for ids in neighbors:
        tempVal = uGrid.GetPointData().GetArray(1).GetTuple(ids)[ELE]
        if tempVal <= currentValue:
            isMinima = False
            break

    if isMaxima:
        maxima.append(pointId)

    if isMinima:
        minima.append(pointId)


# remove duplicate points

# cleanFilter = vtk.vtkCleanUnstructuredGrid()
# type(cleanFilter)
# cleanFilter.SetInputData(inputGrid)
# cleanFilter.Update()
# uGrid = cleanFilter.GetOutput()
# print(uGrid)

# Find the surface points
# surface = ExtractSurface(Input=cleanGrid)
# surface.UpdatePipeline()
# sGrid = surface.GetClientSideObject().GetOutputDataObject(0)

filter = vtk.vtkDataSetSurfaceFilter()
filter.SetInputData(uGrid)
filter.PassThroughPointIdsOn()
filter.Update()
sGrid = filter.GetOutput()

toc = time.perf_counter()

print("number of minima:")
print(len(minima))
print("number of maxima:")
print(len(maxima))


print(f"Execution Time is {toc - tic} seconds")
print("------------------")

polyData = vtk.vtkPolyData()
points = vtk.vtkPoints()
for pid in minima:
    points.InsertNextPoint(uGrid.GetPoint(pid))

polyData.SetPoints(points)
minimaProduce = TrivialProducer()
minimaProduce.GetClientSideObject().SetOutput(polyData)
Show(minimaProduce)

polyData = vtk.vtkPolyData()
points = vtk.vtkPoints()
for pid in maxima:
    points.InsertNextPoint(uGrid.GetPoint(pid))

polyData.SetPoints(points)
maximaProduce = TrivialProducer()
maximaProduce.GetClientSideObject().SetOutput(polyData)
Show(maximaProduce)