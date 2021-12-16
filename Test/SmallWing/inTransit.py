import time
import math
import vtk
import numpy as np
from paraview.simple import *
from paraview import coprocessing
from scipy import stats
from mpi4py import MPI
from vtk.util import numpy_support
from vtkmodules.vtkCommonCore import vtkMath
from vtkmodules.vtkFiltersGeneral import vtkVertexGlyphFilter

NSTEP = 150
BBOX = (-2, 3, -2, 2, 0, 0.05)

# Initialize parameters for KDE
NumScale = 60
xMin, xMax, yMin, yMax, zMin, zMax = BBOX

xMin = -1
yMax = 1
yMin = -1

xNum = (xMax - xMin) * NumScale + 1
yNum = (yMax - yMin) * NumScale + 1
zNum = (zMax - zMin) * NumScale + 1

''' Functions '''

# Save to inviwo file
def saveToInv(fileName, kdeVal):
    dat = """RawFile: {file}
    Resolution: {dimX} {dimY} {dimZ}
    Format: {format}
    DataRange: {dataMin} {dataMax}
    ValueRange: {dataMin} {dataMax}
    BasisVector1: {V1} 0.0 0.0
    BasisVector2: 0.0 {V2} 0.0
    BasisVector3: 0.0 0.0 {V3}
    """

    with open(fileName + ".dat", 'w') as f:
        f.write(dat.format(
            file = fileName.split("/")[-1]+".raw", 
            format = "FLOAT32", 
            dimX = kdeVal.shape[0],
            dimY = kdeVal.shape[1],
            dimZ = kdeVal.shape[2],
            dataMin = kdeVal.min(),
            dataMax = kdeVal.max(),
            V1 = xMax - xMin,
            V2 = yMax - yMin,
            V3 = zMax - zMin
        ))

    arr = np.float32(kdeVal.T)
    arr.tofile(fileName + ".raw")

# Save to csv file
def saveToCsv(fileName, X, Y, Z, kdeVal):
    output = np.vstack([X.T.ravel(), Y.T.ravel(), Z.T.ravel(), kdeVal.T.ravel()])
    np.savetxt(fileName + ".csv", output.T, delimiter = ",", header = "x,y,z,density", comments='')

# Generate uniform grid
X, Y, Z = np.mgrid[xMin:xMax:xNum*1j, yMin:yMax:yNum*1j, zMin:zMax:zNum*1j]
grid = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])    

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

newComm = comm.Split(1, rank)
print(f"[INTS]: **************** This is rank {rank} from inTransit.")
newComm2 = comm.Split(1, rank)

print("[INTS]: **************** InTransit Initialization DONE")

# data = comm.recv(source=2, tag=11)
# print("-*-*-*-*-*-*-*-*-*-*-*- Receive Successful! -*-*-*-*-*-*-*-*-*-*-*-")
# print(data)
# pdata = comm.recv(source=1, tag=0)

# # Receive urid
# pointList = []
# for i in range(rank):
#     tempPoint = comm.recv(source=i, tag=0)
#     pointList = pointList + tempPoint
# print("[INTS]: **************** Grid Receive Successful!")

# # Convert list of coordinates to vtkPolyData
# pointVec = np.array(pointList)
# vtk_data_array = numpy_support.numpy_to_vtk(num_array=pointVec,deep=True,array_type=vtk.VTK_FLOAT)
# points = vtk.vtkPoints()
# points.SetData(vtk_data_array)
# poly = vtk.vtkPolyData()
# poly.SetPoints(points)

# vFilter = vtkVertexGlyphFilter() # Used to convert points to vertices in vtkPolyData
# cFilter = vtk.vtkCleanPolyData() # Used to remove the duplicate points

# vFilter.SetInputData(poly)
# vFilter.Update()
# cFilter.SetInputData(vFilter.GetOutput())
# cFilter.Update()
# newPoly = cFilter.GetOutput()
# # print(newPoly.GetNumberOfPoints())s
# print("[INTS]: **************** Set up Grid Successful!")

# Set up surface points
surfacePoints = []

for i in range(rank):
    tempPoint = comm.recv(source=i, tag=11)
    surfacePoints.append(tempPoint)

print("[INTS]: **************** Receive Surface Points Successful!")

overlapPoints = {}

for i in range(rank):
    for pointCoord in surfacePoints[i]:
        if pointCoord in overlapPoints:
            overlapPoints[pointCoord] += 1
        else:
            overlapPoints[pointCoord] = 1

# count = 0
# for key in overlapPoints:
#     if overlapPoints[key] > 2:
#         count += 1
#         print(f"key:{key}, value:{overlapPoints[key]}")
# print(f"Here is the count: {count}")

print("[INTS]: **************** Extract the Fake Surface Done!")

# Start to iterate
for rounds in range(NSTEP):
    minimaCandi = []
    maximaCandi = []
    for sendRank in range(rank):
        tempMini = comm.recv(source=sendRank, tag=22)
        minimaCandi.append(tempMini)
        tempMaxi = comm.recv(source=sendRank, tag=33)
        maximaCandi.append(tempMaxi)

    print(f"[INTS]: **************** Round[{rounds}] Receive Extrema Successful!")

    minimaCount = {}
    maximaCount = {}

    for i in range(rank):
        for pointCoord in minimaCandi[i]:
            if pointCoord in minimaCount:
                minimaCount[pointCoord] += 1
            else:
                minimaCount[pointCoord] = 1

    for i in range(rank):
        for pointCoord in maximaCandi[i]:
            if pointCoord in maximaCount:
                maximaCount[pointCoord] += 1
            else:
                maximaCount[pointCoord] = 1

    minima = []
    maxima = []

    for key in minimaCount:
        if key in overlapPoints:
            if minimaCount[key] == overlapPoints[key]:
                minima.append(key)
        else:
            minima.append(key)

    for key in maximaCount:
        if key in overlapPoints:
            if maximaCount[key] == overlapPoints[key]:
                maxima.append(key)
        else:
            maxima.append(key)

    print("number of minima:")
    print(len(minima))
    print("number of maxima:")
    print(len(maxima))

    print(f"[INTS]: **************** Round[{rounds}] Extract Extrema Done!")

    # Convert extrema
    exList = np.array(minima + maxima)
    extrema = exList.T

    #calculate density
    kernel = stats.gaussian_kde(extrema)
    # kernel = stats.gaussian_kde(extrema, bw_method=0.2 / np.std(extrema, ddof=1))
    # kernel = stats.gaussian_kde(extrema, bw_method='silverman')
    kdeVal = np.reshape(kernel(grid), X.shape)

    print(f"[INTS]: **************** Round[{rounds}] Calculate KDE Done!")

    # Save to file
    print(kdeVal.shape)
    fileName = "kdeData/kde_" + str(rounds)
    # saveToCsv(fileName, X, Y, Z, kdeVal)
    saveToInv(fileName, kdeVal)
    