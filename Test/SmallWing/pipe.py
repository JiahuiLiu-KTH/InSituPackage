# --------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt = 0
forceOutputAtFirstCall = False

# Global screenshot output options
imageFileNamePadding = 0
rescale_lookuptable = False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays = False

# a root directory under which all Catalyst output goes
rootDirectory = ''

# makes a cinema D index table
make_cinema_table = False

# TODO: Change this number for the real case (to e.g. 10)
global_freq = 2

neighborsBook = []

# --------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.6.0
# --------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing
# import numpy as np
import vtk
import time
from mpi4py import MPI

# def find neighbor function
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

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
    def _CreatePipeline(coprocessor, datadescription):
        class Pipeline:
            # state file generated using paraview version 5.6.0

            # ----------------------------------------------------------------
            # setup the data processing pipelines
            # ----------------------------------------------------------------

            # trace generated using paraview version 5.6.0
            #
            # To ensure correct image size when batch processing, please search
            # for and uncomment the line `# renderView*.ViewSize = [*,*]`

            #### disable automatic camera reset on 'Show'
            paraview.simple._DisableFirstRenderCameraReset()

            # create a new 'PVTrivialProducer'
            # create a producer from a simulation input
            input = coprocessor.CreateProducer(datadescription, 'input')

            ''' ---------------------------------------------------------------------------------------------------- '''

            # remove duplicate points
            cleanGrid = CleantoGrid(Input=input)
            cleanGrid.UpdatePipeline()

            # Get the actual unstructured grid
            uGrid = cleanGrid.GetClientSideObject().GetOutputDataObject(0)

            numPoints = uGrid.GetNumberOfPoints()

            # Find neighbours
            for pointId in range(numPoints):
                neighbors = FindNeighborPoints(pointId, uGrid)
                neighborsBook.append(neighbors)

            # Find the surface points
            filter = vtk.vtkDataSetSurfaceFilter()
            filter.SetInputData(uGrid)
            filter.PassThroughPointIdsOn()
            filter.Update()
            sGrid = filter.GetOutput()

            # Useless Print Out
            worldRank = MPI.COMM_WORLD.Get_rank()
            print(f"[PIPE]: ++++++++++++++++ [Rank {worldRank}] Pipe Initialization DONE")

            inTransitRank = MPI.COMM_WORLD.Get_size() - 1

            # # send grid
            # if MPI.COMM_WORLD.Get_rank() < inTransitRank:
            #     n = uGrid.GetNumberOfPoints()
            #     points = []
            #     for i in range(n):
            #         pointCoord = uGrid.GetPoint(i)
            #         points.append(pointCoord)
            #     MPI.COMM_WORLD.send(points, dest=4, tag=0)
            #     print(f"[PIPE]: ++++++++++++++++ [Rank {worldRank}] Grid Send Successful!")

            surfacePoint = []

            n = sGrid.GetNumberOfPoints()
            for i in range(n):
                pointCoord = sGrid.GetPoint(i)
                surfacePoint.append(pointCoord)

            MPI.COMM_WORLD.send(surfacePoint, dest=inTransitRank, tag=11)

            print(f"[PIPE]: ++++++++++++++++ [Rank {worldRank}] Send Surface Points Successful!")

            ''' --------------------------------------------------------------------------------------------------- '''
            # ----------------------------------------------------------------
            # finally, restore active source
            SetActiveSource(input)
            # SetActiveSource(cleanGrid)
            # SetActiveSource(uGrid)
            # ----------------------------------------------------------------

            # Now any catalyst writers
            # xMLPUnstructuredGridWriter1 = servermanager.writers.XMLPUnstructuredGridWriter(Input=input)
            # coprocessor.RegisterWriter(xMLPUnstructuredGridWriter1, filename='data/input_%t.pvtu', freq=global_freq)

        return Pipeline()

    class CoProcessor(coprocessing.CoProcessor):
        def CreatePipeline(self, datadescription):
            self.Pipeline = _CreatePipeline(self, datadescription)

    coprocessor = CoProcessor()
    # these are the frequencies at which the coprocessor updates.
    freqs = {'input': [global_freq]}
    coprocessor.SetUpdateFrequencies(freqs)
    if requestSpecificArrays:
        arrays = [['pressure', 0], ['temperature', 0], ['velocity', 0]]
        coprocessor.SetRequestedArrays('input', arrays)
    coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt, forceOutputAtFirstCall)

    if rootDirectory:
        coprocessor.SetRootDirectory(rootDirectory)

    if make_cinema_table:
        coprocessor.EnableCinemaDTable()

    return coprocessor


# --------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

# --------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)


# ------------------------ Processing method ------------------------

# TODO: Change this number for the real case (to e.g. 10 and 1)
samples_per_mode = 2
break_between_modes = 2

# Mode: Save data
startTime = 1 + 27 * global_freq * (samples_per_mode - 1 + break_between_modes + 1)
modeStartTimes = [startTime]
modeEndTimes = [startTime + global_freq * (samples_per_mode - 1)]
endTime = modeEndTimes[-1]
numModes = len(modeEndTimes)


def setTimeSettings(set_samples_per_mode, set_break_between_modes, set_startTime):
    global samples_per_mode, break_between_modes
    global modeStartTimes, modeEndTimes
    global startTime, endTime
    samples_per_mode = set_samples_per_mode
    break_between_modes = set_break_between_modes
    startTime = set_startTime
    modeStartTimes = [startTime]
    modeEndTimes = [startTime + global_freq * (samples_per_mode - 1)]
    endTime = modeEndTimes[-1]


def printModes():
    print("Time: [%i, %i],\tData: All, Mode: Save." % (modeStartTimes[0], modeEndTimes[0]))


def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    '''--------------------------------------------------------------------------------------------------------------'''
    # if (MPI.COMM_WORLD.Get_rank() == 2):
    if (MPI.COMM_WORLD.Get_rank() > -1):
        input = coprocessor.Pipeline.input
        cleanGrid = CleantoGrid(Input=input)
        cleanGrid.UpdatePipeline()
        uGrid = cleanGrid.GetClientSideObject().GetOutputDataObject(0)

        numPoints = uGrid.GetNumberOfPoints()

        maxima = []
        minima = []

        for pointId in range(numPoints):
            isMaxima = True
            isMinima = True

            currentValue = uGrid.GetPointData().GetArray(1).GetTuple(pointId)[1]
            neighbors = neighborsBook[pointId]

            for ids in neighbors:
                tempVal = uGrid.GetPointData().GetArray(1).GetTuple(ids)[1]
                if tempVal >= currentValue:
                    isMaxima = False
                    break

            for ids in neighbors:
                tempVal = uGrid.GetPointData().GetArray(1).GetTuple(ids)[1]
                if tempVal <= currentValue:
                    isMinima = False
                    break

            if isMaxima:
                maxima.append(pointId)

            if isMinima:
                minima.append(pointId)

        # print("number of minima:")
        # print(len(minima))
        # print("number of maxima:")
        # print(len(maxima))

        minimaList = []
        maximaList = []

        for pid in minima:
            minimaList.append(uGrid.GetPoint(pid))

        for pid in maxima:
            maximaList.append(uGrid.GetPoint(pid))

        inTransitRank = MPI.COMM_WORLD.Get_size() - 1
        MPI.COMM_WORLD.send(minimaList, dest=inTransitRank, tag=22)
        MPI.COMM_WORLD.send(maximaList, dest=inTransitRank, tag=33)

        # print(uGrid.GetPointData().GetArray(0).GetTuple(0)[0])
        print(f"[PIPE]: ++++++++++++++++ [Rank {MPI.COMM_WORLD.Get_rank()}] Send Extrema Successful!")
    '''--------------------------------------------------------------------------------------------------------------'''
    # Write output data, if appropriate.
    # coprocessor.WriteData(datadescription)

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    # coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
    #                         image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
