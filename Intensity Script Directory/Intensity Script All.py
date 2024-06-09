# Script created using Python version 3.11.0
# Last Edited by: Kyle Hackett on 19SEP2023
# Note: Follow the instructions in the README.txt file

# Imports Used in this Script, and to auto-install all needed
import datetime                 # In Python's standard library
import os                       # In Python's standard library
import time                     # In Python's standard library
import math                     # In Python's standard library
import xlsxwriter               # XlsxWriter version == 3.1.4
import matplotlib.pyplot as plt # matplotlib version == 3.8.0
import numpy as np              # numpy version == 1.26.0
from scipy.optimize import _lsq # scipy version == 1.11.2

# Gets the current, inputs, outputs directories, and variables
directory = os.getcwd()
inputDirectory = directory + '\Inputs'
outputDirectory = directory + '\Outputs'
inputFormat = 'bar.list.txt'

# Defines the output Excel file and gets a sorted array of input file names
inputWorkbook = xlsxwriter.Workbook(outputDirectory + '\data_{}.xlsx'.format(datetime.date.today()))
inputWorksheet = inputWorkbook.add_worksheet()
fileInputList = os.listdir(inputDirectory)
for file in range(0, len(fileInputList)):
    fileInputList[file] = int(fileInputList[file].strip(inputFormat))
fileInputList.sort()
inputList = np.array(fileInputList)

# Creates an array of input contents
assignmentValues = []
intensityValues = []
allValuesByRow = [[0]]
allValuesByCol = []

# Parses through the data of every input file
for filename in range(0, len(inputList)):
    # Opens the current input file
    f = open(os.path.join(inputDirectory, str(inputList[filename]) + inputFormat))
    # Reads the current input file contents
    fileContents = (f.read()).split()
    # Parses the input file contents for Assignment and Intensity values
    for j in range(5, len(fileContents)):
        if (j+3)%4 == 0 and filename == 0:
            assignmentValues.append(int(fileContents[j].strip('N-H')))
        if (j)%4 == 0:
            intensityValues.append(int(fileContents[j]))
    # Formats the output file's Assignment Values
    if inputWorksheet.table.get(1, None) == None:
        allValuesByRow[0] += assignmentValues
        for k in range(0, len(assignmentValues)):
            inputWorksheet.write(1, k+1, assignmentValues[k])
    # Formats the output file's Intensity Values
    inputWorksheet.write(filename+2, 0, inputList[filename])
    tempRow = [inputList[filename]]
    for l in range (0, len(intensityValues)):
        inputWorksheet.write(filename+2, l+1, intensityValues[l])
        tempRow.append(intensityValues[l])
    allValuesByRow.append(tempRow)

# Constructs a 2D array using columns
for col in range(0, len(allValuesByRow[0])):
    tempStorage = []
    for row in range(0, len(allValuesByRow)):
        tempStorage.append(allValuesByRow[row][col])
    allValuesByCol.append(tempStorage)
allValuesByCol = np.array(allValuesByCol)

# print(allValuesByRow) # Formatted correctly
# print(allValuesByCol) # Formatted correctly
# print(assignmentValues) # Formatted correctly
# print(intensityValues) # Formatted correctly




######## This code fits equilibrium 2state pressure data when x(1)=1 ########

# Initialize constants
TEMP = 298
RGC = 0.00197
RGP = 83.1451

# Initialize input variables
inputData = allValuesByRow # All data cells in a 2D array
resName = assignmentValues # Values of the header row only
rowTotal = len(inputList) # Number of rows minus the header row, # of pressures
colTotal = len(assignmentValues) # Number of cols - the header col, measured residue #, curves

# print(inputData) # Correctly formatted
# print(resName) # Correctly formatted
# print(rowTotal) # Correctly formatted
# print(colTotal) # Correctly formatted

# Initialize output variables
IntF = [0 for i in range(colTotal)] # Named from Matlab
sigma_IntF = [0 for i in range(colTotal)] # Named from Matlab
IntU = [0 for i in range(colTotal)] # Named from Matlab
sigma_IntU = [0 for i in range(colTotal)] # Named from Matlab
DeltaG = [0 for i in range(colTotal)] # Named from Matlab
sigma_DG = [0 for i in range(colTotal)] # Named from Matlab
DeltaV = [0 for i in range(colTotal)] # Named from Matlab
sigma_DV = [0 for i in range(colTotal)] # Named from Matlab
fracintfit = [[0]*rowTotal for i in range(colTotal)] # Named from Matlab
pressureval= [[0]*rowTotal for i in range(colTotal)] # Named from Matlab
Ipcalcres = [0 for i in range(colTotal)] # Named from Matlab
Ipcalcresn = [0 for i in range(colTotal)] # Named from Matlab
Ipcalc = [0 for i in range(rowTotal)] # Named from Matlab

# Initialize Excel output files for Parameters and Results
pWorkbook = xlsxwriter.Workbook(outputDirectory + '\data_{}_p.xlsx'.format(datetime.date.today()))
pWorksheet = pWorkbook.add_worksheet()
rWorkbook = xlsxwriter.Workbook(outputDirectory + '\data_{}_r.xlsx'.format(datetime.date.today()))
rWorksheet = rWorkbook.add_worksheet()

#***
#  Options tested:
#   - n, n      => success (no plateau fix)
#   - y, n      => success (no plateau fix)
#   - n, y, y   => success (no plateau fix)
#   - n, y, n   => success (plateau fix every single one)
#   - y, y, y   => success (no plateau fix)
#   - y, y, n   => success (plateau fix every single one)
#***
# Gets the parameter guesses from the user
Kuo = float(input('Enter the initial Ku (0.01 is a good guess): '))
DVu = float(input('Enter the intial deltaVu (50 is a good guess): '))
DVu /= (RGP*TEMP)
finalPlateau = input('Set final plateau to 0 (Y or N)?: ').lower() #str1
initialPlateau = input('Fix initial plateau (Y or N)?: ').lower() #str2
if initialPlateau == 'y':
    allPlateau = input('Fix all initial plateau values at 105 percent of average of first and second values (Y or N)?: ').lower() #str3


# Initialize Functions
# Function for on-hover coordinate display when viewing scatter plots
def motion_hover(event):
    annotation_visibility = annotation.get_visible()
    if event.inaxes == ax:
        is_contained, annotation_index = scatter.contains(event)
        if is_contained:
            data_point_location = scatter.get_offsets()[annotation_index['ind'][0]]
            annotation.xy = data_point_location

            text_label = '({}, {})'.format(data_point_location[0], data_point_location[1])
            annotation.set_text(text_label)

            annotation.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if annotation_visibility:
                annotation.set_visible(False)
                fig.canvas.draw_idle()



# Calculates and plots each point for every entry then assigns to r file
for c in range(1, colTotal):
    # Calculates graph x values, y values, and midpoint
    xdata = allValuesByCol[0][1:]
    ydata = allValuesByCol[c][1:]
    coordinates = np.array([xdata, ydata])
    Int0 = (ydata[0] + ydata[1]) / 2
    Iof = 0
    Int0 = 0
    # Sets an initial guess for low pressue plateau values
    if finalPlateau == 'y':
        IntHP = 0
    else:
        IntHP = ydata[-1]

    # Fits the data for residue
    if initialPlateau == 'y':
        if allPlateau == 'n':
            # Scatter plot of x and y datas with labels:
            fig, ax = plt.subplots()
            scatter = plt.scatter(
                x=xdata,
                y=ydata
            )
            plt.title('High Pressure Dentauration')
            plt.xlabel('Pressure in Bar')
            plt.ylabel('Intensity')

            # Create annotations of data coordinates
            annotation = ax.annotate(
                text = '',
                xy=(0, 0),
                xytext=(-65, 15), # distance from x, y
                textcoords='offset points',
                bbox={'boxstyle': 'round', 'fc': 'w'},
                arrowprops={'arrowstyle': '->'}
            )
            annotation.set_visible(False)
            # Show data annotations on data point hover
            fig.canvas.mpl_connect('motion_notify_event', motion_hover)

            # Display graph and await user initial plateau fix value
            plt.show()
            Iof = int(input('Enter initial plateau value to fix:   '))
            Int0 = Iof
        else:
            Int0 *= 1.05
            Iof = Int0
    else:
        Iof = 0

    # Calculate curve fitting
    if Iof == 0:
        if finalPlateau == 'n':
            p0 = [Int0, IntHP, Kuo, DVu] # Coefficient variables
            totalDiff = lambda x0: (
                (x0[0]+x0[1]*x0[2]*math.e**(x0[3]*coordinates[0])) / (1 + (x0[2]*math.e**(x0[3]*coordinates[0]))) - coordinates[1])
            fitTest = _lsq.least_squares(totalDiff, p0) # only inputs needed
            # Output fitTest.x returns the solutions found (1)
            # Output fitTest.fun returns the residuals found (3)
            # Output fitTest.jac returns the Jacobian matrix fount (7)
            # Output fitTest.status returns why the algorithm stopped (1-4=good)
        else:
            p0 = [Int0, Kuo, DVu]
            totalDiff = lambda x0: ( x0[0]*(1 / (1+ (x0[1]*math.e**(x0[2]*coordinates[0])))) - coordinates[1] )
            fitTest = _lsq.least_squares(totalDiff, p0)
    else:
        if finalPlateau == 'n': # Working?
            p0 = [IntHP, Kuo, DVu]
            totalDiff = lambda x0: ( Iof+x0[0]*x0[1]*math.e**(x0[2]*coordinates[0])) / (1 + (x0[1]*math.e**(x0[2]*coordinates[0])) - coordinates[1] )
            fitTest = _lsq.least_squares(totalDiff, p0)
        else:
            p0 = [Kuo, DVu]
            totalDiff = lambda x0: ( Iof * (1 / (1+ (x0[0]*math.e**(x0[1]*coordinates[0])))) - coordinates[1] )
            fitTest = _lsq.least_squares(totalDiff, p0)

    print(coordinates)
    print(fitTest)
    # Plot the fit and data
        # line 200: plot, scatter, assign legend/title, file name, .tiff

# print(Kuo) # Correctly Formatted
# print(DVu) # Correctly Formatted
# print(finalPlateau) # Correctly Formatted
# print(initialPlateau) # Correctly Formatted
# print(allPlateau) # Correctly Formatted

# print(IntF) # Correctly Formatted
# print(sigma_IntF) # Correctly Formatted
# print(IntU) # Correctly Formatted
# print(sigma_IntU) # Correctly Formatted
# print(DeltaG) # Correctly Formatted
# print(sigma_DG) # Correctly Formatted
# print(DeltaV) # Correctly Formatted
# print(sigma_DV) # Correctly Formatted
# print(fracintfit) # Correctly Formatted
# print(pressureval) # Correctly Formatted
# print(Ipcalcres) # Correctly Formatted
# print(Ipcalcresn) # Correctly Formatted
# print(Ipcalc) # Correctly Formatted

# time.sleep(5)

# output p file formatting
if Iof == 0: # p file
    if finalPlateau == 'n':
        Int_Fold = IntF
        sigma_IntF = sigma_IntF
        Int_Unf = IntU
        sigma_IntU = sigma_IntU
        DeltaGu = DeltaG
        sigma_DG = sigma_DG
        DeltaVu = DeltaV
        sigma_DV = sigma_DV
        pWorksheet.add_table(0, 0, len(assignmentValues), 8, {
            'columns': [
                {'header': 'resname'},
                {'header': 'Int_Fold'},
                {'header': 'sigma_IntF'},
                {'header': 'Int_Unf'},
                {'header': 'sigma_IntU'},
                {'header': 'DeltaGu'},
                {'header': 'sigma_DG'},
                {'header': 'DeltaVu'},
                {'header': 'sigma_DV'}
            ]
        })
        pWorksheet.write_column(1, 0, assignmentValues)
        pWorksheet.write_column(1, 1, Int_Fold)
        pWorksheet.write_column(1, 2, sigma_IntF)
        pWorksheet.write_column(1, 3, Int_Unf)
        pWorksheet.write_column(1, 4, sigma_IntU)
        pWorksheet.write_column(1, 5, DeltaGu)
        pWorksheet.write_column(1, 6, sigma_DG)
        pWorksheet.write_column(1, 7, DeltaVu)
        pWorksheet.write_column(1, 8, sigma_DV)
    else:
        Int_Fold = IntF
        sigma_IntF = sigma_IntF
        DeltaGu = DeltaG
        sigma_DG = sigma_DG
        DeltaVu = DeltaV
        sigma_DV = sigma_DV
        pWorksheet.add_table(0, 0, len(assignmentValues), 6, {
            'columns': [
                {'header': 'resname'},
                {'header': 'Int_Fold'},
                {'header': 'sigma_IntF'},
                {'header': 'DeltaGu'},
                {'header': 'sigma_DG'},
                {'header': 'DeltaVu'},
                {'header': 'sigma_DV'}
            ]
        })
        pWorksheet.write_column(1, 0, assignmentValues)
        pWorksheet.write_column(1, 1, Int_Fold)
        pWorksheet.write_column(1, 2, sigma_IntF)
        pWorksheet.write_column(1, 3, DeltaGu)
        pWorksheet.write_column(1, 4, sigma_DG)
        pWorksheet.write_column(1, 5, DeltaVu)
        pWorksheet.write_column(1, 6, sigma_DV)
else:
    if finalPlateau == 'n':
        Int_Unf = IntU
        sigma_IntU = sigma_IntU
        DeltaGu = DeltaG
        sigma_DG = sigma_DG
        DeltaVu = DeltaV
        sigma_DV = sigma_DV
        pWorksheet.add_table(0, 0, len(assignmentValues), 6, {
            'columns': [
                {'header': 'resname'},
                {'header': 'Int_Unf'},
                {'header': 'sigma_IntU'},
                {'header': 'DeltaGu'},
                {'header': 'sigma_DG'},
                {'header': 'DeltaVu'},
                {'header': 'sigma_DV'}
            ]
        })
        pWorksheet.write_column(1, 0, assignmentValues)
        pWorksheet.write_column(1, 1, Int_Unf)
        pWorksheet.write_column(1, 2, sigma_IntU)
        pWorksheet.write_column(1, 3, DeltaGu)
        pWorksheet.write_column(1, 4, sigma_DG)
        pWorksheet.write_column(1, 5, DeltaVu)
        pWorksheet.write_column(1, 6, sigma_DV)
    else:
        DeltaGu = DeltaG
        sigma_DG = sigma_DG
        DeltaVu = DeltaV
        sigma_DV = sigma_DV
        pWorksheet.add_table(0, 0, len(assignmentValues), 4, {
            'columns': [
                {'header': 'resname'},
                {'header': 'DeltaGu'},
                {'header': 'sigma_DG'},
                {'header': 'DeltaVu'},
                {'header': 'sigma_DV'}
            ]
        })
        pWorksheet.write_column(1, 0, assignmentValues)
        pWorksheet.write_column(1, 1, DeltaGu)
        pWorksheet.write_column(1, 2, sigma_DG)
        pWorksheet.write_column(1, 3, DeltaVu)
        pWorksheet.write_column(1, 4, sigma_DV)

time.sleep(10)

# Sets the column width of all output files to be readable
inputWorksheet.autofit()
pWorksheet.autofit()
rWorksheet.autofit()

# Saves all used Excel files
inputWorkbook.close()
rWorkbook.close()
pWorkbook.close()