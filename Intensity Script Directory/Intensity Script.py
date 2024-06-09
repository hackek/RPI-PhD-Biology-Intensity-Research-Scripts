# Script created using Python version 3.9.0
# Last Edited by: Kyle Hackett on 19SEP2023
# Note: Follow the instructions in the README.txt file

# Imports Used in this Script, and to auto-install all needed
import datetime     # In Python's standard library, so no need to download it
import os           # In Python's standard library, so no need to download it
import time         # In Python's standard library, so no need to download it
import xlsxwriter   # XlsxWriter version == 3.1.4

# Gets the current, inputs, outputs directories, and variables
directory = os.getcwd()
inputDirectory = directory + '\Inputs'
outputDirectory = directory + '\Outputs'
inputFormat = 'bar.list.txt'

# Defines the output Excel file and gets a sorted array of input file names
workbook = xlsxwriter.Workbook(outputDirectory + '\data_{}.xlsx'.format(datetime.date.today()))
worksheet = workbook.add_worksheet()
inputList = os.listdir(inputDirectory)
for file in range(0, len(inputList)):
    inputList[file] = int(inputList[file].strip(inputFormat))
inputList.sort()

# Parses through the data of every input file
for filename in range(0, len(inputList)):
    # Opens the current input file
    f = open(os.path.join(inputDirectory, str(inputList[filename]) + inputFormat))
    # Reads the current input file contents
    fileContents = (f.read()).split()
    assignmentValues = []
    intensityValues = []
    # Parses the input file contents for Assignment and Intensity values
    for j in range(5, len(fileContents)):
        if (j+3)%4 == 0:
            assignmentValues.append(fileContents[j].strip('N-H'))
        if (j)%4 == 0:
            intensityValues.append(fileContents[j])
    # Formats the output file's Assignment Values
    if worksheet.table.get(1, None) == None:
        for k in range(0, len(assignmentValues)):
            worksheet.write(1, k+1, int(assignmentValues[k]))
    # Formats the output file's Intensity Values
    worksheet.write(filename+2, 0, int(str(inputList[filename]).strip(inputFormat)))
    for l in range (0, len(intensityValues)):
        worksheet.write(filename+2, l+1, int(intensityValues[l]))

# Saves the .xlsx file
workbook.close()