import numpy as np


# Base Power
Sb = 100.0

# Method variables
tolerance = 0.00001
iterration = 0


# read data from file and split them properly

readData = open("inputdata.txt", "r").read()

dataSection = readData.split("----\n")
print('dataSection[0]')
print(dataSection[0])
print('dataSection[1]')
print(dataSection[1])
# Create bus objects
buses = dict()
busdata = dataSection[0].split('\n')
print('busdata[0]')
print(busdata[0])
print('busdata[1]')
print(busdata[1])
print('busdata[2]')
print(busdata[2])
