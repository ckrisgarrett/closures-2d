import sys
import numpy

# Check number of arguments.
if len(sys.argv) < 2:
    print 'Usage: python put_data_together.py (<time>)+'
    quit()


# Get data from input.deck
inputDeckFile = open('input.deck')

line = inputDeckFile.readline().split()
solver = line[1]
line = inputDeckFile.readline().split()
NX = int(line[1])
line = inputDeckFile.readline().split()
NY = int(line[1])
line = inputDeckFile.readline().split()
NPartX = int(line[1])
line = inputDeckFile.readline().split()
NPartY = int(line[1])

print "solver = ", solver
print "NX = ", NX
print "NY = ", NY
print "NPartX = ", NPartX
print "NPartY = ", NPartY 

inputDeckFile.close()


# For each time write out data.
times = sys.argv[1:len(sys.argv)]
for time in times:
    print time

    # Get array of input files.
    inFile = []
    for i in range(0, NPartX):
        for j in range(0, NPartY):
            if solver == 'kinetic':
                filename = 'out_{0}_{1}.sn'.format(time, i * NPartY + j)
            else:
                filename = 'out_{0}_{1}.pn'.format(time, i * NPartY + j)
            inFile.append(open(filename, "rb"))

    # Create output file
    if solver == 'kinetic':
        filename = 'outall_{0}.sn'.format(time)
    else:
        filename = 'outall_{0}.pn'.format(time)
    outFile = open(filename, 'wb')
    
    # Read header data.
    sizeX = range(0, NPartX)
    sizeY = range(0, NPartY)
    ax = 0
    ay = 0
    bx = 0
    by = 0
    stride = 0
    header2 = 0
    for i in range(0, NPartX):
        for j in range (0, NPartY):
            sizeX[i] = numpy.fromfile(inFile[i * NPartY + j], numpy.int64, 1)
            sizeY[j] = numpy.fromfile(inFile[i * NPartY + j], numpy.int64, 1)
            ax = numpy.fromfile(inFile[i * NPartY + j], numpy.float64, 1)
            ay = numpy.fromfile(inFile[i * NPartY + j], numpy.float64, 1)
            bx = numpy.fromfile(inFile[i * NPartY + j], numpy.float64, 1)
            by = numpy.fromfile(inFile[i * NPartY + j], numpy.float64, 1)
            stride = numpy.fromfile(inFile[i * NPartY + j], numpy.int64, 1)
            if solver == 'kinetic':
                header2 = inFile[i * NPartY + j].read(stride * 8)
    
    # Write data
    header = numpy.int64([NX, NY])
    outFile.write(header)
    header = numpy.float64([ax, ay, bx, by])
    outFile.write(header)
    header = numpy.int64([stride])
    outFile.write(header)
    if solver == 'kinetic':
        outFile.write(header2);
    for i in range(0, NPartX):
        for l in range(0, sizeX[0]):
            for j in range(0, NPartY):
                data = inFile[i * NPartY + j].read(stride * 8 * sizeY[j])
                outFile.write(data)

    # Close files.
    outFile.close()
    for i in range(0, NPartX):
        for j in range(0, NPartY):
            inFile[i * NPartY + j].close()

