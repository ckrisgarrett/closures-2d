import numpy
import matplotlib.pyplot as plt
import math
import sys

# Get the data
filename = sys.argv[1]
data = numpy.fromfile(filename)
data.dtype = numpy.int64
sizeX = data[0]
sizeY = data[1]
numQuadPoints = data[6]
data.dtype = numpy.double
ax = data[2]
ay = data[3]
bx = data[4]
by = data[5]

quadWeights = data[7:numQuadPoints+7]

offset = numQuadPoints+7
M = numpy.zeros([sizeX, sizeY]);
for i in range(0,sizeX):
    for j in range(0,sizeY):
        k1 = (i + j * sizeY) * numQuadPoints + offset;
        k2 = (i + j * sizeY) * numQuadPoints + offset + numQuadPoints;
        M[i][j] = numpy.dot(data[k1:k2], quadWeights) / (4.0 * math.pi);
    

# Plot heat map
plt.figure()
plt.imshow(M, origin='lower',extent=[ax,bx,ay,by])
plt.colorbar()


# Plot lineout
#plt.figure()
#y = M[math.floor(sizeX/2),:]
#dx = (bx-ax)/(sizeX-1);
#x = numpy.arange(ax, bx+dx/2, dx)
#plt.plot(x, y)


# Show plots
plt.show()
