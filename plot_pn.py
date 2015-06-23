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
numMoments = data[6]
data.dtype = numpy.double
ax = data[2]
ay = data[3]
bx = data[4]
by = data[5]

print sizeX
print sizeY
print numMoments
print ax
print ay
print bx
print by

M = data[7:sizeX*sizeY*numMoments+7:numMoments]
M = M / (2.0 * math.sqrt(math.pi))
M = numpy.reshape(M, (sizeX,sizeY))


# Plot heat map
plt.figure()
plt.imshow(M, origin='lower',extent=[ax,bx,ay,by])
plt.colorbar()


# Plot lineout
#plt.figure()
#y = M[math.floor(sizeX/2),:]
#dx = (ay-ax)/(sizeX-1);
#x = numpy.arange(ax, ay+dx/2, dx)
#plt.plot(x, y)


# Show plots
plt.show()
