import numpy
import copy
import struct
import math
import fractions

import matplotlib.pyplot as plt


class Grid(object):
    def _generalize(self):
        # change class
        self.__class__ = Grid
        # remove format-specific members
        try:
            del self.numQuadPoints
        except AttributeError:
            pass
        try:
            del self.quadWeights
        except AttributeError:
            pass
        try:
            del self.matrix
        except AttributeError:
            pass
        try:
            del self.numMoments
        except AttributeError:
            pass
        return self

    def __getitem__(self, key):
        return self.grid[key]

    def __abs__(self):
        accumulator = 0.0
        for x in range(self.sizeX):
            for y in range(self.sizeY):
                accumulator += self[x, y]**2
        return math.sqrt(accumulator)

    def __neg__(self):
        result = copy.copy(self)
        result.grid = numpy.copy(self.grid)
        result.grid *= -1
        return result._generalize()

    def __add__(self, other):
        if self.sizeX < other.sizeX and self.sizeY < other.sizeY:
            return other + self
        elif self.sizeX >= other.sizeX and self.sizeY >= other.sizeY:
            result = copy.copy(self)
            result.grid = numpy.zeros([self.sizeX, self.sizeY])
            factor = self.sizeX / other.sizeX
            for x in range(other.sizeX):
                for y in range(other.sizeY):
                    sX = (other[(x + 1) % other.sizeX, y] - other[(x - 1) % other.sizeX, y]) / 2.0
                    sY = (other[x, (y + 1) % other.sizeY] - other[x, (y - 1) % other.sizeY]) / 2.0
                    for k in range(int(factor)):
                        for l in range(int(factor)):
                            x2 = x * factor + k
                            y2 = y * factor + l
                            interpolation = other[x, y] + (-0.5 + (0.5 + k) / factor) * sX + (-0.5 + (0.5 + l) / factor) * sY
                            result.grid[x2, y2] = self[x2, y2] + interpolation
            return result._generalize()
        else:
            raise ArithmeticError("Size mismatch")

    def __sub__(self, other):
        return (-other) + self

    def plot(self):
        plt.figure()
        plt.imshow(self.grid, origin="lower", interpolation="nearest",
            extent=[self.ax, self.bx, self.ay, self.by])
        plt.colorbar()
        plt.show()

    def logplot(self):
        plt.figure()
        plt.imshow(numpy.log(self.grid), origin="lower", interpolation="nearest",
            extent=[self.ax, self.bx, self.ay, self.by])
        plt.colorbar()
        plt.show()

    @property
    def mass(self):
        accumulator = 0.0
        for x in range(self.sizeX):
            for y in range(self.sizeY):
                accumulator += self[x, y] * self.dx * self.dy
        return accumulator

class Composite(Grid):
    def __init__(self, parser, filename_pattern, numNodesX, numNodesY):
        self.subgrids = {}
        for x in range(numNodesX):
            for y in range(numNodesY):
                self.subgrids[x,y] = parser(filename_pattern.format(x * numNodesY + y))
        self.sizeX = sum([self.subgrids[i,0].sizeX for i in range(numNodesX)])
        self.sizeY = sum([self.subgrids[0,i].sizeY for i in range(numNodesY)])
        self.ax = self.subgrids[0,0].ax
        self.ay = self.subgrids[0,0].ay
        self.bx = self.subgrids[0,0].bx
        self.by = self.subgrids[0,0].by
        self.dx = (self.bx - self.ax) / self.sizeX
        self.dy = (self.by - self.ay) / self.sizeY
        self.grid = numpy.zeros([self.sizeX, self.sizeY])
        for x in range(numNodesX):
            for y in range(numNodesY):
                xOffset = sum([self.subgrids[i,0].sizeX for i in range(x)])
                yOffset = sum([self.subgrids[0,i].sizeY for i in range(y)])
                self.grid[xOffset:xOffset + self.subgrids[x,y].sizeX,
                    yOffset:yOffset + self.subgrids[x,y].sizeY] = self.subgrids[x,y].grid

class Kinetic(Grid):
    def __init__(self, filename):
        with open(filename, "rb") as f:
            header = struct.Struct("=2q4dq")
            data = f.read()
            (self.sizeX,
             self.sizeY,
             self.ax,
             self.ay,
             self.bx,
             self.by,
             self.numQuadPoints) = header.unpack(data[:header.size])
            self.quadWeights = numpy.fromstring(data[header.size:header.size +
                8 * self.numQuadPoints], dtype=numpy.double)
            self.matrix = numpy.fromstring(data[header.size + 8 * self.numQuadPoints:],
                dtype=numpy.double).reshape(self.sizeX, self.sizeY, self.numQuadPoints)
            self.grid = numpy.zeros([self.sizeX, self.sizeY])
            for x in range(self.sizeX):
                for y in range(self.sizeY):
                    self.grid[x, y] = numpy.dot(self.quadWeights, self.matrix[x, y])
            self.grid /= (4.0 * math.pi)
            self.dx = (self.bx - self.ax) / self.sizeX
            self.dy = (self.by - self.ay) / self.sizeY

class Moment(Grid):
    def __init__(self, filename):
        with open(filename, "rb") as f:
            header = struct.Struct("=2q4dq")
            data = f.read()
            (self.sizeX,
             self.sizeY,
             self.ax,
             self.ay,
             self.bx,
             self.by,
             self.numMoments) = header.unpack(data[:header.size])
            self.matrix = numpy.fromstring(data[header.size:],
                dtype=numpy.double).reshape(self.sizeX, self.sizeY, self.numMoments)
            self.grid = numpy.zeros([self.sizeX, self.sizeY])
            for x in range(self.sizeX):
                for y in range(self.sizeY):
                    self.grid[x, y] = self.matrix[x, y, 0]
            self.grid /= (2.0 * math.sqrt(math.pi))
            self.dx = (self.bx - self.ax) / self.sizeX
            self.dy = (self.by - self.ay) / self.sizeY
