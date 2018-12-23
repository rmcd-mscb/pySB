from .centerline import centerline
import numpy as np
import vtk
from vtk.util import numpy_support
import matplotlib.pyplot as plt
class COGridGen():
    def __init__(self, center_shp, tension, nx, ny, width):
        self.cl = centerline(center_shp)
        self.tension = tension
        self.width = width
        self.nx = nx
        self.ny = ny
        self.xgrid = np.zeros(self.nx*self.ny, dtype=np.double)
        self.ygrid = np.zeros(self.nx*self.ny, dtype=np.double)
        self.sgrid = vtk.vtkStructuredGrid()
        self.sg_points = vtk.vtkPoints()
        self.__buildXYGrid()

    def description(self):
        return r"Uses centerline {} and width {} with nx = {} and ny = {} " \
               "to generate curvilinear orthogonal grid".format(self.cl.getCLShapeFile(), self.width, self.nx, self.ny)

    def getGridDims(self):
        return self.nx, self.ny

    def __buildXYGrid(self):
        clx, cly = self.cl.getinterppts(self.nx, self.tension)
        delt = np.double(self.width/(self.ny - 1))
        nm = np.int((self.ny+1)/2)

        self.sgrid.SetDimensions(self.nx, self.ny, 1)
        sg_points = vtk.vtkPoints()

        for i in range(0,self.nx):
            for j in range(0,self.ny):
                index = (i*self.ny) + j
                self.xgrid[index] = clx[i] + delt * (nm-j-1) * np.sin(self.cl.getphiinterp(i))
                self.ygrid[index] = cly[i] - delt * (nm-j-1) * np.cos(self.cl.getphiinterp(i))
                self.sg_points.InsertNextPoint(self.xgrid[index], self.ygrid[index], 0.0)

        self.sgrid.SetPoints(self.sg_points)
    def getXYGrid(self):
        return self.xgrid, self.ygrid

    def getGrid(self):
        return self.cl, self.sgrid

    def plotmapgrid(self, valname):
        vtknodes = self.sgrid.GetPoints().GetData()
        npnodes = numpy_support.vtk_to_numpy(vtknodes)
        x,y,z = npnodes[:,0], npnodes[:,1], npnodes[:,2]
        nx, ny, nz = self.sgrid.GetDimensions()
        x = np.reshape(x, (nx, ny))
        y = np.reshape(y, (nx, ny))
        val = self.sgrid.GetPointData().GetArray(valname)
        val = np.reshape(val, (nx, ny))
        plt.contourf(x,y,val)
        plt.show()