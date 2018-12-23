from .COGridGen import COGridGen
import geopandas as gpd
import numpy as np
from shapely.geometry import shape
import vtk
from vtk.util import numpy_support
import matplotlib.pyplot as plt
class mapper():
    """

    """
    def __init__(self, grid, valname, ts, tn, wghtpowval, pointset, shp_has_geo='False', XCOL = '', YCOL='', ZCOL=''):
        self.cl, self.grid = grid.getGrid()
        self.valname = valname
        self.ts = ts
        self.tn = tn
        self.wghtpowval = wghtpowval
        self.rmax = self.ts if self.ts > self.tn else self.tn
        self.ptsfile = pointset
        self.gpdata = gpd.read_file(self.ptsfile)
        self.shp_has_geo = shp_has_geo
        self.xcol = XCOL
        self.ycol = YCOL
        self.zcol = ZCOL
        self.xpt = np.empty(0, dtype=np.double)
        self.ypt = np.empty(0, dtype=np.double)
        self.zpt = np.empty(0, dtype=np.double)
        self.numpts = 0
        self.ptlocator = vtk.vtkPointLocator()
        self.interpval = np.empty(0, dtype=np.double)
        self.defval = 0
        self.__initialize_pts()

    def getData(self):
        return self.xpt, self.ypt, self.zpt

    def MapwCLTemplate(self):
        nx, ny, nz = self.grid.GetDimensions()
        gval = vtk.vtkDoubleArray()
        gval.SetNumberOfComponents(1)
        gval.SetName(self.valname)
        coords = self.grid.GetPoints().GetData()
        # print coords
        np_coords = numpy_support.vtk_to_numpy(coords)
        # print np_coords
        xx = np_coords[:, 0]
        yy = np_coords[:, 1]
        print(xx.shape)
        xx = xx.reshape((nx, ny))
        yy = yy.reshape((nx, ny))
        ptsinradius = vtk.vtkIdList()
        self.interpval = np.zeros(shape=(nx,ny), dtype=np.double)
        for i in range(0, nx):
            for j in range(0,ny):
                found = False
                nrpts = 0
                searchRadius = 0.0
                rs = self.ts
                rn = self.tn
                val = 0.0
                rsum = 0.0
                # for radinc in range(0, 3):
                searchRadius += self.rmax
                self.ptlocator.FindPointsWithinRadius(searchRadius, (xx[i,j], yy[i,j], 0.0), ptsinradius)
                numpts = ptsinradius.GetNumberOfIds()
                phi = self.cl.getphiinterp(i)
                nv = 0
                for l in range(0,3):
                    for k in range(0,numpts):
                        xd = self.xpt[ptsinradius.GetId(k)] - xx[i,j]
                        yd = self.ypt[ptsinradius.GetId(k)] - yy[i,j]
                        ds = np.fabs(xd * np.cos(phi) + yd * np.sin(phi))
                        dn = np.fabs(yd * np.cos(phi) - xd * np.sin(phi))
                        if (ds < rs) and (dn < rn):
                            if ds == 0 and dn == 0:
                                self.interpval[i,j] = self.zpt[ptsinradius.GetId(k)]
                                found = True

                            nv += 1
                            r1 = np.power(np.sqrt(np.power(ds,2) + np.power(dn, 2)), self.wghtpowval)
                            rsum = rsum + 1/r1
                            val = val + self.zpt[ptsinradius.GetId(k)]/r1
                    if found: break
                    if nv == 0:
                        rs += self.ts
                        rn += self.tn
                    else:
                        self.interpval[i,j] = val/rsum
                        found = True
                        break
                if found == False:
                    self.interpval[i,j] = self.defval
                gval.InsertNextValue(self.interpval[i,j])
        self.grid.GetPointData().AddArray(gval)

    def plotmapgrid(self):
        vtknodes = self.grid.GetPoints().GetData()
        npnodes = numpy_support.vtk_to_numpy(vtknodes)
        x,y,z = npnodes[:,0], npnodes[:,1], npnodes[:,2]
        nx, ny, nz = self.grid.GetDimensions()
        x = np.reshape(x, (nx, ny))
        y = np.reshape(y, (nx, ny))
        val = self.grid.GetPointData().GetArray(self.valname)
        val = np.reshape(val, (nx, ny))
        plt.contourf(x,y,val)
        plt.show()
        # print(val)

    def __build_point_locator(self):
        points = vtk.vtkPoints()
        points.SetDataTypeToDouble()
        points.SetNumberOfPoints(self.numpts)
        for i in range(0, self.numpts):
            points.SetPoint(i, self.xpt[i], self.ypt[i], 0.0)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        points.ComputeBounds()

        self.ptlocator.SetDataSet(polydata)
        self.ptlocator.SetNumberOfPointsPerBucket(5)
        self.ptlocator.AutomaticOn()
        self.ptlocator.BuildLocator()

    def __initialize_pts(self):
        if self.shp_has_geo:
            tx = []
            ty = []
            tz = []
            for index, row in self.gpdata.iterrows():
                for pt in list(row['geometry'].coords):
                    tx.append(pt[0])
                    ty.append(pt[1])
                    tz.append(pt[2])
            self.xpt = np.array(tx)
            self.ypt = np.array(ty)
            self.zpt = np.array(tz)
            self.numpts = self.xpt.size
        else:
            # print(self.gpdata.columns)
            # xcol, ycol, zcol = self.__input_function()
            self.xpt = self.gpdata[self.xcol].values
            self.ypt = self.gpdata[self.ycol].values
            self.zpt = self. gpdata[self.zcol].values
            self.numpts = self.xpt.size
        self.defval = np.max(self.zpt)
        self.__build_point_locator()



    # def __initialize_pts(self):
    #     tx = []
    #     ty = []
    #     tz = []
    #     for index, row in self.gpdata.iterrows():
    #         if self.gpdata.geometry.geom_type[index] == 'Point':
    #             if self.gpdata.has_z[index] == True:
    #                 for pt in list(row['geometry'].coords):
    #                     # print(pt, type(pt))
    #                     tx.append(pt[0])
    #                     ty.append(pt[1])
    #                     tz.append(pt[2])
    #             else:
    #                 for pt in list(row['geometry'].coords):
    #                     tx.append(pt[0])
    #                     ty.append(pt[1])
    #                     tz.append(row['ELEVATION'])
    #         elif self.gpdata.geometry.geom_type[index] == None:
    #
    #                 tx.append(np.double(row.X))
    #                 ty.append(np.double(row.Y))
    #                 tz.append(np.double(row.ELEVATION))
    #
    #     self.xpt = np.array(tx)
    #     self.ypt = np.array(ty)
    #     self.zpt = np.array(tz)
    #     self.numpts = self.xpt.size
    # def __input_function(self):
    #     xcol = input("Enter column heading for X")
    #     ycol = input("Enter column heading for Y")
    #     zcol = input("Enter column heading for Z")
    #     return xcol, ycol, zcol
