import sys
import numpy as np
from pyFM.Utils.COGridGen import COGridGen
from pyFM.Utils.mapper import mapper
from shapely.geometry import shape
import matplotlib.pyplot as plt
print("The excutable is {}".format(sys.executable))

grid = COGridGen(r"..\data\GreenRiver\centerline.shp", nx=25, ny=11, width=350)

print(grid.description())

x, y = grid.cl.getpoints()
print(x)
xo, yo = grid.cl.getinterppts(grid.nx)
print(xo)

xgrid, ygrid = grid.getXYGrid()
print(xgrid)

# plt.scatter(xgrid,ygrid)
# plt.plot(x,y)
# plt.plot(xo, yo)
# plt.show()

map = mapper(grid, 'Elevation', 60, 5, 1, r"..\data\GreenRiver\channeltopoz.shp", shp_has_geo=True)
map.MapwCLTemplate()
tmp = 1000
print(tmp)
map.plotmapgrid()
grid.plotmapgrid('Elevation')
# map2 = mapper(grid, 'Elevation', 60, 5, 1, r"..\data\GreenRiver\topo_csv.csv",XCOL='X',YCOL='Y',ZCOL='ELEVATION')
# def test1():
#     print("The excutable is {}".format(sys.executable))
#
#     cl = COGrid(r"..\data\GreenRiver\centerline.shp", nx=25, ny=10)
#
#     print(cl.description())
#
#     x,y = cl.getpoints()
#     print(x)
#     cl.getspline(25)
#     xo,yo = cl.getinterppts()
#     print(xo)
#
# if __name__ == '__main__':
#     test1()
