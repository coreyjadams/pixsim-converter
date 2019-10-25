
This script is to calculate the reasonable projections of 3D into X, Y, Z.

The 3D coordinate range is:
X: 0 to 360     [cm] - 900  voxels
Y: -100 to 100  [cm] - 500  voxels
Z: 0 to 500     [cm] - 1250 voxels


The projection is at 0, +/- Theta degrees (3 projections total)
The X point (drift distance) is shared across all projections.

In practice, this means that an X,Y,Z point in 3D gets projected as follows:

Projection 0:
x_2d = x_3d
y_2d = cos(theta)*z_3d - sin(theta)*y_3d


Projection 1:
x_2d = x_3d
y_2d = cos(theta)*z_3d + sin(theta)*z_3d


Projection 2:
x_2d = x_3d
y_2d = z_3d


For the ranges of projections, this means that the projection min/max is (0.4 cm / pixel):

Projection 0:
x: (0, 360) [cm]   - 900  pixels
y: ( cos(theta)*z_min - sin(theta)*y_max to cos(theta)*z_max - sin(theta)*y_min)



Projection 1:

x: (0, 360) [cm]   - 900  pixels
y: ( cos(theta)*z_min + sin(theta)*y_min to cos(theta)*z_max + sin(theta)*y_max)

Projection 2:

x: (0, 360) [cm]   - 900  pixels
y: (0, 500) [cm]   - 1250 pixels
