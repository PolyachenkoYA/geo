-----------------------------------------------------------------------------------------------------------
---------------------------------------------- .prm -------------------------------------------------------
-----------------------------------------------------------------------------------------------------------

                                            file format
---------------------------------------------------------------------------------------------------
comment_line
values               // Nrays, Tmax, dt, Amin, f, tau, eps, use_det, draw_mov, prnt_frames/raw, Xmin.x, Xmin.z, Xmin.z, Xmax.x, Xmax.z, Xmax.z, dX.x, dX.y, dX.z
---------------------------------------------------------------------------------------------------







                                              example
---------------------------------------------------------------------------------------------------
Nrays        Tmax       dt         Amin     f          tau      eps     use_det      draw_mov    prnt_frames/raw    Xmin       Xmax      dX
10000        6          0.01       0.01     0.002      5        0       0            1           1                  0 0 0      0 0 0     1 1 0
---------------------------------------------------------------------------------------------------
// 10000 rays from the start
// 601 timesteps in total
// frequency is 2000 Gz = 2 KGz
// tau is 0.000005 sec == 5 mkSec
// SYS_EPS (10^-10) used
// dont use detectors
// draw a movie of the evolution of the waves
// draw a movie by frames (not by raw data - don't draw each regictered ray separately)
// auto-find Xmin and Xmax
// use 2D grid with step==1 on XoY and with no cells for Z






                                               hints
---------------------------------------------------------------------------------------------------
1)
If prm.eps==0
system eps will be used. SYS_EPS = 1E-10 for now.
2)
if prm.dt==0
you can't draw movie and/or print raw rays for detectors
3) 
prm.f is assumed to be set in   10^6  Gz.
prm.tau is assumed to be set in 10^-6 sec.
4) 
if Xmax == Xmin
Xmax and Xmin will be found automatically
5) if dX.x/y/z == 0
dX.x/y/z = (Xmax.x/y/z - Xmin.x/y/z)*1.01
6)
(with draw_mov==1)
grid 300x400x400x1 uses ~5GB of RAM. 112 bytes per Ray.
7)
prnt_frames/raw = 1/0
---------------------------------------------------------------------------------------------------



















-----------------------------------------------------------------------------------------------------------
---------------------------------------------- .mat -------------------------------------------------------
-----------------------------------------------------------------------------------------------------------

                                            file format
---------------------------------------------------------------------------------------------------
N_of_materials
Cp_1 Cs_1
Cp_2 Cs_2
...
Cp_i Cs_i
...
---------------------------------------------------------------------------------------------------






                                              example
---------------------------------------------------------------------------------------------------
2
999999 999999
800 500
---------------------------------------------------------------------------------------------------
// 2 types of materials:
// background material designed to reflect everything
// usual material








                                               hints
---------------------------------------------------------------------------------------------------
1)
Must be at least 2 materials:
    1st is for background
    2nd is normal material
2)
If material of any side of any polygon is left unset after reading the whole .tet file, then its material is set to background - the very 1st material.
---------------------------------------------------------------------------------------------------





























-----------------------------------------------------------------------------------------------------------
---------------------------------------------- .tet -------------------------------------------------------
-----------------------------------------------------------------------------------------------------------

                                              file format
---------------------------------------------------------------------------------------------------
N_vertices
x_i y_i z_i // coords of the i-th vertex
...
N_tetrahedrons
v1_i v2_i v3_i v4_i mat_i // v1-4 - index of the 1-4th vertex of the i-th tetrahedron; mat_i - material of the i-th tetrahedron 
...
N_of_absorb_polygons N_of_detectors
v1_i v2_i v3_i det_i // v1-3 - index of the 1-3th vertex of the i-th polygon; det_i - index of the detector to assign i-th polygon to
---------------------------------------------------------------------------------------------------








                                                example
---------------------------------------------------------------------------------------------------
8
/* vertices section */
-100 -100 50 // x[0] y[0] z[0]
-100 100 50  // x[1] y[1] z[1]
100 100 50   // ...
100 -100 50
-100 -100 -50
-100 100 -50
100 100 -50
100 -100 -50
/* tetrahedrons section */
6
1 7 2 3 1   // vertices of the 1st tetrahedron - 1(-100 -100 50), 7(100 100 -50), 2(-100 100 50), 3(100 100 50); material of this tetrahedron - n.1 - Cp=800,Cs=500 (look an example in .mat section)
1 7 3 4 1
1 7 2 6 1
1 7 5 6 1
1 7 5 8 1
1 7 4 8 1
/* detectors section */
8 3         // 8 polygons-absorbers, 3 detectors
1 2 6 1     // polygons with vertices 1,2,6 is a part of the 1st detector and an absorber
1 5 6 1
2 3 7 2
2 6 7 2
4 3 7 3
4 8 7 3
1 4 8 0     // polygons with vertices 1,4,8 is an absorber, but not a detector
1 5 8 0
// polygons with vertices 1,2,3 wan't mentioned here, so it can reflect rays (it's not an absorber)
---------------------------------------------------------------------------------------------------










                                                  hints
---------------------------------------------------------------------------------------------------
1)
To make a polygon an absorber but not a detector set i_det=0 (4th number in the detector section)
---------------------------------------------------------------------------------------------------
