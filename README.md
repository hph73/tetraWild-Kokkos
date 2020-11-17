# tetraWild-Kokkos
Static libraries and header files in include are from fTetrawild (https://github.com/wildmeshing/fTetWild)
Libraries are compliled on Linux Centos. You may need to build fTetrawild on your machine to get those static libraries.

To build the project, simply run 
./run_installlibs_cmake.sh
make

To run the code, simply
./main

The code runs fine wiht g++ but fail with kokkos nvcc_wrapper

Here is the output with nvcc_wrapper:
================== loading mesh =====================
     ___________________________________________________________________________________________________________________________________________________________________________________________________________________________ 
    |                                                                                                                                                                                                                           |
    | o-[I/O         ] Loading file ./cube10.stl...                                                                                                                                                                             |
    |                  (FP64) nb_v:36 nb_e:0 nb_f:12 nb_b:36 tri:1 dim:3                                                                                                                                                        |
    |                  Attributes on vertices: point[3]                                                                                                                                                                         |
bbox_diag_length = 17.3205
ideal_edge_length = 0.866025
stage = 2
eps_input = 0.0173205
eps = 0.00958846
eps_simplification = 0.00767077
eps_coplanar = 1.73205e-05
dd = 0.011547
dd_simplification = 0.0092376
================== simplify =====================
collapsing 2.3226
swapping 0.453611
#boundary_e1 = 0
#boundary_e2 = 0
================== tetrahedralize =====================
EXIT_INV


Here is the output with g++:
================== loading mesh =====================
     ___________________________________________________________________________________________________________________________________________________________________________________________________________________________ 
    |                                                                                                                                                                                                                           |
    | o-[I/O         ] Loading file ./cube10.stl...                                                                                                                                                                             |
    |                  (FP64) nb_v:36 nb_e:0 nb_f:12 nb_b:36 tri:1 dim:3                                                                                                                                                        |
    |                  Attributes on vertices: point[3]                                                                                                                                                                         |
bbox_diag_length = 17.3205
ideal_edge_length = 0.866025
stage = 2
eps_input = 0.0173205
eps = 0.00958846
eps_simplification = 0.00767077
eps_coplanar = 1.73205e-05
dd = 0.011547
dd_simplification = 0.0092376
================== simplify =====================
collapsing 2.30477
swapping 0.45167
#boundary_e1 = 0
#boundary_e2 = 0
================== tetrahedralize =====================
================== insert_triangles =====================
#boundary_e1 = 0
#boundary_e2 = 0
known_surface_fs.size = 0
known_not_surface_fs.size = 0
================== optimization =====================
initializing...
edge collapsing...
fixed 0 tangled element
success(env) = 243
success = 592(2888)
success(env) = 25
success = 86(1169)
success(env) = 5
success = 20(416)
success(env) = 1
success = 4(99)
success(env) = 0
success = 2(32)
success(env) = 0
success = 2(17)
success(env) = 0
success = 2(12)
success(env) = 0
success = 1(4)
success(env) = 0
success = 0(4)
edge collapsing done!
time = 1.29745s
#v = 1742
#t = 8437
max_energy = 22.1874
avg_energy = 4.69214
//////////////// pass 0 ////////////////
edge splitting...
fixed 0 tangled element
success = 9522(9522)
edge splitting done!
time = 0.69949s
#v = 11264
#t = 57752
max_energy = 20.4709
avg_energy = 4.40273
edge collapsing...
fixed 0 tangled element
success(env) = 546
success = 6192(33856)
success(env) = 149
success = 720(9828)
success(env) = 18
success = 82(3780)
success(env) = 4
success = 18(959)
success(env) = 1
success = 4(186)
success(env) = 0
success = 1(31)
success(env) = 0
success = 0(8)
edge collapsing done!
time = 5.05729s
#v = 4247
#t = 20551
max_energy = 11.2935
avg_energy = 3.87913
edge swapping...
fixed 0 tangled element
success3 = 131
success4 = 831
success5 = 84
success = 1046(14113)
edge swapping done!
time = 1.28331s
#v = 4247
#t = 20504
max_energy = 10.833
avg_energy = 3.81963
vertex smoothing...
success = 2505(2973)
vertex smoothing done!
time = 5.70559s
#v = 4247
#t = 20504
max_energy = 6.50736
avg_energy = 3.63529
//////////////// postprocessing ////////////////
edge collapsing...
fixed 0 tangled element
success(env) = 72
success = 308(5973)
success(env) = 12
success = 34(2277)
success(env) = 1
success = 8(420)
success(env) = 0
success = 2(101)
success(env) = 1
success = 1(23)
success(env) = 0
success = 0(10)
edge collapsing done!
time = 0.870559s
#v = 3894
#t = 18679
max_energy = 5.80706
avg_energy = 3.60481
================== correct_tracked_surface_orientation =====================
================== filter_outside =====================