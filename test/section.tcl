
# W14x159 
#  A = 46.7 in^2
# Ix = 1900 in^4
# Iy =  748 in^4

set d 15.0
set tw 0.745
set bf 15.6
set tf 1.19

set nfx 40
set nfy 40

set Fy 50.0
set E  29000.0
set Eh [expr $E/1000.0]
set H  [expr $E*$Eh/($E-$Eh)]

uniaxialMaterial Hardening 1 $E $Fy $H 0.0

section Fiber 1 {
    # Flanges
    set nX $nfx
    set nY [expr int(ceil($nfy*$tf/$d))]
    patch rect 1 $nY $nX [expr  0.5*$d-$tf] [expr -0.5*$bf] [expr  0.5*$d]     [expr 0.5*$bf]
    patch rect 1 $nY $nX [expr -0.5*$d]     [expr -0.5*$bf] [expr -0.5*$d+$tf] [expr 0.5*$bf]

    # Web
    set nX [expr int(ceil($nfx*$tw/$bf))]
    set nY [expr int(ceil($nfx*($d-2*$tf)/$bf))]
    patch rect 1 $nY $nX [expr -0.5*$d+$tf] [expr -0.5*$tw] [expr 0.5*$d-$tf] [expr 0.5*$tw]
}