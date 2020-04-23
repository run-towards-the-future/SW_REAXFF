#lattice
#bsub -I -q q_sw_expr -debug -host_stack 256 -share_size 6144 -priv_size 4 -n 4 -b -cgsp 64  ./lmp_sunway  -in in.reaxc.lattice -sf sunway -var t 20 -var S 6 -var SX 1

