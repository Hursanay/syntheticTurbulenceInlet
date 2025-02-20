set t cairolatex pdf standalone size 7,7 linewidth 3.0 
set output "URL_inlet.tex"

# set terminal png size 800,800 enhanced font "Helvetica,12" #linewidth 2 
# set output 'URL_overTime.png'

# set grid
set pointsize 0.5
# set key samplen 1 maxrows 5 #ti "time [s]\n" #spacing 1.2 height 0.5 #width -1.5 # opaque
set multiplot layout 3,2 #title "Outlet area-averaged quantities\n"

file_syn = "syn_init.txt"

#order: 1.dis,2...,3.XX,4.XY,5.XZ,6.YY,7.YZ,8.ZZ,9.UXmean,10.Uymean,11.Uzmean
fileinlet = "../postProcessing/spanAv_inlet/85/patchCutLayerAverage.xy"


#########################################
set xrange [-0.1:2.1]
set xlabel "Position"
# set xtics 0.01

# set key at 0.8,7.0,0
set key center 
#########################################
#dist,U,L,Ruu,Ruv,Ruw,Rvv,Rvw,Rww

set ylabel "Ruu"
plot file_syn u 1:4 w lp pt 6 lc 'black' ti 'DNS',\
    fileinlet u 1:3 w l ti 'OF inlet'

set ylabel "Rvv"      
plot file_syn u 1:7 w lp pt 6 lc 'black' not,\
    fileinlet u 1:6 w l not

set ylabel "Rww"
plot file_syn u 1:9 w lp pt 6 lc 'black' not,\
    fileinlet u 1:8 w l not

set ylabel "Ruv"      
plot file_syn u 1:5 w lp pt 6 lc 'black' not,\
    fileinlet u 1:4 w l not

set ylabel "Ux mean profile"
plot file_syn u 1:2 w lp pt 6 lc 'black' not,\
    fileinlet u 1:9 w l not

set ylabel "L"
nu = 2.532e-3
plot file_syn u 1:3 w lp pt 6 lc 'black' not 


#########################################
# set multiplot out
# set out
