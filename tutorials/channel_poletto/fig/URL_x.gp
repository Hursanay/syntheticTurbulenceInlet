set t cairolatex pdf standalone size 7,7 linewidth 3.0 
set output "URL_x.tex"

# set terminal png size 800,800 enhanced font "Helvetica,12" #linewidth 2 
# set output 'URL_overTime.png'

# set grid
set pointsize 0.5
# set key samplen 1 maxrows 5 #ti "time [s]\n" #spacing 1.2 height 0.5 #width -1.5 # opaque
set multiplot layout 3,2 

file_syn = "syn_init.txt"

#order: 1.dis,2...,3.XX,4.XY,5.XZ,6.YY,7.YZ,8.ZZ,9.UXmean,10.Uymean,11.Uzmean
fileinlet = "../postProcessing/spanAv_inlet/85/patchCutLayerAverage.xy"

#order: 1.dis,2.XX,3.XY,4.XZ,5.YY,6.YZ,7.ZZ
filex1 = "../postProcessing/lineX1/85/line.xy"
filex2 = "../postProcessing/lineX2/85/line.xy"
filex3 = "../postProcessing/lineX3/85/line.xy"
filex4 = "../postProcessing/lineX4/85/line.xy"
filex5 = "../postProcessing/lineX5/85/line.xy"
filex6 = "../postProcessing/lineX6/85/line.xy"
filex7 = "../postProcessing/lineX7/85/line.xy"
filex8 = "../postProcessing/lineX8/85/line.xy"


#########################################
set xrange [-0.1:2.1]
# set xtics 0.01
set xlabel "Position"

# set key at 0.8,7.0,0
set key center top
#########################################

#Ruu,Ruv,Ruw,Rvv,Rvw,Rww
set ylabel "Ruu"
plot file_syn u 1:4 w lp pt 6 lc 'black' ti 'DNS',\
    filex1 u 1:2 w l ti 'x = 0.0',\
    filex2 u 1:2 w l ti 'x = 3.6',\
    filex3 u 1:2 w l ti 'x = 7.5',\
    filex4 u 1:2 w l ti 'x = 11.3',\
    filex5 u 1:2 w l ti 'x = 15.1',\
    filex6 u 1:2 w l ti 'x = 22.6',\
    filex7 u 1:2 w l ti 'x = 30.2',\
    filex8 u 1:2 w l ti 'x = 37.7'
    # fileinlet u 1:3 w l ti 'av. inlet',\

set ylabel "Rvv"      
plot file_syn u 1:7 w lp pt 6 lc 'black' not,\
     filex1 u 1:5 w l not,\
     filex2 u 1:5 w l not,\
     filex3 u 1:5 w l not,\
     filex4 u 1:5 w l not,\
     filex5 u 1:5 w l not,\
     filex6 u 1:5 w l not,\
     filex7 u 1:5 w l not,\
     filex8 u 1:5 w l not
    # fileinlet u 1:6 w l not,\

set ylabel "Rww"
plot file_syn u 1:9 w lp pt 6 lc 'black' not,\
     filex1 u 1:7 w l not,\
     filex2 u 1:7 w l not,\
     filex3 u 1:7 w l not,\
     filex4 u 1:7 w l not,\
     filex5 u 1:7 w l not,\
     filex6 u 1:7 w l not,\
     filex7 u 1:7 w l not,\
     filex8 u 1:7 w l not
    # fileinlet u 1:8 w l not,\

set ylabel "Ruv"      
plot file_syn u 1:5 w lp pt 6 lc 'black' not,\
     filex1 u 1:3 w l not,\
     filex2 u 1:3 w l not,\
     filex3 u 1:3 w l not,\
     filex4 u 1:3 w l not,\
     filex5 u 1:3 w l not,\
     filex6 u 1:3 w l not,\
     filex7 u 1:3 w l not,\
     filex8 u 1:3 w l not
    # fileinlet u 1:4 w l not,\

set ylabel "U mean profile"
plot file_syn u 1:2 w lp pt 6 lc 'black' not,\
     filex1 u 1:8 w l not,\
     filex2 u 1:8 w l not,\
     filex3 u 1:8 w l not,\
     filex4 u 1:8 w l not,\
     filex5 u 1:8 w l not,\
     filex6 u 1:8 w l not,\
     filex7 u 1:8 w l not,\
     filex8 u 1:8 w l not
    # fileinlet u 1:9 w l not

set ylabel "L"
nu = 2.532e-3
plot file_syn u 1:3 w lp pt 6 lc 'black' not #ti 'DNS'



#########################################
# set multiplot out
# set out
