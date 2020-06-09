
###gnuplot commands:
unset key
#set size 0.26,0.5
set size 0.52, 1.0
#set pm3d map  corners2color c1
#set dgrid3d 254,333
set palette defined (1 'white', 2 'blue', 3 'red',4 'yellow' ,5 'pink')
unset xtics
unset ytics
unset colorbox
#set xrange [1:254]
#set yrange [1:333] 


set term pdfcairo
set output "SegSA.pdf"
plot "SegSA.dat"  u 1:(-$2):3 w image t "" 

set output "SegLocal.pdf"
plot "SegLocal.dat"  u 1:(-$2):3 w image t "" 

set term postscript eps enhanced color
set output "SegLocal.ps"
plot "SegLocal.dat"  u 1:(-$2):3 w image t "" 

set output "SegSA.ps"
plot "SegSA.dat"  u 1:(-$2):3 w image t "" 





 