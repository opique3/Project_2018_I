reset
set xlabel 'time'
set ylabel 'Temperature'
set title 'Temperature VS. Time'
unset key

plot 	'data.out' 1:6 w l 

pause -1
reset
