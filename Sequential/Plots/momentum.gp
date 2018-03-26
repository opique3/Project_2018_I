reset
set xlabel 'time'
set ylabel 'Momentum'
set title 'Momentum VS. Time'
unset key

plot 	'data.out' 1:5 w l

pause -1
reset
