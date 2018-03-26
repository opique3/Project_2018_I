reset
set xlabel 'time'
set ylabel 'Energy'
set title 'Different Energies VS. Time'

plot 	'data.out' 1:2 w l title 'Potential Energy', \
	'data.out' 1:3 w l title 'Kinetic Energy', \
	'data.out' 1:4 w l title 'Mechanical Energy'

pause -1
reset
