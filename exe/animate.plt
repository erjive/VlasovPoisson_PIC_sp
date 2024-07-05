######################
###   Parameters   ###
######################

initial = 0
figs = 1600
time_output = 100
r_min = 0.0
r_max = 8.0
p_min =-0.4
p_max = 0.4

##############################
###   Create a directory   ###
##############################
path =  system('pwd')


directory = sprintf('plot_all2D%05.0f',figs)
create_directory = sprintf('mkdir %s', directory)

system(create_directory)

### Change to the directory

cd directory

set loadpath path
##################################
###   Make the contour table   ###
##################################

#set contour base
#set cntrparam levels 10
#set isosample 250,250
#unset surface
#outfile_contour2D = sprintf('contour2D%02.0f.dat',ts)
#set table outfile_contour2D
#splot "vlasov_fdist.2D" u 1:2:3 i ts
#unset table

#################################
###   Final terminal output   ###
#################################

#reset
set encoding locale #Permite poner acentos
#set title 'Evolución PIC Vlasov-Poisson'
#set title 'l=0, sin fondo'
set terminal png nocrop large size 1024,768

#set terminal epslatex
#outfile_plot = sprintf('contourplot%02.0f.tex',ts)
#set output outfile_plot

#Plot range

set xrange [r_min:r_max]
set yrange [p_min:p_max]
#set zrange [0:]
set cbrange [0:]
#Labels and key

set xlabel "r"
set ylabel "pr"
unset key

#Set palette

set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)
#set palette grey

### Plot size
set size square

### Margins

set lmargin 0 #left
set rmargin 0 #right
set bmargin 5 #bottom
set tmargin 7 #top

#set pm3d interpolate 0,0

#plot outfile_table2D with image, outfile_contour2D w l lt -1 lw 1.0
#set pm3d interpolate 0,0
set view map
set pm3d map
#set view 15,30
#unset surface

do for [i=initial:figs]{
 outfile = sprintf('animation%05.0f.png',i)
 set output outfile
 #splot 'vlasov_fdist.2D' u 1:2:3 index i with pm3d notitle
 steptime = sprintf('time step=%i',i*time_output)
 set title steptime
 plot "vlasov_fdist.2D" u 1:2:3 index i w points pt 5 ps 0.2 lc palette z 
 if (i%100 == 0) {print i}

 }

