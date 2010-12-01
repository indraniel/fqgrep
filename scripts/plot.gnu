# specify encoding and output format
#set terminal svg enhanced size 2048 500
#set terminal svg enhanced size 800 1200
set terminal svg enhanced size 1024 2048
set output 'test.fq.svg'
#set size size ratio -1

# setup the datafile format
set datafile separator "\t"

# setup a 3x2 plot layout
set multiplot layout 5, 1

# 1st graph (CDF) #############################################################

# plot labels
set title "Cumulative Filtered Read Distribution"
set ylabel "% of Filtered Reads"
set xlabel "mismatch level"

# axis format
set xtics 1, 1

# enable grid lines
set grid back lt 20 lw 0.5 lc rgbcolor "#dddddd"
show grid

# draw the plot and specify its appearance
plot 'test.fq.cnt' using 1:5 notitle with points 2 6, \
     'test.fq.cnt' using 1:5 notitle with lines

# 2nd graph (PDF) #############################################################

# plot labels
set title "Filtered Read Count at Mismatch Level"
set ylabel "# of Filtered Reads"
set xlabel "mismatch level"

# axis format
set xtics 1, 1

# enable grid lines
set grid back lt 20 lw 0.5 lc rgbcolor "#dddddd"
show grid

# draw the plot and specify its appearance
plot 'test.fq.cnt' using 1:3 notitle with points 3 6, \
     'test.fq.cnt' using 1:3 notitle with lines


# 3rd graph (Read Length Distribution : Mismatch 0) ###########################

# plot labels
set title "Read Length Distribution at Mismatch = 0"
set ylabel "# of Filtered Reads"
set xlabel "Read Length"

# axis format
set xtics 1, 1

# enable grid lines
set grid back lt 20 lw 0.5 lc rgbcolor "#dddddd"
show grid

plot 'test.fq.rlh.dat' index 0 using 1:2 notitle with points 2 6, \
     'test.fq.rlh.dat' index 0 using 1:2 notitle with lines

# 4th graph (Read Length Distribution : Mismatch 1) ###########################

# plot labels
set title "Read Length Distribution at Mismatch = 1"
set ylabel "# of Filtered Reads"
set xlabel "Read Length"

# axis format
set xtics 1, 1

# enable grid lines
set grid back lt 20 lw 0.5 lc rgbcolor "#dddddd"
show grid

plot 'test.fq.rlh.dat' index 1 using 1:2 notitle with points 2 6, \
     'test.fq.rlh.dat' index 1 using 1:2 notitle with lines

# 5th graph (Read Length Distribution : Mismatch 2) ###########################

# plot labels
set title "Read Length Distribution at Mismatch = 2"
set ylabel "# of Filtered Reads"
set xlabel "Read Length"

# axis format
set xtics 1, 1

# enable grid lines
set grid back lt 20 lw 0.5 lc rgbcolor "#dddddd"
show grid

plot 'test.fq.rlh.dat' index 2 using 1:2 notitle with points 2 6, \
     'test.fq.rlh.dat' index 2 using 1:2 notitle with lines


# finish up the overall plot
unset multiplot

# Close the file (so I don't have to close gnuplot to view it)
set output
