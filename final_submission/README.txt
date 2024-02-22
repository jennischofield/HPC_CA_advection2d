To compile the advection program serially, first run: 

gcc -o advection2D -std=c99 advection2D.c -lm

To compile the advection program parallelised, first run: 

gcc -fopenmp -o advection2D -std=c99 advection2D.c -lm

This assumes you are on one of the Lovelace lab computers.

To run the program, use: 

./advection2D

This will generate 3 files - initial.dat, final.dat, and average.dat

To generate the graph for final.dat:

gnuplot plot_final

To generate the graph for initial.dat:

gnuplot plot_initial

To generate the graph for average.dat:

gnuplot plot_average

If testing using a Windows machine, add the following pragma below the #include statements:

#pragma check_stack(off)

This is not necessary for the Lovelace machines.

This folder contains:

advection2D.c - source code for the advection function, with parameters set to the values in 2.3

final_values_2.2 - image displaying the final results for 2.2

final_values_2.3 - image displaying the final results for 2.3

initial_conditons - image displaying starting conditions

plot_average - script to plot the graph of x versus vertically averaged u(x,y) from average.dat, will output image "average.png" 

plot_final - script to plot the final results from final.dat, will output image "final.png"

plot_initial - script to plot the initial results from initial.dat, will output image "initial.png"

README - file containing instructions on how to run code

vertically_averaged - image displaying the graph of x versus vertically averaged u(x,y), for part 2.4
