linsker
=======

A self-adaptive neural network simulator to explore emergent spatio-opponency in cortical cells after application of random noise.

The best way to plot the output data is, in gnuplot:
    set palette model RGB defined ( -1 'red', 1 'green' ); set zrange [-1.0:1.0];
    plot "data.c.99" with points palette

TODO:
- Clean up Network::applyNoise().
- CLean up code in general.
- Add argument parsing routines.
- Write better data output routines.
