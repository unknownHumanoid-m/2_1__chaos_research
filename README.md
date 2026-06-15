Here I will include a brief explination of the files I have included.

The folder labled YAMLs contians the .yml for the miniconda enviroments I used running celmech and heyoka simulations.
There was likley an update in the past few  months to one of celmech's dependancies that causes celmech's integrator to crash, so I use older version of some of the modules.
For the integrator I wrote it works better with the newest updates of everything.

The hk_integrator.py is a python file containing the integrator I wrote along with a couple functions wrote and use often. I used to copy and paste the cel with the functions
to every new jupyter file, but more recently I have been importing them from this file.

The rest of the files are jupyter notebooks files.

tree diagram.ipynb has my code for "tree diagrams" to look at resonance overlaps between different MMRs. Since we are currently only looking at the 2:1  MM and secular resonances, this code is somewhat
outdated, however I felt it might be useful to have.

Term isolation test particle and term isolation massive are where I took away terms from the disturbing function at different orders and looked at their e vs a, FFTs, and e vs t graphs.
The graphs generated from this code are included in the term_isolation_test_test_particle_and_massive_case slideshow.

The expanded grid notebook is when I took my origional chaotic initial conditions (note that these differ from the initital conditions I found that work for both rebound and celmech),
and varried eccentricity and semimajor axis of "Saturn". This was done using rebound considering the massive case.

The super_earths_grid follows the same process expect Jupiter and Saturn have been replaced with super earth sized planets.

pomega_frequency_using_chapter_seven was where I used the equations in chapter seven to analytically find secular resonances and compared them against rebound results.

kirkwood.ipynb and hk_kirkwood.ipynb contain the code where I used both n-body simulations to model the kirkwood gap chaos.

In 5 grid massive I have the code for the slide deck containing massive case time series and ffts. 

In the astroidsNearRessonance.ipynb, and the celmech_with_asteroids.ipynb file I have included code to simulate the motion of asteroids affected by the J-S chaos. For the first file,
I did this in rebound (NOTE: I neglected to divide the time by 2pi). In the second file I used celmech to simulation the motion of Jupiter and Saturn then used their location to 
add an "additional force" in rebound. 
