# QMMM benchmark scripts


## Extract 

**Example usage:**
```
python3 qmmmextract.py extract --function=mp_waitany --notop=5
```
This extracts timing data from CP2K output (.log) files within directories
in order to generate scaling data. The subroutine name can be given as 
arguement (--function) but if this is not supplied the total
run time is given.


For outputs with the same CP2K project name with the same number of
threads, MD steps and processes timing values are averaged. This is
written into a file with the name format:

```
$PROJECT_$function_nsteps_nthreads.out
```


This also writes the top n=notop most costly subroutines and their self
times into a file so they can be plotted as a stacked bar chart. The
files for this will be named as follows:

```
$PROJECT_topcalls_nsteps_nthreads.out
```

To get the normalised timestep difference times for each subrountine
 add the minsteps and maxsteps options, where these are the two 
timesteps to difference, e.g.

```
python3 qmmmextract.py extract --function=total --notop=10 --minsteps=1 --maxsteps=6
```

This will write the normalise times in a file named as follows:

```
$PROJECT_topcalls_diff_maxsteps-minsteps_nthreads.out
```
and also write the total time (6-1)/5 differnce in:

```
$PROJECT_total_diff_maxsteps-minsteps_nthreads.out
```

## Plot

**Example usage:**
```
python3 qmmmextract.py plot --filename=MQAE_topcalls_1steps_1threads.out --plottype=bar
```

This plots a data file (e.g. from generate by extract) given by its filename. The
plot type can be _time_ or _su_ (speed up) or _bar_ (stacked bar for the top calls)

## Compare


