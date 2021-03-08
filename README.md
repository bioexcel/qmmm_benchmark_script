# QMMM benchmark scripts


## extract
```
extract(function)
```
extracts timing data from CP2K output (.log) files within directories 
subroutine name can be given as arguement, empty gives total run time
writes out run times for same thread and step number values
also gets the top 10 subrountines from each run

## plot data

```
plotdata(filename, plottype)
```

plot a data file given by filename
plot type can be time or su or bar

## Get times based on time per step (e.g removing the contribution from the first step)

```
gettimeperstep("MQAE","total",6,1,1)
```
