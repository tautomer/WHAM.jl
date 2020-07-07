# A Julia implemetation of methods to re-construct the free energy profile from 1D umbrella sampling

[Umbrella sampling](https://en.wikipedia.org/wiki/Umbrella_sampling) is one of
the common techniques to overcome the [rare events issue](https://en.wikipedia.org/wiki/Rare_events)
when computing the free energy prfile which is also called the potential of
mean force (PMF).

After umbrella sampling, it is necessary to unbias the biased distribution
obtained from the retrained dynamics. There are two popular ways of doing this,

* [the weighted histogram analysis method](http://membrane.urmc.rochester.edu/sites/default/files/wham/wham_talk.pdf),
* and [umbrella integration](http://www.theochem.uni-stuttgart.de/kaestner/umbrella_integration.html).

Both techniques are supported in this package.

## TODO's

* Expand to multi-dimension, at least 2D.
* Minimize allocations.
* Parallelize WHAM.
* More detailed documentaions.
* Make the package also work as a script.