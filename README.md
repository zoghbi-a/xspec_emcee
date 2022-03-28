# xspec_emcee

### emcee
[emcee](https://emcee.readthedocs.io/en/stable/) is a pure-Python implementation of Goodman & Weare’s [Affine Invariant Markov chain Monte Carlo (MCMC) Ensemble sampler](http://msp.berkeley.edu/camcos/2010/5-1/p04.xhtml).

### xspec
[xspec](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) is an X-Ray Spectral Fitting Package, distributed as part of the high energy astrophysics software package, [HEAsoft](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/) from NASA.

xspec has its own implementation of the GW algorithm, but I find it somewhat difficult to use, so I created my own implementation, which gives more control on the chains.

## INSTALL:
Choose the files relevant to your HEAsoft version by selecting the appropriate branch (master is the latest). Here I will work with v6.28.
- Download all updated 5 files: `Chain.cxx, Chain.h, ChainManager.cxx, ChainManager.h, xsChain.cxx`
- Place them in the right place inside the xspec source structure and recompile the relevant code:
  - `Chain.cxx, Chain.h, ChainManager.cxx, ChainManager.h` inside: `heasoft-6.28/Xspec/src/XSFit/MCMC`, and then inside `heasoft-6.28/Xspec/src/XSFit`, run: `hmake` and `hmake install`. Ensure that HEAsoft is initialized in the standard way before doing this. See their documentaton.
  - `xsChain.cxx` inside: `heasoft-6.30/Xspec/src/XSUser/Handler`, then inside `heasoft-6.30/Xspec/src/XSUser`, run: `hmake` and `hmake install`
- Run the GW chain in the usual way.


## Example:
Assuming the spectra and model have been setup and an initial fit is found, we do:
```tcl
chain len 10000
chain burn 10000
chain walker 100
para walk 30
chain run mcmc.fits
```
This will run the chain, printing progress along the way:
- The chains are initialized using 0.5*sigma from the fit covariance, so a valid fit is needed.
- The progress prints:
  - percentage progress:
  - best statistic in the current run.
  - acceptance fraction. It should be around ~0.2-0.3
  - The last number is the adjustable `a` parameter in the GW algorithm (see the [algorithm paper](https://arxiv.org/abs/1202.3665) for details). It can be adjusted to drive the acceptance fraction towards a desired value. If the acceptance fraction is too small, `a` can be reduced (using `chain temperature 1.5` for example) to increase the acceptance fraction.

```tcl
* Initializing: Using the 0.5* Covariance **

** Done initializaing **
         5%  498.871    0.313333       2
        10%  498.868       0.307       2
        15%  498.869       0.342       2
        20%  498.866       0.342       2
        25%  498.868       0.328       2
        30%  498.865       0.308       2
        35%  498.872       0.327       2
        40%  498.866       0.324       2
        45%  498.866       0.346       2
        50%  498.87       0.331       2
        55%  498.871       0.285       2
        60%  498.868       0.252       2
        65%  498.867       0.301       2
        70%  498.867        0.33       2
        75%  498.869       0.308       2
        80%  498.869       0.317       2
        85%  498.866       0.321       2
        90%  498.873       0.312       2
        95%  498.868       0.318       2
       100%  498.867       0.299       2
  New chain tmp.fits is now loaded.

```
