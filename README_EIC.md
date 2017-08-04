podio for the EIC
=================

The `lcio` contains a version of the LCIO data-model using `podio` called 
`lcio2`.  First you need to build this version of `podio` then build the new 
data-model.

## Installation


### Step 1: Build `podio`

```
git clone git@eicweb.phy.anl.gov:EIC/podio.git
cd podio
mkdir build && cd build
cmake ../. -DCMAKE_INSTALL_PREFIX=/usr/local
make -j4 install
```

### Step 2: Build `lcio2`

Here we assume that the prefix (`/usr/local` in the above case) is already in 
your `PATH` and `LD_LIBRARY_PATH`.

```
cd ../lcio
mkdir build && cd build
cmake ../. -DCMAKE_INSTALL_PREFIX=/usr/local
make -j4 install
```



