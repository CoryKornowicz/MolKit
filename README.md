# MolKit

A Swift-version of a cheminformatics library, __shamelessly__ based off of OpenBabel to make adoption feel intuitive. 
To use openbabel from Swift, it would have needed C++ __<->__ Obj-C++ __<->__ Swift wrappers, so might as well write it in pure Swift then.  

## WARNING  

Currently in very very early development, and not even close to being usable. Once the base molecule class is working, then things will become interesting. Shortly after will arrive the basic file parsers, and then the fun stuff.

One problem looming on the horizon is the issue of Swift's copy-on-write paradigm, which could inflate the sizes of underlying datastructs in memory. I would like to experiment with using true pointers (UnsafeMutablePointer) as that would probably make handling coordinates more intuitive and behave more like the original OpenBabel implementation.


## Installation

### Swift Package Manager
Not included as it should not be used publicly yet.


**Not meant to be used as a cmd line program at the moment** 
