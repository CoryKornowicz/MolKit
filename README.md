# MolKit

A Swift-version of a cheminformatics library, __shamelessly__ based off of OpenBabel to make adoption feel intuitive. 
To use openbabel from Swift, it would have needed C++ __<->__ Obj-C++ __<->__ Swift wrappers, so might as well write it in pure Swift then.  

~ 30,000 lines of code and still does not yet read molecules, but __almost__ there.

### Milestones:  
-  Tests for typers is completed and passing 
-  SmartsParser appears to be fully functional

## WARNINGS AND CRITICAL HINDRANCES   

Currently in very very early development, and not even close to being usable. Once the base molecule class is working, then things will become interesting. Shortly after will arrive the basic file parsers, and then the fun stuff.

One problem looming on the horizon is the issue of Swift's copy-on-write paradigm, which could inflate the sizes of underlying datastructs in memory. I would like to experiment with using true pointers (UnsafeMutablePointer) as that would probably make handling coordinates more intuitive and behave more like the original OpenBabel implementation, however, the first goal is to make it work, then make it fast. 

This project lacks the ability to dynamically switch between Double and Float for the base precision. There are accompantying functions for Float types, but a global type that encompasses Float and Double needs to be inserted. 

There is also a lack of test cases, which are crucial for a production-ready package.  

Plans for CLI applications are scheduled after the package is working and useable for applications, then standalone binaries will be compiled.

File Reader Class is not as expansive as OpenBabel's. Order of importance is SMI, PDB, and then SDF/MOL/MOL2. The goal is to also have trajectory parsing, but that would most likely need other additions as well, and will likely need to wait until the conformer issue is addressed.  

Strong/Weak references (weak var/self) will need to be introduced to reduce the chances of ARC leaving dangling memory. 

If plugins are buggy, static variables may need to be swapped out for class variables to scope the instances better. 

If chains parsing is being difficult, rethink how to implement typeunion for ByteCode.

Huge up and coming debt bill on the the RefValue object. It needs to be refactored into a single protocol + wrapper class(avoid CoW) to supply either a UInt or RefValue 

## Installation

### Swift Package Manager
Not included as it should not be used publicly yet.




Thank you OpenBabel, this project would have never been possible without your commitment to 
being open source. 
