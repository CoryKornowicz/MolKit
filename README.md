# MolKit

A Swift-version of a cheminformatics library, __shamelessly__ based off of OpenBabel to make adoption feel intuitive (got to start somewhere, right?).
To use openbabel from Swift, it would have needed C++ __<->__ Obj-C++ __<->__ Swift wrappers, so might as well write it in pure Swift then.  

~ 30,000 lines of code and still does not yet read molecules, but __almost__ there.

### Milestones:  
[X]  Tests for typers is completed and passing 
[X]  SmartsParser appears to be fully functional
[ ]  MKBuilder Tests
[ ]  SMILES Parser tests
[ ]  MOL2 Parser 
[ ]  PDB/PDBQT Parser

## WARNINGS AND CRITICAL HINDRANCES   

Currently in very very early development, and is not usable in productino. Once the file readers start working, then things will become interesting. Then comes the fun stuff like AutoDock Vina on iOS (but in Swift).

One problem looming on the horizon is the issue of Swift's copy-on-write paradigm, which could inflate the sizes of underlying datastructs in memory. I would like to experiment with using true pointers (UnsafeMutablePointer) as that would probably make handling coordinates more intuitive and behave more like the original OpenBabel implementation, however, the first goal is to make it work, then make it fast. 

This project lacks the ability to dynamically switch between Double and Float for the base precision. A global type that can switch between Float and Double needs to be created especially for devices with little RAM. 

There is also a lack of test cases, which are crucial for a production-ready package.  

There are no plans for CLI applications at the moment.

File Reader Class is not as expansive as OpenBabel's. Order of importance is SMI, PDB, and then SDF/MOL/MOL2. The goal is to also have trajectory parsing, but that would most likely need other additions as well.

Strong/Weak references are not considered yet, maybe lurking chances of ARC leaving dangling memory. 

If chains parsing is being difficult, rethink how to implement typeunion for ByteCode.

Huge up and coming debt bill on the the RefValue object. It needs to be refactored into a single protocol + wrapper class to supply either a UInt or RefValue 

## Installation
🤷‍♂️ spm?

### Swift Package Manager
Not production ready yet.




Thank you OpenBabel, tremendously, this project would have never been possible without your commitment to 
being open source. 
