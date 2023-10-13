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

Currently in very very early development, and is not usable in production. Once the file readers start working, then things will become interesting. Then comes the fun stuff like AutoDock Vina on iOS but in pure Swift.

One problem looming on the horizon is the issue of Swift's copy-on-write paradigm, which could inflate the sizes of underlying datastructs in memory. I would like to experiment with using true pointers (UnsafeMutablePointer) as that would probably make handling coordinates more intuitive and behave more like the original OpenBabel implementation, however, the first goal is to make it work, then make it fast. Making it fast also includes injecting more Swift-esque paradigms. 

This project lacks the ability to dynamically switch between Double and Float for the base precision. A global type that can switch between Float and Double needs to be created especially for devices with little RAM. 

There is also a lack of test cases, which are crucial for a production-ready package.  

There are no plans for CLI applications at the moment.

File Reader Class is not as expansive as OpenBabel's. Order of importance is SMI, PDB, and then SDF/MOL/MOL2. The goal is to also have trajectory parsing, but that would most likely need other additions as well.

Strong/Weak references are not considered yet, maybe lurking chances of ARC leaving dangling memory. 

If chains parsing is being difficult, rethink how to implement typeunion for ByteCode.

Huge up and coming debt bill on the the RefValue object. It needs to be refactored into a single protocol + wrapper class to supply either a UInt or RefValue 

## Installation
### Swift Package Manager
To use [Swift Package Manager](https://swift.org/package-manager/) add MolKit to your Package.swift file: 
```swift
let package = Package(
    name: "your-project",
    dependencies: [
        .package(url: "https://github.com/CoryKornowicz/MolKit.git", .branchItem("main")),
    ],
    targets: [
        .target(
            name: "your-project",
            dependencies: ["MolKit"]),
    ]
)
```


Thank you OpenBabel, tremendously, this project would have never been possible without your commitment to 
being open source. 

## Final Projects ## 
- [ ] Minimize a ligand molecule in a binding pocket.

### Final Project Details ### 
#### Protein Ligand Minimization ####
Essentially, it will look something like this is Swift (adopted from openbabel/src/forcefield.cpp): 
```swift
var mol: MKMol = MKMol()
//
// Read the pocket + ligand (initial guess for position) into mol...
//
var pocket: Bitset // set the bits with atoms indexes for the pocket to 1...
var ligand: Bitset // set the bits with atoms indexes for the ligand to 1...

guard let pFF: MKForceField = MKForceField.findForceField("MMFF94") else { 
    // handle error...
}

// set logging 
pFF.setLogFile(...)
pFF.setLogLevel(.low)

// Fix the binding pocket atoms
var constraints: MKConstraints = MKConstraints()
for a in mol.getAtomIterator() { 
    if pocket.contains(a.getIdx()) { 
        constraints.addAtomConstraint(a.getIdx())
    }
}

// Specify the interacting groups. The pocket atoms are fixed, so there
// is no need to calculate intra- and inter-molecular interactions for
// the binding pocket.

pFF.addIntraGroup(ligand) // bonded interactions in the ligand
pFF.addInterGroup(ligand) // non-bonded between ligand-ligand atoms
pFF.addIntrGroups(ligand, pocket) // non-bonded between ligand and pocket atoms

// We pass the constraints as argument for Setup()
do {
    try pFF.setup(mol, constraints)
    // Perform the actual minimization, maximum 1000 steps
    pFF.conjugateGradients(1000)
} catch { 
    // handle error
}

```
