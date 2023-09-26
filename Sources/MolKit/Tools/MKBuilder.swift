

import Foundation
import Surge
import simd

/** \class OBBuilder builder.h <openbabel/builder.h>
    \brief Class for 3D structure generation
    \since version 2.2

    The OBBuilder class is used for generating 3D structures.

    Below is and example which explain the basics.

    \code
    //
    // code to read molecule from smiles goes here...
    //
    OBBuilder builder;
    builder.Build(mol);
    //
    // code to write molecule to 3D file format goes here...
    //
    \endcode
**/
class MKBuilder {

    var _keeprings: Bool    
    static var _rigid_fragments: [String] = []
    static var _ring_fragments: [Pair<MKSmartsPattern, [Vector<Double>]>] = []
    static var _rigid_fragments_index: [String: Int64] = [:]
    static var _rigid_fragments_cache: [String: [Vector<Double>]] = [:]

    init(_ keepring: Bool) {
        self._keeprings = keepring
    }

    ///@name Setup build parameters
      //@{
      /*! If the molecule already contains 3D coordinates, if you set KeepRings to true it will use
       *  retain the 3D coordinates of the rings. By default KeepRings is false, and ring conformations
       *  are obtained by lookup in a library of ring conformers. However, since the ring conformer library
       *  is not exhaustive, if the ring system is not found in the library, the resulting 3D structure can
       *  be poor, and require geometry optimisation before it is reasonable. If your starting point is
       *  a 3D structure, you can set KeepRings to true, and the conformation will be taken from the input.
       *  The remaining (acyclic) bonds will still all be built by the builder.
       */
    func setKeepRings() {
        self._keeprings = true
    }

    func unsetKeepRings() {
        self._keeprings = false
    }


    //! Used by LoadFragments to check for invalid (all zero coordinates) fragments
    func addRingFragment(_ sp: MKSmartsPattern, _ coords: [Vector<Double>]) {
        var hasAllZeroCoords: Bool = true
        for coord in coords {
            if fabs(coord.x) > 10e-8 || fabs(coord.y) > 10e-8 || fabs(coord.z) > 10e-8 {
                hasAllZeroCoords = false
                break
            }
        }
        if hasAllZeroCoords {
            print("Ring fragment \(sp.getSMARTS()) in ring-fragments.txt has all zero coordinates. Ignoring fragment.")
        } else {
            MKBuilder._ring_fragments.append((sp, coords))
        }
    }

    //! Load fragment info from file, if is it has not already been done
    func loadFragments() {
        // open data/fragments.txt
        guard let filePath = Bundle.module.url(forResource: "rigid-fragments-index", withExtension: "txt", subdirectory: "Data") else {
            fatalError("Could not locate rigid-fragments-index.txt")
        }
        
        filePath.foreachRow { rowContent, lineNum in
            var contents = rowContent.components(separatedBy: .whitespaces)
            contents.removeAll { $0.isEmpty }
            MKBuilder._rigid_fragments.append(contents[0])
            MKBuilder._rigid_fragments_index[contents[0]] = Int64(contents[1].toInt()!)
        }

        guard let filePath = Bundle.module.url(forResource: "ring-fragments", withExtension: "txt", subdirectory: "Data") else {
            fatalError("Could not locate ring-fragments.txt")
        }

        var sp: MKSmartsPattern? = nil
        var coords: [Vector<Double>] = []
        filePath.foreachRow { rowContent, lineNum in
            if !rowContent.starts(with: "#") {
                var contents = rowContent.components(separatedBy: .whitespaces)
                contents.removeAll { $0.isEmpty }
                
                if contents.count == 1 { // SMARTS Pattern
                    if sp != nil {
                        addRingFragment(sp!, coords)
                    }
                    coords = [] 
                    sp = MKSmartsPattern()
                    if !sp!.initialize(contents[0]) {
                        sp = nil 
                        print("Could not parse SMARTS from contribution data file: \(contents[0])")
                    }
                } else if contents.count == 3 { // XYZ Coordinates 
                    let coord = Vector<Double>.init(dimensions: 3, { val in
                        contents[val].toDouble()!
                    })
                    coords.append(coord)
                }
            }
        }
        addRingFragment(sp!, coords)
        
        // return the locale to the original one
        //torsions.txt used to be referenced here too

    }

    func getFragmentCoord(_ smiles: String) -> [Vector<Double>] {
        if MKBuilder._rigid_fragments_cache.filter({$0.key == smiles}).count > 0 {
            return MKBuilder._rigid_fragments_cache[smiles]!
        }
        var coords: [Vector<Double>] = []
        
        if MKBuilder._rigid_fragments_index.filter({$0.key == smiles}).count == 0 {
            return coords
        }
        
        guard let filePath = Bundle.module.url(forResource: "rigid-fragments", withExtension: "txt", subdirectory: "Data") else {
            fatalError("Could not locate rigid-fragments-index.txt")
        }
        
        var hasAllZeroCoords: Bool = true
        
        filePath.foreachRow (offset: MKBuilder._rigid_fragments_index[smiles], { rowContent, lineNum in
            var vs = rowContent.components(separatedBy: .whitespaces)
            vs.removeAll { $0.isEmpty }
            if vs.count == 4 {
                let coord = Vector<Double>.init(dimensions: 3, { val in
                    vs[val+1].toDouble()!
                })
                assert(coord.count == 3, "Coordinates should have 3 dimensions")
                if fabs(coord.x) > 10e-8 || fabs(coord.y) > 10e-8 || fabs(coord.z) > 10e-8 {
                    hasAllZeroCoords = false
                }
                coords.append(coord)
            } else if vs.count == 1 { // SMARTS pattern
                // stop reading and prepare to return
                if hasAllZeroCoords {
                    fatalError("Rigid fragment \(smiles) in rigid-fragments.txt has all zero coordinates.")
                }
                return
            }
        })
        return coords
    }

    /*! Get the position for a new neighbour on atom.  Returns
       * non-finite vector if there is no reasonable location.
       *  \param atom Atom for which we want a new neighbour location.
       *  \returns The position for the new atom.
       */
    static func getNewBondVector(_ atom: MKAtom) -> Vector<Double> {
        getNewBondVector(atom, 1.5)
    }

    static func getNewBondVector(_ atom: MKAtom, _ length: Double) -> Vector<Double> {
        var bond1: Vector<Double> = VZero
        var bond2: Vector<Double> = VZero
        var bond3: Vector<Double> = VZero
        var bond4: Vector<Double> = VZero
        var bond5: Vector<Double> = VZero
        var v1: Vector<Double> = VZero
        var v2: Vector<Double> = VZero
        var newbond: Vector<Double> = VZero


        let dimension = atom.getParent()!.getDimension()

        if dimension != 2 {
            /////////////////
            //   3D or 0D  //
            /////////////////

            //
            //  a   --->   a--*
            //
            if atom.getExplicitDegree() == 0 {
                let newbond = atom.getVector() + VX * length 
                return newbond
            }
            // hyb * = 1
            // ^^^^^^^^^
            //
            //   (a-1)--a   --->   (a-1)--a--*        angle(a-1, a, *) = 180
            //
            // hyb * = 2
            // ^^^^^^^^^
            // make sure we place the new atom trans to a-2 (if there is an a-2 atom)
            //
            //   (a-2)             (a-2)
            //     \                 \                                             //
            //    (a-1)==a   --->   (a-1)==a          angle(a-1, a, *) = 120
            //                              \                                      //
            //                               *
            // hyb * = 3
            // ^^^^^^^^^
            // make sure we place the new atom trans to a-2 (if there is an a-2 atom)
            //
            //   (a-2)             (a-2)
            //     \                 \                                             //
            //    (a-1)--a   --->   (a-1)--a          angle(a-1, a, *) = 109
            //                              \                                      //
            //                               *
            if atom.getExplicitDegree() == 1 {
                let isCarboxylate0 = atom.isCarboxylOxygen()

                for nbr in atom.getNbrAtomIterator()! {
                    if nbr.getHyb() == 1 { continue }
                    bond1 = atom.getVector() - nbr.getVector()
                    for nbr2 in nbr.getNbrAtomIterator()! {
                        if nbr2 != atom {
                            bond2 = nbr.getVector() - nbr2.getVector()
                            if isCarboxylate0 && nbr2.getAtomicNum() == MKElements.Oxygen.atomicNum {
                                break // make sure that the hydrogen is trans to the C=O
                            }
                        }
                    }
                }

                bond1.normalize()
                var v1 = cross3x3(bond1, bond2)
                if (bond2 == VZero || v1 == VZero) {
                    var vrand: Vector<Double> = VZero
                    vrand.randomUnitVector()
                    var angle = fabs(acos(dot(bond1, vrand)).radiansToDegrees)
                    while angle < 45.0 || angle > 135.0 {
                        vrand.randomUnitVector()
                        angle = fabs(acos(dot(bond1, vrand)).radiansToDegrees)
                    }
                    // there is no a-2 atom 
                    v1 = cross3x3(bond1, vrand)
                }
                v2 = cross3x3(bond1, v1)
                v2.normalize()
                // check to see if atom is a square planar in disguise
                if atom.getHyb() == 3 {
                    let stereoFacade: MKStereoFacade = MKStereoFacade(atom.getParent()!, m_perceive: false) // don't perceive stereo
                    if stereoFacade.hasSquarePlanarStereo(atom.getId().rawValue) {
                        atom.setHyb(4) // force sq. planar geometry for sq. planar stereo 
                    }
                }

                if atom.getHyb() == 1 {
                    newbond = bond1 // i.e., in the exact opposite direction
                } else if atom.getHyb() == 2 {
                    newbond = bond1 - v2 * tan(60.0.degreesToRadians)
                } else if atom.getHyb() == 3 {
                    newbond = bond1 - v2 * tan((180.0 - 109.471).degreesToRadians)
                } else if atom.getHyb() == 4 {
                    newbond = bond1 // like 5-coordinate below, we want a 180-degree bond (trans)
                } else if atom.getHyb() == 5 {
                    /* the first two atoms are the axial ones;  the third, fourth, and fifth atom are equatorial */
                    newbond = bond1
                } else if atom.getHyb() == 6 {
                    newbond = bond1 - v2 * tan(90.0.degreesToRadians)
                }

                newbond.normalize()
                newbond *= length
                newbond += atom.getVector()
                return newbond
            }
            //
            //    \	      \                                                     //
            //     X  --->   X--*
            //    /         /
            //
            if atom.getExplicitDegree() == 2 {
                for nbr in atom.getNbrAtomIterator()! {
                    if bond1 == VZero {
                        bond1 = atom.getVector() - nbr.getVector()
                    } else {
                        bond2 = atom.getVector() - nbr.getVector()
                    }
                }
                bond1.normalize()
                bond2.normalize()
                v1 = bond1 + bond2
                v1.normalize()

                if atom.getHyb() == 2 {
                    newbond = v1
                } else if atom.getHyb() == 3 {
                    v2 = cross3x3(bond1, bond2) // find the perpendicular 
                    v2.normalize()
                    // newbond = bond1 - v2 * tan((180.0 - 109.471).degreesToRadians) // old code 
                    newbond = v2 + v1 * (sqrt(2.0) / 2.0) // used to be tan(70.53 degrees/2) which is sqrt(2.0) / 2.0
                } else if atom.getHyb() == 4 || atom.getHyb() == 5 {
                    /* add the first equatorial atom, orthogonally to bond1 (and bond2 = -bond1) */
                    /* is atom order correct?  I don't think it matters, but I might have to ask a chemist
                    * whether PClF4 would be more likely to have an equatorial or axial Cl-P bond */
                    var vrand: Vector<Double> = VZero
                    vrand.randomUnitVector()
                    var angle = fabs(acos(dot(bond1, vrand)).radiansToDegrees)
                    while angle < 45.0 || angle > 135.0 {
                        vrand.randomUnitVector()
                        angle = fabs(acos(dot(bond1, vrand)).radiansToDegrees)
                    }
                    v1 = cross3x3(bond1, vrand)
                    v1.normalize()
                    newbond = v1
                } else if atom.getHyb() == 6 {
                    v2 = cross3x3(bond1, bond2)
                    newbond = v2
                }
                newbond.normalize() //if newbond was not set, it will become non-finite here
                newbond *= length
                newbond += atom.getVector()
                return newbond
            }

            /* UFF:
            *    b lg dg  o  y
            *  b - 45 30 45 30
            * lg    - 45  0 45
            * dg       - 45 30
            *  o          - 45
            *  y             -

            * 94s:
            *    b lg dg  o  y
            *  b - 34 34 34 34
            * lg    - 48 21 48
            * dg       - 48 21
            *  o          - 48
            *  y             -

            //
            //    \	       \
            //   --X  --->  --X--*
            //    /          /
            //
            */
            if atom.getExplicitDegree() == 3 {
                for nbr in atom.getNbrAtomIterator()! {
                    if bond1 == VZero {
                        bond1 = atom.getVector() - nbr.getVector()
                    } else if bond2 == VZero {
                        bond2 = atom.getVector() - nbr.getVector()
                    } else if bond3 == VZero {
                        bond3 = atom.getVector() - nbr.getVector()
                    }
                }
                if atom.getHyb() == 3 {
                    bond1.normalize()
                    bond2.normalize()
                    bond3.normalize()
                    newbond = bond1 + bond2 + bond3
                    newbond.normalize()
                    newbond *= length
                    newbond += atom.getVector()
                    return newbond
                } else if atom.getHyb() == 4 { // OK, we want this at -bond3, since bond1 & bond2 are opposite
                    bond3.normalize() 
                    newbond = bond3 
                    newbond *= length
                    newbond += atom.getVector()
                    return newbond
                } else if atom.getHyb() == 5 {
                    bond1.normalize()
                    bond2.normalize()
                    bond3.normalize()
                    v1 = cross3x3(bond1, bond3) 
                    v1.normalize()
                    newbond = v1 + tan(30.0.degreesToRadians) * bond3
                    newbond.normalize()
                    newbond *= length
                    newbond += atom.getVector()
                    return newbond
                } else if atom.getHyb() == 6 {
                    /* the way things work, newbond is equal to bond1, but will show up at -bond1 next time around */
                    newbond = bond1 
                    newbond.normalize()
                    newbond *= length
                    newbond += atom.getVector()
                    return newbond
                }
            }
            
            if atom.getExplicitDegree() == 4 {
                for nbr in atom.getNbrAtomIterator()! {
                    if bond1 == VZero {
                        bond1 = atom.getVector() - nbr.getVector()
                    } else if bond2 == VZero {
                        bond2 = atom.getVector() - nbr.getVector()
                    } else if bond3 == VZero {
                        bond3 = atom.getVector() - nbr.getVector()
                    } else if bond4 == VZero {
                        bond4 = atom.getVector() - nbr.getVector()
                    }
                }

                if atom.getHyb() == 5 {
                    bond1.normalize()
                    bond2.normalize()
                    bond3.normalize()
                    bond4.normalize()
                    v1 = cross3x3(bond1, bond3)
                    v1.normalize()
                    newbond = (-1 * v1) + tan(30.0.degreesToRadians) * bond3
                    newbond.normalize()
                    newbond *= length
                    newbond += atom.getVector()
                    return newbond
                } 
                if atom.getHyb() == 6 {
                    newbond = bond2 
                    newbond.normalize()
                    newbond *= length
                    newbond += atom.getVector()
                    return newbond
                }
            }

            if atom.getExplicitDegree() == 5 {
                if atom.getHyb() == 6 {
                    for nbr in atom.getNbrAtomIterator()! {
                        if bond1 == VZero {
                            bond1 = atom.getVector() - nbr.getVector()
                        } else if bond2 == VZero {
                            bond2 = atom.getVector() - nbr.getVector()
                        } else if bond3 == VZero {
                            bond3 = atom.getVector() - nbr.getVector()
                        } else if bond4 == VZero {
                            bond4 = atom.getVector() - nbr.getVector()
                        } else if bond5 == VZero {
                            bond5 = atom.getVector() - nbr.getVector()
                        }
                    }
                    newbond = bond3 
                    newbond.normalize()
                    newbond *= length
                    newbond += atom.getVector()
                    return newbond
                }
            }

            // Undefined case -- return a random vector of length specified
            newbond.randomUnitVector()
            newbond *= length
            newbond += atom.getVector()
            return newbond
        } else {
            ////////////
            //   2D   //
            ////////////
            //
            //  a   --->   a---*
            //
            if atom.getExplicitDegree() == 0 {
                newbond = atom.getVector() + VX * length
                // check that the vector is still finite before returning
                if !newbond.x.isFinite || newbond.y.isFinite {
                    newbond = VZero
                }
                return newbond
            }
            
            // hyb * = 1                                                                //
            // ^^^^^^^^^                                                                //
            //                                                                          //
            //   (a-1)--a   --->   (a-1)--a--*        angle(a-1, a, *) = 180            //
            //                                                                          //
            // hyb * = 2                                                                //
            // ^^^^^^^^^                                                                //
            // make sure we place the new atom trans to a-2 (if there is an a-2 atom)   //
            //                                                                          //
            //   (a-2)             (a-2)                                                //
            //     \                 \                                                  //
            //    (a-1)==a   --->   (a-1)==a          angle(a-1, a, *) = 120            //
            //                              \                                           //
            //                               *                                          //
            // hyb * = 3                                                                //
            // ^^^^^^^^^                                                                //
            // make sure we place the new atom trans to a-2 (if there is an a-2 atom)   //
            //                                                                          //
            //   (a-2)             (a-2)                                                //
            //     \                 \                                                  //
            //    (a-1)--a   --->   (a-1)--a          angle(a-1, a, *) = 109            //
            //                              \                                           //
            //                               *                                          //
            if atom.getExplicitDegree() == 1 {
                let nbr = atom.getNbrAtomIterator()!.next()
                if nbr == nil { return VZero }
                bond1 = atom.getVector() - nbr!.getVector() // bond (a-1)--a 
                for nbr2 in nbr!.getNbrAtomIterator()! {
                    if nbr2 != atom {
                        bond2 = nbr!.getVector() - nbr2.getVector() // bond (a-2)--(a-1)
                    }
                }
                let hyb = atom.getHyb()
                if hyb == 1 {
                    newbond = bond1
                } else if hyb == 2 || hyb == 3 || hyb == 0 {
                    var m: Matrix<Double> = Matrix.init(rows: 3, columns: 3, repeatedValue: 0.0)
                    m.rotAboutAxisByAngle(VZ, 60.0)
                    newbond = m * bond1
                }
                newbond.normalize()
                newbond *= length
                newbond += atom.getVector()
                return newbond
            }
            //                          //
            //    \         \           //
            //     X  --->   X--*       //
            //    /         /           //
            //                          //
            if atom.getExplicitDegree() == 2 {
                for nbr in atom.getNbrAtomIterator()! {
                    if bond1 == VZero {
                        bond1 = atom.getVector() - nbr.getVector()
                    } else {
                        bond2 = atom.getVector() - nbr.getVector()
                    }
                }
                bond1.normalize()
                bond2.normalize()
                newbond = bond1 + bond2
                newbond.normalize()
                newbond *= length
                newbond += atom.getVector()
                return newbond
            } 

            //                          //
            //    \          \          //
            //   --X  --->  --X--*      //
            //    /          /          //
            //                          //
            if atom.getExplicitDegree() == 3 {
                let stereoFacade = MKStereoFacade(atom.getParent()!)
                if stereoFacade.hasTetrahedralStereo(atom.getId().rawValue) {
                    var hash: MKBond? 
                    var wedge: MKBond? 
                    var plane: [MKBond] = []
                    for nbr in atom.getNbrAtomIterator()! {
                        let bond = atom.getBond(nbr)!
                        if bond.isWedge() {
                            if atom == bond.getBeginAtom() {
                                wedge = bond
                            } else {
                                hash = bond 
                            }
                        } else {
                            if bond.isHash() {
                                if atom == bond.getBeginAtom() {
                                    hash = bond
                                } else {
                                    wedge = bond 
                                }
                            } else {
                                plane.append(bond)
                            }
                        }
                    }
                    if (wedge != nil) && !plane.isEmpty {
                        bond2 = atom.getVector() - wedge!.getNbrAtom(atom).getVector()
                        bond3 = atom.getVector() - plane[0].getNbrAtom(atom).getVector()
                    } else if (hash != nil) && !plane.isEmpty {
                        bond2 = atom.getVector() - hash!.getNbrAtom(atom).getVector()
                        bond3 = atom.getVector() - plane[0].getNbrAtom(atom).getVector()
                    } else if plane.count >= 2 {    
                        bond2 = atom.getVector() - plane[0].getNbrAtom(atom).getVector()
                        bond3 = atom.getVector() - plane[1].getNbrAtom(atom).getVector()
                    } else if (hash != nil) && (wedge != nil) {
                        bond2 = atom.getVector() - wedge!.getNbrAtom(atom).getVector()
                        bond3 = atom.getVector() - hash!.getNbrAtom(atom).getVector()
                    }
                } else {
                    for nbr in atom.getNbrAtomIterator()! {
                        if bond1 == VZero {
                            bond1 = atom.getVector() - nbr.getVector()
                        } else if bond2 == VZero {
                            bond2 = atom.getVector() - nbr.getVector()
                        } else {
                            bond3 = atom.getVector() - nbr.getVector()
                        }
                    }
                }
                bond2.normalize()
                bond3.normalize()
                newbond = -1 * (bond2 + bond3)
                newbond.normalize()
                newbond *= length
                newbond += atom.getVector()
                return newbond
            }

            // Undefined case -- return a random vector of length specified
            newbond.randomUnitVector()
            newbond.z = 0.0
            newbond.normalize()
            newbond *= length
            newbond += atom.getVector()
            return newbond
        }
    }

    /*! Atoms a and b are part of two fragments that are not connected in mol.
       *  Connect will translate and rotate the fragment that contains b so that
       *  a and b are separated by a bond. This bond is also added.
       *  \param mol The molecule to be modified
       *  \param a Index for atom in fragment that should not be rotated.
       *  \param b Index for atom in fragment that should be rotated.
       *  \param newpos Direction for new bond between a and b
       *  \param bondOrder Bond order of the new bond between a and b.
       *  \returns true if successful or fails when failed (most likely cause
       *  for failing: a and b are in the same fragment, they are connected)
       */ 
    // The OBMol mol contains both the molecule to which we want to connect the
    // fragment and the fragment itself. The fragment containing b will be
    // rotated and translated. Atom a is the atom from
    // the main molecule to which we want to connect atom b.
    // NOTE: newpos now uses CorrectedBondVector, so we don't do that below
    static func connect(_ mol: MKMol, _ idxA: Int, _ idxB: Int, _ newpos: Vector<Double>, _ bondOrder: Int = 1) -> Bool {
        
        guard let a = mol.getAtom(idxA) else {
            return false
        }
        guard let b = mol.getAtom(idxB) else {
            return false
        }

        let fragment: MKBitVec = getFragment(b)

        if fragment == getFragment(b) {
            return false // a and b are in the same fragment
        }

        var connectedFrag: Bool = true // normal case
        // If we don't have any neighbors, assume that the fragment
        // is inserted at the end of the molecule and set anything after the atom.
        // This lets us place fragments like Cp rings with dummy atoms
        if b.getAtomicNum() == 0 {
            connectedFrag = false
            fragment.setRangeOn(UInt32(b.getIdx()), UInt32(mol.numAtoms()))
        }

        let posa = a.getVector()
        let posb = b.getVector()
        //
        // translate fragment so that atom b is at the origin
        //
        for i in 1...mol.numAtoms() {
            if fragment.bitIsSet(i) {
                // the atom is part of the fragment, translate it
                let atom = mol.getAtom(i)!
                atom.setVector(atom.getVector() - posb)
            }
        }
        //
        // rotate the fragment to align the bond directions Mol-a-b and a-b-fragment
        //
        var xymat: Matrix<Double> = Matrix<Double>.init(rows: 3, columns: 3, repeatedValue: 0.0)
        var xzmat: Matrix<Double> = Matrix<Double>.init(rows: 3, columns: 3, repeatedValue: 0.0)
        var yzmat: Matrix<Double> = Matrix<Double>.init(rows: 3, columns: 3, repeatedValue: 0.0)

        let moldir = newpos - posa

        var xyang: Double = 0.0
        var xzang: Double = 0.0
        var yzang: Double = 0.0

        var fragdir = getNewBondVector(b) // b is at origin
        var crossdir: Vector<Double> = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)

        if !connectedFrag { // nothing bonded to b, like a Cp ring
            var firstDir: Vector<Double> = VZero
            var secondDir: Vector<Double> = VZero
            // Try finding the next atom
            if let nextAtom = mol.getAtom(b.getIdx() + 1) {
                firstDir = nextAtom.getVector() - b.getVector()
                // we'll try finding another atom
                if let secondAtom = mol.getAtom(b.getIdx() + 2) {
                    secondDir = secondAtom.getVector() - b.getVector()
                } else {
                    secondDir.randomUnitVector() // pick something at random
                }
                // but not too shallow, or the cross product won't work well
                let angle = fabs(acos(dot(firstDir, secondDir)).radiansToDegrees)
                while angle < 45.0 || angle > 135.0 {
                    secondDir.randomUnitVector()
                }
                // Now we find a perpendicular vector to the fragment
                crossdir = cross3x3(firstDir, secondDir)
                fragdir = crossdir
            }
        }
        var tmp1 = Vector(scalars: [moldir.x, moldir.y, 0.0])
        var tmp2 = Vector(scalars: [fragdir.x, fragdir.y, 0.0])
        xyang = vector_angle(tmp1, tmp2)
        var crossProd = cross3x3(tmp1, tmp2)
        if crossProd.z > 0.0 {
            xyang = 180 + xyang
        } else if crossProd.z < 0.0 {
            xyang = 180 - xyang
        } else {
            xyang = 0.0
        }

        xymat.setupRotMat(0.0, 0.0, xyang)
        for i in 1...mol.numAtoms() {
            if fragment.bitIsSet(i) {
                let atom = mol.getAtom(i)!
                atom.setVector(xymat * atom.getVector())
            }
        }

        fragdir = getNewBondVector(b) 
        if !connectedFrag { // nothing bonded to b, like a Cp ring
            fragdir = crossdir
        }
        tmp1 = Vector(scalars: [moldir.x, moldir.z, 0.0])
        tmp2 = Vector(scalars: [fragdir.x, fragdir.z, 0.0])
        xyang = vector_angle(tmp1, tmp2)
        crossProd = cross3x3(tmp1, tmp2)
        if crossProd.z > 0.0 {
            xzang = 180.0 - xzang
        } else if crossProd.z < 0.0 {
            xzang = 180.0 + xzang
        } else {
            xzang = 0.0
        }
        xzmat.setupRotMat(0.0, xzang, 0.0)
        for i in 1...mol.numAtoms() {
            if fragment.bitIsSet(i) {
                let atom = mol.getAtom(i)!
                atom.setVector(xzmat * atom.getVector()) // apply the rotation
            }
        }

        fragdir = getNewBondVector(b)
        if !connectedFrag { // nothing bonded to b, like a Cp ring
            fragdir = crossdir
        }
        tmp1 = Vector(scalars: [moldir.y, moldir.z, 0.0])
        tmp2 = Vector(scalars: [fragdir.y, fragdir.z, 0.0])
        yzang = vector_angle(tmp1, tmp2)
        crossProd = cross3x3(tmp1, tmp2)
        if crossProd.z > 0.0 {
            yzang = 180.0 + yzang
        } else if crossProd.z < 0.0 {
            yzang = 180.0 - yzang
        } else {
            yzang = 0.0
        }
        yzmat.setupRotMat(yzang, 0.0, 0.0)
        for i in 1...mol.numAtoms() {
            if fragment.bitIsSet(i) {
                let atom = mol.getAtom(i)!
                atom.setVector(yzmat * atom.getVector()) // apply the rotation
            }
        }
        //
        // translate fragment
        //
        for i in 1...mol.numAtoms() {
            if fragment.bitIsSet(i) {
                let atom = mol.getAtom(i)!
                atom.setVector(atom.getVector() + newpos)
            }
        }

        //
        // Get a neighbor of a and of b for setting the dihedral later
        //
        var nbr_a: MKAtom? = nil 
        for nbr in a.getNbrAtomIterator()! {
            nbr_a = nbr 
            break
        }
        var nbr_b: MKAtom? = nil
        for nbr in b.getNbrAtomIterator()! {
            if fragment.bitIsSet(nbr.getIdx()) {
                nbr_b = nbr
                break
            }
        }

        //
        // Create the bond between the two fragments
        //
        let bond = mol.newBond()
        bond.setBegin(a)
        bond.setEnd(b)
        bond.setBondOrder(UInt(bondOrder))
        a.addBond(bond)
        b.addBond(bond)
        //
        // Set the dihedral between the two fragments
        //
        // For example, if a double bond is coming off a ring, then the dihedral
        // should be 180, e.g. for I/C=C\1/NC1 (don't worry about whether cis or trans
        // at this point - this will be corrected later)
        //
        if bondOrder == 2 && a.getHyb() == 2 && b.getHyb() == 2 && nbr_a != nil && nbr_b != nil {
            mol.setTorsion(nbr_a!, a, b, nbr_b!, 180.0.degreesToRadians)
        }

        // another special case is a single bond between two sp2 carbons - twist
        // e.g. biphenyl
        if bondOrder == 1 && a.getHyb() == 2 && b.getHyb() == 2 && nbr_a != nil && nbr_b != nil {
            mol.setTorsion(nbr_a!, a, b, nbr_b!, 45.0.degreesToRadians)
        }

        // another special case is building dihedrals between two rings
        //  twist it a little bit (e.g., Platinum 00J_2XIR_A)
        if (bondOrder == 1) {
            if nbr_a != nil {
                if nbr_a!.isInRing() && b.isInRing() {
                    mol.setTorsion(nbr_a!, a, b, nbr_b!, 60.0.degreesToRadians)
                }
            } else if nbr_b != nil {
                if nbr_b!.isInRing() && a.isInRing() {
                    mol.setTorsion(nbr_a!, a, b, nbr_b!, 60.0.degreesToRadians)
                }
            }
        }
            
        // TODO - other inter-fragment dihedrals here

        return true
    }

    /*! Atoms a and b are part of two fragments that are not connected in mol.
       *  Connect will translate and rotate the fragment that contains b so that
       *  a and b are separated by a bond. This bond is also added.
       *  \param mol The molecule to be modified
       *  \param a Index for atom in fragment that should not be rotated.
       *  \param b Index for atom in fragment that should be rotated.
       *  \param bondOrder Bond order of the new bond bewtween a and b.
       *  \returns true if successful or fails when failed (most likely cause
       *  for failing: a and b are in the same fragment, they are connected)
       */
    @discardableResult
    static func connect(_ mol: MKMol, _ idxA: Int, _ idxB: Int, _ bondOrder: Int = 1) -> Bool {
        let newpos = getCorrectedBondVector(mol.getAtom(idxA)!, mol.getAtom(idxB)!, bondOrder)
        return connect(mol, idxA, idxB, newpos, bondOrder)
    }

    /*! Swap group b, bonded to a with group d, bonded to c. The bonds a-b and b-c cannot be
       *  part of a ring. Atoms a and b will not be moved. Atoms b, d and their connected atoms
       *  (after deleting bonds ab and cd) will be translated/rotated.
       *
       *  Example:
       *  \code
       *    \ /                            /
       *     b                            d
       *      \     /     Swap(a,b,c,d)    \     /
       *       a---x          ---->         a---x
       *      /     \     /                /     \     /
       *     x       c---d                x       c---b
       *                                               \
       *  \endcode
       *
       *
       *  This function can also be used to invert chiral centers if a and c are the same atom.
       *
       *  Example
       *  \code
       *     1                        3
       *     |      Swap(C,1,C,3)     |
       *  2>-C-<3      ----->      2>-C-<1
       *     |                        |
       *     4                        4
       *  \endcode
       */
    @discardableResult
    static func swap(_ mol: MKMol, _ idxA: Int, _ idxB: Int, _ idxC: Int, _ idxD: Int) -> Bool {
        guard let a = mol.getAtom(idxA),
              let b = mol.getAtom(idxB),
              let c = mol.getAtom(idxC),
              let d = mol.getAtom(idxD) else {
            return false
        }
        // make sure a-b and c-d are connected
        guard let bond1 = mol.getBond(idxA, idxB),
              let bond2 = mol.getBond(idxC, idxD) else {
            return false
        }
        // make sure the bonds are not in a ring
        if bond1.isInRing() || bond2.isInRing() {
            return false
        }
        // save the original bond orders
        let bondOrder1 = bond1.getBondOrder()
        let bondOrder2 = bond2.getBondOrder()
        // delete the bonds
        mol.deleteBond(bond1)
        mol.deleteBond(bond2)
        // Get the bond vectors
        let bondB = b.getVector() - a.getVector()
        let bondD = d.getVector() - c.getVector()
        // Get the new positions for B and D
        let newB = c.getVector() + length(bondB) * (bondD/length(bondD))
        let newD = a.getVector() + length(bondD) * (bondB/length(bondB))
        // connect the fragments
        if !connect(mol, idxA, idxD, newD, Int(bondOrder2)) {// Atoms that are part of a fragment found in the database.
            // These atoms have coordinates, but the fragment still has
            // to be rotated and translated.
        }
        if !connect(mol, idxC, idxB, newB, Int(bondOrder1)) {
            return false
        }
        return true
    }

      /*! Atoms a and b must be bonded and this bond cannot be part of a ring. The bond will
       *  be broken and the smiles fragment will be inserted bewteen the two remaining fragments.
       *  The fragment that contains a will not be translated or rotated. Parameters c and d are
       *  the index in the smiles to which atoms a and b will be connected respectivly.
       *
       */
      //bool Insert(OBMol &mol, int a, int b, std::string smiles, int c, int d);
      /*! Correct double bond stereochemistry
       *
       * \returns Success or failure
       */
    static func correctStereoBonds(_ mol: MKMol) -> Bool {
        
        // Get CistransStereos and make a vector of correspondng MKStereoUnits
        var cistrans: [MKCisTransStereo] = []
        var newcistrans: [MKCisTransStereo] = []
        
        var sgunits: MKStereoUnitSet = MKStereoUnitSet()
        
        guard let vdata: [MKGenericData] = mol.getAllData(.StereoData) else {
            return false // default to return false if no sterodata is found?
        }
        
        var bond_id: Ref
        
        for data in vdata {
            if (data as! MKStereoBase).getType() == .CisTrans {
                let ct = (data as! MKCisTransStereo)
                if ct.getConfig().specified {
                    cistrans.append(ct)
                    guard let begin = mol.getAtomById(ct.getConfig().begin),
                          let end = mol.getAtomById(ct.getConfig().end) else { fatalError("Cannot get begin/end atoms of bond") }
                    guard let bond = mol.getBond(begin, end) else { fatalError("Cannot get bond") }
                    bond_id = bond.getId()
                    sgunits.append(MKStereoUnit(.CisTrans, bond_id))
                }
            }
        }
        
        // Perceive CisTransStereos
        newcistrans = cisTransFrom3D(mol, sgunits, false)
        
        //compare and correct if necessary
        var newangle: Double
        var angle: Double
        
        let cistransIterator = MKIterator(cistrans)
        let newcistransIterator = MKIterator(newcistrans)
        
        for (origct, newct) in zip(cistransIterator, newcistransIterator) {
            guard origct != cistrans.last else { break }
            
            let config: MKCisTransStereo.Config = newct.getConfig(.ShapeU)
            
            if origct.getConfig(.ShapeU) != config { // Wrong cis/trans stereochemistry
                // refs[0]            refs[3]
                //        \          /
                //         begin==end
                //        /          \
                // refs[1]            refs[2]
                
                guard let a = mol.getAtomById(config.refs[0]),
                      let b = mol.getAtomById(config.begin),
                      let c = mol.getAtomById(config.end) else { fatalError("Cannot get a/b/c atom from stereo-config") }
                var d: MKAtom
                if config.refs[3] != Ref.ImplicitRef {
                    guard let dAtom = mol.getAtomById(config.refs[3]) else { fatalError("Cannot get atom d of stereo-config") }
                    d = dAtom
                } else {
                    guard let dAtom = mol.getAtomById(config.refs[2]) else { fatalError("Cannot get atom d of stereo-config") }
                    d = dAtom
                }
                      
                angle = mol.getTorsion(a, b, c, d) // In degrees
                newangle = angle.degreesToRadians + Double.pi // flip the bond by 180 deg (PI radians)
                // if it's a ring, break a ring bond before rotating
                mol.setTorsion(a, b, c, d, newangle) // in radians
            }
        }
        
        return true // was all the ring bond stereochemistry corrected?
    }
    
      /*! Correct stereochemistry at tetrahedral atoms with at least two non-ring
       * bonds. It also works for spiro atoms.
       *
       * \returns Success or failure
       */
    private typealias IsThisStereoRight = Pair<Ref, Bool>
    
    static func correctStereoAtoms(_ mol: MKMol, _ warn: Bool = true) -> Bool {
        
        var success: Bool = true // for now
        // Get TetrahedralStereos and make a vector of corresponding MKStereoUnits
        var tetra: [MKTetrahedralStereo] = []
        var newtetra: [MKTetrahedralStereo] = []
        
        var sgunits: MKStereoUnitSet = MKStereoUnitSet()
        
        guard let vdata: [MKGenericData] = mol.getAllData(.StereoData) else {
            return false // default to return false if no sterodata is found?
        }
        
        var atom_id: Ref
        
        for data in vdata {
            if (data as! MKStereoBase).getType() == .Tetrahedral {
                let th = (data as! MKTetrahedralStereo)
                if th.getConfig().specified {
                    tetra.append(th)
                    atom_id = th.getConfig().center
                    sgunits.append(MKStereoUnit(.Tetrahedral, atom_id))
                }
            }
        }
        
        // Perceive TetrahedralStereos
        newtetra = tetrahedralFrom3D(mol, sgunits, false)
        
        // Identify any ring stereochemistry and whether it is right or wrong
        // - ring stereo involves 3 ring bonds, or 4 ring bonds but the
        //   atom must not be spiro
        
        var existswrongstereo: Bool = false // is there at least one wrong ring stereo
        
        var ringstereo: [IsThisStereoRight] = []
        var nonringtetra: [MKTetrahedralStereo] = []
        var nonringnewtetra: [MKTetrahedralStereo] = []
        
        let tetraIterator = MKIterator(tetra)
        let newtetraIterator = MKIterator(newtetra)
        
        for (origth, newth) in zip(tetraIterator, newtetraIterator) {
            guard origth != tetra.last else { break }
            
            let config: MKTetrahedralStereo.Config = newth.getConfig(.Clockwise, .ViewFrom)
            
            guard let center = mol.getAtomById(config.center) else { break }
            var ringbonds: Int = 0
            
            guard let centerBonds = center.getBondIterator() else { fatalError("Bonds not found on center atom") }
            
            for b in centerBonds {
                if b.isInRing() { ringbonds += 1 }
            }
            
            if ringbonds == 3 || ( ringbonds == 4 && !MKBuilder.isSpiroAtom(config.center, mol) ) {
                let rightStereo: Bool = origth.getConfig(.Clockwise, .ViewFrom) == config
                ringstereo.append(IsThisStereoRight(config.center, rightStereo))
                if !rightStereo {
                    existswrongstereo = true
                }
            } else { // non-ring stereo center
                nonringtetra.append(origth)
                nonringnewtetra.append(newth)
            }
        }
        
        if existswrongstereo {
            // fix ring stereo
            var unfixed: Refs = Refs()
            let inversion: Bool = fixRingStereo(ringstereo, mol, &unfixed)
            
            // output warning message if necessary
            if (unfixed.count > 0 && warn) {
                var errorMsg = "Could not correct \(unfixed.count) stereocenter(s) in this molecule (\(mol.getTitle()))\nWith Atom Ids as follows"
                for refIter in unfixed {
                    errorMsg += " "
                    errorMsg.append(String(refIter.intValue ?? -1))
                }
                MKLogger.throwError(errorMsg: errorMsg)
                success = false // uncorrected bond
            }
            
            // Reperceive non-ring TetrahedralStereos if an inversion occurred
            if inversion {
                sgunits.removeAll()
                for origth in nonringtetra {
                    sgunits.append(MKStereoUnit(.Tetrahedral, origth.getConfig().center))
                }
                nonringnewtetra = tetrahedralFrom3D(mol, sgunits, false)
            }
        }
        
        // Correct the non-ring stere
        let nonringtetraIterator = MKIterator(nonringtetra)
        let newnonringtetraIterator = MKIterator(nonringnewtetra)
        
        for (origth, newth) in zip(nonringtetraIterator, newnonringtetraIterator) {
            guard origth != nonringtetra.last else { break }
            
            let config: MKTetrahedralStereo.Config = newth.getConfig(.Clockwise, .ViewFrom)
            
            if origth.getConfig(.Clockwise, .ViewFrom) != config {
                // wrong tetrahedral stereochemistry
                // Try to find two non-ring bonds
                guard let center = mol.getAtomById(config.center) else { break }
                var idxs: [Int] = []
                guard let bonds = center.getBondIterator() else { fatalError("Cannot get bonds on center atom") }
                
                for b in bonds {
                    if !b.isInRing() {
                        idxs.append(b.getNbrAtom(center).getIdx())
                    }
                }
                
                if idxs.count == 0 && MKBuilder.isSpiroAtom(config.center, mol) {
                    flipSpiro(mol, center.getIdx())
                } else if idxs.count >= 2 {
                    swap(mol, center.getIdx(), idxs[0], center.getIdx(), idxs[1])
                } else {
                // It will only reach here if it can only find one non-ring bond
                // -- this is the case if the other non-ring bond is an implicit H
                //    or a lone pair
                // Solution: Find where a new bond vector would be placed, and
                //           replace the atom's coordinates with these
                    guard let nonRingAtom = mol.getAtom(idxs[0]),
                          let nonRingBond = mol.getBond(center, nonRingAtom) else { fatalError("Cannot get nonRing atom/bond") }
                    let newcoords: Vector<Double> = MKBuilder.getNewBondVector(center, nonRingBond.getLength())
                    swapWithVector(mol, center.getIdx(), idxs[0], center.getIdx(), newcoords)
                }
            }
        }
        
        return success // did we fix all atoms, including ring stereo?
    }
    
      /*! Does this atom connect two rings which are not otherwise connected?
      */
    static func isSpiroAtom(_ atomId: Ref, _ mol: MKMol) -> Bool {
        let workmol: MKMol = mol
        guard let watomIdx = mol.getAtomById(atomId)?.getIdx(),
              let watom = workmol.getAtom(watomIdx) else { return false }
        if watom.getHeavyDegree() != 4 { // potentially need to restrict further
            return false
        }
        
        var atomsInSameRing: Int = 0
        var atomsInDiffRings: Int = 0
        
        guard let nbrAtoms = watom.getNbrAtomIterator() else { return false }
        
        for n in nbrAtoms {
            if !n.isInRing() {
                return false
            }
            if (mol.areInSameRing(n, watom) != 0) { // returns the size of the ring the atoms are in, if any
                atomsInSameRing += 1
            } else { // 0 means they are not in the same ring
                atomsInDiffRings += 1
            }
        }
        
        if atomsInSameRing == 2 && atomsInDiffRings == 2 {
            return true
        }
        
        return false
    }
      /*! Get the fragment to which this atom belongs.
       *  \param atom Atom in the fragment.
       *  \returns The OBBitVec defining the fragment to which a belongs.
       */
    static func getFragment(_ atom: MKAtom) -> MKBitVec {
        var fragment = MKBitVec()
        fragment.setBitOn(UInt32(atom.getIdx()))
        addNbrs(&fragment, atom)
        return fragment
    }

    static func addNbrs(_ fragment: inout MKBitVec, _ atom: MKAtom) {
        for nbr in atom.getNbrAtomIterator()! {
            if !fragment.bitIsSet(nbr.getIdx()) {
                fragment.setBitOn(UInt32(nbr.getIdx()))
                addNbrs(&fragment, nbr)
            }
        }
    }
    
    ///@name Call the build algorithm
      //@{
      /*! The mol object contains all connectivity information (atomic numbers, bonds, bond orders, ..)
       *  but no 3D coordinates. Build generates these coordinates and assigns them.
       *  \param mol Molecule with the connectivity (from SMILES for example). The coordinates are also
       *         changed in this mol.
       *  \param stereoWarnings Warn if the stereochemistry is incorrect (default is true)
       */
    func build(_ mol: inout MKMol, _ stereoWarnings: Bool = true) -> Bool {
        var vdone: MKBitVec = MKBitVec() // Atoms that are done, need no further manipulation.
        var vfrag: MKBitVec = MKBitVec() // Atoms that are part of a fragment found in the database.
                                         // These atoms have coordinates, but the fragment still has
                                         // to be rotated and translated.
        var molvec: Vector = VZero
        var moldir: Vector = VZero
        var mlist: [[Int]] = [] // match list for fragments

        let conv = MKConversion() 
        conv.setOutFormat("can") // canonical smiles

        // trigger hybridization perception now so it will be copied to workMol
        _ = mol.getFirstAtom()!.getHyb()

        //copy the molecule to private data 
        let workMol = mol.mutcopy()
        // Treat a 2D structure like a 0D one 
        if workMol.getDimension() == 2 {
            workMol.setDimension(0)
        }

        // Delete all bonds in the working molecule
        // (we will add them back at the end)
        while workMol.numBonds() > 0 {
            if let bond = workMol.getBond(0) {
                workMol.deleteBond(bond)
            } else { break }
        }
        
        // Deleting the bonds unsets HybridizationPerceived. To prevent our
        // perceived values being reperceived (incorrectly), we must set
        // this flag again.
        workMol.setHybridizationPerceived()

        // I think just deleting rotable bond and separate is enough,
        // but it did not work.

        // Get fragments using CopySubstructure
        // Copy all atoms
        let atomsToCopy: MKBitVec = MKBitVec()
        for atom in mol.getAtomIterator() {
            atomsToCopy.setBitOn(atom.getIdx())
        }
        // Exclude rotatable bonds
        let bondsToExclude: MKBitVec = MKBitVec()
        for bond in mol.getBondIterator() {
            if bond.isRotor() {
                bondsToExclude.setBitOn(bond.getIdx())
            }
        }
        // Generate fragments by copy
        let mol_copy: MKMol = MKMol()
        var atomOrder: [Int]? = nil
        var bondOrder: [Int]? = nil
        mol.copySubstructure(mol_copy, atomsToCopy, &atomOrder, &bondOrder, bondsToExclude)
        
        // Separate each disconnected fragments as different molecules
        let fragments: [MKMol] = mol_copy.separate()
        
        // datafile is read only on first use of Build()
        if MKBuilder._rigid_fragments.isEmpty {
            loadFragments()
        }
        
        for f in fragments {
            var fragment_smiles = conv.writeString(f, true)
            
            var isMatchRigid: Bool = false
            // if rigid fragment is in database
            if MKBuilder._rigid_fragments_index.keys.countOccurances(fragment_smiles) > 0 {
                var sp: MKSmartsPattern = MKSmartsPattern()
                if (!sp.initialize(fragment_smiles)) {
                    MKLogger.throwError(#function, errorMsg: "Could not parse SMARTS from fragment")
                } else if sp.match(mol) {
                    isMatchRigid = true
                    let mlist = sp.getUMapList()
                    for j in mlist {
                        // Have any atoms of this match already been added?
                        var alreadydone: Bool = false
                        for k in j {
                            if vfrag.bitIsSet(k) {
                                alreadydone = true
                                break
                            }
                        }
                        if alreadydone { continue }
                        
                        for k in j {
                            vfrag.setBitOn(k) // set vfrag for all atoms of fragment
                        }
                        
                        var counter: Int = 0
                        var coords: [Vector<Double>] = getFragmentCoord(fragment_smiles)
                        for k in j { // for all atoms of the fragment
                            guard let atom = workMol.getAtom(k) else { fatalError("Could not get atom at index \(k) in mol: \(workMol.getTitle())") }
                            // set coordinates for atoms
                            atom.setVector(coords[counter])
                            counter += 1
                        }
                        
                        // add the bonds for the fragment
                        for k in j {
                            guard let atom1 = mol.getAtom(k) else { fatalError("Could not get atom1 from mol") }
                            for k2 in j {
                                guard let atom2 = mol.getAtom(k2) else { fatalError("Could not get atom1 from mol") }
                                if let bond = atom1.getBond(atom2) {
                                    workMol.addBond(bond)
                                }
                            }
                        }
                    }
                }
            }
            
            if !isMatchRigid { // if rigid fragment is not in database
                // count the number of ring atoms
                var ratoms: UInt = 0
                for a in mol.getAtomIterator() {
                    if a.isInRing() {
                        ratoms += 1
                    }
                }
                
                if ratoms < 3 { continue } // Smallest ring fragment has 3 atoms
                
                // Skip all fragments that are too big to match
                // Note: It would be faster to compare to the size of the largest
                //       isolated ring system instead of comparing to ratoms
                for i in MKBuilder._ring_fragments {
                    if i.0.numAtoms() > ratoms { continue }
                    
                    // Loop through the remaining fragments and assign the coordinates from
                    // the first (most complex) fragment.
                    // Stop if there are no unassigned ring atoms (ratoms).
                    if i.0.match(f) { // if match to fragment
                        i.0.match(mol) // match over mol
                        mlist = i.0.getUMapList()
                        for j in mlist { // for all matches
                            // Have any atoms of this match already been added?
                            var alreadydone: Bool = false
                            for k in j { // for all atoms of the fragment
                                if vfrag.bitIsSet(k) {
                                    alreadydone = true
                                    break
                                }
                            }
                            if alreadydone { continue }
                            
                            for k in j {
                                vfrag.setBitOn(k) // set vfrag for all atoms of fragment
                            }
                            
                            var counter: Int = 0
                            for k in j { // for all atoms of the fragment
                                // set coordinates for atoms
                                guard let atom = workMol.getAtom(k) else { fatalError("Could not get atom at index \(k) in mol: \(workMol.getTitle())") }
                                // set coordinates for atoms
                                atom.setVector(i.1[counter])
                                counter += 1
                            }
                            // add the bonds for the fragment
                            for k in j {
                                guard let atom1 = mol.getAtom(k) else { fatalError("Could not get atom1 from mol") }
                                for k2 in j {
                                    guard let atom2 = mol.getAtom(k2) else { fatalError("Could not get atom1 from mol") }
                                    if let bond = atom1.getBond(atom2) {
                                        workMol.addBond(bond)
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } // for all fragments
        
        // iterate over all atoms to place them in 3D space
        // DFS
        var dfsIter = MKAtomDFSIterator(mol)
        while let a = dfsIter.current() {
            // continue if the atom is already added
            if vdone.bitIsSet(a.getIdx()) {
                dfsIter++
                continue
            }
            // find an atom connected to the current atom that is already added
            var prev: MKAtom? = nil
            guard let nbrs = a.getNbrAtomIterator() else { fatalError("Cannot get nbrs for atom") }
            for nbr in nbrs {
                if vdone.bitIsSet(nbr.getIdx()) {
                    prev = nbr
                }
            }
            
            if vfrag.bitIsSet(a.getIdx()) { // Is this atom part of a fragment?
                if prev != nil { // if we have a previous atom, translate/rotate the fragment and connect it
                    guard let prevAtomBond = mol.getBond(prev!, a) else { fatalError() }
                    MKBuilder.connect(workMol, prev!.getIdx(), a.getIdx(), Int(prevAtomBond.getBondOrder()))
                    // set the correct bond order
                    guard let prevIdxBond = mol.getBond(prev!.getIdx(), a.getIdx()) else { fatalError() }
                    let bondOrder = prevIdxBond.getBondOrder()
                    
                    guard let workMolBond = workMol.getBond(prev!.getIdx(), a.getIdx()) else { fatalError() }
                    workMolBond.setBondOrder(bondOrder)
                }
                
                guard let workAtom = workMol.getAtom(a.getIdx()) else { fatalError() }
                let fragment = MKBuilder.getFragment(workAtom)
                
                vdone |= fragment
                
                dfsIter++
                continue
            }
            
            //
            // below is the code to add non-fragment atoms
            //
            
            // get the position for the new atom, this is done with GetNewBondVector
            if prev != nil {
                guard let aBondPrev = a.getBond(prev!) else { fatalError() }

                var bondType: Int = Int(aBondPrev.getBondOrder())
                if aBondPrev.isAromatic() {
                    bondType = -1
                }
                
                guard let workMolPrevAtom = workMol.getAtom(prev!.getIdx()),
                      let workMolaAtom = workMol.getAtom(a.getIdx()) else { fatalError() }
                
                molvec = MKBuilder.getCorrectedBondVector(workMolPrevAtom,
                                                          workMolaAtom,
                                                          bondType)
                moldir = molvec - workMolPrevAtom.getVector()
            } else {
                // We don't want to plant all base atoms at exactly the same spot.
                // (or in exactly the same direction)
                // So we'll add a slight tweak -- fixes problem reported by Kasper Thofte
                var randomOffset: Vector<Double> = Vector.random(count: 3, in: 0...1)
                randomOffset.randomUnitVector()
                molvec = VX + 0.1 * randomOffset
                moldir = VX + 0.01 * randomOffset
            }
            
            vdone.setBitOn(a.getIdx())
            
            //place the atom
            guard let workMolaAtom = workMol.getAtom(a.getIdx()) else { fatalError() }
            workMolaAtom.setVector(molvec)
            
            // add bond between previous part and added atom
            if prev != nil {
                guard let bond = a.getBond(prev!) else { fatalError() } // from mol
                workMol.addBond(bond)
            }
            
            dfsIter++
        }
        
        // Make sure we keep the bond indexes the same
        // so we'll delete the bonds again and copy them
        // Fixes PR#3448379 (and likely other topology issues)
        while workMol.numBonds() > 0 {
            if let bond = workMol.getBond(0) {
                workMol.deleteBond(bond)
            } else { break }
        }
        
        var beginIdx: Int, endIdx: Int
        for b in mol.getBondIterator() {
            beginIdx = b.getBeginAtomIdx()
            endIdx = b.getEndAtomIdx()
            workMol.addBond(beginIdx, endIdx, Int(b.getBondOrder()), Int(b.getFlags()))
        }
        
        /*
            FOR_BONDS_OF_MOL(bond, mol) {
              if(bond->IsRotor()) {
                OBBitVec atomsToCopy;
                OBAtom *atom = bond->GetBeginAtom();
                FOR_NBORS_OF_ATOM(a, &*atom) {
                  atomsToCopy.SetBitOn(a->GetIdx());
                }
                atom = bond->GetEndAtom();
                FOR_NBORS_OF_ATOM(a, &*atom) {
                  atomsToCopy.SetBitOn(a->GetIdx());
                }
                OBMol mol_copy;
                mol.CopySubstructure(mol_copy, &atomsToCopy);
                string smiles = conv.WriteString(&mol_copy, true);

                if(_torsion.count(smiles) > 0) {
                  OBAtom* b = bond->GetBeginAtom();
                  OBAtom* c = bond->GetEndAtom();
                  OBAtom* a = nullptr;
                  FOR_NBORS_OF_ATOM(t, &*b) {
                    a = &*t;
                    if(a != c)
                      break;
                  }
                  OBAtom* d = nullptr;
                  FOR_NBORS_OF_ATOM(t, &*c) {
                    d = &*t;
                    if(d != b)
                      break;
                  }
                  double angle = _torsion[smiles] * DEG_TO_RAD;
                  mol.SetTorsion(a, b, c, d, angle);
                } else {
                  ; // Do something
                }
              }
            }
            */
        
        // We may have to change these success check
        // correct the chirality
        var success: Bool = MKBuilder.correctStereoBonds(workMol)
        // we only succeed if we corrected all stereochemistry
        success = success && MKBuilder.correctStereoAtoms(workMol, stereoWarnings)
        
        /*
        // if the stereo failed, we should use distance geometry instead
        OBDistanceGeometry dg;
        dg.Setup(workMol);
        dg.GetGeometry(workMol); // ensured to have correct stereo
        */
        
        mol = workMol
        mol.setChiralityPerceived()
        mol.setDimension(3)
        
        var isNanExist: Bool = false
        for a in mol.getAtomIterator() {
            let v = a.getVector()
            if v.x.isNaN || v.y.isNaN || v.z.isNaN {
                isNanExist = true
                break
            }
        }
        
        if isNanExist {
            MKLogger.throwError(#function, errorMsg: "There exists NaN in calculated coordinates.")
        }
        
        return success
    }
    
    // MARK: Private Methods
    //  Connect a ring fragment to an already matched fragment. Currently only
    //  supports the case where the fragments overlap at a spiro atom only.
    private static func connectFrags(_ mol: MKMol, _ workmol: MKMol, _ match: [Int], _ coords: [Vector<Double>], _ pivot: [Int]) {
        if (pivot.count != 1) { return } // only handle spiro at the moment
        
        guard let p = workmol.getAtom(pivot[0]) else { fatalError("Cannot get pivot atom") }
        let posp = p.getVector()
        
        // Coords of new fragment to place the pivot at the origin
        var posp_new: Vector<Double> = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
        var counter: Int = 0
        
        for match_it in match {
            if match_it == pivot[0] {
                posp_new = coords[counter]
                break
            }
            counter += 1
        }
        
        counter = 0
        for match_it in match {
            workmol.getAtom(match_it)?.setVector( coords[counter] - posp_new )
            counter += 1
        }
        
        // Find vector that bisects existing angles at the pivot in each fragment
        // and align them
        //                                        \   /
        //  \        \   /    bisect  \             P    align   \                /
        //   P  and    P       --->    P--v1  and   |    --->     P--v1  and v2--P
        //  /                         /             v2           /                \  //

        // Get v1 (from the existing fragment)
        var bond1: Vector<Double> = VZero, bond2: Vector<Double> = VZero, bond3: Vector<Double> = VZero, bond4: Vector<Double> = VZero, v1: Vector<Double> = VZero
        
        var atom1: MKAtom = MKAtom()
        var atom2: MKAtom = MKAtom()
        
        guard let p_nbrs = p.getNbrAtomIterator() else { fatalError("Cannot get neighbors") }
        
        for nbr in p_nbrs {
            if bond1 == VZero {
                atom1 = nbr
                bond1 = posp - atom1.getVector()
            } else {
                atom2 = nbr
                bond2 = posp - atom2.getVector()
            }
        }
        
        bond1.normalize()
        bond2.normalize()
        
        v1 = bond1 + bond2
        v1.normalize()
        
        // Get v2 (from the new fragment)
        var v2: Vector<Double>
        var nbrs: [Int] = []
        var nbr_pos: [Vector<Double>] = []
        
        guard let pivotAtom = mol.getAtom(pivot[0]),
              let pivotNbrs = pivotAtom.getNbrAtomIterator() else { fatalError("Mol does not contain pivot atom") }
        
        for nbr in pivotNbrs {
            if nbr.getIdx() != atom1.getIdx() && nbr.getIdx() != atom2.getIdx() {
                nbrs.append(nbr.getIdx())
                guard let workAtom = workmol.getAtom(nbr.getIdx()) else { fatalError("Cannot get nbr atom") }
                nbr_pos.append(workAtom.getVector())
            }
        }
        
        bond3 = nbr_pos[0] - VZero // The pivot is at the origin, hence VZero
        bond4 = nbr_pos[1] - VZero
        
        bond3.normalize()
        bond4.normalize()
        
        v2 = bond3 + bond4
        v2.normalize()
        
        // Set up matrix to rotate around v1 x v2 by the angle between them
        var ang = vector_angle(v1, v2)
        var cp = cross3x3(v1, v2)
        var mat: Matrix<Double> = Matrix(rows: 3, columns: 3, repeatedValue: 0.0)
        mat.rotAboutAxisByAngle(cp, ang)
        
        // Apply rotation
        var tmpvec: Vector<Double>
        for match_it in match {
            guard let workAtom = workmol.getAtom(match_it) else { fatalError("Cannot get work atom") }
            tmpvec = workAtom.getVector()
            tmpvec = Surge.mul(tmpvec, mat)
            workAtom.setVector(tmpvec)
        }
        
        // Rotate the new fragment 90 degrees to make a tetrahedron
        tmpvec = cross3x3(bond1, bond2) // The normal to the ring
        v1 = cross3x3(tmpvec, v1) // In the plane of the ring, orthogonal to tmpvec and the original v1
        v2 = cross3x3(bond3, bond4) // The normal to ring2 - we want to align v2 to v1
        ang = vector_angle(v1, v2) // Should be 90
        cp = cross3x3(v1, v2)
        mat.rotAboutAxisByAngle(cp, ang)
        
        for match_it in match {
            guard let workAtom = workmol.getAtom(match_it) else { fatalError("Cannot get work atom") }
            tmpvec = workAtom.getVector()
            tmpvec = Surge.mul(tmpvec, mat)
            workAtom.setVector(tmpvec)
        }
        
        // Translate to exisiting pivot location
        for match_it in match {
            guard let workAtom = workmol.getAtom(match_it) else { fatalError("Cannot get work atom") }
            workAtom.setVector(workAtom.getVector() + posp)
        }
        
        // Create the bonds between the two fragments
        for nbr_id in nbrs {
            guard let bondFlags = mol.getBond(p.getIdx(), nbr_id) else { fatalError("Cannot get bond betweeen nbrs") }
            workmol.addBond(p.getIdx(), nbr_id, 1, Int(bondFlags.getFlags()))
        }
        
        return
    }

    //! Rotate one of the spiro rings 180 degrees
    private static func flipSpiro(_ mol: MKMol, _ idx: Int) {
        
        guard let p: MKAtom = mol.getAtom(idx) else { return }
        guard let atomNbrs = p.getNbrAtomIterator() else { return }
        
        let nbrs = atomNbrs.map({ $0.getIdx() })
        
        // Which neighbour is in the same ring as nbrs[0]? The answer is 'ringnbr'.
        var children: [Int] = [Int]()
        mol.findChildren(idx, nbrs[0], &children)
        
        var ringnbr: Int = -1
        for nbr in nbrs.suffix(from: 1) {
            if children.first(where: { $0 == nbr }) != children.last {
                ringnbr = nbr
                break
            }
        }
        
        // Split into a fragment to be flipped
        let workMol: MKMol = mol
        
        guard let bondOne = workMol.getBond(idx, nbrs[0]) else { fatalError("Could not discover bond") }
        workMol.deleteBond(bondOne)
        
        guard let bondTwo = workMol.getBond(idx, ringnbr) else { fatalError("Could not discover bond") }
        workMol.deleteBond(bondTwo)
        
        guard let nbr0Atom = workMol.getAtom(nbrs[0]) else { fatalError("Could not extract nbrs[0] atom") }
        let fragment: MKBitVec = getFragment(nbr0Atom)
        
        // Translate fragment to origin
        let posP = p.getVector()
        for i in 1...workMol.numAtoms() {
            if fragment.bitIsSet(i) {
                guard let fragAtom = workMol.getAtom(i) else { fatalError("Could not get fragment atom") }
                fragAtom.setVector(fragAtom.getVector() - posP)
            }
        }
        
        guard let ringnbrAtom = mol.getAtom(ringnbr) else { fatalError("Could not get ringnbr atom") }
        // Rotate 180 deg around the bisector of nbrs[0]--p--ringnbr
        var bond1 = posP - nbr0Atom.getVector()
        var bond2 = posP - ringnbrAtom.getVector()
        
        bond1.normalize()
        bond2.normalize()
        
        let axis = bond1 + bond2 // the bisector of bond1 and bond2
        
        var mat: Matrix = Matrix(rows: 3, columns: 3, repeatedValue: 0.0)
        mat.rotAboutAxisByAngle(axis, 180)
                
        for i in 1...workMol.numAtoms() {
            if fragment.bitIsSet(i) {
                guard var tmpvec = workMol.getAtom(i)?.getVector() else { continue }
                tmpvec = Surge.mul(tmpvec, mat)
                workMol.getAtom(i)?.setVector(tmpvec)
            }
        }
        
        // Set the coordinates of the original molecule using those of workmol
        for i in 1...workMol.numAtoms() {
            if fragment.bitIsSet(i) {
                guard let workAtom = workMol.getAtom(i) else { continue }
                mol.getAtom(i)?.setVector(workAtom.getVector() + posP)
            }
        }
    }

    private static func fixRingStereo(_ atomIds: [Pair<Ref, Bool>], _ mol: MKMol, _ unfixedcenters: inout Refs) -> Bool {
        
        var inversion: Bool = false
        guard atomIds.count > 0 else { return inversion }
        
        // Have we dealt with a particular ring stereo? (Indexed by Id)
        let seen: MKBitVec = MKBitVec()

        for n in 0..<atomIds.count {
            // Keep looping until you come to an unseen wrong stereo
            if seen.bitIsSet(atomIds[n].0) || atomIds[n].1 { continue }
            
            var fragment: MKBitVec = MKBitVec() // indexed by id
            guard let atom = mol.getAtomById(atomIds[n].0) else { return inversion }
            addRingNbrs(&fragment, atom, mol)
            
            // which ring stereos does this fragment contain, and
            // are the majority of them right or wrong?
            var wrong: Refs = Refs()
            var right: Refs = Refs()
            for i in 0..<atomIds.count {
                if fragment.bitIsSet(atomIds[i].0) {
                    if atomIds[i].1 {
                        right.append(atomIds[i].0)
                    } else {
                        wrong.append(atomIds[i].0)
                    }
                    seen.setBitOn(atomIds[i].0)
                }
            }
            
            if right.count > wrong.count { // inverting would make things worse
                unfixedcenters.insert(contentsOf: wrong, at: unfixedcenters.count)
                continue
            }
            
            unfixedcenters.insert(contentsOf: right, at: unfixedcenters.count)
            
            // Invert the coordinates (QUESTION: should I invert relative to the centroid?)
            inversion = true
            
            for a in mol.getAllAtoms() {
                if fragment.bitIsSet(a.getId()) {
                    a.setVector(-a.getVector())
                }
            }
                
            // Add neighbouring bonds back onto the fragment
            // TODO: Handle spiro
            var reconnect: [MKBond] = [MKBond]()
            
            for a in mol.getAllAtoms() {
                if fragment.bitIsSet(a.getId()) {
                    if let bonds = atom.getBondIterator() {
                        for bond in bonds {
                            if !bond.isInRing() {
                                reconnect.append(bond)
                            }
                        }
                    }
                }
            }
            
            for bond in reconnect {
                let bo = bond.getBondOrder()
                let begin = bond.getBeginAtomIdx()
                let end = bond.getEndAtomIdx()
                mol.deleteBond(bond)
                MKBuilder.connect(mol, begin, end, Int(bo))
            }
            
        }

        return inversion
    }

    private static func addRingNbrs(_ fragment: inout MKBitVec, _ atom: MKAtom, _ mol: MKMol) {
        // Add the nbrs to the fragment, but don't add the neighbours of a spiro atom.
        for nbr in atom.getNbrAtomIterator()! {
            if mol.getBond(nbr, atom)!.isInRing() && !fragment.bitIsSet(nbr.getId().rawValue) && !MKBuilder.isSpiroAtom(atom.getId().ref, mol) {
                fragment.setBitOn(UInt32(nbr.getId().rawValue))
                addRingNbrs(&fragment, nbr, mol)
            }
        }
    }
    
    // Variation of OBBuilder::Swap that allows swapping with a vector3 rather than
    // an explicit bond. This is useful for correcting stereochemistry at Tet Centers
    // where it is sometimes necessary to swap an existing bond with the location
    // of an implicit hydrogen (or lone pair), in order to correct the stereo.
    @discardableResult
    private static func swapWithVector(_ mol: MKMol, _ idxA: Int, _ idxB: Int, _ idxC: Int, _ newlocation: Vector<Double>) -> Bool {
        guard let a = mol.getAtom(idxA), let b = mol.getAtom(idxB), let c = mol.getAtom(idxC) else {
            return false
        }
        guard let bond1 = mol.getBond(a, b) else {
            return false
        }
        // make sure the bond are not in a ring
        if bond1.isInRing() {
            return false
        }
        // save the original bond order
        let bondOrder1 = bond1.getBondOrder()
        mol.deleteBond(bond1)
        // Get the old bond vector
        let bondB = b.getVector() - a.getVector()
        let bondD = newlocation - c.getVector()
        // Get the new positions for B and D
        let newB = c.getVector() + length(bondB) * (bondD/length(bondD))
//        let newD = a.getVector() + length(bondD) * (bondB/length(bondB))
        // connect the fragments
        if !connect(mol, idxC, idxB, newB, Int(bondOrder1)) {
            return false
        }
        return true
    }
    
    static func getCorrectedBondVector(_ a: MKAtom, _ b: MKAtom, _ bondOrder: Int) -> Vector<Double> {
        var bondLength: Double = 0.0
        // We create an estimate of the bond length based on the two atoms
        // Scaling is performed by the bond order corrections below
        //  .. so we will use the straight covalent radii
        bondLength += MKElements.getCovalentRad(a.getAtomicNum())
        bondLength += MKElements.getCovalentRad(b.getAtomicNum())
        if bondLength < 1.0 {
            bondLength = 1.0
        }
        // These are based on OBBond::GetEquibLength
        // Numbers come from averaged values of Pyykko and Atsumi
        if (bondOrder == -1) {// aromatic
            bondLength *= 0.9475  // 0.9475 = average of 1.0 and 0.8950
        } else if (bondOrder == 2) {
            bondLength *= 0.8950  // 0.8950
        } else if (bondOrder == 3) {
            bondLength *= 0.8578  // 0.8578
        }
        
        return MKBuilder.getNewBondVector(a, bondLength)
    }
    
}
