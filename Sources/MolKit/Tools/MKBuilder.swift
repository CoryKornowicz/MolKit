

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
        fatalError()
    }
    
      /*! Correct stereochemistry at tetrahedral atoms with at least two non-ring
       * bonds. It also works for spiro atoms.
       *
       * \returns Success or failure
       */
    static func correctStereoAtoms(_ mol: MKMol, _ warn: Bool = true) -> Bool {
        
        var success: Bool = true // for now
        // Get TetrahedralStereos and make a vector of corresponding MKStereoUnits
        
        fatalError()
    }
    
      /*! Does this atom connect two rings which are not otherwise connected?
      */
    static func isSpiroAtom(_ atomId: Int, _ mol: MKMol) -> Bool {
        fatalError()
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
    func build(_ mol: MKMol, _ stereWarnings: Bool = true) -> Bool {
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

        
        return false
    }
    
    // Private Methods
    //! Connect a ring fragment to an already matched fragment. Currently only
    //  supports the case where the fragments overlap at a spiro atom only.
    private static func connectFrags(_ mol: MKMol, _ workmol: MKMol, _ match: [Int], _ coords: [Vector<Double>], _ pivot: [Int]) {
        
    }

    //! Rotate one of the spiro rings 180 degrees
    private static func flipSpiro(_ mol: MKMol, _ idx: Int) {
        
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
            if mol.getBond(nbr, atom)!.isInRing() && !fragment.bitIsSet(nbr.getId().rawValue) && !MKBuilder.isSpiroAtom(atom.getId().rawValue, mol) {
                fragment.setBitOn(UInt32(nbr.getId().rawValue))
                addRingNbrs(&fragment, nbr, mol)
            }
        }
    }
    
    // Variation of OBBuilder::Swap that allows swapping with a vector3 rather than
    // an explicit bond. This is useful for correcting stereochemistry at Tet Centers
    // where it is sometimes necessary to swap an existing bond with the location
    // of an implicit hydrogen (or lone pair), in order to correct the stereo.
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
