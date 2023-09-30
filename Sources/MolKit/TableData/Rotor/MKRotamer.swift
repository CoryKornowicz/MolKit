//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 3/7/23.
//

import Foundation
/** \class OBRotamerList rotamer.h <openbabel/rotamer.h>
 
 A high-level class for rotamer / conformer generation, intended mainly
 for use with the related class OBRotorList and the OBRotorRules database
 
 Rotamers represent conformational isomers formed simply by rotation of
 dihedral angles. On the other hand, conformers may include geometric
 relaxation (i.e., slight modification of bond lengths, bond angles, etc.)
 
 The following shows an example of generating 2 conformers using different
 rotor states. Similar code could be used for systematic or Monte Carlo
 conformer sampling when combined with energy evaluation (molecular
 mechanics or otherwise).
 
 \code
 OBRotorList rl; // used to sample all rotatable bonds via the OBRotorRules data
 // If you want to "fix" any particular atoms (i.e., freeze them in space)
 // then set up an OBBitVec of the fixed atoms and call
 // rl.SetFixAtoms(bitvec);
 rl.Setup(mol);
 
 // How many rotatable bonds are there?
 cerr << " Number of rotors: " << rl.Size() << endl;
 
 // indexed from 1, rotorKey[0] = 0
 std::vector<int> rotorKey(rl.Size() + 1, 0);
 
 // each entry represents the configuration of a rotor
 // e.g. indexes into OBRotor::GetResolution() -- the different angles
 //   to sample for a rotamer search
 for (unsigned int i = 0; i < rl.Size() + 1; ++i)
 rotorKey[i] = 0; // could be anything from 0 .. OBRotor->GetResolution().size()
 // -1 is for no rotation
 
 // The OBRotamerList can generate conformations (i.e., coordinate sets)
 OBRotamerList rotamers;
 rotamers.SetBaseCoordinateSets(mol);
 rotamers.Setup(mol, rl);
 
 rotamers.AddRotamer(rotorKey);
 rotorKey[1] = 2; // switch one rotor
 rotamers.AddRotamer(rotorKey);
 
 rotamers.ExpandConformerList(mol, mol.GetConformers());
 
 // change the molecule conformation
 mol.SetConformer(0); // rotorKey 0, 0, ...
 conv.Write(&mol);
 
 mol.SetConformer(1); // rotorKey 0, 2, ...
 
 \endcode
 
 **/
class MKRotamerList: MKGenericData {
    private var _NBaseCoords: UInt = 0                 //! Number of atoms in the base coordinate set (i.e., OBMol::NumAtoms())
    private var _c: [[Double]]  = []                    //! Base coordinate sets (i.e., existing conformers to be modified)
    private var _vrotor: [Pair<[MKAtom], [Int]>] = [] //! Individual bond rotors (from an OBRotor object or other)
    private var _vres: [[Double]] = []                //! \brief Index of each rotor's different sampling states ("resolution")
                                              //! Usually from OBRotor::GetResolution()
    private var _vrotamer: [[Int]] = []               //! Individual rotamer states (i.e., the array of rotor settings)
    private var _vrings: [[Int]] = []                 //! Rotors in rings
    private var _vringTors: [[Double]] = []           //! Dihedral angles of ring bonds

    init() {
        super.init("RotamerList", .RotamerList, .any)
        _NBaseCoords = 0 
    }
    
    //! Set up a rotamer list based on an already created OBRotorList
    func setup(_ mol: MKMol, _ rlist: MKRotorList) {
        fatalError()
    }
    //! Set up a rotamer list based on the supplied reference atoms and the number of rotors
    //! \param mol The molecule to evaluate
    //! \param ref An array of the 4 dihedral atoms for each rotor
    //! \param nrotors The number of rotors (i.e., the size of ref / 4)
    func setup(_ mol: MKMol, _ ref: [Int], _ nrotors: Int) {
        fatalError()
    }
    
    //! \return the number of rotatable bonds considered
    func numRotors() -> Int {
        return _vrotor.count
    }
    
    //! \return the number of rotamer (conformation) coordinate sets
    func numRotamers() -> Int {
        return _vrotamer.count
    }
    
    func numAtom() -> UInt {
        return _NBaseCoords
    }
    
    //! \return The number of "base" coordinate sets (i.e., the number of conformers in the base OBMol)
    func numBaseCoordinateSets() -> Int {
        return _c.count
    }
    
    //! Get specific conformer
    func getBaseCoordinateSet(_ i: UInt) -> [Double]? {
        return i < _c.count ? _c[Int(i)] : nil
    }
    
    //! Add a rotamer to the list based on the supplied coordinate set as a double*
//    void AddRotamer(double*);
    
//    //! Add a rotamer to the list based on @p key as a configuration of the individual rotor bonds
//    void AddRotamer(int *key);
    
//    //! Add a rotamer to the list based on @p key as a configuration of the individual rotor bonds
//    void AddRotamer(std::vector<int> key);
    
//    //! Add a rotamer to the list based on @p key as a configuration of the individual rotor bonds
//    void AddRotamer(unsigned char *key);
    
//    //! Add @p nconf rotamers based on @p as an array of configurations much like AddRotamer()
//    void AddRotamers(unsigned char *arr,int nconf);
    
//    //! \return A reference array (as used by AddRotamer() as a configuration of the individual rotor bonds
//    void GetReferenceArray(unsigned char*) const;
    func getReferenceArray(_ ref: inout [Int]) {
        var j: Int = 0
        for i in _vrotor {
            ref[j++] = i.0[0].getIdx()
            ref[j++] = i.0[1].getIdx()
            ref[j++] = i.0[2].getIdx()
            ref[j++] = i.0[3].getIdx()
        }
    }
    
    //! \brief Create a conformer list using the internal base set of coordinates
    //! \return The set of coordinates by rotating the bonds in each rotamer
    func createConformerList(_ mol: MKMol) -> [[Double]] {
        fatalError()
    }
    
    //! \brief Create a conformer list using the internal base set of coordinates
    //! \return The set of coordinates as a reference in @p confs
    func expandConformerList(_ mol: MKMol, _ confs: inout [[Double]]) {
        fatalError()
    }
    
    func setCurrentCoordinates(_ mol: MKMol, _ arr: [Int]) -> Bool {
        var angle: Double
        if arr.count != _vrotor.count + 1 {
            return false // wrong size key
        }
        
        // gotta check for weird ring torsion combinatins
        if !_vrings.isEmpty {
            // go through each ring and update teh possible torsions
            for j in 0..<_vrings.count {
                var path = _vrings[j]
                var torsionSum: Double = 0.0
                
                // go around the loop and add up the torsions
                for i in 0..<path.count {
                    if path[i] == -1 {
                        // not a rotor, use the fixed value
                        torsionSum += _vringTors[j][i]
                        continue
                    }
                    // what angles are we trying to use with this key?
                    angle = _vres[ path[i] ][arr[path[i] + 1]]
                    while angle < 0.0 {
                        angle += 360.0
                    }
                    while angle > 360.0 {
                        angle -= 360.0
                    }
                    
                    //updaate the ring torsion for this setting
                    _vringTors[j][i] = angle
                    torsionSum += angle
                }
                
                // if the sum of the ring torsions is not ~0, bad move
                if fabs(torsionSum) > 45.0 {
                    return false // bad move, do not take
                }
            }
        }
        
        for i in 0..<_vrotor.count {
            if arr[i+1] == -1 { // skip this rotor
                continue
            } else {
                angle = _vres[i][arr[i+1]]
                while angle < 0.0 {
                    angle += 360.0
                }
                while angle > 360.0 {
                    angle -= 360.0
                }
                setRotorToAngle(&mol.coordinates, _vrotor[i].0, angle, _vrotor[i].1)
            }//set an angle
        }// for rotors
        
        return true
    }
    
    //! \brief Copies the mol's conformers (the coordinates, NOT the pointers)
    //! into the object as base coordinates
    func setBaseCoordinates(_ mol: MKMol) {
        setBaseCoordinateSet(mol.getConformers(), N: UInt(mol.numAtoms()))
    }
    
    //! Copies the coordinates in bc, NOT the pointers, into this object
    /** 
     \param bc The conformer set for the molecule
     \param N  The number of atoms in the molecule
     **/
    func setBaseCoordinateSet(_ bc: [[Double]], N: UInt) {
        _c.removeAll()
        _c = bc
        _NBaseCoords = N
    }

//    override func clone<T>(_ mol: T) -> MKGenericData? where T : MKBase {
//
//    }
}

//! Rotate the coordinates of 'atoms'
//! such that tor == ang.
//! Atoms in 'tor' should be ordered such that the 3rd atom is
//! the pivot around which atoms rotate (ang is in degrees)
func setRotorToAngle(_ c: inout [Double], _ ref: [MKAtom], _ ang: Double, _ atoms: [Int]) {
    var v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z: Double
    var c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z: Double
    var c1mag, c2mag, radang, costheta: Double
    var x, y, z, mag, rotang, sn, cs, t, tx, ty, tz: Double
    var m: [Double] = Array(repeating: 0.0, count: 9)
    
    let tor: [Int] = [
        Int(ref[0].getCoordinateIdx()),
        Int(ref[1].getCoordinateIdx()),
        Int(ref[2].getCoordinateIdx()),
        Int(ref[3].getCoordinateIdx())
    ]
    
    // calculate the torsion angle
    v1x = c[tor[0]] - c[tor[1]]
    v2x = c[tor[1]] - c[tor[2]]
    v1y = c[tor[0] + 1] - c[tor[1] + 1]
    v2y = c[tor[1] + 1] - c[tor[2] + 1]
    v1z = c[tor[0] + 2] - c[tor[1] + 2]
    v2z = c[tor[1] + 2] - c[tor[2] + 2]
    v3x = c[tor[2]] - c[tor[3]]
    v3y = c[tor[2] + 1] - c[tor[3] + 1]
    v3z = c[tor[2] + 2] - c[tor[3] + 2]
    
    c1x = v1y * v2z - v1z * v2y
    c2x = v2y * v3z - v2z * v3y
    c1y = -v1x * v2z + v1z * v2x
    c2y = -v2x * v3z + v2z * v3x
    c1z = v1x * v2y - v1y * v2x
    c2z = v2x * v3y - v2y * v3x
    c3x = c1y * c2z - c1z * c2y
    c3y = -c1x * c2z + c1z * c2x
    c3z = c1x * c2y - c1y * c2x
    
    c1mag = c1x * c1x + c1y * c1y + c1z * c1z
    c2mag = c2x * c2x + c2y * c2y + c2z * c2z
    if c1mag * c2mag < 0.01 {
        costheta = 1.0
    } else {
        costheta = (c1x * c2x + c1y * c2y + c1z * c2z) / sqrt(c1mag * c2mag)
    }
    
    if costheta < -0.999999 {
        costheta = -0.999999
    }
    if costheta > 0.999999 {
        costheta = 0.999999
    }
    
    if (v2x * c3x + v2y * c3y + v2z * c3z) > 0.0 {
        radang = -acos(costheta)
    } else {
        radang = acos(costheta)
    }
    
    // now we have the torsion angle (radang) - set up the rot matrix
    rotang = ang.degreesToRadians - radang
    
    sn = sin(rotang)
    cs = cos(rotang)
    t = 1 - cs
    mag = sqrt(v2x * v2x + v2y * v2y + v2z * v2z)
    if mag < 0.1 {
        mag = 0.1
    }
    x = v2x / mag
    y = v2y / mag
    z = v2z / mag
    
    // set up the rotation matrix
    m[0] = t * x * x + cs
    m[1] = t * x * y + sn * z
    m[2] = t * x * z - sn * y
    m[3] = t * x * y - sn * z
    m[4] = t * y * y + cs
    m[5] = t * y * z + sn * x
    m[6] = t * x * z + sn * y
    m[7] = t * y * z - sn * x
    m[8] = t * z * z + cs
    
    // now the matrix is set - time to rotate the atoms
    tx = c[tor[1]]
    ty = c[tor[1] + 1]
    tz = c[tor[1] + 2]
    
    for atom in atoms {
        let j = (atom - 1) * 3
        c[j] -= tx
        c[j + 1] -= ty
        c[j + 2] -= tz
        x = c[j] * m[0] + c[j + 1] * m[1] + c[j + 2] * m[2]
        y = c[j] * m[3] + c[j + 1] * m[4] + c[j + 2] * m[5]
        z = c[j] * m[6] + c[j + 1] * m[7] + c[j + 2] * m[8]
        c[j] = x
        c[j + 1] = y
        c[j + 2] = z
        c[j] += tx
        c[j + 1] += ty
        c[j + 2] += tz
    }
}
