//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 3/7/23.
//

import Foundation
import Bitset


//! \class OBRotorRule rotor.h <openbabel/rotor.h>
//! \brief A rule for torsional conformer searching, defined by a SMARTS pattern
//!
//! Rules define a SMARTS pattern to match and a set of 4 reference atoms
//! defining the dihedral angle. The rule can either define a set of possible
//! dihedral angles in degrees and/or a "delta" (i.e., all multiples of delta will
//! be considered)

private let MK_DEFAULT_DELTA = 15.0

class MKRotorRule {
    private var _ref: [Int] = [Int]()                        //!< Reference atoms specifying the dihedral angle (as integers), numbered from 1 inside the SMARTS pattern
    private var _delta: Double = 0.0                         //!< (optional) the resolution of a dihedral step in degrees
    private var _s: String = ""                              //!< Text of the SMARTS pattern
    private var _sp: MKSmartsPattern                         //!< The SMARTS pattern for the rotation rule
    private var _vals: [Double] = [Double]()                 //!< At least one torsion angle (in radians) to evaluate

    init(_ buffer: String, _ ref: [Int], _ vals: [Double], _ delta: Double) {
        _delta = delta
        _s = buffer
        _vals = vals
        _ref = ref
        _sp = MKSmartsPattern()
        _sp.initialize(_s)
    }

    //! \return whether this rotor rule is valid (i.e., is the SMARTS pattern valid)
    func isValid() -> Bool {
        return(_sp.isValid())
    }
    //! \return a copy of the reference atom indexes inside the SMARTS pattern
    func getReferenceAtoms(_ ref: inout [Int]) {
        ref = _ref
    }
    //! Set the resolution (delta) of a torsional step in degrees
    func setDelta(_ d: Double) {
        _delta = d
    }
    //! \return the resolution (delta) of a torsional step in degrees
    func getDelta() -> Double {
        return(_delta)
    }
    //! \return a reference to the dihedral angles to evaluate (in radians)
    func getTorsionVals() -> [Double] {
        return(_vals)
    }
    //! \return the text of the SMARTS pattern for this rule
    func getSmartsString() -> String {
        return(_s)
    }
    //! \return the exact OBSmartsPattern object for this rule
    func getSmartsPattern() -> MKSmartsPattern {
        return(_sp)
    }
}

//! \class OBRotorRules rotor.h <openbabel/rotor.h>
//! \brief Database of default hybridization torsional rules and SMARTS-defined OBRotorRule objects
//!
//! Use to automatically evaluate potentially rotatable bonds to generate
//! lists of dihedral angles to consider.
//! e.g., rotamer/conformer energy calculations
class MKRotorRules: MKGlobalDataBase {

    var _quiet: Bool = false                 //!< Control debugging output from GetRotorIncrements()
    var _vr: [MKRotorRule] = [MKRotorRule]() //!< Database of specific OBRotorRules defined by SMARTS patterns
    var _sp3sp3: [Double] = [Double]()       //!< Default dihedral angles to check for generic sp3 - sp3 hybridized rotatable bonds (in radians)
    var _sp3sp2: [Double] = [Double]()       //!< Default dihedral angles to check for generic sp3 - sp2 hybridized rotatable bonds (in radians)
    var _sp2sp2: [Double] = [Double]()       //!< Default dihedral angles to check for generic sp2 - sp2 hybridized rotatable bonds (in radians)

    init() {
        super.init(fileName: "torlib", subDir: "Data")
        self.readFile()
    }

    /// Reinitializes the class
    func initialize() {
        _vr.removeAll()
        _sp3sp3.removeAll()
        _sp3sp2.removeAll()
        _sp2sp2.removeAll()
        self.readFile()
    }
    
    override func readFile() {
        guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }
        
        filePath.foreachRow { rowContent, lineNum in
            if !rowContent.starts(with: "#") {
                let vs = rowContent.components(separatedBy: .whitespaces)
                
                if rowContent.starts(with: "SP3-SP3") {
                    _sp3sp3.removeAll()
                    for j in vs[1...] { // skip the tag
                        if let val = Double(j) {
                            _sp3sp3.append(val.degreesToRadians)
                        }
                    }
                } else if rowContent.starts(with: "SP3-SP2") {
                    _sp3sp2.removeAll()
                    for j in vs[1...] { // skip the tag
                        if let val = Double(j) {
                            _sp3sp2.append(val.degreesToRadians)
                        }
                    }
                } else if rowContent.starts(with: "SP2-SP2") {
                    _sp2sp2.removeAll()
                    for j in vs[1...] { // skip the tag
                        if let val = Double(j) {
                            _sp2sp2.append(val.degreesToRadians)
                        }
                    }
                } else if vs.count > 5 {
                    var vals: [Double] = []
                    var delta: Double = MK_DEFAULT_DELTA
                    var refs: [Int] = [Int].init(repeating: 0, count: 4)
                    let temp_buffer = vs[0]
                    //reference atoms
                    for i in 0..<4 {
                        guard let ref_i = Int(vs[i+1]) else { fatalError("could not read Int from \(vs[i+1])") }
                        refs[i] = ref_i - 1
                    }
                    //possible torsions
                    var i: Int = 5
                    repeat {
                        if (i == vs.count - 2) && vs[i] == "Delta" {
                            guard let newDelta = Double(vs[i+1]) else { fatalError("Could not read Double from \(vs[i+1])") }
                            delta = newDelta
                            i += 2
                        } else {
                            guard let newDelta = Double(vs[i]) else { fatalError("Could not read Double from \(vs[i])") }
                            vals.append(newDelta.degreesToRadians)
                        }
                        i += 1
                    } while i < vs.count
                    
                    if vals.isEmpty {
                        var err = "The following rule has no associated torsions: "
                        err += vs[0]
                        MKLogger.throwError(#function, errorMsg: err)
                    }
                    
                    let rr: MKRotorRule = MKRotorRule(temp_buffer, refs, vals, delta)
                    if rr.isValid() {
                        _vr.append(rr)
                    }
                }
            }
        }
    }
    
    override func getSize() -> Int {
        return _vr.count
    }
    
    func setFilename(_ s: String) {
        _filename = s
    }

    //! Determine the torsional angles to evaluate based on the database
    //! \param mol molecule to evaluate
    //! \param bond rotatable bond to evaluate
    //! \param refs set to be the atom indexes (in mol) of the dihedral angle
    //! \param vals set to be the list of angles to evaluate (in radians)
    //! \param delta potential dihedral angle steps (in degrees)
    func getRotorIncrements(_ mol: MKMol, _ bond: MKBond, refs: inout [Int], vals: inout [Double], delta: inout Double) {
        if !_init {
            initialize()
        }
        
        vals.removeAll()
        var vpr: [Pair<Int, Int>] = []
        vpr.append(Pair<Int, Int>(0, bond.getBeginAtomIdx()))
        vpr.append(Pair<Int, Int>(0, bond.getEndAtomIdx()))
        
        delta = MK_DEFAULT_DELTA
        
        var sp: MKSmartsPattern
        var map: [[Int]] = []
        
        for i in _vr {
            sp = i.getSmartsPattern()
            i.getReferenceAtoms(&refs)
            vpr[0].0 = refs[1]
            vpr[1].0 = refs[2]
            
            if !sp.restrictedMatch(mol, vpr, true) {
                // perform swap with local variables
                var a = vpr[0].0
                var b = vpr[1].0
                a <-> b
                vpr[0].0 = a
                vpr[1].0 = b
                if !sp.restrictedMatch(mol, vpr, true) {
                    continue
                }
            }
            
            map = sp.getUMapList()
            for j in 0..<4 {
                refs[j] = map[0][refs[j]]
            }
            
            vals = i.getTorsionVals()
            delta = i.getDelta()
            
            guard var a1: MKAtom = mol.getAtom(refs[0]),
                  var a4: MKAtom = mol.getAtom(refs[3]) else { fatalError("Cannot get end atoms") }
            if a1.getAtomicNum() == MKElements.Hydrogen.atomicNum && a4.getAtomicNum() == MKElements.Hydrogen.atomicNum {
                continue // don't allow hydrogens at both ends
            }
            
            if a1.getAtomicNum() == MKElements.Hydrogen.atomicNum || a4.getAtomicNum() == MKElements.Hydrogen.atomicNum { //need a heavy atom reference - can use hydrogen
                
                var swapped: Bool = false
                guard var a2 = mol.getAtom(refs[1]),
                      var a3 = mol.getAtom(refs[2]) else { fatalError("Cannot get inside atoms") }
                
                if a4.getAtomicNum() == MKElements.Hydrogen.atomicNum {
                    swap(&a1, &a4)
                    swap(&a2, &a3)
                    swapped = true
                }
                
                var r: MKAtom? = nil
                guard let a2Nbrs = a2.getNbrAtomIterator() else { fatalError() }
                
                for a2Nbr in a2Nbrs {
                    if a2Nbr.getAtomicNum() != MKElements.Hydrogen.atomicNum && a2Nbr != a3 {
                        r = a2Nbr
                        break
                    }
                }
                
                if r == nil {
                    continue // unable to find reference heavy atom
                }
                
                let t1 = mol.getTorsion(a1, a2, a3, a4)
                let t2 = mol.getTorsion(r!, a2, a3, a4)
                var diff = t2 - t1
                
                if diff > 180.0 {
                    diff -= 360.0
                }
                
                if diff < -180.0 {
                    diff += 360.0
                }
                
                diff = diff.degreesToRadians
                
                for var m in vals {
                    m += diff
                    if m < Double.pi {
                        m += 2*Double.pi
                    }
                    if m > Double.pi {
                        m -= 2*Double.pi
                    }
                }
                
                if swapped {
                    refs[3] = r!.getIdx()
                } else {
                    refs[0] = r!.getIdx()
                }
            }
            
            if !_quiet {
                let errorMsg = "".appendingFormat("%3d%3d%3d%3d %s", refs[0], refs[1], refs[2], refs[3], i.getSmartsString())
                MKLogger.throwError(#function, errorMsg: errorMsg)
            }
            return
        }
        
        //***didn't match any rules - assign based on hybridization***
        let a2 = bond.getBeginAtom()
        let a3 = bond.getEndAtom()
        
        guard let a2Nbrs = a2.getNbrAtomIterator(),
              let a3Nbrs = a3.getNbrAtomIterator() else { fatalError("Cannot get nbrs for a2, a3 atoms") }
        
        guard let a1 = a2Nbrs.first(where: { $0.getAtomicNum() != MKElements.Hydrogen.atomicNum && $0 != a3 }) else { fatalError("cannot get a1 anchor") }
        guard let a4 = a3Nbrs.first(where: { $0.getAtomicNum() != MKElements.Hydrogen.atomicNum && $0 != a2 }) else { fatalError("cannot get a4 anchor") }
        
        refs[0] = a1.getIdx()
        refs[1] = a2.getIdx()
        refs[2] = a3.getIdx()
        refs[3] = a4.getIdx()
        
        if a2.getHyb() == 3 && a3.getHyb() == 3 { // sp3-sp3
            vals = _sp3sp3
            if !_quiet {
                let errorMsg = "".appendingFormat("%3d%3d%3d%3d %s", refs[0], refs[1], refs[2], refs[3], "sp3-sp3")
                MKLogger.throwError(#function, errorMsg: errorMsg)
            }
        } else {
            if a2.getHyb() == 2 && a3.getHyb() == 2 { // sp2-sp2
                vals = _sp2sp2
                if !_quiet {
                    let errorMsg = "".appendingFormat("%3d%3d%3d%3d %s", refs[0], refs[1], refs[2], refs[3], "sp2-sp2")
                    MKLogger.throwError(#function, errorMsg: errorMsg)
                }
            } else { // must be sp2-sp3
                vals = _sp3sp2
                if !_quiet {
                    let errorMsg = "".appendingFormat("%3d%3d%3d%3d %s", refs[0], refs[1], refs[2], refs[3], "sp3-sp2")
                    MKLogger.throwError(#function, errorMsg: errorMsg)
                }
            }
        }
    }
    
    //! Turn off debugging output from GetRotorIncrements()
    func quiet() { _quiet = true }
    
}

/**
 * @class OBRotor rotor.h <openbabel/rotor.h>
 * @brief A single rotatable OBBond as part of rotamer searching
 */
class MKRotor {
    var _idx: Int = 0 //!< the index in an OBRotorList
    var _rotatoms: [Int] = [Int]() //!< the atoms to rotate
    var _imag: Double = 0.0, _refang: Double = 0.0 //!< inverse magnitude and reference angle (see Precompute())
    var _bond: MKBond? //!< the bond associated with this rotor
    var _ref: [Int] = [Int]() //!< indexes for atom coordinates (from 0, multiplied by 3)
    var _torsion: [Int] = [Int]() //!< indexes for atom coordinates (from 0, multiplied by 3)
    var _fixedatoms: Bitset = Bitset() //!< fixed atoms/bonds
    var _fixedbonds: Bitset = Bitset() //!< fixed atoms/bonds
    var _evalatoms: Bitset = Bitset() //!< fixed atoms/bonds
    var _torsionAngles: [Double] = [Double]()  //!< torsion resolution
    var _invmag: [Double] = [Double]() //!< the inverse magnitudes (see Precalc)
    var _sn: [[Double]] = [[Double]]() //!< the rotation matrix (see Precalc())
    var _cs: [[Double]] = [[Double]]() //!< the rotation matrix (see Precalc())
    var _t: [[Double]] = [[Double]]() //!< the rotation matrix (see Precalc())
    var _rings: [MKRing] = [MKRing]() //!< the parent ring (if this is a rotor in a ring)
    
    
    init() {}

    ///@name Setup
    ///@{
    /**
     * Set the OBBond associated with this OBRotor.
     */
    func setBond(_ bond: MKBond) {
        _bond = bond
        setRings()
    }
    /**
     * Set the rings associated with this bond (if it's a ring bond)
     * \since Version 2.4
     */
    func setRings() {
        _rings.removeAll()
        if _bond == nil { return } // do nothing
                
        guard let mol = _bond?.getParent() else { return } // nothing to do
        
        let rlist = mol.getSSSR()
        for i in rlist {
            if i.isMember(_bond!) {
                _rings.append(i)
            }
        }
    }
    /**
     * Set the index for this rotor. Used by OBRotorList
     */
    func setIdx(_ idx: Int) {
        _idx = idx
    }
    /**
     * Set the dihedral atoms.
     * @param ref The dihedral atom indexes. These indexes start from 1.
     */
    func setDihedralAtoms(_ ref: [Int]) {
        _ref = ref
        // convert the indexes (start from 0, multiplied by 3) for easy access to coordinates
        _torsion[0] = (ref[0] - 1)*3
        _torsion[1] = (ref[1] - 1)*3
        _torsion[2] = (ref[2] - 1)*3
        _torsion[3] = (ref[3] - 1)*3
    }
    /**
     * Set the atom indexes that will be displaced when this rotor
     * changes torsion angle. These indexes start from 0 and are multiplied
     * by 3 for easy coordinate access.
     */
    func setRotAtoms(_ atoms: [Int]) {
        _rotatoms = atoms
    }
    /**
     * Set the possible torsion values or angles.
     */
    func setTorsionValues(_ angles: [Double]) {
        _torsionAngles = angles
    }
    /**
     * Set the bonds that will be fixed.
     */
    func setFixedBonds(_ bv: Bitset) {
        _fixedbonds = bv
    }
    
    func setEvalAtoms(_ bv: Bitset) {
        _evalatoms = bv
    }
    
    ///@}
    
    ///@name Performing rotations
    ///@{
    /**
     * Rotate the atoms in the specified @p coordinates to the specified angle.
     * @param coordinates The coordinates to rotate.
     * @param setang The new torsion angle in radians.
     */
    func setToAngle(_ coordinates: inout [Double], _ setang: Double) {
        var sn: Double = 0.0, cs: Double = 0.0, t: Double = 0.0, ang: Double = 0.0, mag: Double = 0.0
        // compute the angle to rotate (radians)
        ang = setang - calcTorsion(coordinates)
        // if the angle to rotate is too small, we're done
        if (fabs(ang) < 1e-5) {
            return
        }
        // compute the bond length
        mag = calcBondLength(coordinates)
        // compute some rotation matrix elements
        sn = sin(ang)
        cs = cos(ang)
        t = 1 - cs
        // perform rotation
        set(&coordinates, sn, cs, t, 1.0 / mag)
    }

    /**
     * Rotate the atoms in the specified @p coordinates. This function does not
     * require any precomputation and will compute all needed information when
     * needed.
     * @param coordinates The coordinates for the molecules as pointer to double.
     * @param next The index of the new rotor angle. This is an index for the
     * GetTorsionValues() list.
     * @param prev If specified, the torsion current torsion angle can be
     * looked up and does not have to be calculated again.
     */
    func setRotor(_ c: inout [Double], _ idx: Int, _ prev: Int = -1) {
        var ang, sn, cs, t, dx, dy, dz, mag: Double

        if prev == -1 {
            ang = _torsionAngles[idx] - calcTorsion(c)
        } else {
            ang = _torsionAngles[idx] - _torsionAngles[prev]
        }

        sn = sin(ang)
        cs = cos(ang)
        t = 1 - cs

        dx = c[_torsion[1]] - c[_torsion[2]]
        dy = c[_torsion[1] + 1] - c[_torsion[2] + 1]
        dz = c[_torsion[1] + 2] - c[_torsion[2] + 2]
        mag = sqrt(dx * dx + dy * dy + dz * dz)
        
        set(&c, sn, cs, t, 1.0 / mag)
    }
    /**
     * Rotate the specified @p coordinates by using the specified rotation matrix.
     *
     */
    func set(_ c: inout [Double], _ sn: Double, _ cs: Double, _ t: Double, _ invmag: Double) {
        var x, y, z, tx, ty, tz: Double
        var m: [Double] = Array(repeating: 0.0, count: 9)
        
        x = c[_torsion[1]] - c[_torsion[2]]
        y = c[_torsion[1] + 1] - c[_torsion[2] + 1]
        z = c[_torsion[1] + 2] - c[_torsion[2] + 2]
        
        // normalize the rotation vector
        x *= invmag
        y *= invmag
        z *= invmag
        
        // set up the rotation matrix
        tx = t * x
        ty = t * y
        tz = t * z
        m[0] = tx * x + cs
        m[1] = tx * y + sn * z
        m[2] = tx * z - sn * y
        m[3] = tx * y - sn * z
        m[4] = ty * y + cs
        m[5] = ty * z + sn * x
        m[6] = tx * z + sn * y
        m[7] = ty * z - sn * x
        m[8] = tz * z + cs
        
        // now the matrix is set - time to rotate the atoms
        tx = c[_torsion[1]]
        ty = c[_torsion[1] + 1]
        tz = c[_torsion[1] + 2]
        
        for i in _rotatoms {
            let j = i
            c[j] -= tx
            c[j + 1] -= ty
            c[j + 2] -= tz
            x = c[j] * m[0] + c[j + 1] * m[1] + c[j + 2] * m[2]
            y = c[j] * m[3] + c[j + 1] * m[4] + c[j + 2] * m[5]
            z = c[j] * m[6] + c[j + 1] * m[7] + c[j + 2] * m[8]
            c[j] = x + tx
            c[j + 1] = y + ty
            c[j + 2] = z + tz
        }
    }
    /**
     * Precompute the reference angle and inverse bond length of this rotor for
     * a single conformer. This function should be used in combination with
     * Set(double *coordinates, int idx).
     * @param coordinates The coordinates to use in the computation.
     *
     * @code
     * OBMol mol;
     * ...
     *
     * unsigned int numCoords = mol.NumAtoms() * 3;
     * double *coords = mol.GetCoordinates();
     * OBRotor rotor;
     * rotor.SetBond(mol.GetBond(3));
     *
     * // set the possible torsion values
     * std::vector<double> angles;
     * angles.push_back(0.0);
     * angles.push_back(3.1415);
     * rotor.SetTorsionValues(angles);
     *
     * // precompute inverse bond length (i.e. the bond length of bond with index 3
     * // using the specified coordinates) and reference angle (i.e. torsion angle
     * //in coords)
     * rotor.Precompute(coords);
     *
     * // copy coordinates to coords_1
     * double *coords_1 = new double[numCoords];
     * for (unsigned int i = 0; i < numCoords; ++i)
     *   coords_1[i] = coords[i];
     * // rotate the atoms in coords_1 to angle with index 0 (i.e. 0.0 degrees)
     * // note: on input, the coordinates should be the same as the coordinates used
     * //       to precompute the inverse bond length and reference angle (in other
     * //       words, the inverse magnitude and reference angle in the specfied
     * //       coordinates should be the same as the one used for Precompute)
     * rotor.Set(coords_1, 0)
     *
     * // copy coordinates to coords_2
     * double *coords_2 = new double[numCoords];
     * for (unsigned int i = 0; i < numCoords; ++i)
     *   coords_2[i] = coords[i];
     * // rotate the atoms in coords_2 to angle with index 1 (i.e. 180.0 degrees)
     * rotor.Set(coords_2, 1)
     *
     * delete coords_1;
     * delete coords_2;
     * @endcode
     */
    func precompute(_ coordinates: [Double]) {
        _imag = 1.0 / calcBondLength(coordinates)
        _refang = calcTorsion(coordinates)
    }
    /**
     * Rotate the @p coordinates to set the torsion angle of this rotor to the angle
     * specified by the index @p idx. Make sure to call Precompute before calling
     * this function.
     * @param coordinates The coordinates to rotate.
     * @param idx The index of the torsion angle in the GetTorsionValues() list.
     */
    func set(_ coordinates: inout [Double], _ idx: Int) {
        var ang, sn, cs, t: Double
        // compute the rotation angle
        ang = _torsionAngles[idx] - _refang
        // compute some values for the rotation matrix
        sn = sin(ang)
        cs = cos(ang)
        t = 1 - cs
        
        var x, y, z, tx, ty, tz: Double
        var m: [Double] = Array(repeating: 0.0, count: 9)
        
        // compute the bond vector
        x = coordinates[_torsion[1]] - coordinates[_torsion[2]]
        y = coordinates[_torsion[1] + 1] - coordinates[_torsion[2] + 1]
        z = coordinates[_torsion[1] + 2] - coordinates[_torsion[2] + 2]
        
        // normalize the bond vector
        x *= _imag
        y *= _imag
        z *= _imag
        
        // set up the rotation matrix
        tx = t * x
        ty = t * y
        tz = t * z
        m[0] = tx * x + cs
        m[1] = tx * y + sn * z
        m[2] = tx * z - sn * y
        m[3] = tx * y - sn * z
        m[4] = ty * y + cs
        m[5] = ty * z + sn * x
        m[6] = tx * z + sn * y
        m[7] = ty * z - sn * x
        m[8] = tz * z + cs
        
        // now the matrix is set - time to rotate the atoms
        tx = coordinates[_torsion[1]]
        ty = coordinates[_torsion[1] + 1]
        tz = coordinates[_torsion[1] + 2]
        
        for i in _rotatoms {
            let j = i
            coordinates[j] -= tx
            coordinates[j + 1] -= ty
            coordinates[j + 2] -= tz
            x = coordinates[j] * m[0] + coordinates[j + 1] * m[1] + coordinates[j + 2] * m[2]
            y = coordinates[j] * m[3] + coordinates[j + 1] * m[4] + coordinates[j + 2] * m[5]
            z = coordinates[j] * m[6] + coordinates[j + 1] * m[7] + coordinates[j + 2] * m[8]
            coordinates[j] = x + tx
            coordinates[j + 1] = y + ty
            coordinates[j + 2] = z + tz
        }
    }
    /**
     * Precompute the inverse bond lengths, rotation matrices for all
     * specified conformers and all possible torsion values. This method is
     * used in combination with Set(double *coordinates, int conformer, int idx).
     * @param conformers The pointers to the conformer coordinates
     */
    func precalc(_ cv: inout [[Double]]) {
        for c in cv {
            var cs: [Double] = []
            var sn: [Double] = []
            var t: [Double] = []
            let ang = calcTorsion(c)
            
            for angle in _torsionAngles {
                cs.append(cos(angle - ang))
                sn.append(sin(angle - ang))
                t.append(1 - cos(angle - ang))
            }
            
            _cs.append(cs)
            _sn.append(sn)
            _t.append(t)
            _invmag.append(1.0 / calcBondLength(c))
        }
    }
    /**
     * Rotate the @p coordinates to set the torsion to the torsion value with the
     * specified @p index. The coordinates should be the same as the conformer used
     * for calling Precalc (i.e. conformers[conformer] == coordinates). Make sure
     * to call Precalc before calling this method.
     * @param coordinates The conformer coordinates.
     * @param conformer The conformer index in the conformer list given to Precalc().
     * @param idx The torsion value index in the GetTorsionValues() list.
     */
    func set(_ coordinates: inout [Double], _ conformer: Int, _ idx: Int) {
        set(&coordinates, _sn[conformer][idx], _cs[conformer][idx], _t[conformer][idx], _invmag[conformer])
    }

    ///@name Methods to retrieve information
    ///@{
    /**
     * Get the OBBond object associated with this OBRotor.
     */
    func getBond() -> MKBond? {
        return _bond
    }
    /**
     * Get the number of possible torsion angles for this OBRotor. This
     * is the length of the GetTorsionValues() list.
     */
    func size() -> Int {
        return _torsionAngles.count
    }
    /**
     * Get the index for this rotor (index in an OBRotorList).
     */
    func getIdx() -> Int {
        return _idx
    }
    /**
     * Get the dihedral atom indexes. These indexes start from 1.
     */
    func getDihedralAtoms(_ ref: inout [Int]){
        for i in 0..<4 {
            ref[i] = _ref[i]
        }
    }
    /**
     * Get the dihedral atom indexes. These indexes start from 1.
     */
    func getDihedralAtoms() -> [Int] {
        return _ref
    }
    /**
     * Get the atom indexes that will be displaced when this rotor changes
     * torsion angle. These indexes start from 1.
     */
    func getRotAtoms() -> [Int] {
        return _rotatoms
    }
    /**
     * Get the possible torsion angles for this OBRotor.
     */
    func getTorsionValues() -> [Double] {
        return _torsionAngles
    }
    /**
     * Get an OBBitVec objects with bits set for all bonds that are fixed.
     * Bonds are indexed from 0.
     */
    func getFixedBonds() -> Bitset {
        return _fixedbonds
    }
    /**
     * Calculate the torsion for this OBRotor using the specified coordinates.
     * @param coordinates The coordinates (e.g. OBMol::GetCoordinates()).
     * @return The torsion angle in radians.
     */
    func calcTorsion(_ c: [Double]) -> Double {
        var v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z: Double
        var c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z: Double
        var c1mag, c2mag, ang, costheta: Double

        // calculate the torsion angle
        v1x = c[_torsion[0]] - c[_torsion[1]]
        v1y = c[_torsion[0] + 1] - c[_torsion[1] + 1]
        v1z = c[_torsion[0] + 2] - c[_torsion[1] + 2]
        v2x = c[_torsion[1]] - c[_torsion[2]]
        v2y = c[_torsion[1] + 1] - c[_torsion[2] + 1]
        v2z = c[_torsion[1] + 2] - c[_torsion[2] + 2]
        v3x = c[_torsion[2]] - c[_torsion[3]]
        v3y = c[_torsion[2] + 1] - c[_torsion[3] + 1]
        v3z = c[_torsion[2] + 2] - c[_torsion[3] + 2]

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
            costheta = 1.0 // avoid div by zero error
        } else {
            costheta = (c1x * c2x + c1y * c2y + c1z * c2z) / (sqrt(c1mag * c2mag))
        }

        if costheta < -0.9999999 {
            costheta = -0.9999999
        }
        if costheta > 0.9999999 {
            costheta = 0.9999999
        }

        if (v2x * c3x + v2y * c3y + v2z * c3z) > 0.0 {
            ang = -acos(costheta)
        } else {
            ang = acos(costheta)
        }

        return ang
    }
    /**
     * Calculate the bond length for this OBRotor using the specified coordinates.
     * @param coordinates The coordinates (e.g. OBMol::GetCoordinates()).
     */
    func calcBondLength(_ c: [Double]) -> Double {
        // compute the difference
        let dx = c[_torsion[1]] - c[_torsion[2]]
        let dy = c[_torsion[1] + 1] - c[_torsion[2] + 1]
        let dz = c[_torsion[1] + 2] - c[_torsion[2] + 2]
        // compute the length
        return sqrt(dx * dx + dy * dy + dz * dz)
    }
    ///@}
    
    func removeSymTorsionValues(_ fold: Int) {
        if _torsionAngles.count == 1 { return }
        var tv: [Double] = []
        for i in _torsionAngles {
            if i >= 0.0 && i < (2*Double.pi / Double(fold)) {
                tv.append(i)
            }
        }
        if tv.isEmpty { return }
        _torsionAngles = tv
    }
}

@discardableResult
func getDFFVector(_ mol: MKMol, _ dffv: inout [Int], _ bv: Bitset) -> Bool {
    dffv.removeAll()
    dffv.reserveCapacity(mol.numAtoms())
    
    var dffcount: Int, natom: Int
    var used: Bitset = Bitset()
    var curr: Bitset = Bitset()
    let next: Bitset = Bitset()
    
    for atom in mol.getAtomIterator() {
        if bv[atom.getIdx()] {
            dffv[atom.getIdx() - 1] = 0
            continue
        }
        
        dffcount = 0
        used.removeAll()
        curr.removeAll()
        used.add(atom.getIdx())
        curr.add(atom.getIdx())
        
        while !curr.isEmpty() && (bv&curr).isEmpty() {
            next.removeAll()
            for natom in curr {
                guard let atom1 = mol.getAtom(natom),
                      let bonds = atom1.getBondIterator() else { fatalError() }
                for bond in bonds {
                    if (!used.contains(bond.getNbrAtomIdx(atom1)) &&
                        !curr.contains(bond.getNbrAtomIdx(atom1))) {
                        if bond.getNbrAtom(atom1).getAtomicNum() != MKElements.Hydrogen.atomicNum {
                            next.add(bond.getNbrAtomIdx(atom1))
                        }
                    }
                }
            }
            
            used |= next
            curr = next
            dffcount += 1
        }
        
        dffv[atom.getIdx() - 1] = dffcount
    }
    
    return true
}

func compareRotor(_ a: Pair<MKBond,Int>,_ b: Pair<MKBond,Int>) -> Bool {
//    return a.1 > a.1 // outside->in
    return a.1 < b.1 // inside->out
}

 /**
   * @class OBRotorList rotor.h <openbabel/rotor.h>
   * @brief Given an OBMol, set up a list of possibly rotatable torsions,
   */
class MKRotorList {
    private var _quiet: Bool = false                         // Control debugging output
    private var _removesym: Bool = false                     // Control removal of symmetric rotations
    private var _ringRotors: Bool = false                    // Are there ring rotors
    private var _fixedatoms: Bitset = Bitset()           // Bit vector of fixed (i.e., invariant) atoms
    private var _fixedbonds: Bitset = Bitset()           // Bit vector of fixed (i.e., invariant) atoms
    private var _rr: MKRotorRules = MKRotorRules()           // Database of rotatable bonds and dihedral angles to test
    private var _dffv: [Int] = []                            // Distance from fixed
    private var _rotor: [MKRotor] = []                       // List of individual OBRotor torsions
    private var _vsym2: [(MKSmartsPattern, (Int, Int))] = [] // !
    private var _vsym3: [(MKSmartsPattern, (Int, Int))] = [] // !


    init() {
        _rotor.removeAll()
        _quiet = true
        _removesym = true
        _ringRotors = false
    }

    /**
     * Clear the internal list of rotors and reset.
     */
    func clear() {
        _rotor.removeAll()
        _ringRotors = false
    }

    /**
     * @return the number of rotors in this list
     */
    func size() -> Int {
        return _rotor.count
    }

    /**
     * When no atoms/bonds are fixed or when bonds are fixed, this function will
     * return true if the bond is fixed. When using the deprecated fixed atoms,
     * this function will return true if the bond and at least one neighboring
     * bond has fixed atoms.
     */
    func isFixedBond(_ bond: MKBond) -> Bool {
        if _fixedatoms.isEmpty() && _fixedbonds.isEmpty() {
            return false
        }

        // new fixed bonds
        if !_fixedbonds.isEmpty() {
            return _fixedbonds.contains(Int(bond.getIdx()))
        }

        if _fixedatoms.isEmpty() {
            return false
        }

        let a1 = bond.getBeginAtom()
        let a2 = bond.getEndAtom()
        if !_fixedatoms[a1.getIdx()] || !_fixedatoms[a2.getIdx()] {
            return false
        }
        var isfixed: Bool = false 
        for a3 in a1.getNbrAtomIterator()! {
            if a3 != a2 && _fixedatoms[a3.getIdx()] {
                isfixed = true
                break
            }
        }
        if !isfixed {
            return false 
        }
        isfixed = false
        for a3 in a2.getNbrAtomIterator()! {
            if a3 != a1 && _fixedatoms[a3.getIdx()] {
                isfixed = true
                break
            }
        }
        return isfixed
    }

    /**
     * @return True if this rotor list has any fixed bonds.
     */
    func hasFixedBonds() -> Bool {
        return !_fixedbonds.isEmpty()
    }

    /**
     * @return whether this rotor list has any fixed (invariant) atoms
     * @deprecated See HasFixedBonds()
     */
    func hasFixedAtoms() -> Bool {
        return !_fixedatoms.isEmpty()
    }
    
    //! Rotates each bond to zero and 180 degrees and tests
    //! if the 2 conformers are duplicates.  if so - the symmetric torsion
    //! values are removed from consideration during a search
    func removeSymVals(_ mol: MKMol) {
        
        let gs: MKGraphSym = MKGraphSym(mol)
        var sym_classes: [Ref] = []
        gs.getSymmetry(&sym_classes)
        
        var syms: Set<Ref> = Set<Ref>()
        
        for rotor in _rotor {
            guard let bond = rotor.getBond() else { fatalError() }
            let end = bond.getEndAtom()
            let begin = bond.getBeginAtom()
            var N_fold_symmetry: Int = 1
            for here in 0...1 { // try each side of the bond in turn
                var this_side: MKAtom, other_side: MKAtom
                
                if here == 0 {
                    this_side = begin
                    other_side = end
                } else {
                    this_side = end
                    other_side = begin
                }
                
                for hyb in 2...3 { // sp2 and sp3 carbons, with explicit Hs
                    if this_side.getAtomicNum() == 6 && this_side.getHyb() == hyb && this_side.getExplicitDegree() == (hyb + 1) {
                        syms.removeAll()
                        guard let nbrs = this_side.getNbrAtomIterator() else { fatalError() }
                        for nbr in nbrs {
                            if nbr == other_side { continue }
                            syms.insert(sym_classes[nbr.getIdx() - 1])
                        }
                        
                        if syms.count == 1 { // All of the rotated atoms have the same symmetry class
                            N_fold_symmetry *= hyb
                        }
                    }
                }
                
            }
            
            if N_fold_symmetry > 1 {
                let old_size = rotor.size()
                rotor.removeSymTorsionValues(N_fold_symmetry)
                if !_quiet {
                    var printMsg = "....\(N_fold_symmetry)-fold symmetry at rotor between \(begin.getIdx()) and \(end.getIdx())"
                    printMsg += " - reduced from \(old_size) to \(rotor.size())\n"
                    print(printMsg)
                }
            }
        }
    }

    /**
     * @return True if this rotor list has any ring bonds.
     * @since version 2.4
     */
    func hasRingRotors() -> Bool {
        return _ringRotors
    }

    ///@name Setup
    /**
     * Setup this rotor list for the supplied molecule. This method calls
     * FindRotors(), SetEvalAtoms(), and AssignTorVals().
     * @param mol The molecule.
     * @param sampleRings Whether to sample ring conformers - default = false
     * @return True if rotatable bonds were found.
     */
    func setup(_ mol: MKMol, _ sampleRingBonds: Bool = false) -> Bool {
        clear()
        // find the roatable bonds 
        findRotors(mol, sampleRingBonds)
        if size() == 0 {
            return false
        }
        // set the atoms that should be evaluated when this rotor changes
        setEvalAtoms(mol)
        assignTorVals(mol)

        for rotor in _rotor {
            if rotor.size() == 0 {
                var ref: [Int] = [0, 0, 0, 0]
                rotor.getDihedralAtoms(&ref)
                print("ERROR: The rotor has no associated torsion values -> \(ref[0]) \(ref[1]) \(ref[2]) \(ref[3])")
            }
        }
        // reduce the number of torsions to be checked through symmetry considerations
        if _removesym {
            removeSymVals(mol)
        }
        return true
    }

    /**
     * Set the bonds that will be fixed.
     */
    func setFixedBonds(_ fixedbonds: Bitset) {
        _fixedbonds = fixedbonds
        _fixedatoms.removeAll() // why is this here? 
    }

    /**
     * Initialize the private OBRotorRules database from a specific file.
     */
    func initialize(_ filename: String) {
        _rr.setFilename(filename)
        _rr.initialize()
    }
    /**
     * Turn off debugging output.
     */
    func setQuiet(_ quiet: Bool) {
        _quiet = quiet
        _rr.quiet()
    }
    
    /**
     * Find all potentially rotatable bonds in the molecule. This method uses
     * OBBond::IsRotor() for initial evaluation which depends on ring perception
     * (i.e. ring bonds are considered rotatable). Fixed bonds, specified using the
     * deprecated fixed atoms or the new fixed bonds methods are not added to
     * the list. All rotatable bonds will be sorted by their graph theoretical
     * distance (GTD) score (see OBMol::GetGTDVector()). This results in the
     * the rotors going from the inside to the outside of the mol.
     * @param mol The molecule.
     * @param sampleRingBonds whether to sample ring bonds from analysis (default = false)
     * @return True.
     */
    @discardableResult
    func findRotors(_ mol: MKMol, _ sampleRingBonds: Bool = false) -> Bool {
        // Find ring atoms & bonds
        // This function will set OBBond::IsRotor().
        mol.findRingAtomsAndBonds()

        //
        // Score the bonds using the graph theoretical distance (GTD).
        // The GTD is the distance from atom i to every other atom j.
        // Atoms on the "inside" of the molecule will have a lower GTD
        // value than atoms on the "outside"
        //
        // The scoring will rank "inside" bonds first.
        //
        var gtd: [Int] = []
        mol.getGTDVector(&gtd)
        // compute the scores
        var vtmp: [Pair<MKBond, Int>] = []
        for bond in mol.getBondIterator() {
            // check if the bond is "rotatable"
            if bond.isRotor(sampleRingBonds) {
                // check if the bond is fixed (using deprecated fixed atoms or new fixed bonds)
                if (hasFixedAtoms() || hasFixedBonds()) && isFixedBond(bond) {
                    continue
                }

                if bond.isInRing() {
                    //otherwise mark that we have them and add it to the pile
                    _ringRotors = true
                }

                let score: Int = gtd[bond.getBeginAtomIdx() - 1] + gtd[bond.getEndAtomIdx() - 1]
                // compute the GTD bond score as sum of atom GTD scores
                vtmp.append(Pair<MKBond, Int>(bond, score))
            }
        }

        // sort the rotatable bonds by GTD score
        vtmp.sort { $0.1 < $1.1 }
        var count = 0
        for (bond, _) in vtmp {
            let rotor = MKRotor()
            rotor.setBond(bond)
            rotor.setIdx(count)
//            rotor.setNumCoords(mol.numAtoms() * 3)  @DEPRECATED
            _rotor.append(rotor)
            count += 1
        }

        return true
    }
   
    //! Determines which atoms should be used to calculate the internal energy
    //! if the dihedral angle of the rotor is modified
    //! \return True
    @discardableResult
    func setEvalAtoms(_ mol: MKMol) -> Bool {
        
        var eval: Bitset = Bitset()
        var curr: Bitset = Bitset()
        let next: Bitset = Bitset()

        for rotor in _rotor {
            guard let bond = rotor.getBond() else { fatalError("Cannot get rotor bond") }
            curr.removeAll()
            eval.removeAll()
            curr.add(bond.getBeginAtomIdx())
            curr.add(bond.getEndAtomIdx())
            eval |= curr
            //follow all non-rotor bonds and add atoms to eval list
            repeat {
                next.removeAll()
                for j in curr {
                    guard let a1 = mol.getAtom(j),
                          let a1Nbrs = a1.getNbrAtomBondIterator() else { fatalError("Cannot get atom \(j)") }
                    for (a2, k) in a1Nbrs {
                        if !eval[a2.getIdx()] { // warning: maybe this needs to be subtracted by 1
                            if !k.isRotor(_ringRotors) || ((hasFixedAtoms() || hasFixedBonds()) && isFixedBond(k)) {
                                next.add(a2.getIdx())
                                eval.add(a2.getIdx())
                            }
                        }
                    }
                    
                }
                curr = next
            } while !curr.isEmpty()
            
            //add atoms alpha to eval list
            next.removeAll()
            for j in eval {
                guard let a1 = mol.getAtom(j),
                      let a1Nbrs = a1.getNbrAtomIterator() else { fatalError("Cannot get atom j") }
                for a2 in a1Nbrs {
                    next.add(a2.getIdx())
                }
            }
            
            eval |= next
            rotor.setEvalAtoms(eval)
        }
        return true
    }
    
    /**
     * Using the OBRotorRules database, set the torsion values (and delta)
     * to be evaluated and tested. This method also sets the rotatable atoms
     * for the rotors to the smallest possible set. For each bond there are
     * two candidate sets, one on either side. The smallest of these is used
     * and the torsion indexes for the rotor are inverted if needed
     * (i.e. a-b-c-d -> d-c-b-a).
     */
    @discardableResult
    func assignTorVals(_ mol: MKMol) -> Bool {
        for rotor in _rotor {
            guard let bond = rotor.getBond() else { fatalError("Cannot retrieve rotor bond") }
            // query the rotor database
            var refs: [Int] = []
            var angles: [Double] = []
            var delta: Double = 0.0
            
            _rr.getRotorIncrements(mol, bond, refs: &refs, vals: &angles, delta: &delta)
            rotor.setTorsionValues(angles)
//            rotor.setDelta() deprecated
            // Find the smallest set of atoms to rotate. There are two candidate sets,
            // one on either side of the bond. If the first tried set size plus one is
            // larger than half of the number of atoms, the other set is smaller and the
            // rotor atom indexes are inverted (.i.e. a-b-c-d -> d-c-b-a).
            var atoms: [Int] = []
            // find all atoms for which there is a path to ref[2] without goinf through ref[1]
            mol.findChildren(refs[1], refs[2], &atoms)
            
            if atoms.count + 1 > mol.numAtoms() / 2 {
                atoms.removeAll()
                // select the other smaller set
                mol.findChildren(refs[2], refs[1], &atoms)
                // invert the rotor
                refs.swapAt(0, 3)
                refs.swapAt(1, 2)
            }
            
            // translate the rotate atom indexes to coordinate indexes (i.e. from 0, multiplied by 3)
            for var j in atoms {
                j = (j - 1) * 3
            }
            // set the rotate atoms and dihedral atom indices
            rotor.setRotAtoms(atoms)
            rotor.setDihedralAtoms(refs)
        }
        return true
    }
    
    //! \brief Set the atoms to rotate from the dihedral atoms for each rotor
    //! Uses OBRotor->GetDihedralAtoms() to call OBRotor->SetRotAtoms()
    //! and standarizes the dihedral angles via OBRotor->SetDihedralAtoms()
    //! \return True
    @discardableResult
    func setRotAtoms(_ mol: MKMol) -> Bool {
        var refs: [Int] = []
        var rotatoms: [Int] = []
        for rotor in _rotor {
            rotor.getDihedralAtoms(&refs)
            
            mol.findChildren(refs[1], refs[2], &rotatoms)
            if rotatoms.count + 1 > mol.numAtoms() / 2 {
                rotatoms.removeAll()
                mol.findChildren(refs[2], refs[1], &rotatoms)
                refs.swapAt(0, 3)
                refs.swapAt(1, 2)
            }
            
            for var j in rotatoms {
                j = (j - 1) * 3
            }
            rotor.setRotAtoms(rotatoms)
            rotor.setDihedralAtoms(refs)
        }
        return true
    }

    //! \brief Set the atoms to rotate from the dihedral atoms for each rotor
    //! Insures the fixed atoms are respected, but otherwise functions like
    //! SetRotAtoms()
    func setRotAtomsByFix(_ mol: MKMol) {
        
        var refs: [Int] = []
        var rotatoms: [Int] = []
        
        getDFFVector(mol, &_dffv, _fixedatoms)
        
        for rotor in _rotor {
            rotatoms.removeAll()
            rotor.getDihedralAtoms(&refs)
            
            if _fixedatoms[refs[1]] && _fixedatoms[refs[2]] {
                refs.swapAt(0, 3)
                refs.swapAt(1, 2)
                mol.findChildren(refs[1], refs[2], &rotatoms)
                for var j in rotatoms {
                    j = (j - 1) * 3
                }
                rotor.setRotAtoms(rotatoms)
                rotor.setDihedralAtoms(refs)
            } else {
                if _dffv[refs[1] - 1] > _dffv[refs[2] - 1] {
                    refs.swapAt(0, 3)
                    refs.swapAt(1, 2)
                    mol.findChildren(refs[1], refs[2], &rotatoms)
                    for var j in rotatoms {
                        j = (j - 1) * 3
                    }
                    rotor.setRotAtoms(rotatoms)
                    rotor.setDihedralAtoms(refs)
                }
            }
        }
    }

}


//! \class OBRotorKeys rotor.h <openbabel/rotor.h>
//! \brief A class to generate all possible rotorKeys
class MKRotorKeys {
    /**
      \brief A class to generate all possible rotorKeys

      This class can generate all possible rotor keys for a set of OBRotors
      which can all have their own resolution. Thanks to Yongjin Xu for this
      patch.

      the code blow is taken from  OBForceField::SystematicRotorSearch():
      \code
      #include <openbabel/rotor.h>
      #include <openbabel/mol.h>

      // See OBConversion class to fill the mol object.
      OBMol mol;
      OBRotorList rl;
      OBRotamerList rotamers;

      rl.Setup(_mol);
      rotamers.SetBaseCoordinateSets(_mol);
      rotamers.Setup(_mol, rl);

      cout << "number of rotatable bonds: " <<  rl.Size() << endl;

      if (!rl.Size()) { // only one conformer
        cout << "generated only one conformer" << endl;
        // exit here
      }

      OBRotorKeys rotorKeys;
      OBRotorIterator ri;
      OBRotor *rotor = rl.BeginRotor(ri);
      for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) { // foreach rotor
        rotorKeys.AddRotor(rotor->GetResolution().size());
      }

      while (rotorKeys.Next()) {
        std::vector<int> rotorKey = rotorKeys.GetKey();
        cout << "rotorKey = " << rotorKey[1] << " " << rotorKey[2] << endl;
        rotamers.AddRotamer(rotorKey);
      }

      rotamers.ExpandConformerList(_mol, _mol.GetConformers());
      \endcode
      **/

    private var _vr: [rotor_digit] = []

    init() { }

    //! Clear all rotors
    func clear() {
        _vr.removeAll()
    }
    //! Number of rotor keys (= number of possible conformers)
    func numKeys() -> Int {
        var keys: Int = 0
        while next() {
            keys += 1
        }
        return keys
    }
    //! Add a rotor
    //! \param size the rotor resolution
    func addRotor(_ size: UInt) {
        _vr.append(rotor_digit(size))
    }
    //! Select the next rotor key
    //! \return true if there are more rotor keys
    func next() -> Bool {
        if _vr.count == 0 {
            return false
        }
        var carry = _vr[0].next()
        var i = 1
        while carry {
            if i == _vr.count {
                return false
            }
            carry = _vr[i].next()
            i += 1
        }
        return true
    }
    //! Get the currently selected rotor key
    //! \return current rotor key
    func getKey() -> [Int] {
        var rt: [Int] = []
        rt.append(0)
        for i in 0..<_vr.count {
            rt.append(_vr[i].get_state())
        }
        return rt
    }

}

private class rotor_digit {
    private var resolution_size: UInt = 0
    private var state: Int = 0

    init() {} 
    init(_ rs: UInt) {
        resolution_size = rs
    }

    func setSize(_ rs: UInt) {
        resolution_size = rs
    }

    func setState(_ s: Int) {
        state = s
    }

    func get_state() -> Int {
        return state
    }

    func size() -> UInt {
        return resolution_size
    }

    func next() -> Bool {
        if state < Int(resolution_size - 1) {
            state += 1
            return false
        } else {
            state = 0
        }
        return true
    }
}
