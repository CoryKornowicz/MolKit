//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation
import Surge
import simd

public let OB_SSSR_MOL: UInt = 1<<1
public let OB_RINGFLAGS_MOL: UInt = 1<<2
public let OB_AROMATIC_MOL: UInt = 1<<3
public let OB_ATOMTYPES_MOL: UInt = 1<<4
public let OB_CHIRALITY_MOL: UInt = 1<<5
public let OB_PCHARGE_MOL: UInt = 1<<6
public let OB_HYBRID_MOL: UInt = 1<<8
public let OB_CLOSURE_MOL: UInt = 1<<11
public let OB_H_ADDED_MOL: UInt = 1<<12
public let OB_PH_CORRECTED_MOL: UInt = 1<<13
public let OB_CHAINS_MOL: UInt = 1<<15
public let OB_TCHARGE_MOL: UInt = 1<<16
public let OB_TSPIN_MOL: UInt = 1<<17
public let OB_RINGTYPES_MOL: UInt = 1<<18
public let OB_PATTERN_STRUCTURE: UInt = 1<<19
public let OB_LSSR_MOL: UInt = 1<<20
public let OB_ATOMSPIN_MOL: UInt = 1<<21
public let OB_REACTION_MOL: UInt = 1<<22
public let OB_PERIODIC_MOL: UInt = 1<<23

// Flags 24-32 are unspecified 

public let OB_CURRENT_CONFORMER	= -1

enum HydrogenType: Int {
    case AllHydrogen = 0
    case PolarHydrogen = 1
    case NonPolarHydrogen = 2
}

let NumElements = 118 + 2
let alphabetical: [Int] = [ 89, 47, 13, 95, 18, 33, 85, 79, 5, 56, 4, 107, 83, 97, 35, 6, 20, 48,
                            58, 98, 17, 96, 112, 27, 24, 55, 29, NumElements-1,
                            105, 110, 66, 68, 99, 63, 9, 26, 114, 100, 87, 31,
                            64, 32, 1, 2, 72, 80, 67, 108, 53, 49, 77, 19, 36, 57, 3, 103, 71, 116, 115, 101,
                            12, 25, 42, 109, 7, 11, 41, 60, 10, 113, 28, 102, 93, 8, 118, 76, 15, 91, 82, 46,
                            61, 84, 59, 78, 94, 88, 37, 75, 104, 111, 45, 86, 44, 16, 51, 21, 34, 106, 14,
                            62, 50, 38, NumElements, 73, 65, 43, 52, 90, 22, 81, 69, 117, 92, 23, 74, 54, 39, 70,
                            30, 40 ]

class MKMol: MKBase {
    
    private var _autoPartialCharge: Bool = false
    private var _autoFormalCharge: Bool = false

    private var _title: String = ""
    private var _vatom: [MKAtom] = []
    private var _vatomIds: [MKAtom] = []
    private var _vbond: [MKBond] = []
    private var _vbondIds: [MKBond] = []
    private var _dimension: UInt16 = 0
    private var _totalCharge: Int = 0
    private var _totalSpin: UInt = 0
    private var _c: Array<Double> = []
    private var _vconf: [[[Double]]] = [] // MARK: Not sure how I feel about this, would like to simplify into one unifying structure
    private var _energy: Double = 0.0
    private var _residue: [MKResidue] = []
    private var _internals: [MKInternalCoord]? = []
    private var _mod: UInt16 = 0 


    private var _natoms: Int {
        return _vatom.count
    }
    
    private var _nbonds: Int {
        return _vbond.count
    }

    public override init() {
        super.init()
        _mod = 0;
        _totalCharge = 0
        _dimension = 3
        _vatom = []
        _vatomIds = []
        _vbond = []
        _vbondIds = []
        _title = ""
        _c = []
        _vconf = []
        _autoPartialCharge = true
        _autoFormalCharge = true
        _energy = 0.0
    }
    
    public override init(_ base: MKBase) {
        super.init(base)
        _mod = 0;
        _totalCharge = 0
        _dimension = 3
        _vatom = []
        _vatomIds = []
        _vbond = []
        _vbondIds = []
        _title = ""
        _c = []
        _vconf = []
        _autoPartialCharge = true
        _autoFormalCharge = true
        _energy = 0.0
    }
    
    public func reserveAtoms(_ n: Int) {
        if n > 0 && (self._mod != 0) {
            self._vatom.reserveCapacity(n)
            self._vatomIds.reserveCapacity(n)
        }
    }
    
    func getTitle() -> String {
        return self._title
    }

    func setTitle(_ title: String) {
        self._title = title
    }

    //! Returns a pointer to the atom after a safety check
    //! 0 < idx <= NumAtoms
    func getAtom(_ idx: Int) -> MKAtom? {
        if idx > self._natoms || idx < 1 { return nil }
        return self._vatom[idx - 1]
    }
    
    func getAtomById(_ id: Int) -> MKAtom? {
        if id >= self._vatomIds.count { return nil }
        return self._vatomIds[id]
    }

    func getFirstAtom() -> MKAtom? {
        if self._natoms == 0 { return nil }
        return self._vatom[0]
    }

    func getAllAtoms() -> [MKAtom] {
        return self._vatom
    }

    func getAtomIterator() -> MKIterator<MKAtom> {
        return MKIterator<MKAtom>(self._vatom)
    }
    
    func getBondIterator() -> MKIterator<MKBond> {
        return MKIterator<MKBond>(self._vbond)
    }

    //! Returns a pointer to the bond after a safety check
    //! 0 <= idx < NumBonds
    func getBond(_ idx: Int) -> MKBond? {
        if idx >= self._nbonds || idx < 0 { return nil }
        return self._vbond[idx]
    }

    func getBondById(_ id: Int) -> MKBond? {
        if id >= self._vbondIds.count { return nil }
        return self._vbondIds[id]
    }

    func getBond(_ bgn: MKAtom, _ end: MKAtom) -> MKBond? {
        if bgn == end { return nil }
        guard let nbrBonds = bgn.getBondIterator() else { return nil }
        for bond in nbrBonds {
            if bond.getNbrAtom(bgn) == end {
                return bond
            }
        }
        return nil
    }

    func getBond(_ bgn: Int, _ end: Int) -> MKBond? {
        if bgn == end { return nil }
        guard let bgnAtom = self.getAtom(bgn) else { return nil }
        guard let endAtom = self.getAtom(end) else { return nil }
        return self.getBond(bgnAtom, endAtom)
    }

    func numHeavyAtoms() -> Int {
        var count: Int = 0
        for atom in self.getAtomIterator() {
            if atom.getAtomicNum() != 1 {
                count += 1
            }
        }
        return count
    }

    func getResidue(_ idx: Int) -> MKResidue? {
        if idx >= self._residue.count || idx < 0 { return nil }
        return self._residue[idx]
    }
    
    func getInternalCoord() -> [MKInternalCoord]? {
        return nil
    }

    //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#findSmallestSetOfSmallestRings">blue-obelisk:findSmallestSetOfSmallestRings</a>.
    func getSSSR() -> [MKRing] {
        if !self.hasSSSRPerceived() {
            self.findSSSR()
        }

        var ringData = MKRingData() 
        if !self.hasData("SSSR") {
            ringData.setAttribute("SSSR")
            self.setData(ringData)
        }

        ringData = self.getData("SSSR") as! MKRingData
        ringData.setOrigin(.perceived)
        return ringData.getData()
    }

    func getLSSR() -> [MKRing] {
        if !self.hasLSSRPerceived() {
            self.findLSSR()
        }

        var ringData = MKRingData() 
        if !self.hasData("LSSR") {
            ringData.setAttribute("LSSR")
            self.setData(ringData)
        }

        ringData = self.getData("LSSR") as! MKRingData
        ringData.setOrigin(.perceived)
        return ringData.getData()
    }


    // TODO: Merge these two into one base function with a boolean parameter
    func getMolWt(_ implicitH: Bool) -> Double {
        var molWt: Double = 0.0
        let hMass = MKElements.getMass(1)

        for atom in self.getAtomIterator() {
            molWt += atom.getAtomicMass()
            if implicitH {
                molWt += Double(atom.getImplicitHCount()) * hMass
            }
        }
        return molWt
    }

    func getExactMass(_ implicitH: Bool) -> Double {
        var molWt: Double = 0.0
        let hMass = MKElements.getMass(1)

        for atom in self.getAtomIterator() {
            molWt += atom.getExactMass()
            if implicitH {
                molWt += Double(atom.getImplicitHCount()) * hMass
            }
        }
        return molWt
    }

    //! Stochoimetric formula in spaced format e.g. C 4 H 6 O 1
    //! No pair data is stored. Normally use without parameters: GetSpacedFormula()
    //! \since version 2.1
    func getSpacedFormula(_ ones: Int, _ sp: String, _ implicitH: Bool) -> String {
        //Default ones=0, sp=" ".
        //Using ones=1 and sp="" will give unspaced formula (and no pair data entry)
        // These are the atomic numbers of the elements in alphabetical order, plus
        // pseudo atomic numbers for D, T isotopes.
        var atomicCount = [Int](repeating: 0, count: NumElements)
        var formula = ""

        var useImplicitH = (self.numBonds() != 0 || self.numAtoms() == 1)
        // Do not use implicit hydrogens if explicitly required not to
        if !implicitH {
            useImplicitH = false
        }

        let hasHvyAtoms = self.numHeavyAtoms() > 0

        for atom: MKAtom in self.getAtomIterator() {
            var anum: Int = atom.getAtomicNum()
            if anum == 0 {
                continue
            }

            if anum > (NumElements - 2) {
                print("Skipping unknown element with atomic number \(anum)")
                continue
            }

            let isHiso: Bool = anum == 1 && atom.getIsotope() >= 2

            if useImplicitH {
                if anum == 1 && !isHiso && hasHvyAtoms {
                    continue // skip explicit hydrogens except D,T
                }

                if anum == 1 {
                    if isHiso && hasHvyAtoms {
                        atomicCount[0] -= 1 //one of the implicit hydrogens is now explicit
                    }
                } else {
                    atomicCount[0] += Int(atom.getImplicitHCount() + atom.explicitHydrogenCount())
                }
            }
            
            if isHiso {
                anum = NumElements + Int(atom.getIsotope()) - 3 //pseudo AtNo for D, T
            }
            atomicCount[anum - 1] += 1
        }

        if atomicCount[5] != 0 { // Carbon (i.e. 6 - 1 = 5)
            if atomicCount[5] > ones {
                formula += "C" + sp + String(atomicCount[5]) + sp
            } else if atomicCount[5] == 1 {
                formula += "C"
            }

            atomicCount[5] = 0 // So we don't output C twice

            // only output H if there's also carbon -- otherwise do it alphabetical
            if atomicCount[0] != 0 { // Hydrogen (i.e., 1 - 1 = 0)
                if atomicCount[0] > ones {
                    formula += "H" + sp + String(atomicCount[0]) + sp
                } else if atomicCount[0] == 1 {
                    formula += "H"
                }

                atomicCount[0] = 0
            }
        }

        for j in 0..<NumElements {
            let alph = alphabetical[j] - 1
            if atomicCount[alph] != 0 {
                var symb: String
                if alph == NumElements - 1 {
                    symb = "T" //T
                } else if alph == NumElements - 2 {
                    symb = "D" //D
                } else {
                    symb = MKElements.getSymbol(alphabetical[j])
                }

                formula += symb + sp
                if atomicCount[alph] > ones {
                    formula += String(atomicCount[alph]) + sp
                }
            }
        }

        var chg = self.getTotalCharge()
        let ch = chg > 0 ? "+" : "-"
        chg = abs(chg)
        while chg > 0 {
            formula += ch + sp
            chg -= 1
        }
        return formula.trimmingCharacters(in: .newlines)
    }

    //! Stochoimetric formula (e.g., C4H6O).
    //!   This is either set by OBMol::SetFormula() or generated on-the-fly
    //!   using the "Hill order" -- i.e., C first if present, then H if present
    //!   all other elements in alphabetical order.
    func getFormula() -> String {
        if let dp = self.getData("Formula") as! MKPairData<String>? {
            // TODO: Potential Crash here is dp value is nil, but it shouldn't be if it was created through this method...
            return dp.getValue()!
        }

        let formula = self.getSpacedFormula(1, "", false)
        let dp = MKPairData<String>()
        dp.setAttribute("Formula")
        dp.setValue(formula)
        dp.setOrigin(.perceived)
        self.setData(dp)
        return formula
    }

    func setFormula(_ formula: String) {
        if let dp = self.getData("Formula") as! MKPairData<String>? {
            dp.setValue(formula)
            dp.setOrigin(.fileformatInput)
        } else {
            let dp = MKPairData<String>()
            dp.setAttribute("Formula")
            dp.setValue(formula)
            dp.setOrigin(.fileformatInput)
            self.setData(dp)
        } 
    }


    func setTotalCharge(_ charge: Int) {
        self.setFlag(OB_TCHARGE_MOL)
        self._totalCharge = charge
    }

    //! Returns the total molecular charge -- if it has not previously been set
    //!  it is calculated from the atomic formal charge information.
    //!  (This may or may not be correct!)
    //!  If you set atomic charges with OBAtom::SetFormalCharge()
    //!   you really should set the molecular charge with OBMol::SetTotalCharge()
    func getTotalCharge() -> Int {
        if self.hasFlag(OB_TCHARGE_MOL) {
            return self._totalCharge
        } else {
            print("GetTotalCharge -- calculated from formal charges")
            var charge = 0
            for atom: MKAtom in self.getAtomIterator() {
                charge += atom.getFormalCharge()
            }
            // Addition from original code, to avoid repeating this calculation
            self.setFlag(OB_TCHARGE_MOL)
            self._totalCharge = charge
            //
            return charge
        }
    }

    func setTotalSpinMultiplicity(_ spin: UInt) {
        self.setFlag(OB_TSPIN_MOL)
        self._totalSpin = spin
    }

    func setInteralCoord(_ int_coord: [MKInternalCoord]) {
        // The original implementation adds a nullptr to the start of the internal coordinate array
        // TODO: Let's see what happens if we only stick to number of atoms 

        if int_coord.count != self.numAtoms() {
            print("Error: Internal coordinate array size does not match number of atoms")
            return
        }

        self._internals = int_coord
    }


    //! Returns the total spin multiplicity -- if it has not previously been set
    //!  It is calculated from the atomic spin multiplicity information
    //!  assuming the high-spin case (i.e. it simply sums the number of unpaired
    //!  electrons assuming no further pairing of spins.
    //!  if it fails (gives singlet for odd number of electronic systems),
    //!  then assign wrt parity of the total electrons.
    func getTotalSpinMultiplicity() -> UInt {
        if self.hasFlag(OB_TSPIN_MOL) {
            return self._totalSpin
        } else {
            print("GetTotalSpinMultiplicity -- calculating from atomic spins assuming high spin case") 
            var unpairedElectrons = 0
            var charge = self.getTotalCharge()
            for atom: MKAtom in self.getAtomIterator() {
                if atom.getSpinMultiplicity() > 1 {
                    unpairedElectrons += Int(atom.getSpinMultiplicity() - 1)
                }
                charge += atom.getAtomicNum()
            }
            // Deviate from original by setting the flag here as well
            if charge % 2 != unpairedElectrons % 2 {
                self.setFlag(OB_TSPIN_MOL)
                self._totalSpin = UInt(abs(charge) % 2 + 1)
                return self._totalSpin
            } else {
                self.setFlag(OB_TSPIN_MOL)
                self._totalSpin = UInt(unpairedElectrons + 1)
                return self._totalSpin
            }

        }
    }

    // Cannot override default assignment operator in Swift, maybe that is not bad thing?

    static public func += (lhs: inout MKMol, rhs: MKMol) {
        // TODO: IMPLEMENT
        // After implementing StereoBase and Peception class 
    } 

    override func clear() {
        // Destroy Atom list 
        self._vatom.removeAll()
        // Destroy Bond list
        self._vbond.removeAll()
        // Destroy AtomID list
        self._vatomIds.removeAll()
        // Destroy BondID list
        self._vbondIds.removeAll()
        // Delete Residues
        self._residue.removeAll()
        // Clear multi-conformer data
        self._vconf.removeAll()
        //Clear flags except OB_PATTERN_STRUCTURE which is left the same
        self._flags &= OB_PATTERN_STRUCTURE

        self._c = Array(repeating: 0.0, count: 3)
        self._mod = 0
    }

    // Swift has ARC and does not need a manual dealloc for these
    // private func destroyAtom(_ atom: MKAtom) { }
    // private func destroyBond(_ bond: MKBond) { }
    // private func destroyResidue(_ residue: MKResidue) { }

    func deleteAtom(_ atom: MKAtom) {
        
    }

    func deleteBond(_ bond: MKBond) {
        
    }


    func numAtoms() -> Int {
        return self._vatom.count
    }
    
    func numBonds() -> Int {
        return self._vbond.count
    }

    func addBond(_ start_idx: Int, _ end_idx: Int, _ type: Int) {}


    func automaticPartialCharge() -> Bool {
        return false
    }

    func beginModify() {
        // This should probably be retired as a more Swifty approach is adopted
        //suck coordinates from _c into _v for each atom
        if (self._mod == 0) && !self.isEmpty() {
            for atom in self._vatom {
                atom.setVector()
                atom.clearCoordPtr()
            }
            
            self._c = []
            self._vconf.removeAll()
            
            //Destroy rotamer list if necessary
//            if ((OBRotamerList *)GetData(OBGenericDataType::RotamerList))
//              {
//                delete (OBRotamerList *)GetData(OBGenericDataType::RotamerList);
//                DeleteData(OBGenericDataType::RotamerList);
//              }
        }
        
        self._mod += 1 
    }
    
    func endModify(_ nukePercievedData: Bool) {
        if self._mod == 0 {
            print("_mod is negative - EndModify() called too many times")
            return 
        }

        self._mod -= 1

        if nukePercievedData {
        // wipe all but whether it has aromaticity perceived, is a reaction, or has periodic boundaries enabled
            self._flags &= (OB_AROMATIC_MOL|OB_REACTION_MOL|OB_PERIODIC_MOL)
        }

        // TODO: Does it really need to have an extra _c array? Why couldn't the normal _v array suffice? 
        // Does this verbosity matter with Swift in a performance sense, and how can we streamline this?

        self._c = []
        if self.isEmpty() {
            return
        }

        //if atoms present convert coords into array
        var idx = 0
        var c: [[Double]] = []
        for atom in self._vatom {
            atom.setIdx(idx+1)
            c.append(atom.getVector().scalars)
            atom.setCoordPtr(c[idx])
            idx += 1
        }

        self._vconf.append(c)

        // Always remove angle and torsion data, since they will interfere with the iterators
        // PR#2812013
        self.deleteData(.AngleData)
        self.deleteData(.TorsionData)
    }

    func newAtom() -> MKAtom {
        return self.newAtom(UInt(self._vatomIds.count))
    }

    //! \brief Instantiate a New Atom and add it to the molecule
    //!
    //! Checks bond_queue for any bonds that should be made to the new atom
    //! and updates atom indexes.
    func newAtom(_ id: UInt) -> MKAtom {
        
        // Begin Modify is commented out here in original code, taken as should be called beforehand 

        if id >= self._vatomIds.count {
            self._vatomIds.reserveCapacity(self._vatomIds.count+1)
            //  Do not need to prefill the array with nils, as it is already done
        }
        
        let atom = MKAtom()
        atom.setIdx(self._natoms+1)
        atom.setParent(self)

        self._vatomIds[Int(id)] = atom
        atom.setId(Int(id))
        
        self._vatom.append(atom)

        if self.hasData(.VirtualBondData) {
            /*add bonds that have been queued*/

            guard let bondQueue: [MKVirtualBond] = self.getDataVector(.VirtualBondData)! as? [MKVirtualBond] else { return atom }

            for bond: MKVirtualBond in bondQueue {
                if bond.getBgn() > self._natoms || bond.getEnd() > self._natoms {
                    continue
                }
                if atom.getIdx() == bond.getBgn() || atom.getIdx() == bond.getEnd() {
                    self.addBond(Int(bond.getBgn()), Int(bond.getEnd()), Int(bond.getOrder()))
                }
            }

            self.deleteData(bondQueue)
        }   

        // End Modify 
        return atom 
    }

    func newResidue() -> MKResidue {
        let newRes = MKResidue()
        newRes.setIdx(self._residue.count)
        self._residue.append(newRes)
        return newRes
    }

    func newBond() -> MKBond {
        return self.newBond(UInt(self._vbondIds.count))
    }

  //! \since version 2.1
  //! \brief Instantiate a New Bond and add it to the molecule
  //!
  //! Sets the proper Bond index and insures this molecule is set as the parent.
    func newBond(_ id: UInt) -> MKBond {
        // Begin Modify is commented out here in original code, taken as should be called beforehand 

        if id >= self._vbondIds.count {
            self._vbondIds.reserveCapacity(self._vbondIds.count+1)
            self._vbondIds.removeAll()
        }

        let bond = MKBond()
        bond.setIdx(self._nbonds+1)
        bond.setParent(self)

        self._vbondIds[Int(id)] = bond
        bond.setId(Int(id))

        self._vbond.append(bond)

        // End Modify 
        return bond 
    }

//! \brief Add an atom to a molecule
//!
//! Also checks bond_queue for any bonds that should be made to the new atom
    func addAtom(_ atom: MKAtom, _ forceNewId: Bool) -> Bool {
        //    BeginModify();
        // Use the existing atom Id unless either it's invalid or forceNewId has been specified
        var id: MKBaseID
        if forceNewId {
            id = MKBaseID._id(self._vatomIds.count)
        } else {
            id = atom.getId()
            if id == MKBaseID.NoId {
                id = MKBaseID._id(self._vatomIds.count)
            }
        }

        atom.setIdx(self._natoms+1)
        atom.setParent(self)
        if id.rawValue >= self._vatomIds.count {
            self._vatomIds.reserveCapacity(id.rawValue+1)
        }

        atom.setId(id.rawValue)
        self._vatomIds[id.rawValue] = atom

        if self._natoms+1 >= self._vatom.count {
            self._vatom.reserveCapacity(self._natoms+1)
        }

        self._vatom[self._natoms] = atom

        if self.hasData(.VirtualBondData) {
            /*add bonds that have been queued*/

            guard let bondQueue: [MKVirtualBond] = self.getDataVector(.VirtualBondData)! as? [MKVirtualBond] else { return true }

            for bond: MKVirtualBond in bondQueue {
                if bond.getBgn() > self._natoms || bond.getEnd() > self._natoms {
                    continue
                }
                if atom.getIdx() == bond.getBgn() || atom.getIdx() == bond.getEnd() {
                    self.addBond(Int(bond.getBgn()), Int(bond.getEnd()), Int(bond.getOrder()))
                }
            }

            self.deleteData(bondQueue)
        }

        // End Modify

        return true
    }

    func insertAtom(_ atom: MKAtom) -> Bool {
        self.beginModify()
        _ = self.addAtom(atom, false)
        self.endModify(false)

        return true
    }

    func addResidue(_ residue: MKResidue) -> Bool{
        self.beginModify()
        residue.setIdx(self._residue.count)
        self._residue.append(residue)
        self.endModify(false)

        return true
    }

    func stripSalts(_ threshold: UInt) -> Bool {
        var cfl = Array<Array<Int>>()
        self.contigFragList(&cfl)
        return false 
    }

    //    MARK: HAS/IS Functions

    func has2D(_ not3D: Bool = false) -> Bool {
        return false
    }
    
    func has3D() -> Bool {
        return false
    }
    
    func hasNonZeroCoords() -> Bool {
        return false
    }
    
    func hasAromaticPerceived() -> Bool {
        return self.hasFlag(OB_AROMATIC_MOL)
    }
    
    func hasSSSRPerceived() -> Bool {
        return self.hasFlag(OB_SSSR_MOL)
    }
    
    func hasLSSRPerceived() -> Bool {
        return self.hasFlag(OB_LSSR_MOL)
    }
    
    func hasRingAtomsAndBondsPerceived() -> Bool {
        return self.hasFlag(OB_RINGFLAGS_MOL)
    }
    
    func hasAtomTypesPerceived() -> Bool {
        return self.hasFlag(OB_ATOMTYPES_MOL)
    }
    
    func hasRingTypesPerceived() -> Bool {
        return self.hasFlag(OB_RINGTYPES_MOL)
    }
    
    func hasChiralityPerceived() -> Bool {
        return self.hasFlag(OB_CHIRALITY_MOL)
    }
    
    func hasPartialChargesPerceived() -> Bool {
        return self.hasFlag(OB_PCHARGE_MOL)
    }
    
    func hasHybridizationPerceived() -> Bool {
        return self.hasFlag(OB_HYBRID_MOL)
    }
    
    func hasClosureBondsPerceived() -> Bool {
        return self.hasFlag(OB_CLOSURE_MOL)
    }

    func hasChainsPerceived() -> Bool {
        return self.hasFlag(OB_CHAINS_MOL)
    }
    
    func hasHydrogensAdded() -> Bool {
        return self.hasFlag(OB_H_ADDED_MOL)
    }
    
    func isReaction() -> Bool {
        return self.hasFlag(OB_REACTION_MOL)
    }
    
    func isEmpty() -> Bool {
        return self._natoms == 0
    }
    
    func isPeriodic() -> Bool {
        return self.hasFlag(OB_PERIODIC_MOL)
    }
    
//    MARK: FIND Functions
    
    //! locates all atoms for which there exists a path to 'end'
    //! without going through 'bgn'
    //! children must not include 'end'
    func findChildren(_ bgn: MKAtom, _ end: MKAtom, _ children: inout [MKAtom]) {
        
        var used: MKBitVec = MKBitVec()
        var curr: MKBitVec = MKBitVec()
        var next: MKBitVec = MKBitVec()

        used |= bgn.getIdx()
        used |= end.getIdx()
        curr |= end.getIdx()
        children.removeAll()
        
        // MARK: Potentially bad code but unsure how to convert for (;;) loop
        while true {
            next.clear()
            for i in curr.nextBit(-1)..<curr.endBit() {
                guard let atom = self.getAtom(i) else { break }
                guard let neighA = atom.getNbrAtomIterator() else { continue }
                for nbr in neighA {
                    if !used[nbr.getIdx()] {
                        children.append(nbr)
                        next |= nbr.getIdx()
                        used |= nbr.getIdx()
                    }
                }
                
            }
            if next.isEmpty() { break }
            curr = next
        }
    }

    //! locates all atoms for which there exists a path to 'second'
    //! without going through 'first'
    //! children must not include 'second'
    func findChildren(_ first: Int, _ second: Int, _ children: inout [Int]) {

        var used: MKBitVec = MKBitVec()
        var curr: MKBitVec = MKBitVec()
        let next: MKBitVec = MKBitVec()
        
        used.setBitOn(UInt32(first))
        used.setBitOn(UInt32(second))
        curr.setBitOn(UInt32(second))

        while !curr.isEmpty() {
            next.clear()
            for i in curr.nextBit(-1)..<curr.endBit() {
                guard let atom = self.getAtom(i) else { break }
                guard let bondA = atom.getBondIterator() else { continue }
                for bond in bondA {
                    if !used[bond.getNbrAtomIdx(atom)] {
                        next.setBitOn(UInt32(bond.getNbrAtomIdx(atom)))
                    }
                }
            }
            used |= next
            curr = next
        }
        
        used.setBitOff(UInt32(first))
        used.setBitOff(UInt32(second))
        used.toVecInt(&children)
    }
    
    func findAngles() {
        if self.hasData(.AngleData) {
            return
        }
        
        let newAngleData: MKAngleData = MKAngleData()
        newAngleData.setOrigin(.perceived)
        self.setData(newAngleData)
                
        for atom in self.getAtomIterator() {
            
            if atom.getAtomicNum() == 1 { // Hydrogen
                continue
            }
            
            guard let neigh = atom.getNbrAtomIterator() else { continue }
            
            for neighA in neigh {
                for neighB in neigh {
                    if neighA != neighB {
//                        MARK: Why does this not fill in the real angle? Maybe add a member function \
//                        to MKAngle to automagically calculate the angle?
                        let angle = MKAngle()
                        angle.setAtoms(atom, neighA, neighB)
                        newAngleData.setData(angle)
                    }
                }
            }
            
        }
    }
    
    func findTorsions() {
        if self.hasData(.TorsionData) {
            return
        }
        
        let newTorsionData: MKTorsionData = MKTorsionData()
        newTorsionData.setOrigin(.perceived)
        self.setData(newTorsionData)

        for bond in self.getBondIterator() {
            let b = bond.getBeginAtom()
            let c = bond.getEndAtom()
            let torsion = MKTorsion()
            
            if b.getAtomicNum() == 1 || c.getAtomicNum() == 1 {
                continue
            }
            
            guard let bNeigh = b.getNbrAtomIterator() else { continue }
            
            for a in bNeigh {
                if a == c {
                    continue
                }
                
                guard let cNeigh = c.getNbrAtomIterator() else { continue }
                
                for d in cNeigh {
                    if d == b || d == a {
                        continue
                    }
                    
                    
                    _ = torsion.addTorsion(a, b, c, d)
                }
            }
            //add torsion to torsionData
            if torsion.getSize() > 0{
                newTorsionData.setData(torsion)
            }
        }
    }
    
    func findLargestFragment(_ lf: inout MKBitVec) {
        var used: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        var curr: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        let next: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        var frag: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        
        lf.clear()
        while used.countBits() < numAtoms() {
            for atom in self.getAtomIterator() {
                if !used.bitIsSet(atom.getIdx()) {
                    curr.setBitOn(UInt32(atom.getIdx()))
                    break
                }
            }
            frag |= curr
            while !curr.isEmpty() {
                for j in curr.nextBit(-1)..<curr.endBit() {
                    guard let atom = self.getAtom(j) else { continue }
                    guard let bonds = atom.getBondIterator() else { continue }
                    for bond in bonds {
                        if !used.bitIsSet(bond.getNbrAtomIdx(atom)) {
                            next.setBitOn(UInt32(bond.getNbrAtomIdx(atom)))
                        }
                    }
                }
                used |= curr
                used |= next
                frag |= next
                curr = next
            }
            if lf.isEmpty() || lf.countBits() < frag.countBits() {
                lf = frag
            }
        }

    }

    func findRingAtomsAndBonds() {
        
    }

    func findSSSR() {
        
    }
    
    func findLSSR() {
        
    }
    
    
    
//    MARK: Island of Misfit Toys

    func getAngle(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom) -> Double {
        return a.getAngle(b, c)
    }

    func getTorsion(_ a: Int, _ b: Int, _ c: Int, _ d: Int) -> Double {
        let i = self._vatom[a-1]
        let j = self._vatom[b-1]
        let k = self._vatom[c-1]
        let m = self._vatom[d-1]
        return self.getTorsion(i,j,k,m)
    }
    
    func getTorsion(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom, _ d: MKAtom) -> Double {
        
        if !self.isPeriodic() {
            return calculateTorsionAngle(a.getVector(), b.getVector(), c.getVector(), d.getVector())
        }else {
            let v1 = a.getVector()
            var v2 = b.getVector()
            var v3 = c.getVector()
            var v4 = d.getVector()

            guard let box = self.getData(.UnitCell) as? MKUnitCell else { return 0.0 }
            v2 = box.unwrapCartesianNear(v2, v1)
            v3 = box.unwrapCartesianNear(v3, v2)
            v4 = box.unwrapCartesianNear(v4, v3)
            return calculateTorsionAngle(v1, v2, v3, v4)
        }
        
    }

//    Angle needs to be in radians
    func setTorsion(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom, _ d: MKAtom, _ angle: Double) {
        
        var tors: [Int] = []
        var m: [Double] = []
        m.reserveCapacity(9)
        
        
        tors.append(Int(a.getCoordinateIdx()))
        tors.append(Int(b.getCoordinateIdx()))
        tors.append(Int(c.getCoordinateIdx()))
        tors.append(Int(d.getCoordinateIdx()))
        
        var children: [Int] = []
        self.findChildren(b.getIdx(), c.getIdx(), &children)
        
        for j in 0..<children.count {
            children[j] = (children[j] - 1) * 3
        }
        
//        in radians!
        //calculate the torsion angle
        // TODO: fix this calculation for periodic systems
        let radang = calculateTorsionAngle(a.getVector(), b.getVector(), c.getVector(), d.getVector()).degreesToRadians
        // now we have the torsion angle (radang) - set up the rot matrix
        
        //find the difference between current and requested
        let rotang = angle - radang
        
        let sn = sin(rotang)
        let cs = cos(rotang)
        
        let t = 1 - cs
        
        let v2x = self._c[tors[1]] - self._c[tors[2]]
        let v2y = self._c[tors[1]+1] - self._c[tors[2]+1]
        let v2z = self._c[tors[1]+2] - self._c[tors[2]+2]

        let mag = sqrt(v2x * v2x + v2y * v2y + v2z * v2z)
        let x = v2x / mag
        let y = v2y / mag
        let z = v2z / mag

        m[0] = t * x * x + cs
        m[1] = t * x * y + sn * z
        m[2] = t * x * z - sn * y
        m[3] = t * x * y - sn * z
        m[4] = t * y * y + cs
        m[5] = t * y * z + sn * x
        m[6] = t * x * z + sn * y
        m[7] = t * y * z - sn * x
        m[8] = t * z * z + cs

        //
        //now the matrix is set - time to rotate the atoms
        //

        let tx = self._c[tors[1]]
        let ty = self._c[tors[1]+1]
        let tz = self._c[tors[1]+2]

        for j in children {
            self._c[j] -= tx
            self._c[j+1] -= ty
            self._c[j+2] -= tz
            let x = self._c[j] * m[0] + self._c[j+1] * m[1] + self._c[j+2] * m[2]
            let y = self._c[j] * m[3] + self._c[j+1] * m[4] + self._c[j+2] * m[5]
            let z = self._c[j] * m[6] + self._c[j+1] * m[7] + self._c[j+2] * m[8]
            self._c[j] = x
            self._c[j+1] = y
            self._c[j+2] = z
            self._c[j] += tx
            self._c[j+1] += ty
            self._c[j+2] += tz
        }
        
    }

    func setPeriodicMol(_ value: Bool = true) {
        if value {
            self.setFlag(OB_PERIODIC_MOL)
        } else {
            self.unsetFlag(OB_PERIODIC_MOL)
        }
    }
    
    func contigFragList(_ cfl: inout Array<Array<Int>>) {
       
        var used: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        var curr: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        let next: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        var frag: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        

        while used.countBits() < self.numAtoms() {
            curr.clear()
            frag.clear()
            for atom in self.getAtomIterator() {
                if !used.bitIsSet(atom.getIdx()) {
                    curr.setBitOn(UInt32(atom.getIdx()))
                    break
                }
            }
            frag |= curr
            while !curr.isEmpty() {
                next.clear()
                for j in curr.nextBit(-1)..<curr.endBit() {
                    guard let atom = self.getAtom(j) else { continue }
                    for bond in atom.getBondIterator()! {
                        if !used.bitIsSet(bond.getNbrAtomIdx(atom)) {
                            next.setBitOn(UInt32(bond.getNbrAtomIdx(atom)))
                        }
                    }
                }
                used |= curr
                used |= next
                frag |= next
                curr = next
            }
            var tmp: [Int] = []
            frag.toVecInt(&tmp)
            cfl.append(tmp)
        }

        cfl.sort(by: { (a, b) -> Bool in
            return a.count < b.count
        })

    }
    
    /*! \brief Calculates the graph theoretical distance (GTD) of each atom.
    *
    * Creates a vector (indexed from zero) containing, for each atom
    * in the molecule, the number of bonds plus one to the most
    * distant non-H atom.
    *
    * For example, for the molecule H3CC(=O)Cl the GTD value for C1
    * would be 3, as the most distant non-H atom (either Cl or O) is
    * 2 bonds away.
    *
    * Since the GTD measures the distance to non-H atoms, the GTD values
    * for terminal H atoms tend to be larger than for non-H terminal atoms.
    * In the example above, the GTD values for the H atoms are all 4.
    */
    func getGTDVector( _ gtd: inout [Double]) {
        //calculates the graph theoretical distance for every atom
        //and puts it into gtd
        gtd = [Double](repeating: 0.0, count: self.numAtoms())

        var gtdcount: Int = 0

        var used: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        var curr: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))
        let next: MKBitVec = MKBitVec(UInt32(self.numAtoms() + 1))

        next.clear()

        for atom in self.getAtomIterator() {
            gtdcount = 0
            used.clear()
            curr.clear()
            used.setBitOn(UInt32(atom.getIdx()))
            curr.setBitOn(UInt32(atom.getIdx()))

            while !curr.isEmpty() {
                next.clear()
                for natom in curr.nextBit(-1)..<curr.endBit() {
                    guard let atom1 = self.getAtom(natom) else { continue }
                    for bond in atom1.getBondIterator()! {
                        if !used.bitIsSet(bond.getNbrAtomIdx(atom1)) && !curr.bitIsSet(bond.getNbrAtomIdx(atom1)) {
                            if bond.getNbrAtom(atom1).getAtomicNum() != 1 { // Hydrogen
                                next.setBitOn(UInt32(bond.getNbrAtomIdx(atom1)))
                            }
                        }
                    }
                }
                used |= next
                curr = next
                gtdcount += 1
            }
            gtd[atom.getIdx() - 1] = Double(gtdcount)
        }
    }

    /*!
    **\brief Calculates a set of graph invariant indexes using
    ** the graph theoretical distance, number of connected heavy atoms,
    ** aromatic boolean, ring boolean, atomic number, and
    ** summation of bond orders connected to the atom.
    ** Vector is indexed from zero
    */    
    func getGIVector(_ vid: inout [UInt]) {
        
        vid = [UInt](repeating: 0, count: self.numAtoms() + 1)

        var v: [Double] = []
        self.getGTDVector(&v)

        var i: Int = 0

        for atom in self.getAtomIterator() {
            vid[i] = UInt(v[i])
            vid[i] += UInt(atom.getHeavyDegree() * 100)
            vid[i] += UInt((atom.isAromatic() ? 1 : 0) * 1000)
            vid[i] += UInt((atom.isInRing() ? 1 : 0) * 10000)
            vid[i] += UInt(atom.getAtomicNum() * 100000)
            vid[i] += UInt(atom.getImplicitHCount() * 10000000)
            i += 1
        }
    }
    
    



    
}
