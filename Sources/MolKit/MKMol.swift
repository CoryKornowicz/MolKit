//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation
import Collections
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
    private var _vconf: Array<Array<Double>> = [] // MARK: Not sure how I feel about this, would like to simplify into one unifying structure
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
    
    private func incrementMod() { self._mod += 1 }
    private func decrementMod() { self._mod -= 1 }
    
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
    
    func getConformer(_ i: Int) -> [Double] {
        return self._vconf[i]
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
    
    func numConformers() -> Int {
        return self._vconf.isEmpty ? 0 : self._vconf.count
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

        self._c = []
        self._mod = 0
    }

    // Swift has ARC and does not need a manual dealloc for these
    // private func destroyAtom(_ atom: MKAtom) { }
    // private func destroyBond(_ bond: MKBond) { }
    // private func destroyResidue(_ residue: MKResidue) { }

    func deleteAtom(_ atom: MKAtom, _ destroyAtom: Bool = false) -> Bool {
        if atom.getAtomicNum() == 1 {
            return self.deleteHydrogen(atom)
        }
        
        self.beginModify()
        //don't need to do anything with coordinates b/c
        //BeginModify() blows away coordinates
        
        //delete bonds to neighbors
        if var bonds = atom.getBondIterator() {
            for bond in bonds {
                self.deleteBond(bond)
            }
        }
        
        self._vatomIds.remove(at: atom.getId().rawValue)
        self._vatom.remove(at: atom.getIdx())
        
        //reset all the indices to the atoms
        var _idx = 0
        for atomi in self.getAtomIterator() {
            atom.setIdx(_idx)
            _idx += 1
        }
        
        self.endModify(false)
        
        // Delete any stereo objects involving this atom
        let id: Ref = .Ref(atom.getId().rawValue)
        MKMol.deleteStereoOnAtom(self, id)
        
//        if destroyAtom {}
        self.setSSSRPerceived(false)
        self.setLSSRPerceived(false)
        return true
    }

    func deleteResidue(_ residue: MKResidue, _ destroyResidue: Bool = false) -> Bool {
        var idx = residue.getIdx()
        self._residue.remove(at: Int(idx))
        var _idx = 0

        for res in self._residue {
            res.setIdx(_idx)
            _idx += 1
        }
        // if destroyResidue { }
        self.setSSSRPerceived(false)
        self.setLSSRPerceived(false)
        return true
    }
    
    func deleteBond(_ bond: MKBond, _ destroyBond: Bool = false) -> Bool {
        self.beginModify()

        bond.getBeginAtom().deleteBond(bond)
        bond.getEndAtom().deleteBond(bond)
        self._vbond.remove(at: Int(bond.getIdx()))
        self._vbondIds.remove(at: bond.getId())

        //reset all the indices to the atoms
        var _idx = 0
        for bond in self.getBondIterator() {
            bond.setIdx(_idx)
            _idx += 1
        }

        self.endModify(false)

        // if destroyBond { }
        self.setSSSRPerceived(false)
        self.setLSSRPerceived(false)
        return true
    }

    func numAtoms() -> Int {
        return self._vatom.count
    }
    
    func numBonds() -> Int {
        return self._vbond.count
    }

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
    
    func endModify(_ nukePercievedData: Bool = false) {
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
        var c: Array<Double> = []
        for atom in self._vatom {
            atom.setIdx(idx+1)
            c.append(contentsOf: atom.getVector().scalars)
            atom.setCoordPtr(atom.getVector().scalars)
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
    
    func addBond(_ first: Int, _ second: Int, _ order: Int, _ flags: Int = 0, insertpos: Int = -1) -> Bool {
        
//       Does not invoke beginModify/endModify
        
        if first == second || (self.getBond(first, second) != nil) {
            return false
        }
        
        if first <= self.numAtoms() && second <= self.numAtoms() {
            let bond = MKBond()
            guard let bgn = self.getAtom(first) else { return false }
            guard let end = self.getAtom(second) else { return false }
            
            bond.set(_nbonds, bgn, end, UInt(order), UInt(flags))
            bond.setParent(self)
            bond.setId(self._vbondIds.count)
            self._vbondIds.append(bond)
            
            self._vbond.append(bond)
            if insertpos == -1 {
                bgn.addBond(bond)
                end.addBond(bond)
            } else {
                if insertpos >= bgn.getExplicitDegree() {
                    bgn.addBond(bond)
                } else//need to insert the bond for the connectivity order to be preserved
                 {    //otherwise stereochemistry gets screwed up
                    bgn.insertBond(insertpos, bond)
                }
                end.addBond(bond)
            }
        } else {
            //at least one atom doesn't exist yet - add to bond_q
            self.setData(MKVirtualBond(UInt(first), UInt(second), UInt(order), flags))
        }
        
        return true
    }

    func addBond(_ bond: MKBond) -> Bool {
        if !self.addBond(bond.getBeginAtomIdx(), bond.getEndAtomIdx(), Int(bond.getBondOrder()), Int(bond.getFlags())) {
            return false
        }
        
        //copy the bond's generic data
        for di in MKIterator(bond.getDataVector() ?? []){
            guard let newBond = self.getBond(self.numBonds()-1) else { return false }
            newBond.cloneData(di)
        }
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

        if cfl.isEmpty || cfl.count == 1 {
            return false
        }

        let maxElem = cfl.max { $0.count < $1.count }
        
        var delatoms: [MKAtom] = []
        var atomIndices: Set<Int> = []
        
        for i in cfl {
            if (i.count < threshold) || (threshold == 0 && i != maxElem) {
                for j in i {
                    if atomIndices.firstIndex(of: j) == atomIndices.endIndex { // MARK: Does this transfer well from C++
                        delatoms.append(self.getAtom(j)!) // MARK: Force unwrap here could be detrimental
                        atomIndices.insert(j)
                    }
                }
            }
        }
        
        if !delatoms.isEmpty {
            //      int tmpflags = _flags & (~(OB_SSSR_MOL));
            self.beginModify()
            for atom in delatoms {
                self.deleteAtom(atom)
            }
            self.endModify(false)
            //      _flags = tmpflags;  // Gave crash when SmartsPattern::Match()
            // was called susequently
            // Hans De Winter; 03-nov-2010
        }
        
        return true
    }
    
    // Convenience function used by the DeleteHydrogens methods
    private func isSuppressibleHydrogen(_ atom: MKAtom) -> Bool {
        return atom.getIsotope() == 0 && atom.getHeavyDegree() == 1 && atom.getFormalCharge() == 0 && (atom.getData("Atom Class") == nil)
    }
    
    func deletePolarHydrogens() -> Bool {
        var delatoms: [MKAtom] = []
        
        for atom in self.getAtomIterator() {
            if atom.isPolarHydrogen() && self.isSuppressibleHydrogen(atom) {
                delatoms.append(atom)
            }
        }
        
        if delatoms.isEmpty {
            return true
        }
        
        self.incrementMod()
        
        for atom in delatoms {
            self.deleteAtom(atom)
        }
        
        self.decrementMod()
        
        self.setSSSRPerceived(false)
        self.setLSSRPerceived(false)
        return true
    }
    
    func deleteNonPolarHydrogens() -> Bool {
        var delatoms: [MKAtom] = []
        
        for atom in self.getAtomIterator() {
            if atom.isNonPolarHydrogen() && self.isSuppressibleHydrogen(atom) {
                delatoms.append(atom)
            }
        }
        
        if delatoms.isEmpty {
            return true
        }
        
        self.incrementMod()
        
        for atom in delatoms {
            self.deleteAtom(atom)
        }
        
        self.decrementMod()
        
        self.setSSSRPerceived(false)
        self.setLSSRPerceived(false)
        return true
    }
    
    func deleteHydrogens() -> Bool {
        var delatoms: [MKAtom] = []
        
        for atom in self.getAtomIterator() {
            if atom.getAtomicNum() == 1 && self.isSuppressibleHydrogen(atom) {
                delatoms.append(atom)
            }
        }
        
        self.setHydrogensAdded(false)
        
        if delatoms.isEmpty {
            return true
        }
        
        /* decide whether these flags need to be reset
           _flags &= (~(OB_ATOMTYPES_MOL));
           _flags &= (~(OB_HYBRID_MOL));
           _flags &= (~(OB_PCHARGE_MOL));
           _flags &= (~(OB_IMPVAL_MOL));
        */
        
        self.incrementMod()
        
        // This has the potential to be slow -- we still need methods to delete a set of atoms
        //  and to delete a set of bonds
        for atom in delatoms {
//            Get neighbor atom, should only have one since it is a hydrogen atom
            guard let nbr = atom.getNbrAtomIterator() else {
                self.deleteAtom(atom)
                continue
            }
            
            for nb in nbr {
                nb.setImplicitHCount(nb.getImplicitHCount() + 1)
            }
            
            self.deleteAtom(atom)
        }
        
        self.decrementMod()
        
        self.setSSSRPerceived(false)
        self.setLSSRPerceived(false)
        return true

    }

    func deleteHydrogens(_ atom: MKAtom) -> Bool {
        var delatoms: [MKAtom] = []
        
        guard let nbrs = atom.getNbrAtomIterator() else { return false }
        
        for nb in nbrs {
            if nb.getAtomicNum() == 1 && self.isSuppressibleHydrogen(nb) {
                delatoms.append(nb)
            }
        }
        
        if delatoms.isEmpty {
            return true
        }
        
        self.incrementMod()
        
        for atom in delatoms {
            _ = self.deleteHydrogen(atom)
        }
        
        self.decrementMod()
        self.setHydrogensAdded(false)
        self.setSSSRPerceived(false)
        self.setLSSRPerceived(false)
        return true
    }
    
    //deletes the hydrogen atom passed to the function
    func deleteHydrogen(_ atom: MKAtom) -> Bool {
        if atom.getAtomicNum() != 1 { return false }
        
        let atomIdx = atom.getIdx()
        
        //find bonds to delete
        guard let bonds = atom.getBondIterator() else { return false } // Probably want to continue to delete atom
        self.incrementMod()
        for bond in bonds {
            self.deleteBond(bond)
        }
        self.decrementMod()
        
        if atomIdx != self.numAtoms() {
            let idx = Int(atom.getCoordinateIdx())
//            let size = self.numAtoms() - atom.getIdx()
//            Resize the _vconf arrays to not include this atom 
            for i in 0..<self._vconf.count {
                self._vconf[i].removeSubrange(idx..<idx+3)
            }
        }
        
        // Deleting hydrogens does not invalidate the stereo objects
        // - however, any explicit refs to the hydrogen atom must be
        //   converted to implicit refs
        let ref = Ref(rawValue: atom.getId().rawValue)
        stereoRefToImplicit(self, ref!)
//      Cast Ref back to int
        var id = 0
        switch ref {
        case .Ref(let value):
            id = value
        default:
            id = 0
//            MARK: Throw error here since this should never error
        }
        
        self._vatomIds.remove(at: id)
        self._vatom.remove(at: id)
        var _idx = 1
        //reset all the indices to the atoms
        for atom in self._vatom {
            atom.setIdx(_idx)
            _idx += 1
        }
        
        self.setHydrogensAdded(false)
//        destroy atom??
        self.setSSSRPerceived(false)
        self.setLSSRPerceived(false)
        return true
        
    }
    
    /*
      this has become a wrapper for backward compatibility
    */
    
    func addHydrogens(_ polaronly: Bool, _ correctForPH: Bool, _ pH: Double) -> Bool {
        return self.addNewHydrogens(polaronly ? .PolarHydrogen : .AllHydrogen, correctForPH, pH)
    }
    
    private func atomIsNSOP(_ atom: MKAtom) -> Bool {
        switch atom.getAtomicNum() {
        case MKElements.getAtomicNum("N"), MKElements.getAtomicNum("S"), MKElements.getAtomicNum("O") ,MKElements.getAtomicNum("P"):
            return true
        default:
            return false
        }
    }
    
    private func isNotCorH(_ atom: MKAtom) -> Bool {
        switch atom.getAtomicNum() {
        case MKElements.getAtomicNum("C"), MKElements.getAtomicNum("H"):
            return false
        default: break
        }
        return true
    }
    
    //! \return a "corrected" bonding radius based on the hybridization.
    //! Scales the covalent radius by 0.95 for sp2 and 0.90 for sp hybrids
    func correctedBondRad(_ ele: Int, _ hyb: Int) -> Double {
        let rad: Double = MKElements.getCovalentRad(ele)
        switch (hyb) {
            case 2:
                return rad * 0.95
            case 1:
                return rad * 0.90
            default:
                return rad
        }
    }
    
    func addNewHydrogens(_ whichHydrogen: HydrogenType, _ correctForPH: Bool, _ pH: Double = 0.0) -> Bool {
        
        if !self.isCorrectedForPH() && correctForPH {
            return self.correctForPH(pH)
        }
        
        if self.hasHydrogensAdded() {
            return true
        }
        
        var hasChiralityPerceived = self.hasChiralityPerceived()
        
        /*
        //
        // This was causing bug #1892844 in avogadro. We also want to add hydrogens if the molecule has no bonds.
        //
        if(NumBonds()==0 && NumAtoms()!=1)
        {
        obErrorLog.ThrowError(__FUNCTION__,
        "Did not run OpenBabel::AddHydrogens on molecule with no bonds", obAuditMsg);
        return true;
        }
        */
//        TODO: Add logging here of which hydrogen was used
        
        // Make sure we have conformers (PR#1665519)
        if !self._vconf.isEmpty && !self.isEmpty() {
            for atom in self.getAtomIterator() {
                atom.setVector()
            }
        }
        
        self.setHydrogensAdded(true) // This must come after EndModify() as EndModify() wipes the flags
        // If chirality was already perceived, remember this (to avoid wiping information
        if hasChiralityPerceived {
            self.setChiralityPerceived(true)
        }
        
        //count up number of hydrogens to add
        
        var count = 0
        var hcount = 0
        
        var vhadd: [Pair<MKAtom, Int>] = []
        
        for atom in self.getAtomIterator() {
            if whichHydrogen == .PolarHydrogen && !self.atomIsNSOP(atom) { continue }
            if whichHydrogen == .NonPolarHydrogen && self.atomIsNSOP(atom) { continue }
            
            hcount = Int(atom.getImplicitHCount())
            atom.setImplicitHCount(0)
            
            if hcount != 0 {
                vhadd.append((atom, hcount))
                count += hcount
            }
        }
        
        if count == 0 {
            // Make sure to clear SSSR and aromatic flags we may have tripped above
            self._flags = (~(OB_SSSR_MOL | OB_AROMATIC_MOL))
            return true
        }
                
        //realloc memory in coordinate arrays for new hydrogens
//        MARK: Not really a thing in Swift....
        
        self.incrementMod()
        
        let hbrad = self.correctedBondRad(1, 0)
        
        for k in vhadd {
            let atom = k.0
            let bondlen = hbrad + self.correctedBondRad(atom.getAtomicNum(), atom.getHyb())
            
            for _ in 0..<k.1 {
                var badh = 0
                for n in 0..<self.numConformers() {
                    self.setConformer(n)
                    if self.hasNonZeroCoords() { // MARK: Move check for coords to be per loop                                       since conformers is changing
                        // Ensure that add hydrogens only returns finite coords
                        //atom->GetNewBondVector(v,bondlen);
                        let v = MKBuilder.sharedInstance.getNewBondVector(atom, bondlen)
                        if v.x.isFinite || v.y.isFinite || v.z.isFinite { // MARK: Why are we only ensuring one??
                            self._c[self.numAtoms()*3]       = v.x
                            self._c[(self.numAtoms()*3) + 1] = v.y
                            self._c[(self.numAtoms()*3) + 2] = v.z
                        } else {
                            self._c[self.numAtoms()*3]       = 0.0
                            self._c[(self.numAtoms()*3) + 1] = 0.0
                            self._c[(self.numAtoms()*3) + 2] = 0.0
//                          MARK: Throw error here, or log error
                            print("AddHydrogens -- no reasonable bond geometry for desired hydrogen")
                            badh+=1
                        }
                    } else {
                        self._c.reserveCapacity(self.numAtoms()+3)
                    }
                }
                if badh == 0 || badh < self.numConformers() {
                    // Add the new H atom to the appropriate residue list
                    //but avoid doing perception by checking for existence of residue
                    //just in case perception is trigger, make sure GetResidue is called
                    //before adding the hydrogen to the molecule
                    let h = self.newAtom()
                    h.setType("H")
                    h.setAtomicNum(1)
                    let aname = "H"
                    
                    if let res = atom.getResidue() {
                        res.addAtom(h)
                        res.setAtomID(h, aname)
                        //hydrogen should inherit hetatm status of heteroatom (default is false)
                        if res.isHetAtom(atom) {
                            res.setHtmAtom(h, true)
                        }
                    }
                    
                    let bondFlags = 0
//                    MARK: Error catch here
                    _ = self.addBond(atom.getIdx(), h.getIdx(), 1, bondFlags)
                    h.setCoordPtr(self._c)
                    implicitRefToStereo(self, .Ref(atom.getId().rawValue), .Ref(h.getId().rawValue))
                }
            }
        }
        
        self.decrementMod()
        
        //reset atom type and partial charge flags
        self._flags &= (~(OB_PCHARGE_MOL|OB_ATOMTYPES_MOL|OB_SSSR_MOL|OB_AROMATIC_MOL|OB_HYBRID_MOL))

        return true
    }
    
    func addPolarHydrogens() -> Bool {
        return self.addNewHydrogens(.PolarHydrogen, false)
    }
    
    func addNonPolarHydrogens() -> Bool {
        return self.addNewHydrogens(.NonPolarHydrogen, false)
    }
    
    func addHydrogens(_ atom: MKAtom) -> Bool {
        let hcount = atom.getImplicitHCount()
        if hcount == 0 { return true }
        
        atom.setImplicitHCount(0)
        
        let vhadd: [Pair<MKAtom, UInt>] = [(atom, hcount)]
        
        //realloc memory in coordinate arrays for new hydroges
//        MARK: Not applicable to Swift...
        self.incrementMod()
        let hbrad = self.correctedBondRad(1, 0)
        for k in vhadd {
            let atom = k.0
            let bondlen = hbrad + self.correctedBondRad(atom.getAtomicNum(), atom.getHyb())
            for _ in 0..<k.1 {
                for n in 0..<self.numConformers() {
                    self.setConformer(n)
                    let v = MKBuilder.sharedInstance.getNewBondVector(atom, bondlen)
                    self._c[self.numAtoms()*3]       = v.x
                    self._c[(self.numAtoms()*3) + 1] = v.y
                    self._c[(self.numAtoms()*3) + 2] = v.z
                }
                let h = self.newAtom()
                h.setType("H")
                h.setAtomicNum(1)
                
                let bondFlags = 0
//                    MARK: Error catch here
                _ = self.addBond(atom.getIdx(), h.getIdx(), 1, bondFlags)
                h.setCoordPtr(self._c)
                implicitRefToStereo(self, .Ref(atom.getId().rawValue), .Ref(h.getId().rawValue))
            }
        }

        self.decrementMod()
        self.setConformer(0)
        //reset atom type and partial charge flags
        //_flags &= (~(OB_PCHARGE_MOL|OB_ATOMTYPES_MOL));
        return true
    }
    
    // Used by DeleteAtom below. Code based on StereoRefToImplicit
    static func deleteStereoOnAtom(_ mol: MKMol, _ atomId: Ref) {
        guard let vdata = mol.getDataVector(.StereoData) else { return }
        for dat in vdata {
            let dataType: MKStereo.TType = (dat as! MKStereoBase).getType()
            if dataType != .CisTrans && dataType != .Tetrahedral {
                print("This function should be updated to handle additional stereo types.\nSome stereochemistry objects may contain explicit refs to hydrogens which have been removed.")
                continue
            }
            
            if dataType == .CisTrans {
                let ct = (dat as? MKCisTransStereo)
                let ct_cfg: MKCisTransStereo.Config? = ct?.getConfig()
//                let ct_cfg = ct.getConfig()
//                TODO: IMPL
            } else if dataType == .Tetrahedral {
//                TODO: IMPL
            }
        }
        
    }
    
    //    MARK: HAS/IS Functions

    func has2D(_ not3D: Bool = false) -> Bool {
        var hasX, hasY: Bool
        hasX = false
        hasY = false
        for atom in self.getAtomIterator() {
            if (!hasX && !isNearZero(atom.getX())) {
                hasX = true
            }
            if (!hasY && !isNearZero(atom.getY())) {
                hasY = true
            }
            if (not3D && atom.getZ() != 0.0) {
                return false
            }
        }
        if hasX || hasY {
            return true
        }
        return false
    }
    
    func has3D() -> Bool {
        var hasX, hasY, hasZ: Bool
        hasX = false
        hasY = false
        hasZ = false
        for atom in self.getAtomIterator() {
            if (!hasX && !isNearZero(atom.getX())) {
                hasX = true
            }
            if (!hasY && !isNearZero(atom.getY())) {
                hasY = true
            }
            if (!hasZ && !isNearZero(atom.getZ())) {
                hasZ = true
            }
            if hasX && hasY && hasZ {
                return true
            }
        }
        return false
    }
    
    func hasNonZeroCoords() -> Bool {
        for atom in self.getAtomIterator() {
            if atom.getVector() == VZero {
                return false
            }
        }
        return true
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
    
    func isCorrectedForPH() -> Bool {
        return self.hasFlag(OB_PH_CORRECTED_MOL)
    }
    
//    MARK: Set Functions

    //! Set the heat of formation for this molecule (in kcal/mol)
    func setEnergy(_ energy: Double) {
        self._energy = energy
    }

    //! Set the dimension of this molecule (i.e., 0, 1 , 2, 3)
    func setDimension(_ dim: UInt16) {
        self._dimension = dim
    }

    //! Set the flag for determining automatic formal charges with pH (default=true)
    func setAutomaticFormalCharge(_ b: Bool) {
        self._autoFormalCharge = b
    }

    //! Set the flag for determining partial charges automatically (default=true)
    func setAutomaticPartialCharge(_ b: Bool) {
        self._autoPartialCharge = b
    }

    //! Mark that aromaticity has been perceived for this molecule (see OBAromaticTyper)
    func setAromaticPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_AROMATIC_MOL, value)
    }
    
    //! Mark that Smallest Set of Smallest Rings has been run (see OBRing class)
    func setSSSRPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_SSSR_MOL, value)
    }
    
    //! Mark that Largest Set of Smallest Rings has been run (see OBRing class)
    func setLSSRPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_LSSR_MOL, value)
    }
    
    //! Mark that rings have been perceived (see OBRing class for details)
    func setRingAtomsAndBondsPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_RINGFLAGS_MOL, value)
    }

    //! Mark that atom types have been perceived (see OBAtomTyper for details)
    func setAtomTypesPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_ATOMTYPES_MOL, value)
    }

    //! Mark that ring types have been perceived (see OBRingTyper for details)
    func setRingTypesPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_RINGTYPES_MOL, value)
    }

    //! Mark that chains and residues have been perceived (see OBChainsParser)
    func setChainsPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_CHAINS_MOL, value)
    }

    //! Mark that chirality has been perceived
    func setChiralityPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_CHIRALITY_MOL, value)
    }
    
    //! Mark that partial charges have been assigned
    func setPartialChargesPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_PCHARGE_MOL, value)
    }

    //! Mark that hybridization of all atoms has been assigned
    func setHybridizationPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_HYBRID_MOL, value)
    }

    //! Mark that ring closure bonds have been assigned by graph traversal
    func setClosureBondsPerceived(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_CLOSURE_MOL, value)
    }

    func setHydrogensAdded(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_H_ADDED_MOL, value)
    }

    func setCorrectedForPH(_ value: Bool = true) {
        self.set_or_unsetFlag(OB_PH_CORRECTED_MOL, value)
    }

    //! The OBMol is a pattern, not a complete molecule. Left unchanged by Clear()
    func setSpinMultiplicityAssigned(_ value: Bool) {
        self.set_or_unsetFlag(OB_ATOMSPIN_MOL, value)
    }

    func setIsPattenStructure(_ value: Bool) {
        self.set_or_unsetFlag(OB_PATTERN_STRUCTURE, value)
    }

    func setIsReaction(_ value: Bool) {
        self.set_or_unsetFlag(OB_REACTION_MOL, value)
    }

    //! Mark that distance calculations, etc., should apply periodic boundary conditions through the minimimum image convention.
    //! Does not automatically recalculate bonding.
    func setPeriodicMol(_ value: Bool) {
        self.set_or_unsetFlag(OB_PERIODIC_MOL, value)
    }
    
    func setConformer(_ i: Int) {
        if i < self._vconf.count {
            self._c = self._vconf[i].compactMap({ $0 })
        }
    }

    func setCoordinates(_ newCoords: Array<Double>) {
        var noCptr = self._c.isEmpty
        if noCptr {
            self._c.reserveCapacity(3)
        }
        
        // copy from external to internal, swift makes this easy...
        self._c = newCoords
        
        if noCptr {
            for atom in self.getAtomIterator() {
                atom.setCoordPtr(self._c)
            }
            self._vconf.append(newCoords)
        }
    }
    
    //! Renumber the atoms according to the order of indexes in the supplied vector
    //! This with assemble an atom vector and call RenumberAtoms(vector<OBAtom*>)
    //! It will return without action if the supplied vector is empty or does not
    //! have the same number of atoms as the molecule.
    //!
    //! \since version 2.3
    func renumberAtoms(_ v: [Int]) {
        
        if self.isEmpty() || v.count != self.numAtoms() {
            return
        }
        var va: [MKAtom] = []
        va.reserveCapacity(self.numAtoms())
        for i in v {
            guard let at = self.getAtom(i) else { continue }
            va.append(at)
        }
        self.renumberAtoms(va)
    }
    
    //! Renumber the atoms in this molecule according to the order in the supplied
    //! vector. This will return without action if the supplied vector is empty or
    //! does not have the same number of atoms as the molecule.
    func renumberAtoms(_ v: [MKAtom]) {
        if self.isEmpty() { return }
        var va = v
        if va.isEmpty || va.count != self.numAtoms() { return }
        
        var bv = MKBitVec()
        
        for i in va {
            bv |= i.getIdx()
        }
        
        for atom in self.getAtomIterator() {
            if !bv[atom.getIdx()] {
                va.append(atom)
            }
        }
        
        // double *c;
        // double *ctmp = new double [NumAtoms()*3];

        // for (j = 0;j < NumConformers();++j)
        // {
        //     c = GetConformer(j);
        //     for (k=0,i = va.begin();i != va.end(); ++i,++k)
        //     memcpy((char*)&ctmp[k*3],(char*)&c[((OBAtom*)*i)->GetCoordinateIdx()],sizeof(double)*3);
        //     memcpy((char*)c,(char*)ctmp,sizeof(double)*3*NumAtoms());
        // }
        // This is the C++ way of only keeping the coordinates from atoms in the va vector 
        // in swift this is not necessary, and we can make this much shorter and easier to read

        for j in 0..<self.numConformers() {
            var c = self.getConformer(j)
            var k = 0
            for i in va {
                c[k*3] = c[Int(i.getCoordinateIdx())]
                c[k*3+1] = c[Int(i.getCoordinateIdx())+1]
                c[k*3+2] = c[Int(i.getCoordinateIdx())+2]
                k += 1
            }
        }
        // this reads in swift as for each conformer, extract the coordinates 
        // for each atom in the va vector, and copy their coordinates back into the c array
        // in the necessary index, whether that is their original index or not, could be quite
        // expensive and does need to be replaced in the future 

        var k = 0
        for i in va {
            i.setIdx(k)
            k += 1
        }
        
        self._vatom.removeAll()
        for i in va {
            self._vatom.append(i)
        }

        self.deleteData(.RingData)
        self.deleteData("OpenBabel Symmetry Classes")
        self.deleteData("LSSR")
        self.deleteData("SSSR")
        self.unsetFlag(OB_LSSR_MOL)
        self.unsetFlag(OB_SSSR_MOL)
    }
    
    //check that unreasonable bonds aren't being added
    func validAdditionalBond(_ a: MKAtom, _ n: MKAtom) -> Bool {
        if a.getExplicitValence() == 5 && a.getAtomicNum() == 15 {
            //only allow octhedral bonding for F and Cl
            if n.getAtomicNum() == 9 || n.getAtomicNum() == 17 {
                return true
            } else {
                return false
            }
        }
        //other things to check?
        return true
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
        if self.hasFlag(OB_RINGFLAGS_MOL) { return }
//        TODO: Catch this int response
        findRingAtomsAndBonds2(self)
    }

    func findSSSR() {
        if self.hasSSSRPerceived() { return }
        self.setSSSRPerceived()
        
        // Delete any old data before we start finding new rings
        // The following procedure is slow
        // So if client code is multi-threaded, we don't want to make them wait
        
        if self.hasData("SSSR") {
            self.deleteData("SSSR")
        }
        
        //get Fr�rejacque taking int account multiple possible spanning graphs
        let frj = determineFRJ(self)
        if frj != 0 {
            var vr: [MKRing] = []
            self.findRingAtomsAndBonds()
            
            //restrict search for rings around closure bonds
            var cbonds: [MKBond] = self.getBondIterator().compactMap { bond in
                bond.isClosure() ? bond : nil
            }
            
            if !cbonds.isEmpty {
                let rs = MKRingSearch()
                //search for all rings about closures
                cbonds.forEach { bond in
                    rs.addRingFromClosure(self, bond)
                }
                
                rs.sortRings()
                rs.removeRedundant(frj)
                //store the SSSR set
                for j in rs._rlist {
                    let ring = MKRing(j._path, numAtoms()+1)
                    ring.setParent(self)
                    vr.append(ring)
                }
            }
            
            let rd = MKRingData()
            rd.setOrigin(.perceived)
            rd.setAttribute("SSSR")
            rd.setData(vr)
            self.setData(rd)
        }
    }
    
    func findLSSR() {
        if self.hasLSSRPerceived() { return }
        self.setLSSRPerceived()
        
        // Delete any old data before we start finding new rings
        // The following procedure is slow
        // So if client code is multi-threaded, we don't want to make them wait
        
        if self.hasData("LSSR") {
            self.deleteData("LSSR")
        }
        
        //get frerejaque taking int account multiple possible spanning graphs
        let frj = determineFRJ(self)
        if frj != 0 {
            var vr: [MKRing] = []
            self.findRingAtomsAndBonds()
            
            //restrict search for rings around closure bonds
            var cbonds: [MKBond] = self.getBondIterator().compactMap { bond in
                bond.isClosure() ? bond : nil
            }
            
            if !cbonds.isEmpty {
                let rs = MKRingSearch()
                //search for all rings about closures
                cbonds.forEach { bond in
                    rs.addRingFromClosure(self, bond)
                }
                
                rs.sortRings()
                rs.removeRedundant(-1) // -1 means LSSR
                //store the LSSR set
                for j in rs._rlist {
                    let ring = MKRing(j._path, numAtoms()+1)
                    ring.setParent(self)
                    vr.append(ring)
                }
            }
            
            let rd = MKRingData()
            rd.setOrigin(.perceived)
            rd.setAttribute("LSSR")
            rd.setData(vr)
            self.setData(rd)
        }
    }
    
//    MARK: Island of Misfit Toys

    /*! This method adds single bonds between all atoms
    closer than their combined atomic covalent radii,
    then "cleans up" making sure bonded atoms are not
    closer than 0.4A and the atom does not exceed its valence.
    It implements blue-obelisk:rebondFrom3DCoordinates.
    */
    func connectTheDots() {
        
        if self.isEmpty() { return }

        if _dimension != 3 {
            // TODO: throw exception
            print("connectTheDots() only works in 3D")
            return 
        }

        var atom : MKAtom
        var maxrad: Double = 0.0
        var zsortedAtoms: [Pair<MKAtom, Double>] = []
        var rad: [Double] = []
        var zsorted: [Int] = []
        var bondCount: [Int] = []
//        var unset = false

        var c: [Double] = []
        c.reserveCapacity(numAtoms()*3)
        rad.reserveCapacity(numAtoms())
        var j = 0
        for atom in self.getAtomIterator() {
            bondCount.append(atom.getExplicitDegree())
            //don't consider atoms with a full valance already
            //this is both for correctness (trust existing bonds) and performance
            if atom.getExplicitValence() >= MKElements.getMaxBonds(atom.getAtomicNum()) { continue }
            if atom.getAtomicNum() == 7 && atom.getFormalCharge() == 0 && atom.getExplicitValence() >= 3 { continue }
            c[j] = atom.getX()
            c[j+1] = atom.getY()
            c[j+2] = atom.getZ()
            zsortedAtoms.append((atom, atom.getZ()))
            j+=1
        }

        zsortedAtoms.sort(by: { $0.1 < $1.1 })
        
        let maxx = zsortedAtoms.count
        for j in 0..<maxx {
            atom = zsortedAtoms[j].0
            rad[j] = MKElements.getCovalentRad(atom.getAtomicNum())
            maxrad = max(rad[j], maxrad)
            zsorted.append(atom.getIdx()-1)
        }
        
        var idx1, idx2 : Int
        var d2, cutoff, zd : Double
        var atom1, atom2, wrapped_coords: Vector<Double>
        
        for j in 0..<maxx {
            let maxcutoff = square(rad[j] + maxrad + 0.45)
            idx1 = zsorted[j]
            for k in j+1..<maxx {
                idx2 = zsorted[k]
                // bonded if closer than elemental Rcov + tolerance
                cutoff = square(rad[j] + rad[k] + 0.45)
                // Use minimum image convention if the unit cell is periodic
                // Otherwise, use a simpler (faster) distance calculation based on raw coordinates
                if self.isPeriodic() {
                    atom1 = Vector<Double>(arrayLiteral: c[idx1*3], c[idx1*3+1], c[idx1*3+2])
                    atom2 = Vector<Double>(arrayLiteral: c[idx2*3], c[idx2*3+1], c[idx2*3+2])
                    let unitCell = self.getData(.UnitCell) as! MKUnitCell
                    wrapped_coords = unitCell.minimumImageCartesian(atom1 - atom2)
                    d2 = length_2(wrapped_coords)
                } else {
                    zd = square(c[idx1*3+2] - c[idx2*3+2])
                    // bigger than max cutoff, which is determined using largest radius,
                    // not the radius of k (which might be small, ie H, and cause an early  termination)
                    // since we sort by z, anything beyond k will also fail
                    if zd > maxcutoff { break }
                    
                    d2 = square(c[idx1*3] - c[idx2*3])
                    if d2 > cutoff {
                        continue // x's bigger than cutoff
                    }
                    d2 += square(c[idx1*3+1] - c[idx2*3+1])
                    if d2 > cutoff {
                        continue // x^2 + y^2 bigger than cutoff
                    }
                    d2 += zd
                }
                
                if d2 > cutoff { continue }
                if d2 < 0.16 { // 0.4 * 0.4 = 0.16
                    continue
                }
                
                guard let atom = self.getAtom(idx1+1) else { continue }
                guard let nbr = self.getAtom(idx2+1) else { continue }
                
                if atom.isConnected(nbr) { continue }
                
                if !self.validAdditionalBond(atom, nbr) || !validAdditionalBond(nbr, atom) { continue }
                
                self.addBond(idx1+1, idx2+1, 1)
            }
        }
        
        // If between BeginModify and EndModify, coord pointers are NULL
        // setup molecule to handle current coordinates
        if self._c.isEmpty {
            self._c = c
            for atom in self.getAtomIterator() {
                atom.setCoordPtr(self._c)
            }
            self._vconf.append(c)
//            unset = true
        }
        
        // Cleanup -- delete long bonds that exceed max valence
        var bond, maxbond : MKBond
        var maxlength: Double
        var valCount: Int
        var changed: Bool
        
        self.beginModify() //prevent needless re-perception in DeleteBond
        
        for atom in self.getAtomIterator() {
            
            while atom.getExplicitValence() > MKElements.getMaxBonds(atom.getAtomicNum()) || atom.smallestBondAngle() < 45.0 {
                // no new bonds added for this atom, just skip it
                guard let bondIterator = atom.getBondIterator() else { break }
                
                guard var bond = bondIterator.next() else { break }
                maxbond = bond
                // Fix from Liu Zhiguo 2008-01-26
                // loop past any bonds
                // which existed before ConnectTheDots was called
                // (e.g., from PDB resdata.txt)
                valCount = 0
                while valCount < bondCount[atom.getIdx()-1] {
                    guard let bond = bondIterator.next() else { break }
                    maxbond = bond
                    valCount += 1
                }
                
                // delete bonds between hydrogens when over max valence
                if atom.getAtomicNum() == MKElements.getAtomicNum("H") {
                    let m = bondIterator
                    changed = false
                    for bond in m {
                        if bond.getNbrAtom(atom).getAtomicNum() == MKElements.getAtomicNum("H") {
                            self.deleteBond(bond)
                            changed = true
                            break
                        }
                    }
                    if changed {
                        continue
                    } else {
                        bond = maxbond
                    }
                }
                
                maxlength = maxbond.getLength()
                for bond in bondIterator {
                    if bond.getLength() > maxlength {
                        maxbond = bond
                        maxlength = bond.getLength()
                    }
                }
                self.deleteBond(maxbond)
            }
        }
        self.endModify()
        
//        apparently, we are supposed to clear _c is the unset flag is set and resize conformers to one less
//        TODO: find out if this matters
//        if unset { }
//        otherwise we just garbage clean the c tmp array, but swift does that for us
    }
    
    /*! This method uses bond angles and geometries from current
      connectivity to guess atom types and then filling empty valences
      with multiple bonds. It currently has a pass to detect some
      frequent functional groups. It still needs a pass to detect aromatic
      rings to "clean up."
      AssignSpinMultiplicity(true) is called at the end of the function. The true
      states that there are no implict hydrogens in the molecule.
    */
    func percieveBondOrders() {
        
        if self.isEmpty() { return }
        
        if _dimension != 3 { return }
        
        var atom, b, c: MKAtom
        var v1, v2: Vector<Double>
        var angle: Double
        
        //  BeginModify();

        // Pass 1: Assign estimated hybridization based on avg. angles
        for atom in self.getAtomIterator() {
            angle = atom.averageBondAngle()
            if angle > 155.0 {
                atom.setHyb(1)
            } else if angle <= 155.0 && angle > 115.0 {
                atom.setHyb(2)
            }
            
            // special case for imines
            if atom.getAtomicNum() == MKElements.getAtomicNum("N") &&
                atom.explicitHydrogenCount() == 1 &&
                atom.getExplicitDegree() == 2 &&
                angle > 109.5 {
                atom.setHyb(2)
            } else if atom.getAtomicNum() == MKElements.getAtomicNum("N") &&
                        atom.getExplicitDegree() == 2 &&
                        atom.isInRing() { // azete
                atom.setHyb(2)
            }
        } // pass 1
        
        // Make sure upcoming calls to GetHyb() don't kill these temporary values
        self.setHybridizationPerceived()
        
        // Pass 2: look for 5-member rings with torsions <= 7.5 degrees
        //         and 6-member rings with torsions <= 12 degrees
        //         (set all atoms with at least two bonds to sp2)
        
        var path: [Int]
        var torsions: Double = 0.0
        
        if !self.hasSSSRPerceived() {
            findSSSR()
        }
        
        var rlist: [MKRing] = getSSSR()
        for ring in rlist {
            if ring.size() == 5 {
                path = ring._path
                torsions = (abs(self.getTorsion(path[0], path[1], path[2], path[3])) +
                            abs(self.getTorsion(path[1], path[2], path[3], path[4])) +
                            abs(self.getTorsion(path[2], path[3], path[4], path[0])) +
                            abs(self.getTorsion(path[3], path[4], path[0], path[1])) +
                            abs(self.getTorsion(path[4], path[0], path[1], path[2]))
                            ) / 5.0
                if torsions <= 7.5 {
                    for ringAtom in 0..<path.count {
//                      TODO: No way this could possibly error...
                        guard let b = self.getAtom(path[ringAtom]) else { break }
                        if b.getExplicitDegree() == 2 || b.getExplicitDegree() == 3 {
                            b.setHyb(2)
                        }
                    }
                }
            } else if ring.size() == 6 {
                path = ring._path
                torsions = (abs(self.getTorsion(path[0], path[1], path[2], path[3])) +
                            abs(self.getTorsion(path[1], path[2], path[3], path[4])) +
                            abs(self.getTorsion(path[2], path[3], path[4], path[5])) +
                            abs(self.getTorsion(path[3], path[4], path[5], path[0])) +
                            abs(self.getTorsion(path[4], path[5], path[0], path[1])) +
                            abs(self.getTorsion(path[5], path[0], path[1], path[2]))
                            ) / 6.0
                if torsions <= 12.0 {
                    for ringAtom in 0..<path.count {
//                      TODO: No way this could possibly error...
                        guard let b = self.getAtom(path[ringAtom]) else { break }
                        if b.getExplicitDegree() == 2 || b.getExplicitDegree() == 3 {
                            b.setHyb(2)
                        }
                    }
                }
            }
        }
        
        // Pass 3: "Antialiasing" If an atom marked as sp hybrid isn't
        //          bonded to another or an sp2 hybrid isn't bonded
        //          to another (or terminal atoms in both cases)
        //          mark them to a lower hybridization for now
        
        var openNbr: Bool
        
        for atom in self.getAtomIterator() {
            if atom.getHyb() == 2 || atom.getHyb() == 1 {
                openNbr = false
                guard let nbrs = atom.getNbrAtomIterator() else { continue }
                for b in nbrs {
                    if b.getHyb() < 3 || b.getExplicitDegree() == 1 {
                        openNbr = true
                        break
                    }
                }
                if !openNbr && atom.getHyb() == 2 { atom.setHyb(3) }
                else if !openNbr && atom.getHyb() == 1 { atom.setHyb(2) }
            }
        } // pass 3
        
        // Pass 4: Check for known functional group patterns and assign bonds
        //         to the canonical form
        //      Currently we have explicit code to do this, but a "bond typer"
        //      is in progress to make it simpler to test and debug.
        MolKit._BondTyper.assignFunctionalGroupBonds(self)
        
        // Pass 5: Check for aromatic rings and assign bonds as appropriate
        // This is just a quick and dirty approximation that marks everything
        //  as potentially aromatic

        // This doesn't work perfectly, but it's pretty decent.
        //  Need to have a list of SMARTS patterns for common rings
        //  which would "break ties" on complicated multi-ring systems
        // (Most of the current problems lie in the interface with the
        //   Kekulize code anyway, not in marking everything as potentially aromatic)
        
        var needs_kekulization: Bool = false
        var typed: Bool
        var loop, loopSize: Int
        
        for ringit in rlist {
            typed = false
            loopSize = ringit.size()
            if loopSize == 5 || loopSize == 6 || loopSize == 7 {
                path = ringit._path
                for loop in 0..<loopSize {
                    guard let atom = self.getAtom(path[loop]) else { break }
                    if atom.hasBondOfOrder(2) || atom.hasBondOfOrder(3) || atom.getHyb() != 2 {
                        typed = true
                        break
                    }
                }
                
                if !typed {
                    for loop in 0..<loopSize {
                        self.getBond(path[loop], path[(loop+1) % loopSize])?.setAromatic()
                        needs_kekulization = true
                    }
                }
            }
        }
        
        // Kekulization is necessary if an aromatic bond is present
        if needs_kekulization {
            self.setAromaticPerceived()
            // First of all, set the atoms at the ends of the aromatic bonds to also
            // be aromatic. This information is required for OBKekulize.
            for bond in self.getBondIterator() {
                if bond.isAromatic() {
                    bond.getBeginAtom().setAromatic()
                    bond.getEndAtom().setAromatic()
                }
            }
            var ok: Bool = MKKekulize(self)
            if !ok {
//                TODO: Handle this error with more indepth error codes for sure
                print("Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders")
                print(self.getTitle())
            }
            self.setAromaticPerceived(false)
        }
        
        // Quick pass.. eliminate inter-ring sulfur atom multiple bonds
        // dkoes - I have removed this code - if double bonds are set,
        // we should trust them.  See pdb_ligands_sdf/4iph_1fj.sdf for
        // a case where the charge isn't set, but we break the molecule
        // if we remove the double bond.  Also, the previous code was
        // fragile - relying on the total mol charge being set.  If we
        // are going to do anything, we should "perceive" a formal charge
        // in the case of a ring sulfur with a double bond (thiopyrylium)
        
        // Pass 6: Assign remaining bond types, ordered by atom electronegativity
        
        var bondLength, testLength: Double
        var maxElNeg, shortestBond, currentElNeg: Double
        var maxx: Int
        var sortedAtoms: [Pair<MKAtom, Double>] = []
        
        for atom in self.getAtomIterator() {
            // if atoms have the same electronegativity, make sure those with shorter bonds
            // are handled first (helps with assignment of conjugated single/double bonds)
            shortestBond = 1.0e5
            guard let nbratoms = atom.getNbrAtomIterator() else { continue }
            for b in nbratoms {
                if b.getAtomicNum() != 1 {
                    shortestBond = min(shortestBond, atom.getBond(b)!.getLength())
                }
            }
            let entry = Pair<MKAtom, Double>(atom, MKElements.getElectroNeg(atom.getAtomicNum()) * 1.0e6 + shortestBond)
            sortedAtoms.append(entry)
        }
        
        sortedAtoms.sort { $0.1 < $1.1 }

        maxx = sortedAtoms.count
        
        for iter in 0..<maxx {
            let atom = sortedAtoms[iter].0
            
            // Possible sp-hybrids
            if (atom.getHyb() == 1 || atom.getExplicitDegree() == 1) && atom.getExplicitValence() + 2 <= MKElements.getMaxBonds(atom.getAtomicNum()) {
                // loop through the neighbors looking for a hybrid or terminal atom
                // (and pick the one with highest electronegativity first)
                // *or* pick a neighbor that's a terminal atom
                if atom.hasNonSingleBond() || (atom.getAtomicNum() == 7 && atom.getExplicitValence() + 2 > 3) { continue }
                
                maxElNeg = 0.0
                shortestBond = 5000.0
                var c: MKAtom?
                guard let nbratoms = atom.getNbrAtomIterator() else { continue }
                for b in nbratoms {
                    currentElNeg = MKElements.getElectroNeg(b.getAtomicNum())
                    if (b.getHyb() == 1 || b.getExplicitDegree() == 1) && (b.getExplicitValence() + 2 <= MKElements.getMaxBonds(b.getAtomicNum())) &&
                        (currentElNeg > maxElNeg || (isApprox(currentElNeg, to: maxElNeg, epsilon: 1.0e-6) && atom.getBond(b)!.getLength() < shortestBond)) {
                        if b.hasNonSingleBond() || (b.getAtomicNum() == 7 && b.getExplicitValence() + 2 > 3) { continue }
                        // Test terminal bonds against expected triple bond lengths
                        bondLength = atom.getBond(b)!.getLength()
                        if atom.getExplicitDegree() == 1 || b.getExplicitDegree() == 1 {
                            testLength = correctedBondRad(atom.getAtomicNum(), atom.getHyb()) +
                                         correctedBondRad(b.getAtomicNum(), b.getHyb())
                            if (bondLength > 0.9 * testLength) { continue } // too long, ignore it
                        }
                        shortestBond = bondLength
                        maxElNeg = MKElements.getElectroNeg(b.getAtomicNum())
                        c = b
                    }
                }
                if c != nil {
                    atom.getBond(c!)!.setBondOrder(3)
                }
            }
            // Possible sp2-hybrid atoms
            else if (atom.getHyb() == 2 || atom.getExplicitDegree() == 1) && atom.getExplicitValence() + 1 <= MKElements.getMaxBonds(atom.getAtomicNum()) {
                
                if atom.hasNonSingleBond() || (atom.getAtomicNum() == 7 && atom.getExplicitValence() + 1 > 3) { continue }
                // Don't build multiple bonds to ring sulfurs
                // except thiopyrylium
                if atom.isInRing() && atom.getAtomicNum() == 16 {
                    if _totalCharge > 1 && atom.getFormalCharge() == 0 {
                        atom.setFormalCharge(+1)
                    } else {
                        continue
                    }
                }
                
                maxElNeg = 0.0
                shortestBond = 5000.0
                var c: MKAtom?
                
                guard let nbratoms = atom.getNbrAtomIterator() else { continue }
                for b in nbratoms {
                    currentElNeg = MKElements.getElectroNeg(b.getAtomicNum())
                    
                    if (b.getHyb() == 2 || b.getExplicitDegree() == 1) && (b.getExplicitValence() + 1 <= MKElements.getMaxBonds(b.getAtomicNum())) &&
                        (currentElNeg > maxElNeg || (isApprox(currentElNeg, to: maxElNeg, epsilon: 1.0e-6) && self.getBond(atom, b)!.isDoubleBondGeometry())) {
                        
                        if b.hasNonSingleBond() || (b.getAtomicNum() == 7 && b.getExplicitValence() + 1 > 3) { continue }
                        
                        if b.isInRing() && b.getAtomicNum() == 16 {
                            if _totalCharge > 1 && b.getFormalCharge() == 0 {
                                b.setFormalCharge(+1)
                            } else {
                                continue
                            }
                        }
                        
                        // Test terminal bonds against expected double bond lengths
                        bondLength = atom.getBond(b)!.getLength()
                        if atom.getExplicitDegree() == 1 || b.getExplicitDegree() == 1 {
                            testLength = correctedBondRad(atom.getAtomicNum(), atom.getHyb()) +
                                         correctedBondRad(b.getAtomicNum(), b.getHyb())
                            if (bondLength > 0.93 * testLength) { continue } // too long, ignore it
                        }
                        
                        // OK, see if this is better than the previous choice
                        // If it's much shorter, pick it (e.g., fulvene)
                        // If they're close (0.1A) then prefer the bond in the ring
                        let difference = shortestBond - bondLength
                        if ((difference > 0.1) || ((difference > -0.01) && (!atom.isInRing() || (c == nil) || !c!.isInRing() || b.isInRing())
                                                   || (atom.isInRing() && (c != nil) && !c!.isInRing() && b.isInRing()))) {
                            shortestBond = bondLength
                            maxElNeg = MKElements.getElectroNeg(b.getAtomicNum())
                            c = b // save this atom for later use
                        }  // is this bond better than previous choices
                    }
                }// loop through neighbors
                
                if c != nil {
                    atom.getBond(c!)!.setBondOrder(2)
                }
                
            }
            
        } // pass 6
        
        // Now let the atom typer go to work again
        self._flags &= (~OB_HYBRID_MOL)
        self._flags &= (~OB_AROMATIC_MOL)
        self._flags &= (~OB_ATOMTYPES_MOL)
        //  EndModify(true); // "nuke" perceived data
        
        //Set _spinMultiplicity other than zero for atoms which are hydrogen
        //deficient and which have implicit valency definitions (essentially the
        //organic subset in SMILES). There are assumed to no implicit hydrogens.
        //AssignSpinMultiplicity(true); // TODO: sort out radicals
    }
    

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
        
        let mag = sqrt(pow(v2x, 2.0) + pow(v2y, 2.0) + pow(v2z, 2.0))
        
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
    
    func align(_ a1: MKAtom, _ a2: MKAtom, _ p1: Vector<Double>, _ p2: Vector<Double>) {
        
        var _children: [Int] = []
        //find which atoms to rotate
        self.findChildren(a1.getIdx(), a2.getIdx(), &_children)
        _children.append(a2.getIdx())
        
        //find the rotation vector and angle
        let v1 = p2 - p1
        let v2 = a2.getVector() - a1.getVector()
        let v3 = cross3x3(v1, v2)
        let angle = vector_angle(v1, v2)
        
        //find the rotation matrix
        var m: Matrix<Double> = Matrix(rows: 3, columns: 3, repeatedValue: 0.0)
        m.rotAboutAxisByAngle(v3, angle)
        
        //rotate atoms
        for _child in _children {
            guard let atom = self.getAtom(_child) else { continue }
            var v = atom.getVector()
            v -= a1.getVector()
            v = v * m //rotate the point
            v += p1   //translate the vector
            atom.setVector(v)
        }
        //set a1 = p1
        a1.setVector(p1)
    }
    
    func toInertialFrame() {
        for conf in 0..<self.numConformers() {
            self.toInertialFrame(conf)
        }
    }
    
    func toInertialFrame(_ conf: Int) {
        
        var center = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
        var m: Matrix<Double> = Matrix.init(rows: 3, columns: 3, repeatedValue: 0.0)
        var mass: Double = 0
        var x,y,z : Double
        
        self.setConformer(conf)
        
        //find center of mass
        for atom in self.getAtomIterator() {
            let mi = atom.getAtomicMass()
            center[0] += mi * atom.getX()
            center[1] += mi * atom.getY()
            center[2] += mi * atom.getZ()
            mass += mi
        }
        
        center[0] /= mass
        center[1] /= mass
        center[2] /= mass
        
        //calculate inertial tensor
        for atom in self.getAtomIterator() {
            let x = atom.getX() - center[0]
            let y = atom.getY() - center[1]
            let z = atom.getZ() - center[2]
            let mi = atom.getAtomicMass()

            m[0, 0] += mi * (y * y + z * z)
            m[0, 1] -= mi * x * y
            m[0, 2] -= mi * x * z
            // m[1][0] -= mi*x*y
            m[1, 1] += mi * (x * x + z * z)
            m[1, 2] -= mi * y * z
            // m[2][0] -= mi*x*z
            // m[2][1] -= mi*y*z
            m[2, 2] += mi * (x * x + y * y)
        }

        // Fill in the lower triangle using symmetry across the diagonal
        m[1, 0] = m[0, 1]
        m[2, 0] = m[0, 2]
        m[2, 1] = m[1, 2]
        
        /* find rotation matrix for moment of inertia */
        var rmat: Matrix<Double> = Matrix(rows: 3, columns: 3, repeatedValue: 0.0)
        mk_make_rmat(&m, &rmat)
        
        /* rotate all coordinates */
        var c = self.getConformer(conf)
        for i in 0..<self.numAtoms() {
            x = c[i*3]     - center[0]
            y = c[i*3 + 1] - center[1]
            z = c[i*3 + 2] - center[2]
            c[i*3]     = x * rmat[0, 0] + y * rmat[0, 1] + z * rmat[0, 2]
            c[i*3 + 1] = x * rmat[1, 0] + y * rmat[1, 1] + z * rmat[1, 2]
            c[i*3 + 2] = x * rmat[2, 0] + y * rmat[2, 1] + z * rmat[2, 2]
        }
    }
    
    
    func correctForPH(_ pH: Double) -> Bool {
        if self.isCorrectedForPH() {
            return true
        }
//      MARK: Should probably wrap this in a try/catch block
        MKPhModel.sharedInstance.correctForPH(self, pH)
        
        return true
    }

    //! \brief set spin multiplicity for H-deficient atoms
      /**
         If NoImplicitH is true then the molecule has no implicit hydrogens. Individual atoms
         on which ForceNoH() has been called also have no implicit hydrogens.
         If NoImplicitH is false (the default), then if there are any explicit hydrogens
         on an atom then they constitute all the hydrogen on that atom. However, a hydrogen
         atom with its _isotope!=0 is not considered explicit hydrogen for this purpose.
         In addition, an atom which has had ForceImplH()called for it is never considered
         hydrogen deficient, e.g. unbracketed atoms in SMILES.
         Any discrepancy with the expected atom valency is interpreted as the atom being a
         radical of some sort and iits _spinMultiplicity is set to 2 when it is one hydrogen short
         and 3 when it is two hydrogens short and similarly for greater hydrogen deficiency.

         So SMILES C[CH] is interpreted as methyl carbene, CC[H][H] as ethane, and CC[2H] as CH3CH2D.
      **/
    func assignSpinMultiplicity(_ noImplicitH: Bool) -> Bool {
        // TODO: The following functions simply returns true, as it has been made
        // redundant by changes to the handling of implicit hydrogens, and spin.
        // This needs to be sorted out properly at some point.
        return true
    }
    
    // Put the specified molecular charge on a single atom (which is expected for InChIFormat).
    // Assumes all the hydrogen is explicitly included in the molecule,
    // and that SetTotalCharge() has not been called. (This function is an alternative.)
    // Returns false if cannot assign all the charge.
    // Not robust in the general case, but see below for the more common simpler cases.
    func assignTotalChargeToAtoms(_ charge: Int) -> Bool {
        
        var extraCharge = charge - self.getTotalCharge()
        
        for atom in self.getAtomIterator() {
            let atomicNum = atom.getAtomicNum()
            if atomicNum == 1 { continue }
            let charge = atom.getFormalCharge()
            let bosum = atom.getExplicitValence()
            let totalValence = bosum + atom.getImplicitHCount()
            let typicalValence = getTypicalValence(UInt(atomicNum), bosum, charge)
            let diff = typicalValence - totalValence
            if diff != 0 {
                var c: Int
                if extraCharge == 0 {
                    c = diff > 0 ? -1 : +1 //e.g. CH3C(=O)O, NH4 respectively
                } else {
                    c = extraCharge < 0 ? -1 : 1
                }
                if totalValence == getTypicalValence(UInt(atomicNum), bosum, charge + c) {
                    atom.setFormalCharge(charge + c)
                    extraCharge -= c
                }
            }
        }
        
        if extraCharge != 0 {
//            MARK: Should return error here
            print("Unable to assign all the charge to atoms")
            return false
        }
        return true
    }
    /* These cases work ok
       original      charge  result
      [NH4]             +1   [NH4+]
      -C(=O)[O]         -1   -C(=O)[O-]
      -[CH2]            +1   -C[CH2+]
      -[CH2]            -1   -C[CH2-]
      [NH3]CC(=O)[O]     0   [NH3+]CC(=O)[O-]
      S(=O)(=O)([O])[O] -2   S(=O)(=O)([O-])[O-]
      [NH4].[Cl]         0   [NH4+].[Cl-]
      */
    
}
