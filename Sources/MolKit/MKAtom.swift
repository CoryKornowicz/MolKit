//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 1/31/23.
//

import Foundation
import simd

//ATOM Property Macros (flags)
//! Atom is in a 4-membered ring
public let OB_4RING_ATOM: UInt =     1<<1  
//! Atom is in a 3-membered ring
public let OB_3RING_ATOM: UInt =     1<<2
//! Atom is aromatic
public let OB_AROMATIC_ATOM: UInt =  1<<3
//! Atom is in a ring
public let OB_RING_ATOM: UInt =      1<<4
//! Atom is an electron donor
public let OB_DONOR_ATOM: UInt =     1<<7
//! Atom is an electron acceptor
public let OB_ACCEPTOR_ATOM: UInt =  1<<8

public let OBATOM_TYPE_LEN: UInt = 6

class MKAtom: MKBase {
    
    private var _parent: MKBase? = nil
    private var _residue: MKResidue? = nil                         //!< parent residue (if applicable)
    
    private var _id: MKBaseID = .NoId                              //!< unique id
    private var _idx: Int = 0                                      //!< unique node index (GetIdx(), SetIdx())
    
    private var _ele: Int = 0                                      //!< atomic number (type unsigned char to minimize space -- allows for 0..255 elements)
    private var _imph: UInt = 0                                     //!< number of implicit hydrogens
    
    private var _type: String = ""                                 //!< atomic type
    
    private var _c: SIMD3<Double> = simd_double3(0.0,0.0,0.0)     //!< coordinate array in double*
    
    private var _flags: UInt = 0                                   //!< bitwise flags (e.g. aromaticity)
    private var _hyb: Int = 0                                     //!< hybridization
    
    private var _isotope: UInt = 0                                 //!< isotope (0 = most abundant)
    private var _spinmultiplicity: Int = 0                         //!< atomic spin, e.g., 2 for radical  1 or 3 for carbene
    private var _fcharge: Int = 0                                  //!< formal charge
    private var _pcharge: Double = 0                               //!< partial charge
    
    private var _vbond: [MKBond]? = nil                            //!< bond array
    
    private var Visit: Bool? = nil
    
    override init() {
        super.init()
        self._parent = nil
        self.clear()
    }
    
    func getFlag() -> UInt {
        return self._flags
    }
    
    func setFlag(_ flag: UInt) {
        self._flags |= flag
    }
    
    func unsetFlag(_ flag: UInt) {
        self._flags &= ~flag
    }
    
    func hasFlag(_ flag: UInt) -> Bool {
        return (self._flags & flag) != 0
    }
    
    func setIdx(_ idx: Int) {
        self._idx = idx
    }
    
    func getIdx() -> Int {
        return self._idx
    }
    
    func setId(_ id: Int) {
        self._id = ._id(id)
    }
    
    func getId() -> MKBaseID {
        return self._id
    }
    
    func getIndex() -> Int {
        return self._idx - 1
    }
    
    func getAtomicMass() -> Double {
        if self._isotope == 0 {
            return MKElements.getMass(self._ele)
        } else {
            return self.getExactMass()
        }
    }
    
    func getExactMass() -> Double {
        return MKElements.getExactMass(self._ele, self._isotope)
    }
    
    func setType(_ type: String) {
        self._type = type
        if self._ele == 1 && type == "D" {
            self._isotope = 2
        }
    }
    
    func getType() -> String {
        guard let mol: MKMol = (self._parent as? MKMol) else { return self._type }
        if !mol.hasAtomTypesPerceived() {
            MKAtomTyper.sharedInstance.assignTypes(mol)
        }
        return self._type
    }
    
    func setHyb(_ hyb: Int) {
        self._hyb = hyb
    }
    
    func getHyb() -> Int {
        //hybridization is assigned when atoms are typed
        guard let mol: MKMol = (self._parent as? MKMol) else { return _hyb }
        if !mol.hasHybridizationPerceived() {
            MKAtomTyper.sharedInstance.assignHyb(mol)
        }

        return _hyb
    }
    
    func setAtomicNum(_ atomicNum: Int) {
        self._ele = atomicNum
    }
    
    func getAtomicNum() -> Int {
        return self._ele
    }
    
    func getHeavyDegree() -> UInt {
        // map reduce over neighbor atoms if they are not Hydrogens
        guard let bonds = self._vbond else { return 0 }
        return bonds.map({ $0.getNbrAtom(self).getAtomicNum() != MKElements.getAtomicNum("H") ? 1 : 0 }).reduce(0, +)
    }
    
    func getHeteroDegree() -> UInt {
        guard let bonds = self._vbond else { return 0 }
        return bonds.map({ $0.getNbrAtom(self).isHeteroatom() ? 1 : 0 }).reduce(0, +)
    }
    
    func setIsotope(_ isotope: UInt) {
        self._isotope = isotope
    }
    
    func getIsotope() -> UInt {
        return self._isotope
    }
    
    func setImplicitHCount(_ imphCount: UInt) {
        self._imph = imphCount
    }
    
    func getImplicitHCount() -> UInt {
        return self._imph
    }
    
    func getExplicitDegree() -> Int {
        return (self._vbond != nil) ? self._vbond!.count : 0
    }
    
    func getTotalDegree() -> UInt {
        guard let vbond_count = self._vbond?.count else { return self._imph }
        return UInt(vbond_count) + self._imph
    }
    
    func getExplicitValence() -> UInt {
        var bosum: UInt = 0
        guard let bonds: [MKBond] = self._vbond else { return 0 }
        for bond in bonds {
            bosum += bond.getBondOrder()
        }
        return bosum
    }
    
    func getTotalValence() -> UInt {
        return self.getExplicitValence() + self._imph
    }

    func setFormalCharge(_ fcharge: Int) {
        self._fcharge = fcharge
    }

    func getFormalCharge() -> Int {
        return self._fcharge
    }
    
    func setSpinMultiplicity(_ spin: Int) {
        self._spinmultiplicity = spin
    }
    
    func getSpinMultiplicity() -> Int {
        return self._spinmultiplicity
    }
    
    func setPartialCharge(_ pcharge: Double) {
        self._pcharge = pcharge
    }
    
    func getPartialCharge() -> Double {
        if let mol = self.getParent() {
            if !mol.automaticPartialCharge() {
                return self._pcharge
            }
            if !mol.hasPartialChargesPerceived() {
                for atom in mol.atoms {
                    atom.setPartialCharge(0.0)
                }
                MKPhModel.sharedInstance.assignSeedPartialCharge(mol)
                let gc: MKGastChrg = MKGastChrg()
                gc.assignPartialCharges(mol)
            }
            return self._pcharge
        } else {
            return self._pcharge
        }
    }

    func setVector(_ vec3: SIMD3<Double>) {
        self._c = vec3
    }

    func setVector(_ x: Double, _ y: Double, _ z: Double) {
        self._c = simd_double3(x,y,z)
    }

    // Keeping both for redundancy and semantics : possible remove later
    func getVector() -> SIMD3<Double> { return self._c }

    func getCoordinates() -> SIMD3<Double> { return self._c}
    
    func clearCoordPtr() {
        self._c = simd_double3(0.0,0.0,0.0)
    }

    func getX() -> Double {
        return self._c[0]
    }
    
    func getY() -> Double {
        return self._c[1]
    }
    
    func getZ() -> Double {
        return self._c[2]
    }
    
    func setResidue(_ res: MKResidue) {
        self._residue = res
    }

    func getResidue() -> MKResidue {
        guard let mol: MKMol = (self._parent as? MKMol) else { return MKResidue() }
        if !mol.hasChainsPerceived() {
            MKChainsParser.sharedInstance.perceiveChains(mol)
        }
        return self._residue!
    }

    func setParent(_ mol: MKMol) {
        self._parent = mol
    }
    
    func getParent() -> MKMol? {
        guard let parent = (self._parent as? MKMol) else {return nil}
        return parent
    }
    
    func setAromatic(_ value: Bool) {
        value ? self.setFlag(OB_AROMATIC_ATOM) : self.unsetFlag(OB_AROMATIC_ATOM)
    }

    func setInRing(_ value: Bool) {
        value ? self.setFlag(OB_RING_ATOM) : self.unsetFlag(OB_RING_ATOM)
    }
    
    // func getNewBondVector(vector3 &v,double length) {}

    func getBond(_ atom: MKAtom) -> MKBond? {
        guard let bonds = self._vbond else { return nil }
        for bond in bonds {
            if bond.getNbrAtom(self) === atom {
                return bond
            }
        }
        return nil
    }
    
    //! \return the distance to the atom defined by OBMol::GetAtom()
    // double GetDistance(int index);
    // //! \return the distance to the supplied OBAtom
    // double GetDistance(OBAtom*);
    // //! \return the distance to the coordinates of the supplied vector3
    // //! \since version 2.4
    // double GetDistance(vector3* v);
    // //! \return the angle defined by this atom -> b (vertex) -> c
    // double GetAngle(int b, int c);
    // //! \return the angle defined by this atom -> b (vertex) -> c
    // double GetAngle(OBAtom *b, OBAtom *c);
    
    //! \name Addition of residue/bond info. for an atom
    //@{
    
    //! If no residue has been set for this atom, create a new one
    // void NewResidue()
    // {
    //     if (!_residue)
    //     _residue = new OBResidue;
    // }
    // //! Add (set) the residue for this atom
    // void AddResidue(OBResidue *res) { SetResidue(res); }
    // //! Delete any residue associated with this atom
    // void DeleteResidue(){
    // if (_residue) {
    //     delete _residue;
    //     _residue = nullptr; // Make sure to clear that a residue existed
    // }
    // }
    // //! Add a bond to the internal list. Does not update the bond.
    // void AddBond(OBBond *bond) { _vbond.push_back(bond); }
    // //! \brief Insert @p bond into the internal list at the position from @p i
    // //! Does not modify the bond
    // void InsertBond(OBBondIterator &i, OBBond *bond)
    // {
    //     _vbond.insert(i, bond);
    // }
    // //! Find @p bond and remove it from the internal list. Does not update the bond.
    // bool DeleteBond(OBBond* bond);
    // //! Clear all bonding information in this atom (does not delete them)
    // void ClearBond() {_vbond.clear();}
    // //@}
    
    // //! \name Builder utilities
    // //@{
    
    // //! \brief If this is a hydrogen atom, transform into a methyl group
    // //! \return success or failure
    // bool HtoMethyl();
    // //! Change the hybridization of this atom and modify the geometry accordingly
    // //! \return success or failure
    // //! \deprecated This will be removed in future versions of Open Babel
    // OB_DEPRECATED bool SetHybAndGeom(int);
    // //@}
    
    //! \name Property information
    //! \return The number of oxygen atoms connected that only have one heavy valence
    func countFreeOxygens() -> UInt {
        guard let bonds = self._vbond else { return 0 }
        return bonds.map { bond in
            let neatom = bond.getNbrAtom(self)
            return (neatom.getAtomicNum() == MKElements.getAtomicNum("O") && neatom.getHeavyDegree() == 1) ? 1 : 0
        }.reduce(0, +)
    }

    //! \return The number of sulfur atoms connected that only have one heavy valence
    //! since version 2.4
    func countFreeSulfurs() -> UInt {
        guard let bonds = self._vbond else { return 0 }
        return bonds.map { bond in
            let neatom = bond.getNbrAtom(self)
            return (neatom.getAtomicNum() == MKElements.getAtomicNum("S") && neatom.getHeavyDegree() == 1) ? 1 : 0
        }.reduce(0, +)
    }
    
    //! \return The number of hydrogens explicitly bound to this atom, optionally excluding D,T and isotope explicitly set to 1
    func explicitHydrogenCount(_ excludeIsotopes: Bool = false) -> UInt {
        var numH: UInt = 0
        guard let neigh = self._vbond?.map({ $0.getNbrAtom(self) }) else { return numH }
        for ne in neigh {
            if (ne.getAtomicNum() == MKElements.getAtomicNum("H") && !(excludeIsotopes && ne.getIsotope() != 0)) {
                numH+=1
            }
        }
        return numH
    }
    // //! \return The number of rings that contain this atom
    // unsigned int  MemberOfRingCount()     const;
    // //! \return The size of the smallest ring that contains this atom (0 if not in a ring)
    // unsigned int  MemberOfRingSize()	  const;
    // //! \return The number of explicit ring connections to this atom
    // unsigned int  CountRingBonds() const;
    // //! \return The smallest angle of bonds to this atom
    // double	  SmallestBondAngle();
    // //! \return The average angle of bonds to this atom
    // double	  AverageBondAngle();
    // /** Lewis acid/base vacancies for this atom
    // *  \return A pair of integers, where first is acid count and second is base count
    // *  \since version 2.3
    // */
    // std::pair<int, int> LewisAcidBaseCounts() const;
    
    // //! \return Is there any residue information?
    // bool HasResidue()    { return(_residue != nullptr);    }
    
    // //! \return Is this a HETATM in a residue (returns false if not in a residue)
    // //! \since version 2.4
    func isHetAtom() -> Bool {
        guard let residue = self._residue else { return false }
        return residue.isHetAtom(self)
    }
    
    //! \return Is the specified element, as specified by atom number (see OBElement namespace)?
    func isElement(_ e: Int) -> Bool {
        return self._ele == e;
    }
    
    // //! \return Is the atom aromatic?
    func isAromatic() -> Bool {
        guard let mol = self.getParent() else { return false }
        if !mol.hasAromaticPerceived() {
            MKAromaticTyper.sharedInstance.assignAromaticFlags(mol)
        }
        if self.hasFlag(OB_AROMATIC_ATOM) {
            return true
        }
        return false
    }

    // //! \return Is the atom in a ring?
    func isInRing() -> Bool {
        guard let mol: MKMol = (self._parent as? MKMol) else {
            return false;
        }
        if !mol.hasRingAtomsAndBondsPerceived(){
            mol.findRingAtomsAndBonds()
        }

        if self.hasFlag(OB_RING_ATOM) {
            return true
        }

        return false
    }

    // //! \return Is the atom in a ring of a given size?
    // bool IsInRingSize(int) const;
    
    //! \return Is this atom an element in the 15th or 16th main groups
    //!  (i.e., N, O, P, S ...) ?
    func isHeteroatom() -> Bool {
        switch(self.getAtomicNum()) {
        case 7,8,15,16,33,34,51,52,83,84 :
            return true
        default :
            return false
        }
    }
    
    //! \return Is this atom directly connected to the supplied OBAtom?
    func isConnected(_ atom: MKAtom) -> Bool {
        if let bonds = self._vbond {
            for bond in bonds {
                if bond.getBeginAtom() == atom || bond.getEndAtom() == atom {
                    return true
                }
            }
        }
        return false
    }
    
    //! \return Is this atom related to the supplied OBAtom in
    //!  a 1,3 bonding pattern?
    func isOneThree(_ atom: MKAtom) -> Bool {
        guard let bonds = self._vbond else { return false }
        guard let a2Bonds = atom._vbond else { return false }
        
        for bond in bonds {
            for bond2 in a2Bonds {
                if bond.getNbrAtom(self) == bond2.getNbrAtom(atom) {
                    return true
                }
            }
        }
        return false
    }
    
    //! \return Is this atom related to the supplied OBAtom in
    //!  a 1,4 bonding pattern?
    func isOneFour(_ atom: MKAtom) -> Bool {
        guard let bonds = self._vbond else { return false }
        guard let a2Bonds = atom._vbond else { return false }
        
        for bond in bonds {
            for bond2 in a2Bonds {
                if bond.getNbrAtom(self).isConnected(bond2.getNbrAtom(atom)) {
                    return true
                }
            }
        }
        return false
    }
    
    // //! \return Is this atom an oxygen in a carboxyl (-CO2 or CO2H) group?
    func isCarboxylOxygen() -> Bool {
        if self.getAtomicNum() != MKElements.getAtomicNum("O") || self.getHeavyDegree() != 1 {
            return false
        }
        guard let bonds = self._vbond else { return false }
        let atom: MKAtom? = bonds.first { bond in
            bond.getNbrAtom(self).getAtomicNum() == MKElements.getAtomicNum("C")
        }?.getNbrAtom(self)
        
        if atom != nil {
            if (!(atom!.countFreeOxygens() == 2) && !(atom!.countFreeOxygens() == 1 && atom!.countFreeSulfurs() == 1)) {
                return false
            } else {
                return true
            }
        } else {
            return false
        }
    }
    
    // //! \return Is this atom an oxygen in a phosphate (R-PO3) group?
    func isPhosphateOxygen() -> Bool {
        if self.getAtomicNum() != MKElements.getAtomicNum("O") || self.getHeavyDegree() != 1 {
            return false
        }
        guard let bonds = self._vbond else { return false }
        let atom: MKAtom? = bonds.first { bond in
            bond.getNbrAtom(self).getAtomicNum() == MKElements.getAtomicNum("P")
        }?.getNbrAtom(self)
        
        if atom != nil {
            if atom!.countFreeOxygens() > 2 {
                return true
            } else {
                return false
            }
        } else {
            return false
        }
    }
    
    // //! \return Is this atom an oxygen in a sulfate (-SO3) group?
    func isSulfateOxygen() -> Bool {
        if self.getAtomicNum() != MKElements.getAtomicNum("O") || self.getHeavyDegree() != 1 {
            return false
        }
        guard let bonds = self._vbond else { return false }
        let atom: MKAtom? = bonds.first { bond in
            bond.getNbrAtom(self).getAtomicNum() == MKElements.getAtomicNum("S")
        }?.getNbrAtom(self)
        
        if atom != nil {
            if atom!.countFreeOxygens() < 3 {
                return false
            } else {
                return true
            }
        } else {
            return false
        }
    }
    
    // Helper function for IsHBondAcceptor
    static func isSulfoneOxygen(_ atom: MKAtom) -> Bool {
        if atom.getAtomicNum() != MKElements.getAtomicNum("O") || atom.getHeavyDegree() != 1 {
            return false
        }
        guard let bonds = atom._vbond else { return false }
        let atom2: MKAtom? = bonds.first { bond in
            bond.getNbrAtom(atom).getAtomicNum() == MKElements.getAtomicNum("S")
        }?.getNbrAtom(atom)
        
        if atom2 != nil {
            // check for sulfate
            if atom2!.countFreeOxygens() != 2 {
                return false
            } else {
                // check for sulfonamide
                guard let bonds2 = atom2!._vbond else { return true }
                let atom3: MKAtom? = bonds2.first { bond in
                    bond.getNbrAtom(atom2!).getAtomicNum() == MKElements.getAtomicNum("N")
                }?.getNbrAtom(atom2!)
                return (atom3 == nil) ? true : false
            }
        } else {
            return false
        }
    }

    // //! \return Is this atom an oxygen in a nitro (-NO2) group?
    func isNitroOxygen() -> Bool {
        if self.getAtomicNum() != MKElements.getAtomicNum("O") || self.getHeavyDegree() != 1 {
            return false
        }
        guard let bonds = self._vbond else { return false }
        let atom: MKAtom? = bonds.first { bond in
            bond.getNbrAtom(self).getAtomicNum() == MKElements.getAtomicNum("N")
        }?.getNbrAtom(self)
        
        if atom != nil {
            if atom!.countFreeOxygens() != 2 {
                return false
            } else {
                return true
            }
        } else {
            return false
        }
    }
    
    // //! \return Is this atom a nitrogen in an amide (-C(=O)NR2) group?
    func isAmideNitrogen() -> Bool {
        
        if self.getAtomicNum() != MKElements.getAtomicNum("N") {
            return false
        }
        
        guard let bonds = self._vbond else { return false }
        
        for bond in bonds {
            var nbatom = bond.getNbrAtom(self)
            guard let nb_bonds = nbatom._vbond else { continue }
            for nbbond in nb_bonds {
                if (nbbond.getBondOrder() == 2 &&
                    ((nbbond.getNbrAtom(nbatom).getAtomicNum() == 8) ||
                     (nbbond.getNbrAtom(nbatom).getAtomicNum() == 16))) {
                    return true
                }
            }
        }
        return false
    }
    
    // //! \return Is this atom a hydrogen connected to a polar atom
    // //!  (i.e., N, O, P, S)
    func isPolarHydrogen() -> Bool {
        if self._ele != 1 {
            return false
        }
        if let bonds = self._vbond {
            for bond in bonds {
                let at_num = bond.getNbrAtom(self).getAtomicNum()
                if (at_num == 7 || at_num == 8 || at_num == 15 || at_num == 16) {
                    return true
                }
            }
        }
        return false
    } 
    
    // //! \return Is this atom a hydrogen connected to a non-polar atom
    // //!  (i.e., C)
    func isNonPolarHydrogen() -> Bool {
        if self._ele != 1 {
            return false
        }
        if let bonds = self._vbond {
            for bond in bonds {
                let at_num = bond.getNbrAtom(self).getAtomicNum()
                if (at_num == 6) {
                    return true
                }
            }
        }
        return false
    }
    
    // //! \return Is this atom an aromatic nitrogen with at least one
    // //!  double bond to an oxygen atom
    func isAromaticNOxide() -> Bool {
        
        if (self.getAtomicNum() != MKElements.getAtomicNum("N")) || !self.isAromatic() {
            return false
        }

        guard let bonds = self._vbond else { return false }

        for bond in bonds {
            if bond.getNbrAtom(self).getAtomicNum() == MKElements.getAtomicNum("O") && !bond.isInRing() && bond.getBondOrder() == 2 {
                return true
            }
        }

        return false
    }
    
    // //! \return Is this atom chiral?
    func isChiral() -> Bool {
        guard let mol = self.getParent() else { return false }
        
    }
    
    // //! \return Is the atom part of a periodic unit cell?
    // bool IsPeriodic() const;
    
    //! \return Is this atom an axial atom in a ring
    func isAxial() -> Bool {
        var tors: Double = 0.0

        guard let mol: MKMol = (self._parent as? MKMol) else { return false }
        guard let neigh = self._vbond?.map({ $0.getNbrAtom(self) }) else { return false }

        for a in neigh {
            if a.getHyb() == 3 && a.isInRing() && !self.isInRing() {
                guard let neighB: [MKAtom] = a._vbond?.map({ $0.getNbrAtom(a) }) else { continue }
                for b in neighB {
                    if b != self && b.isInRing() && b.getHyb() == 3 {
                        guard let neighC: [MKAtom] = b._vbond?.map({ $0.getNbrAtom(b) }) else { continue }
                        for c in neighC {
                            if c != a && c.isInRing() {
                                tors = fabs(mol.getTorsion(self, a, b, c))
                                return (tors > 55.0 && tors < 75.0)
                            }
                        }
                    }
                }
            }
        }

        return false
    }
    
    // //! \return Is this atom a hydrogen-bond acceptor  (considering also atom surrounding)
    // new function, Stefano Forli
    // Incorporate ideas and data from Kubyni and others.
    // [1] Kubinyi, H. "Changing paradigms in drug discovery.
    //    "SPECIAL PUBLICATION-ROYAL SOCIETY OF CHEMISTRY 304.1 (2006): 219-232.
    //
    // [2] Kingsbury, Charles A. "Why are the Nitro and Sulfone
    //     Groups Poor Hydrogen Bonders?." (2015).
    //
    // [3] Per Restorp, Orion B. Berryman, Aaron C. Sather, Dariush Ajami
    //     and Julius Rebek Jr., Chem. Commun., 2009, 5692 DOI: 10.1039/b914171e
    //
    // [4] Dunitz, Taylor. "Organic fluorine hardly ever accepts
    //     hydrogen bonds." Chemistry-A European Journal 3.1 (1997): 83-92.
    //
    // This function has a finer grain than the original
    // implementation, checking also the neighbor atoms.
    // Accordingly to these rules, the function will return:
    //
    //    aliph-O-aliph ether   -> true   [1]
    //    hydroxy O-sp3         -> true   [1]
    //    aro-O-aliph ether     -> true   [1]
    //    ester O-sp2           -> true   [1]
    //    sulfate O (R-SO3)     -> true   [2]
    //    sulfoxyde O (R-SO-R)  -> true   [2]
    //    organoboron-F (R-BF3) -> true   [3]
    //    ester O-sp3           -> false  [1]
    //    sulfone (R1-SO2-R2 )  -> false  [2]
    //    aro-O-aro             -> false  [1]
    //    aromatic O            -> false  [1]
    //    O-nitro               -> false  [2]
    //    organic F (R-F)       -> false  [4]
    //
    func isHbondAcceptor() -> Bool {
        // oxygen
        if self._ele == 8 {
            // oxygen; this should likely be a separate function
            // something like IsHbondAcceptorOxygen() -- Possibly Look into this 
            var aroCount: UInt = 0
            // maybe could be a bool option in the function?
            // aromatic oxygen (furan) (NO)
            // sulfone (NO)
            if self.isNitroOxygen() || self.isAromatic() || MKAtom.isSulfoneOxygen(self) {
                return false
            }
            guard let neighs = self._vbond?.map({ $0.getNbrAtom(self) }) else { return true }
            for neigh in neighs {
                if neigh.isAromatic() {
                    aroCount += 1
                    if aroCount == 2 { // aromatic ether (aro-O-aro) (NO)
                        return false
                    }
                } else {
                    if neigh.getAtomicNum() == MKElements.getAtomicNum("H") {
                        return true // hydroxyl (YES)
                    } else {
                        guard let bond = neigh.getBond(self) else { continue }
                        if bond.isEster() && (!self.isCarboxylOxygen()) {
                            return false
                        }
                    }
                }
            }
            return true
        }
        // fluorine
        else if self._ele == 9 {
            guard let neighs = self._vbond?.map({ $0.getNbrAtom(self) }) else { return true }
            // organic fluorine (NO)
            for n in neighs {
                if n.getAtomicNum() != MKElements.getAtomicNum("C") {
                    return false
                }
            }
            return true
        }
        // N+ ions and sp2 hybrid N with 3 valences should not be Hbond acceptors
        else if self._ele == 7 {
            if (!((self.getExplicitDegree() == 4 && self.getHyb() == 3) || (self.getExplicitDegree() == 3 && self.getHyb() == 2))) {
                return true
            } else {
                return false
            }
        }
        // Changes from Paolo Tosco
        else if self._ele == 16 && self.getFormalCharge() == -1 {
            return true
        }
        // everything else 
        else { return false }
    }
    
    //! \return Is this atom a hydrogen-bond acceptor (old function)?
    func isHbondAcceptorSimple() -> Bool {
        if self._ele == 8 || self._ele == 9 {
            return true
        }
        
        if self._ele == 7 {
            // N+ ions and sp2 hybrid N with 3 valences should not be Hbond acceptors
            if !((self.getExplicitDegree() == 4 && self.getHyb() == 3) || (self.getExplicitDegree() == 3 && self.getHyb() == 2)) {
                return true
            }
        }
        
        if self._ele == 16 && self.getFormalCharge() == -1 {
            return true
        }
        return false
    }
    
    //! \return Is this atom a hydrogen-bond donor?
    func isHbondDonor() -> Bool {
        if ![7,8,9].contains(self._ele) {
            return false
        }
        // Reduce bond array to neighbors and loop
        guard let bonds = self._vbond else { return false }
        
        for bond in bonds {
            let neigh = bond.getNbrAtom(self)
            if neigh.getAtomicNum() == MKElements.getAtomicNum("H") {
                return true
            }
        }
        
        return false
    }
    
    //! \return Is this a hydrogen atom attached to a hydrogen-bond donor?
    func isHbondDonorH() -> Bool {
        guard let bonds = self._vbond else { return false }
        for bond in bonds {
            if bond.getNbrAtom(self).isHbondDonor() {
                return true
            }
        }
        return false
    }
    
    //! \return Is this atom a metal?
    //! \since version 2.4
    private let NMETALS = 78
    private let metals: Array<Int> = [ 3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,
                                       30,31,37,38,39,40,41,42,43,44,45,46,47,48,49,50,55,56,57,58,59,60,61,62,63,
                                       64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,87,88,89,90,91,
                                       92,93,94,95,96,97,98,99,100,101,102,103 ]
    func isMetal() -> Bool {
        return metals.contains(self._ele)
    }
    
    //! \return Whether a neighboring atom (alpha) has an unsaturated bond
    //!   to a third atom (beta).
    //! \param includePandS Whether to include phosphorus and sulfur neighbors
    //! in this determination (or to exclude them)
    /**     This can be sketched as follows
             '*'
               \
                a=b
        where a and b are the 'apha' and 'beta' atoms, respectively and '*'
        indicates the current atom.
    **/
    func hasAlphaBetaUnsat(_ includePandS: Bool = true) -> Bool {
        guard let neighs = self._vbond?.map({ $0.getNbrAtom(self) }) else { return false }
        
        for neigh in neighs {
            let ne_atomnum = neigh.getAtomicNum()
            if (includePandS || (ne_atomnum != MKElements.getAtomicNum("S") && ne_atomnum != MKElements.getAtomicNum("P"))) {
                guard let ne_bonds = neigh._vbond else { continue }
                for ne_bond in ne_bonds {
                    let ne_neigh = ne_bond.getNbrAtom(neigh)
                    if ne_neigh != self && (ne_bond.getBondOrder() == 2 || ne_bond.getBondOrder() == 3 || ne_bond.getBondOrder() == 5) {
                        return true
                    }
                }
            }
        }
        return false
    }

    // //! \return Whether this atom is connected to any bond with order == @p bo
    func hasBondOfOrder(_ bo: UInt) -> Bool {
        guard let bonds = self._vbond else { return false }
        for bond in bonds {
            if bond.getBondOrder() == bo {
                return true
            }
        }
        return false
    }

    // //! \return The count of bonds connected to this atom with order == @p bo
    func countBondsOfOrder(_ bo: UInt) -> Int {
        guard let bonds = self._vbond else { return 0 }
        var count = 0
        for bond in bonds {
            if bond.getBondOrder() == bo {
                count += 1
            }
        }
        return count
    }
    // //! \return The maximum bond order for this atom
    func highestBondOrder() -> UInt {
        guard let bonds = self._vbond else { return 0 }
        var max: UInt = 0
        for bond in bonds {
            if bond.getBondOrder() > max {
                max = bond.getBondOrder()
            }
        }
        return max
    }

    //! \return Whether this atom is connected to any bond with order != 1
    func hasNonSingleBond() -> Bool {
        guard let bonds = self._vbond else { return false }
        for bond in bonds {
            if bond.getBondOrder() != 1 {
                return true
            }
        }
        return false
    }

    //! \return Does this atom have a single bond
    func hasSingleBond() -> Bool { 
        return self.hasBondOfOrder(1)
    }

    // //! \return Does this atom have a double bond
    func hasDoubleBond() -> Bool {
        return self.hasBondOfOrder(2)
    }
    
    // //! \return Does this atom have an aromatic bond
    func hasAromaticBond() -> Bool {
        return self.hasBondOfOrder(5)
    }
    
    // //! \return Whether this atom matches the first atom in a given SMARTS pattern
    // bool MatchesSMARTS(const char *);
    
    
    static func == (lhs: MKAtom, rhs: MKAtom) -> Bool {
        return lhs._idx == rhs._idx
    }
    
    override func clear() {
        super.clear()
        self._flags = 0
        self._idx = 0
        self._hyb = 0
        self._ele = 0
        self._isotope = 0
        self._spinmultiplicity = 0 // CM 18 Sept 2003
        self._imph = 0
        self._fcharge = 0
        self._type = ""
        self._pcharge = 0.0
        self._vbond?.removeAll()
        self._residue = nil;
        self._id = .NoId
        self._c = simd_double3(0.0,0.0,0.0)
    }
    
    deinit {
        if let res = self._residue {
            _ = res.removeAtom(self)
        }
    }
    
}
