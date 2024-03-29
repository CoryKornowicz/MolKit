//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 1/31/23.
//

import Foundation
import Algorithms
import Surge

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

// The number of valence electrons in a free atom
public let VALENCE: [Int] = [0,1,2,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,
                             11,12,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,1,2,
                             4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,
                             8,1,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,5,6,7,8,9,10,11,12]
// The number of electrons required to make up a full valence shell
public let SHELL: [Int] = [0,2,2,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,18,18,18,18,18,18,
                           18,18,18,18,8,8,8,8,8,8,8,8,18,18,18,18,18,18,18,18,18,18,8,
                           8,8,8,8,8,8,8,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,
                           18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,8,8,18,18,18,18,
                           18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18]


//! \return Is this atom a metal?
//! \since version 2.4
public let NMETALS = 78
public let METALS: Array<Int> = [ 3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,
                                   30,31,37,38,39,40,41,42,43,44,45,46,47,48,49,50,55,56,57,58,59,60,61,62,63,
                                   64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,87,88,89,90,91,
                                   92,93,94,95,96,97,98,99,100,101,102,103 ]

public class MKAtom: MKBase {
    
    private var _parent: MKBase? = nil
    private var _residue: MKResidue? = nil                         //!< parent residue (if applicable)
    
    private var _id: Ref = .Ref(generateUUID())               //!< unique id
    private var _idx: Int = 0                                      //!< unique node index (GetIdx(), SetIdx())
    
    private var _ele: Int = 0                                      //!< atomic number (type unsigned char to minimize space -- allows for 0..255 elements)
    private var _imph: UInt = 0                                    //!< number of implicit hydrogens
    
    private var _type: String = ""                                 //!< atomic type
    
    // TODO: private var _v: [Vector<Double>] a true vector of 3 coordinates while _c could __technically__ handle more than 3 dimensions
    private var _c: Array<Double> = []                             //!< coordinate array in double*
    private var _v: Vector<Double> = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)     //!< coordinate array in double*
    private var _cidx: UInt = 0                                    //!< coordinate index
    private var _hyb: Int = 0                                      //!< hybridization
    
    private var _isotope: UInt = 0                                 //!< isotope (0 = most abundant)
    private var _spinmultiplicity: Int = 0                         //!< atomic spin, e.g., 2 for radical  1 or 3 for carbene
    private var _fcharge: Int = 0                                  //!< formal charge
    private var _pcharge: Double = 0                               //!< partial charge
    
    private var _vbond: [MKBond]? = nil                            //!< bond array
    
    private var Visit: Bool? = nil
    
    public override init() {
        super.init()
        self._parent = nil
        self.clear()
    }
    
    func duplicate(_ other: MKAtom) {
        self._hyb = other._hyb
        self._ele = other._ele
        self._imph = other._imph
        self._isotope = other._isotope
        self._fcharge = other._fcharge
        self._pcharge = other._pcharge
        self._spinmultiplicity = other._spinmultiplicity
        self._type = other._type
        self._v = other.getVector()
        self._flags = other._flags
        self._residue = nil
        self._id = other._id

        self._vdata?.removeAll()
        
        guard let otherData = other._vdata else { return }
        
        for data in MKIterator<MKGenericData>(otherData) {
            self.setData(data)
        }
        
    }

    public func setIdx(_ idx: Int) {
        self._idx = idx
        self._cidx = UInt(max((idx - 1), 0) * 3)
    }
    
    public func getIdx() -> Int {
        return self._idx
    }
    
    public func setId(_ id: Int) {
        self._id = .Ref(id)
    }
    
    public func getId() -> Ref {
        return self._id
    }
    
    public func getIndex() -> Int {
        return self._idx - 1
    }
    
    public func getAtomicMass() -> Double {
        if self._isotope == 0 {
            return MKElements.getMass(self._ele)
        } else {
            return self.getExactMass()
        }
    }
    
    public func getExactMass() -> Double {
        return MKElements.getExactMass(self._ele, self._isotope)
    }
    
    public func setType(_ type: String) {
        self._type = type
        if self._ele == 1 && type == "D" {
            self._isotope = 2
        }
    }
    
    public func getType() -> String {
        guard let mol: MKMol = (self._parent as? MKMol) else { return self._type }
        if !mol.hasAtomTypesPerceived() {
            MolKit._AtomTyper.assignTypes(mol)
        }
        return self._type
    }
    
    public func setHyb(_ hyb: Int) {
        self._hyb = hyb
    }
    
    public func getHyb() -> Int {
        //hybridization is assigned when atoms are typed
        guard let mol: MKMol = (self._parent as? MKMol) else { return _hyb }
        if !mol.hasHybridizationPerceived() {
            MolKit._AtomTyper.assignHyb(mol)
        }

        return _hyb
    }
    
    public func setAtomicNum(_ atomicNum: Int) {
        self._ele = atomicNum
    }
    
    public func getAtomicNum() -> Int {
        return self._ele
    }
    
    func getHeavyDegree() -> UInt {
        // map reduce over neighbor atoms if they are not Hydrogens
        guard let bonds = self._vbond else { return 0 }
        return bonds.map({ $0.getNbrAtom(self).getAtomicNum() != MKElements.Hydrogen.atomicNum ? 1 : 0 }).reduce(0, +)
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
                for atom in mol.getAtomIterator() {
                    atom.setPartialCharge(0.0)
                }
                MolKit._PhModel.assignSeedPartialCharge(mol)
                let gc: MKGastChrg = MKGastChrg()
                gc.assignPartialCharges(mol)
            }
            return self._pcharge
        } else {
            return self._pcharge
        }
    }
    
    func setVector() {
        if !self._c.isEmpty {
            self._v[0] = self._c[Int(self._cidx)]
            self._v[1] = self._c[Int(self._cidx)+1]
            self._v[2] = self._c[Int(self._cidx)+2]
        }
    }

    func setVector(_ vec3: Vector<Double>) {
        if self._c.isEmpty {
            self._v = vec3
        } else {
            self._c[Int(self._cidx)] = vec3[0]
            self._c[Int(self._cidx)+1] = vec3[1]
            self._c[Int(self._cidx)+2] = vec3[2]
        }
    }

    func setVector(_ x: Double, _ y: Double, _ z: Double) {
        if self._c.isEmpty {
            self._v = Vector<Double>.init(scalars: [x,y,z])
        } else {
            self._c[Int(self._cidx)] = x
            self._c[Int(self._cidx)+1] = y
            self._c[Int(self._cidx)+2] = z
        }
    }

    func getVector() -> Vector<Double> {
        if self._c.isEmpty {
            return self._v
        }
        return Vector.init(dimensions: 3) { dim in self._c[Int(self._cidx) + dim] }
    }

    public func getCoordinate() -> Double? {
        if !self._c.isEmpty { return self._c[Int(self._cidx)] }
        return nil
    }

    func getCoordinateIdx() -> UInt { return self._cidx }

    func setCoordPtr( _ ptr: Array<Double>) {
        self._c = ptr
        self._cidx = UInt((self.getIdx() - 1) * 3)
    }
    
    func getNbrAtomIterator() -> MKIterator<MKAtom>? {
        guard let neigh = self._vbond?.map({ $0.getNbrAtom(self) }) else { return nil }
        return MKIterator<MKAtom>(neigh)
    }
    
    func getNbrAtomBondIterator() -> MKIterator<(MKAtom, MKBond)>? {
        guard let neigh = self._vbond?.map({ ($0.getNbrAtom(self), $0) }) else { return nil }
        return MKIterator<(MKAtom, MKBond)>(neigh)
    }
    
    func getBondIterator() -> MKIterator<MKBond>? {
        guard let bonds = self._vbond else { return nil }
        return MKIterator<MKBond>(bonds)
    }
    
    func getNumBonds() -> Int? {
        return self._vbond?.count
    }
    
    func clearCoordPtr() {
        self._c = []
        self._cidx = 0
    }

    func getX() -> Double {
        if !self._c.isEmpty { return self._c[0] }
        return self._v[0]
    }
    
    func getY() -> Double {
        if !self._c.isEmpty { return self._c[1] }
        return self._v[1]
    }
    
    func getZ() -> Double {
        if !self._c.isEmpty { return self._c[2] }
        return self._v[2]
    }
    
    func setResidue(_ res: MKResidue?) {
        self._residue = res
    }

    func getResidue() -> MKResidue? {
        if let mol: MKMol = (self._parent as? MKMol) {
            if !mol.hasChainsPerceived() {
                MolKit._ChainsParser.perceiveChains(mol)
            }
            return self._residue!
        } else if self.hasResidue() { 
            return self._residue!
        } else {
            return nil
        }
    }

    func setParent(_ mol: MKMol) {
        self._parent = mol
    }
    
    func getParent() -> MKMol? {
        guard let parent = (self._parent as? MKMol) else {return nil}
        return parent
    }
    
    func setAromatic(_ value: Bool = true) {
        value ? self.setFlag(OB_AROMATIC_ATOM) : self.unsetFlag(OB_AROMATIC_ATOM)
    }

    func setInRing(_ value: Bool = true) {
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
    func getDistance(_ index: Int) -> Double {
        guard let mol = self.getParent() else { return 0.0 }
        return self.getDistance(mol.getAtom(index)!)
    }

    // //! \return the distance to the supplied OBAtom
    func getDistance(_ atom: MKAtom) -> Double {
        if !self.isPeriodic() {
            return length((self.getVector() - atom.getVector()))
        } else {
            guard let data: MKUnitCell = self.getParent()?.getData(.UnitCell) as? MKUnitCell else { return 0.0 }
            return length(data.minimumImageCartesian((self.getVector() - atom.getVector())))
        }
    }

    // //! \return the distance to the coordinates of the supplied vector3
    // //! \since version 2.4
    func getDistance(_ v: Vector<Double>) -> Double {
        return length((self.getVector() - v))
    }

    // //! \return the angle defined by this atom -> b (vertex) -> c
    func getAngle(_ b: Int, _ c: Int) -> Double {
        guard let mol = self.getParent() else { return 0.0 }
        return self.getAngle(mol.getAtom(b)!, mol.getAtom(c)!)
    }

    // //! \return the angle defined by this atom -> b (vertex) -> c
    func getAngle(_ b: MKAtom, _ c: MKAtom) -> Double {
        
        var v1: Vector<Double> = self.getVector() - b.getVector()
        var v2: Vector<Double> = c.getVector() - b.getVector()
        
        if self.isPeriodic() {
            guard let data: MKUnitCell = self.getParent()?.getData(.UnitCell) as? MKUnitCell else { return 0.0 }
            v1 = data.minimumImageCartesian(v1)
            v2 = data.minimumImageCartesian(v2)
        }

        if (isNearZero(length(v1), 1.0e-3) || isNearZero(length(v2), 1.0e-3)) {
            return(0.0)
        }

        return vector_angle(v1, v2)
    }
    
    //! \name Addition of residue/bond info. for an atom
    //@{
    //! If no residue has been set for this atom, create a new one
    func newResidue() {
        if (self._residue == nil) {
            self._residue = MKResidue()
        }
    }

    // //! Add (set) the residue for this atom
    func addResidue(_ res: MKResidue) { 
        self.setResidue(res)
    }

    // //! Delete any residue associated with this atom
    func deleteResidue() {
        if (self._residue != nil) {
            self._residue = nil
        }
    }

    //! Add a bond to the internal list. Does not update the bond.
    func addBond(_ bond: MKBond) {
        if self._vbond != nil {
            self._vbond!.append(bond)
        } else {
            self._vbond = [bond]
        }
    }

    //! \brief Insert @p bond into the internal list at the position from @p i
    //! Does not modify the bond
    func insertBond(_ i: Int, _ bond: MKBond) {
        if self._vbond != nil {
            self._vbond!.insert(bond, at: i)
        } else {
            self._vbond = [bond]
        }
    }

    //! Find @p bond and remove it from the internal list. Does not update the bond.
    @discardableResult
    func deleteBond(_ bond: MKBond) -> Bool {
        guard var bonds = self._vbond else { return false }
        if let index = bonds.firstIndex(of: bond) {
            bonds.remove(at: index)
            return true
        }
        return false
    }

    //! Clear all bonding information in this atom (does not delete them)
    func clearBond() {
        guard var bonds = self._vbond else { self._vbond = []; return }
        bonds.removeAll()
    }
    //@}
    
    // //! \name Builder utilities
    // //@{
    // MARK: Apparently not used?? 
    func getNewBondVector(_ length: Double) -> Vector<Double> {
        return MKBuilder.getNewBondVector(self, length)
    } 

    static func correctedBondRad(_ ele: Int, _ hyb: Int) -> Double {
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

    // //! \brief If this is a hydrogen atom, transform into a methyl group
    // //! \return success or failure
    func hToMethyl() -> Bool {

        if self.getAtomicNum() != MKElements.Hydrogen.atomicNum {
            return false
        }

        var mol: MKMol? = self.getParent()
        
        if mol == nil {
            mol = MKMol(self)
            self.setParent(mol!)
        }

        mol!.beginModify()

        self.setAtomicNum(6)
        self.setType("C3")
        self.setHyb(3)

        let neigh: MKAtom? = self._vbond?.first?.getNbrAtom(self)
        let bond: MKBond? = self._vbond?.first

        if neigh == nil {
            mol!.endModify(false)
            return false
        }

        let br1: Double = MKAtom.correctedBondRad(6,3)
        var br2: Double = MKAtom.correctedBondRad(neigh!.getAtomicNum(), neigh!.getHyb())

        bond!.setLength(neigh!, br1 + br2)

        br2 = MKAtom.correctedBondRad(1,0)

        for _ in 0..<3 {
            let hatom = mol!.newAtom()
            hatom.setAtomicNum(1)
            hatom.setType("H")

            let v = self.getNewBondVector(br1+br2)
            hatom.setVector(v)
            mol!.addBond(self.getIdx(), mol!.numAtoms(), 1)
        }

        mol!.endModify(false)
        return true
    }

    static func applyRotMatToBond(_ mol: MKMol, _ m: Matrix<Double>, _ atom1: MKAtom, _ atom2: MKAtom) {
        
        var children: [Int] = []
        mol.findChildren(atom1.getIdx(), atom2.getIdx(), &children)
        children.append(atom2.getIdx())
        
        for i in children {
            var v: Vector<Double> = mol.getAtom(i)!.getVector()
            v -= atom1.getVector()
            v = v * m
            v += atom1.getVector()
            mol.getAtom(i)!.setVector(v)
        }
        
    }

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
            return (neatom.getAtomicNum() == MKElements.Oxygen.atomicNum && neatom.getHeavyDegree() == 1) ? 1 : 0
        }.reduce(0, +)
    }

    //! \return The number of sulfur atoms connected that only have one heavy valence
    //! since version 2.4
    func countFreeSulfurs() -> UInt {
        guard let bonds = self._vbond else { return 0 }
        return bonds.map { bond in
            let neatom = bond.getNbrAtom(self)
            return (neatom.getAtomicNum() == MKElements.Sulfur.atomicNum && neatom.getHeavyDegree() == 1) ? 1 : 0
        }.reduce(0, +)
    }
    
    //! \return The number of hydrogens explicitly bound to this atom, optionally excluding D,T and isotope explicitly set to 1
    func explicitHydrogenCount(_ excludeIsotopes: Bool = false) -> UInt {
        var numH: UInt = 0
        guard let neigh = self._vbond?.map({ $0.getNbrAtom(self) }) else { return numH }
        for ne in neigh {
            if (ne.getAtomicNum() == MKElements.Hydrogen.atomicNum && !(excludeIsotopes && ne.getIsotope() != 0)) {
                numH+=1
            }
        }
        return numH
    }

    //! \return The number of rings that contain this atom
    func memberOfRingCount() -> UInt {
        var count: UInt = 0
        guard var mol = self.getParent() else { return count }
        
        if (!mol.hasSSSRPerceived()) {
            mol.findSSSR()
        }

        if !self.isInRing() {
            return count
        }
        
        let rings = mol.getSSSR()
        
        for ring in rings {
            if ring.isInRing(self.getIdx()) {
                count+=1
            }
        }
        return count
    }

    //! \return The size of the smallest ring that contains this atom (0 if not in a ring)
    func memberOfRingSize()	-> Int {
        guard var mol = self.getParent() else { return 0 }
        
        if (!mol.hasSSSRPerceived()) {
            mol.findSSSR()
        }

        if !self.isInRing() {
            return 0
        }
        
        let rings = mol.getSSSR()
        
        for ring in rings {
            if ring.isInRing(self.getIdx()) {
                return ring.size()
            }
        }
        return 0
    }

    //! \return The number of explicit ring connections to this atom
    func countRingBonds() -> Int {
        guard let bonds = self._vbond else { return 0 }
        return bonds.filter({ $0.isInRing() }).count
    }

    //! \return The smallest angle of bonds to this atom
    func smallestBondAngle() -> Double {
        var minAngle: Double = 360.0

        guard let neigh_pairs = self._vbond?.map({ $0.getNbrAtom(self) }).adjacentPairs() else { return 0.0 }

        for pair in neigh_pairs {
            let degrees = pair.0.getAngle(self, pair.1)
            if degrees < minAngle {
                minAngle = degrees
            }
        }

        return minAngle
    }

    // //! \return The average angle of bonds to this atom
    func averageBondAngle() -> Double {
        var sumAngle: Double = 0.0
        var count: UInt = 0

        guard let neigh_pairs = self._vbond?.map({ $0.getNbrAtom(self) }).adjacentPairs() else { return 0.0 }

        for pair in neigh_pairs {
            sumAngle += pair.0.getAngle(self, pair.1)
            count+=1
        }

        return sumAngle / Double(count)
    }

    // /** Lewis acid/base vacancies for this atom
    // *  \return A pair of integers, where first is acid count and second is base count
    // *  \since version 2.3
    // */
    func lewisAcidBaseCounts() -> Pair<Int, Int> {
        var counts: Pair<Int, Int> = (0, 0)
        
        let N = self.getAtomicNum()

        if (N == 0 || N > 122) {
            return counts
        } else {
            let S = SHELL[N]
            let V = VALENCE[N]
            let C = self.getFormalCharge()
            let B = Int(self.getTotalValence())

            counts.0 = (S - V - B + C) / 2
            counts.1 = (V - B - C) / 2

            return counts
        }
    }
    
    // //! \return Is there any residue information?
    func hasResidue() -> Bool { 
        return self._residue != nil
    }
    
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
            MolKit._AromTyper.assignAromaticFlags(mol)
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
    func isInRingSize(_ size: Int) -> Bool {
        guard var mol = self.getParent() else { return false }
        if !mol.hasSSSRPerceived() {
            mol.findSSSR()
        }

        if !self.hasFlag(OB_RING_ATOM) {
            return false
        }

        let rings = mol.getSSSR()
        
        for ring in rings {
            if ring.isInRing(self.getIdx()) && ring.pathSize() == size {
                return true
            }
        }
        return false
    }
    
    //! \return Is this atom an element in the 15th or 16th main groups
    //!  (i.e., N, O, P, S ...) ?
    func isHeteroatom() -> Bool {
        switch(self.getAtomicNum()) {
        case 7, 8, 15, 16, 33, 34, 51, 52, 83, 84:
            return true
        default :
            return false
        }
    }
    
    //! \return Is this atom directly connected to the supplied OBAtom?
    func isConnected(_ atom: MKAtom) -> Bool {
        if let bonds = self._vbond {
            for bond in bonds {
                if bond.getBeginAtom().getIdx() == atom.getIdx() || bond.getEndAtom().getIdx() == atom.getIdx() {
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
                if bond.getNbrAtom(self).getIdx() == bond2.getNbrAtom(atom).getIdx() {
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
        if self.getAtomicNum() != MKElements.Oxygen.atomicNum {
            return false
        }
        
        if self.getHeavyDegree() != 1 {
            return false
        }
        
        guard let bonds = self._vbond else { return false }
        
        var atom: MKAtom?
        
        for bond in bonds {
            if bond.getNbrAtom(self).getAtomicNum() == MKElements.Carbon.atomicNum {
                atom = bond.getNbrAtom(self)
                break
            }
        }
        
        if atom != nil {
            if ( !(atom!.countFreeOxygens() == 2) &&
                 !(atom!.countFreeOxygens() == 1 && atom!.countFreeSulfurs() == 1) ) {
                return false
            }
            return true
        } else {
            return false
        }
    }
    
    // //! \return Is this atom an oxygen in a phosphate (R-PO3) group?
    func isPhosphateOxygen() -> Bool {
        if self.getAtomicNum() != MKElements.Oxygen.atomicNum {
            return false
        }
        if self.getHeavyDegree() != 1 {
            return false
        }
        
        guard let bonds = self._vbond else { return false }
        var atom: MKAtom?
        for bond in bonds {
            if bond.getNbrAtom(self).getAtomicNum() == MKElements.Phosphorus.atomicNum {
                atom = bond.getNbrAtom(self)
                break
            }
        }
            
        if atom != nil {
            if atom!.countFreeOxygens() > 2 {
                return true
            }
            return false
        } else {
            return false
        }
    }
    
    // //! \return Is this atom an oxygen in a sulfate (-SO3) group?
    func isSulfateOxygen() -> Bool {
        if self.getAtomicNum() != MKElements.Oxygen.atomicNum  {
            return false
        }
        
        if self.getHeavyDegree() != 1 {
            return false
        }
        
        guard let bonds = self._vbond else { return false }
        var atom: MKAtom?
        
        for bond in bonds {
            if bond.getNbrAtom(self).getAtomicNum() == MKElements.Sulfur.atomicNum {
                atom = bond.getNbrAtom(self)
                break
            }
        }
        
        if atom != nil {
            if atom!.countFreeOxygens() < 3 {
                return false
            }
            return true
        } else {
            return false
        }
    }
    
    // Helper function for IsHBondAcceptor
    func isSulfoneOxygen() -> Bool {
        if self.getAtomicNum() != MKElements.Oxygen.atomicNum {
            return false
        }
        
        if self.getHeavyDegree() != 1 {
            return false
        }
        
        guard let bonds = self._vbond else { return false }
        var nbr: MKAtom?
        
        for bond in bonds {
            if bond.getNbrAtom(self).getAtomicNum() == MKElements.Sulfur.atomicNum {
                nbr = bond.getNbrAtom(self)
                break
            }
        }
        
        if nbr != nil {
            // check for sulfate
            if nbr!.countFreeOxygens() != 2 {
                return false
            }
            // check for sulfonamide
            guard let bonds2 = nbr!._vbond else { return true }
            
            for bond in bonds2 {
                if bond.getNbrAtom(nbr!).getAtomicNum() == MKElements.Nitrogen.atomicNum {
                    return false
                }
            }
            return true
        } else {
            return false
        }
    }

    // //! \return Is this atom an oxygen in a nitro (-NO2) group?
    func isNitroOxygen() -> Bool {
        if self.getAtomicNum() != MKElements.Oxygen.atomicNum || self.getHeavyDegree() != 1 {
            return false
        }
        guard let bonds = self._vbond else { return false }
        let atom: MKAtom? = bonds.first { bond in
            bond.getNbrAtom(self).getAtomicNum() == MKElements.Nitrogen.atomicNum
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
        
        if self.getAtomicNum() != MKElements.Nitrogen.atomicNum {
            return false
        }
        
        guard let bonds = self._vbond else { return false }
        
        for bond in bonds {
            let nbatom = bond.getNbrAtom(self)
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
        
        if (self.getAtomicNum() != MKElements.Nitrogen.atomicNum) || !self.isAromatic() {
            return false
        }

        guard let bonds = self._vbond else { return false }

        for bond in bonds {
            if bond.getNbrAtom(self).getAtomicNum() == MKElements.Oxygen.atomicNum && !bond.isInRing() && bond.getBondOrder() == 2 {
                return true
            }
        }

        return false
    }
    
    // //! \return Is this atom chiral?
    func isChiral() -> Bool {
        guard let mol = self.getParent() else { return false }
        let stereoFacade = MKStereoFacade(mol)
        return stereoFacade.hasTetrahedralStereo(self._id.intValue)
    }
    
    // //! \return Is the atom part of a periodic unit cell?
    func isPeriodic() -> Bool {
        guard let mol = self.getParent() else { return false }
        return mol.isPeriodic()
    }
    
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
            if self.isNitroOxygen() || self.isAromatic() || self.isSulfoneOxygen() {
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
                    if neigh.getAtomicNum() == MKElements.Hydrogen.atomicNum {
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
            if neigh.getAtomicNum() == MKElements.Hydrogen.atomicNum {
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
    
    func isMetal() -> Bool {
        return METALS.contains(self._ele)
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
            if (includePandS || (ne_atomnum != MKElements.Sulfur.atomicNum && ne_atomnum != MKElements.Phosphorus.atomicNum)) {
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
        self._id = .Ref(generateUUID())
        self._c = []
        self._v = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
    }
    
    
    deinit {
        if let res = self._residue {
            res.removeAtom(self)
        }
    }
    
}
