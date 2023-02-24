

import Foundation
import Surge
import simd

enum MKAtomFlag: UInt {
    case Aromatic = 2    // 1 << 1
    case Ring = 16       // 1 << 4
    case Closure = 1024  // 1 << 10
}

enum MKAtomStereoFlag: UInt {
    case Wedge = 4           // 1 << 2
    case Hash = 8            // 1 << 3
    case WedgeOrHash =  2048 // 1 << 11
    case CisOrTrans = 4096   // 1 << 12
}

public let OB_AROMATIC_BOND: UInt = 1<<1
public let OB_WEDGE_BOND: UInt = 1<<2
public let OB_HASH_BOND: UInt = 1<<3
public let OB_RING_BOND: UInt = 1<<4
public let OB_CLOSURE_BOND: UInt = 1<<10
public let OB_WEDGE_OR_HASH_BOND: UInt = 1<<11

class MKBond: MKBase {

    private var _idx: UInt = 0
    private var _parent: MKMol? = nil
    private var _begin: MKAtom? = nil
    private var _end: MKAtom? = nil
    private var _order: UInt = 0  // 0 = single, 1 = double, 2 = triple, 3 = aromatic
    private var _id: MKBaseID = ._id(generateUUID())

    override init() {
        super.init()
        self._idx = 0
        self._order = 0
        self._parent = nil
        self._begin = nil
        self._end = nil
    }

    func set(_ idx: Int, _ begin: MKAtom, _ end: MKAtom, _ order: UInt, _ flags: UInt) {
        self.setIdx(idx)
        self.setBegin(begin)
        self.setEnd(end)
        self.setBondOrder(order)
        self.setFlag(flags)
    }

    // TODO: Need to write implementations
    func getId() -> Int {
        return self._id.rawValue
    }

    func setId(_ id: Int) {
        self._id = ._id(id)
    }

    func getIdx() -> UInt {
        return self._idx
    }
    func setIdx(_ idx: Int) {
        self._idx = UInt(idx)
    }

    func getBondOrder() -> UInt { 
        return self._order
    }

    func setBondOrder(_ order: UInt) {
        self._order = order
    }

    func getBeginAtom() -> MKAtom {
        return self._begin!
    }

    func getBeginAtomIdx() -> Int {
        return self._begin!.getIdx()
    }

    func setBegin(_ begin: MKAtom) {
        self._begin = begin
    }

    func getEndAtom() -> MKAtom {
        return self._end!
    }

    func getEndAtomIdx() -> Int {
        return self._end!.getIdx()
    }

    func setEnd(_ end: MKAtom) {
        self._end = end
    }

    func getNbrAtom(_ atom: MKAtom) -> MKAtom {
        return self._begin! != atom ? self._begin! : self._end!
    }

    //! \return The index to the neighboring atom of @p ptr (i.e., the end if @p ptr is the start)
    /** \warning If @p ptr is not part of the bond, the beginning atom
        index will always be returned **/
    func getNbrAtomIdx(_ atom: MKAtom) -> Int {
        return self._begin! != atom ? self._begin!.getIdx() : self._end!.getIdx()
    }

    func setLength(_ fixed: MKAtom, _ len: Double) {

        guard let mol: MKMol = fixed.getParent() else {
            return
        }

        let a: Int = fixed.getIdx()
        let b: Int = self.getNbrAtom(fixed).getIdx()

        if a == b { return }

        var children: [Int] = []
        mol.findChildren(a, b, &children)
        children.append(b)

        var v1: Vector<Double> = self.getNbrAtom(fixed).getVector()
        let v2: Vector<Double> = fixed.getVector()
        var v3: Vector<Double> = v1 - v2

        if (isNearZero(length(v3), 1e-6)) {
            // set v3 to a random unit vector
            print("Warning: bond length is near zero")
            v3 = Vector<Double>.randomNormal(count: 3)
        } else {
            v3 = normalize(v3)
        }

        v3 *= len
        v3 += v2
        let v4: Vector<Double> = v3 - v1

        for i in children {
            v1 = mol.getAtom(i)!.getVector()
            v1 += v4 
            mol.getAtom(i)!.setVector(v1)
        }
    }

    func setLength(_ length: Double) {

        guard let begin = self._begin, let end = self._end else {
            return
        }

        let firstLength: Double = length + ((self.getLength() - length) / 2)

        self.setLength(begin, firstLength)
        self.setLength(end, length)
    }

    func getParent() -> MKMol? {
        return self._parent
    }

    func setParent(_ parent: MKMol) {
        self._parent = parent
    }

    func setAromatic(_ value: Bool = true) {
        value ? self.setFlag(OB_AROMATIC_BOND) : self.unsetFlag(OB_AROMATIC_BOND)
    }

    func setWedge(_ value: Bool = true) {
        value ? self.setFlag(MKAtomStereoFlag.Wedge.rawValue) : self.unsetFlag(MKAtomStereoFlag.Wedge.rawValue)
    }

    func setHash(_ value: Bool = true) {
        value ? self.setFlag(MKAtomStereoFlag.Hash.rawValue) : self.unsetFlag(MKAtomStereoFlag.Hash.rawValue)
    }

    func setWedgeOrHash(_ value: Bool = true) {
        value ? self.setFlag(MKAtomStereoFlag.WedgeOrHash.rawValue) : self.unsetFlag(MKAtomStereoFlag.WedgeOrHash.rawValue)
    }

    func setCisOrTrans(_ value: Bool = true) {
        value ? self.setFlag(MKAtomStereoFlag.CisOrTrans.rawValue) : self.unsetFlag(MKAtomStereoFlag.CisOrTrans.rawValue)
    }

    func setInRing(_ value: Bool = true) {
        value ? self.setFlag(OB_RING_BOND) : self.unsetFlag(OB_RING_BOND)
    }

    func setClosure(_ value: Bool = true) {
        value ? self.setFlag(OB_CLOSURE_BOND) : self.unsetFlag(OB_CLOSURE_BOND)
    }

      //! \return The expected "equilibrium" length based on the covalent radii and bond order
      /** Length is given in Angstroms **/
    func getEquibLength() -> Double {
        let beg = self._begin!
        let end = self._end!
        let length: Double = MKAtom.correctedBondRad(beg.getAtomicNum(), beg.getHyb()) +
                             MKAtom.correctedBondRad(end.getAtomicNum(), end.getHyb())
        if self.isAromatic() { return length * 0.93 }
        
        switch (self._order) {
        case 3:
            return length * 0.87
        case 2:
            return length * 0.91
        default:
            return length
        }
    }
      //! \return The current length of this bond in Angstroms
    func getLength() -> Double {
        if (!self.isPeriodic()) {
            return dist(self._begin!.getVector(), self._end!.getVector())
        } else {
            return 0.0
//            TODO: implement
//            OBMol *mol = (OBMol*)((OBBond*)this)->GetParent();
//            OBUnitCell *box = (OBUnitCell*)mol->GetData(OBGenericDataType::UnitCell);
//            return (box->MinimumImageCartesian(begin->GetVector() - end->GetVector())).length();
        }
    }
  
    func findSmallestRing() -> MKRing? {
        
        guard let mol = self.getParent() else { return nil }
        let ring_list: [MKRing] = mol.getSSSR()
        
        var min_size: Int = .max
        var min_ring: MKRing? = nil
        for ring in MKIterator<MKRing>(ring_list) {
            if ring.isMember(self) && ring.size() < min_size {
                min_size = ring.size()
                min_ring = ring
            }
        }
        return min_ring
    }
      //@}

      //! \name property request methods
      //@{
      //! \return Is the bond aromatic?
      //!  (Note that the two atoms of the bond may be aromatic,
      //!   but not the bond)
    func isAromatic() -> Bool {
        guard let mol = self.getParent() else { return false }
        if (!mol.hasAromaticPerceived()) {
            MolKit._AromTyper.assignAromaticFlags(mol)
        }
        return self.hasFlag(OB_AROMATIC_BOND)
    }
    
    //! \return Is the bond part of a ring?
    func isInRing() -> Bool {
        guard let mol = self.getParent() else { return false }
        if (!mol.hasClosureBondsPerceived()) {
            mol.findRingAtomsAndBonds()
        }
        return self.hasFlag(OB_RING_BOND)
    }

      //! \return Is the bond a rotatable bond?
      /**  Currently, this function classifies any bond with at least one heavy
           atom, no sp-hybrid atoms (e.g., a triple bond somewhere) not in a ring
           as a potential rotor if includeRingsBonds is false.  If true, rotors in
           rings with more than 3 atoms may be included. No other bond typing is attempted.
           For more detailed rotor detection, check the OBRotorList and
           OBRotorRules classes **/
    func isRotor(_ includeRingBonds: Bool = false) -> Bool {
        if self.getBondOrder() != 1 { return false }
        guard let begin = self._begin, let end = self._end else { return false }
        // not in a ring, or in a large ring
        // and if it's a ring, not sp2
        if let ring = self.findSmallestRing() {
            if !includeRingBonds { return false }
            if ring.size() <= 3 { return false }
            if (begin.getHyb() == 2 || end.getHyb() == 2) { return false }
        }
        
        if (begin.getHyb() == 1 || end.getHyb() == 1) { return false }
        // not just an -OH or -NH2, etc.
        // maybe we want to add this as an option
//        TODO: add option of adding rotatable flag
        //    rotatable = rotatable && ((_bgn->IsHeteroatom() || _bgn->GetHvyDegree() > 1)
        //                               && (_end->IsHeteroatom() || _end->GetHvyDegree() > 1) );
        return (begin.getHeavyDegree() > 1 && end.getHeavyDegree() > 1)
    }

    //! \return Is the bond within a periodic unit cell?
    func isPeriodic() -> Bool {
        guard let mol = self._parent else { return false }
        return mol.isPeriodic()
    }
    
    /** \return Is the bond an amide link (i.e., between a carbonyl C and a N)?
     No distinction is made between primary, secondary, and tertiary amides. **/
    func isAmide() -> Bool {

        guard let begin = self._begin, let end = self._end else {
            return false
        }
        var c,n : MKAtom
        if (begin.getAtomicNum() == 6 && end.getAtomicNum() == 7) {
            c = begin
            n = end
        } else if (begin.getAtomicNum() == 7 && end.getAtomicNum() == 6) {
            c = end
            n = begin
        } else {
            return false
        }

        if self.getBondOrder() != 1 { return false }
        if n.getTotalDegree() != 3 { return false }

        // Make sure c isCarbonyl
        guard let c_bonds = c.getBondIterator() else { return false }
        
        for bond in c_bonds {
            if bond.isCarbonyl() { return true }
        }
        
        return false
    }
    
    /** \return Is the bond a primary amide (i.e., between carbonyl C and a NH2)?
     In versions prior to 2.3, this function incorrectly identified secondary amides. **/
    func isPrimaryAmide() -> Bool {
        guard let begin = self._begin, let end = self._end else {
            return false
        }
        var c,n : MKAtom
        if (begin.getAtomicNum() == 6 && end.getAtomicNum() == 7) {
            c = begin
            n = end
        } else if (begin.getAtomicNum() == 7 && end.getAtomicNum() == 6) {
            c = end
            n = begin
        } else {
            return false
        }

        if self.getBondOrder() != 1 { return false }
        if n.getTotalDegree() != 3 { return false }
        
        if n.getHeavyDegree() != 1 { return false }
        
        // Make sure c isCarbonyl
        guard let c_bonds = c.getBondIterator() else { return false }
        
        for bond in c_bonds {
            if bond.isCarbonyl() { return true }
        }
        
        return false
    }
    
    /** \return Is the bond a secondary amide (i.e., between a carbonyl C and a NH1)?
     In versions prior to 2.3, this function incorrectly identified tertiary amides. **/
    func isSecondaryAmide() -> Bool {
        guard let begin = self._begin, let end = self._end else {
            return false
        }
        var c,n : MKAtom
        if (begin.getAtomicNum() == 6 && end.getAtomicNum() == 7) {
            c = begin
            n = end
        } else if (begin.getAtomicNum() == 7 && end.getAtomicNum() == 6) {
            c = end
            n = begin
        } else {
            return false
        }

        if self.getBondOrder() != 1 { return false }
        if n.getTotalDegree() != 3 { return false }
        
        if n.getHeavyDegree() != 2 { return false }
        
        // Make sure c isCarbonyl
        guard let c_bonds = c.getBondIterator() else { return false }
        
        for bond in c_bonds {
            if bond.isCarbonyl() { return true }
        }
        
        return false
    }
    
    //! \return Is the bond a teriary amide (i.e., between a carbonyl C and a NH0)?
    //!  \since version 2.3.
    func IsTertiaryAmide() -> Bool {
        guard let begin = self._begin, let end = self._end else {
            return false
        }
        var c,n : MKAtom
        if (begin.getAtomicNum() == 6 && end.getAtomicNum() == 7) {
            c = begin
            n = end
        } else if (begin.getAtomicNum() == 7 && end.getAtomicNum() == 6) {
            c = end
            n = begin
        } else {
            return false
        }

        if self.getBondOrder() != 1 { return false }
        if n.getTotalDegree() != 3 { return false }
        
        if n.getHeavyDegree() != 3 { return false }
        
        // Make sure c isCarbonyl
        guard let c_bonds = c.getBondIterator() else { return false }
        
        for bond in c_bonds {
            if bond.isCarbonyl() { return true }
        }
        
        return false
    }

    func isEster() -> Bool {
        guard let begin = self._begin, let end = self._end else {
            return false
        }
        var a : MKAtom
        if (begin.getAtomicNum() == 6 && end.getAtomicNum() == 8) {
            a = begin
        } else if (begin.getAtomicNum() == 8 && end.getAtomicNum() == 6) {
            a = end
        } else {
            return false
        }

        if self.getBondOrder() != 1 { return false }
        
        if let bonds = a.getBondIterator() {
            for bond in bonds {
                if bond.isCarbonyl() {
                    return true
                }
            }
        }
        return false
    }

    func isCarbonyl() -> Bool {
        if self.getBondOrder() != 2 { return false }
        guard let begin = self._begin, let end = self._end else { return false }
        if (begin.getAtomicNum() == 6 && end.getAtomicNum() == 8) ||
            (begin.getAtomicNum() == 8 && end.getAtomicNum() == 6) {
            return true
        }
        return false
    }

    func isClosure() -> Bool {
        guard let mol = self.getParent() else { return false }
        if (!mol.hasClosureBondsPerceived()) {
            mol.findRingAtomsAndBonds()
        }
        return self.hasFlag(OB_CLOSURE_BOND)
    }

    func isWedge() -> Bool {
        return self.hasFlag(OB_WEDGE_BOND)
    }

    func isHash() -> Bool {
        return self.hasFlag(OB_HASH_BOND)
    }

    func isWedgeOrHash() -> Bool {
        return self.hasFlag(MKAtomStereoFlag.WedgeOrHash.rawValue)
    }

    func isCisOrTrans() -> Bool { 
        return self.hasFlag(MKAtomStereoFlag.CisOrTrans.rawValue)
    }

    func isDoubleBondGeometry() -> Bool {
        guard let begin = self._begin, let end = self._end else { return false }
        guard let mol = self.getParent() else { return false }
        // We concentrate on sp2 atoms with valence up to 3 and ignore the rest (like sp1 or S,P)
        // As this is called from PerceiveBondOrders, GetHyb() may still be undefined. 
        if (begin.getHyb() == 1 && begin.getExplicitDegree()>3) || (end.getHyb() == 1 && end.getExplicitDegree()>3) {
            return true
        }

        if let beginNeighbors = begin.getNbrAtomIterator(), let endNeighbors = end.getNbrAtomIterator() {
            for bgnNeighbor in beginNeighbors {
                if bgnNeighbor != end {
                    for endNeighbor in endNeighbors {
                        if endNeighbor != begin {
                            let torsion = abs(mol.getTorsion(bgnNeighbor, begin, end, endNeighbor))
                            // >12&&<168 not enough
                            if ( torsion > 15.0  && torsion < 160.0 ) {
                                    // Geometry does not match a double bond
                                    return false
                                }
                        }
                    }
                }
            }
        }
        return true
    }

}
