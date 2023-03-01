

import Foundation
import Surge
import Collections

class MKRTree {
    
    var _atom: MKAtom
    var _prv: MKRTree? = nil
    
    init(_ atom: MKAtom, _ prv: MKRTree? = nil) {
        self._atom = atom
        self._prv = prv
    }
    
    func getAtomIdx() -> Int {
        return self._atom.getIdx()
    }
    
    func pathToRoot(_ path: inout [MKAtom]) {
        path.append(self._atom)
        if (self._prv != nil) {
            self._prv?.pathToRoot(&path)
        }
    }
}



class MKRing: MKBase {

    private var _parent: MKMol? = nil
    private var _type: String = ""

    var ring_id: Int = 0
    var _path: [Int] = []
    var _pathSet: MKBitVec = MKBitVec()
    
    override init() {
        super.init()
    }
    
    init(_ path: [Int], _ size: Int) {
        self._path = path
        self._pathSet.fromVecInt(_path)
        self._pathSet.resize(UInt32(size))
        self._type = ""
        
        super.init()
    }
    
    init(_ path: [Int], _ set: MKBitVec) {
        self._path = path
        self._pathSet = set
        self._type = ""
        super.init()
    }
    
    func getParent() -> MKMol? {
        return self._parent
    }
    
    func setParent(_ parent: MKMol) {
        self._parent = parent
    }
    
    func setType(_ type: String) {
        self._type = type
    }
    
    func getType() -> String {
        guard let mol = self._parent else { return "" }
        if !mol.hasRingTypesPerceived() {
            MolKit._RingTyper.assignTypes(mol)
        }
        return _type
    }
    
    func getRootAtom() -> Int {
        guard let mol = self._parent else { return 0 }
        
        if size() == 6 {
            for i in _path {
                guard let atom = mol.getAtom(i) else { return 0 }
                if atom.getAtomicNum() != MKElements.getAtomicNum("C") {
                    return i
                }
            }
        }
        
        if size() == 5 {
            for i in _path {
                guard let atom = mol.getAtom(i) else { return 0 }
                switch atom.getAtomicNum() {
                case MKElements.getAtomicNum("S"), MKElements.getAtomicNum("O"):
                    if atom.getExplicitDegree() == 2 {
                        return i
                    }
                case MKElements.getAtomicNum("N"):
                    if atom.getExplicitValence() == atom.getExplicitDegree() {
                        return i
                    }
                default: break
                }
            }
        }
        
        return 0
    }

    func size() -> Int {
        return _path.count
    }
    
    func pathSize() -> Int {
        return _path.count
    }

    //! \return Whether @p i as an atom index is in this ring
    func isInRing(_ i: Int) -> Bool {
       return _pathSet.bitIsSet(i)
    }
    
    func isMember(_ atom: MKAtom) -> Bool {
        return _pathSet.bitIsSet(atom.getIdx())
        
    }
    
    func isMember(_ bond: MKBond) -> Bool {
        return _pathSet.bitIsSet(bond.getBeginAtomIdx()) && _pathSet.bitIsSet(bond.getEndAtomIdx())
    }
    
    func isAromatic() -> Bool {
        guard let mol = self._parent else { return false }
        for i in _path {
            guard let chAtom = mol.getAtom(i) else { return false }
            if !chAtom.isAromatic() { return false }
        }
        return true
    }
    
    func findCenterAndNormal(center: inout Vector<Double>, norm1: inout Vector<Double>, norm2: inout Vector<Double>) -> Bool {
        
        guard let mol = self._parent else { return false }
        let nA = self._path.count
        var tmp: Vector<Double>
        
        center = [0.0,0.0,0.0]
        norm1 = [0.0,0.0,0.0]
        norm2 = [0.0,0.0,0.0]

        for j in 0..<nA {
            guard let atom = mol.getAtom(self._path[j]) else { return false }
            center += atom.getVector()
        }

        center /= Double(nA)

        for j in 0..<nA {
            guard let atom = mol.getAtom(self._path[j]) else { return false }
            guard let atom2 = mol.getAtom(self._path[j+1==nA ? 0 : j+1]) else { continue }
            let v1 = atom.getVector() - center
            let v2 = atom2.getVector() - center
            tmp = cross3x3(v1, v2)
            norm1 += tmp
        }
        
        norm1 /= Double(nA)
        norm1 = normalize(norm1)
        norm2 = norm1 
        norm2 *= -1.0
        return true
    }

}

class MKRingSearch {
    
//    var _bonds: [MKBond] ... Apparently deprecated
    var _rlist: [MKRing] = []
    
    init() {}
    
    func sortRings() {
        var ring_id = 0
        for ring in self._rlist {
            ring.ring_id = ring_id
            ring_id += 1
        }
        self._rlist = self._rlist.sorted { lring, rring in
            lring.size() == rring.size() ? lring.ring_id < rring.ring_id : lring.size() < rring.size()
        }
    }
    
    //! Starting with a full ring set - reduce to SSSR set
    func removeRedundant(_ frj: Int) {        
        //remove identical rings
        for i in (1..._rlist.count-1).reversed() {
            for j in (0...i-1).reversed() {
                if _rlist[i]._pathSet == _rlist[j]._pathSet {
                    _rlist.remove(at: i)
                    break
                }
            }
        }
        
        if _rlist.count == 0 { return } // nothing to do
        
        // handle LSSR
        if frj < 0 {
            guard let mol = _rlist[0].getParent() else { return }
            var rlist: [MKRing] = []
            var rignored: [MKRing] = []
            for i in 0..<_rlist.count {
                visitRing(mol, _rlist[i], &rlist, &rignored)
            }
            for i in 0..<rignored.count {
                // delete rignored[i]
                // in c++ this would free the pointer which arises from _rlist
                // in swift we will simply delete the ring from _rlist
                // TODO: this may have unforseen consequences
                let delRing = rignored[i]
                if let delindex = rlist.firstIndex(of: delRing) {
                    rlist.remove(at: delindex)
                }
            }
            self._rlist = rlist
            return
        }
        
        // exit if we already have frj rings
        if _rlist.count == frj { return }
        
        //make sure tmp is the same size as the rings
        var tmp: MKBitVec = MKBitVec()
        for j in 0..<_rlist.count {
            tmp = _rlist[j]._pathSet
        }
        
        //remove larger rings that cover the same atoms as smaller rings
        for i in (0..._rlist.count-1).reversed() {
            tmp.clear()
            for j in 0..<_rlist.count {
                if (_rlist[j]._path.count <= _rlist[i]._path.count) && (i != j) {
                    tmp |= _rlist[j]._pathSet
                }
            }
            
            tmp = tmp & _rlist[i]._pathSet
            
            if tmp == _rlist[i]._pathSet {
                _rlist.remove(at: i)
            }
            
            if _rlist.count == frj { break }
        }
    }
    
    //! Add a new ring from a "closure" bond: See OBBond::IsClosure()
    func addRingFromClosure(_ mol: MKMol, _ cbond: MKBond) {
        var t1: [MKRTree] = []
        var t2: [MKRTree] = []
        let bv1: MKBitVec = MKBitVec()
        let bv2: MKBitVec = MKBitVec()
        
        var pathok: Bool = false
        
        bv1.setBitOn(UInt32(cbond.getEndAtomIdx()))
        bv2.setBitOn(UInt32(cbond.getBeginAtomIdx()))
        
        buildMKRTreeVector(cbond.getBeginAtom(), nil, &t1, bv1)
        buildMKRTreeVector(cbond.getEndAtom(), nil, &t2, bv2)
        
        var path1: [MKAtom] = []
        var p1: Deque<Int> = []
        var path2: [MKAtom] = []
        var p2: Deque<Int> = []
        
        for tree1 in t1 {
            path1.removeAll()
            tree1.pathToRoot(&path1)
//            Check if the index is greater than the size of t2, if not, the element will be in there
            if t2.count > tree1.getAtomIdx() {
                pathok = true
                path2.removeAll()
                t2[tree1.getAtomIdx()].pathToRoot(&path2)
                
                p1.removeAll()
//              m = path1.begin();
//              if (m != path1.end())
//              p1.push_back((*m)->GetIdx());
                // Not sure what exactly the logic is doing here since it seems to add the first node
                // twice, as it enters the loop over path1
                // we will add the first node and then start the loop beginning at the second node
                if path1.count > 1 {
                    p1.append(path1[0].getIdx())
                }
                for m in path1[1...] {
                    p1.append(m.getIdx())
                    p2.removeAll()
                    for n in path2 {
                        p2.insert(n.getIdx(), at: 0)
                        if n == m {
                            p2.popFirst()
                            if p1.count + p2.count > 2 {
                                saveUniqueRings(p1, p2)
                            }
                            pathok = false
                            break
                        }
                        if n.isConnected(m) && (p1.count + p2.count > 2) {
                            saveUniqueRings(p1, p2)
                        }
                    }
                    if !pathok {
                        break
                    }
                }
            }
        }
        //clean up OBRTree vectors???
        defer {
            t1 = []
            t2 = []
        }
        // set parent for all rings
        for j in 0..<_rlist.count {
            _rlist[j].setParent(mol)
        }
    }
    
    @discardableResult
    func saveUniqueRings(_ d1: Deque<Int>, _ d2: Deque<Int>) -> Bool {
        var path: [Int] = []
        let bv = MKBitVec()
        for i in d1 {
            path.append(i)
            bv.setBitOn(UInt32(i))
        }
        for i in d2 {
            path.append(i)
            bv.setBitOn(UInt32(i))
        }
        for j in _rlist {
            if bv == j._pathSet {
                return false
            }
        }
        let ring = MKRing(path, bv)
        _rlist.append(ring)
        return true
    }
    
}

func atomRingToBondRing(_ mol: MKMol, _ atoms: [Int]) -> [UInt] {
    var bonds: [UInt] = []
    for i in 0..<atoms.count-1 {
        let beginIndex = atoms[i]
        let endIndex = atoms[i+1]
        guard let index = mol.getBond(beginIndex, endIndex)?.getIdx() else { continue }
        bonds.append(index)
    }
    guard let lastBond = mol.getBond(atoms[0], atoms[atoms.count-1]) else { return bonds }
    bonds.append(lastBond.getIdx())
    return bonds
}

/**
   * This function finds the LSSR containing all relevant cycles. A cycle is
   * relevant if it belongs to at least one minimum cycle basis. Another
   * description is more useful though:
   *
   * A cycle (C) is relevant if:
   * - no smaller cycles C_i, ..., C_k exist such that C = C_1 + ... + C_k
   * - both bonds & atoms are checked
   *
   * This is based on lemma 1 from:
   *
   * P. Vismara, Union of all the minimum cycle bases of a graph, The electronic
   * journal of combinatorics, Vol. 4, 1997
   * http://www.emis.de/journals/EJC/Volume_4/PostScriptfiles/v4i1r9.ps
   */
func visitRing(_ mol: MKMol, _ ring: MKRing, _ rlist: inout [MKRing],_ rignored: inout [MKRing]) {
    
    var mask = MKBitVec()
    // Make sure mask is the same size as the maximum ring atom/bond index.
    mask.setBitOn(UInt32(mol.numAtoms()))
    mask.setBitOn(UInt32(mol.numBonds()))
    
    //
    // Remove larger rings that cover the same atoms as smaller rings.
    //
    
    mask.clear()
    for j in 0..<rlist.count {
        // Here we select only smaller rings.
        if rlist[j]._path.count < ring._path.count {
            mask |= rlist[j]._pathSet
        }
    }
    
    mask = mask & ring._pathSet
    
    let containsSmallerAtomRing = mask == ring._pathSet ? true : false
    
    // Translate ring atom indexes to ring bond indexes.
    let bonds = atomRingToBondRing(mol, ring._path)
    let bondset : MKBitVec = MKBitVec()
    bonds.forEach { bd in
        bondset.setBitOn(UInt32(bd))
    }
    
    //
    // Remove larger rings that cover the same bonds as smaller rings.
    //
    mask.clear()
    for j in 0..<rlist.count {
        let otherBonds = atomRingToBondRing(mol, rlist[j]._path)
        let bs = MKBitVec()
        otherBonds.forEach { bd in
            bs.setBitOn(UInt32(bd))
        }
        // Here we select only smaller rings.
        if otherBonds.count < bonds.count {
            mask |= bs
        }
    }
    
    mask = mask & bondset
    
    let containsSmallerBondRing = mask == bondset ? true : false
    
    // The ring is part of the LSSR if all it's atoms and bonds are not
    // found in smaller rings.
    if !containsSmallerAtomRing || !containsSmallerBondRing {
        rlist.append(ring)
    } else {
        rignored.append(ring)
    }
}

func determineFRJ(_ mol: MKMol) -> Int {
    
    if !mol.hasClosureBondsPerceived() {
        return Int(findRingAtomsAndBonds2(mol))
    }
    
    var frj: Int = 0
    for bond in mol.getBondIterator() {
        if bond.isClosure() { // bond->HasFlag(OB_CLOSURE_BOND)?
            frj += 1
        }
    }
    return frj
}

/* A recursive O(N) traversal of the molecule */
@discardableResult
func findRings(_ atom: MKAtom, _ avisit: inout [Int], _ bvisit: inout [Int], _ frj: inout UInt, _ depth: Int) -> Int {
    var result = -1
    
    guard let bonds = atom.getBondIterator() else { return result }
    
    for bond in bonds {
        let bidx = bond.getIdx()
        
        if bvisit[Int(bidx)] == 0 {
            bvisit[Int(bidx)] = 1
            let nbor = bond.getNbrAtom(atom)
            let nidx = nbor.getIdx()
            var nvisit = avisit[nidx]
            if nvisit == 0 {
                avisit[nidx] = depth+1
                nvisit = findRings(nbor, &avisit, &bvisit, &frj, depth+1)
                if nvisit > 0 {
                    if nvisit <= depth {
                        bond.setInRing()
                        if result < 0 || nvisit < result { result = nvisit }
                    }
                }
            } else {
                if result < 0 || nvisit < result { result = nvisit }
                bond.setClosure()
                bond.setInRing()
                frj+=1
            }
        }
    }
    
    if result > 0 && result <= depth { atom.setInRing() }
    return result
}

@discardableResult
func findRingAtomsAndBonds2(_ mol: MKMol) -> UInt {
    
    mol.setRingAtomsAndBondsPerceived()  // mol.SetFlag(OB_RINGFLAGS_MOL);
    mol.setClosureBondsPerceived()       // mol.SetFlag(OB_CLOSURE_MOL);
        
    for atom in mol.getAtomIterator() {
        atom.setInRing(false)
    }
    for bond in mol.getBondIterator() {
        bond.setInRing(false)
        bond.setClosure(false)
    }
    
    let bsize = mol.numBonds() + 1
    var bvisit: [Int] = Array<Int>.init(repeating: 0, count: bsize)
    
    let acount = mol.numAtoms()
    let asize = acount+1
    var avisit: [Int] = Array<Int>.init(repeating: 0, count: asize)
    
    var frj: UInt = 0
    for i in 1...acount {
        if avisit[i] == 0 {
            avisit[i] = 1
            guard let atom = mol.getAtom(i) else { break }
            findRings(atom, &avisit, &bvisit, &frj, 1)
        }
    }
    
    return frj
}

private let MK_RTREE_CUTOFF = 20

func buildMKRTreeVector(_ atom: MKAtom, _ prv: MKRTree?, _ vt: inout [MKRTree], _ bv: MKBitVec) {
    
    vt[atom.getIdx()] = MKRTree(atom, prv)
    
    guard let mol = atom.getParent() else { return }
    var curr: MKBitVec = MKBitVec()
    var used: MKBitVec = MKBitVec()
    var next: MKBitVec = MKBitVec()
    var i: Int = 0
    curr |= atom.getIdx()
    used = bv|curr
    
    var level = 0
    
    repeat {
        next.clear()
        repeat {
            i = curr.nextBit(i)
            guard let aom = mol.getAtom(i) else { break }
            guard let neighs = aom.getNbrAtomIterator() else { break }
            for nbr in neighs {
                if !used[nbr.getIdx()] {
                    next |= nbr.getIdx()
                    used |= nbr.getIdx()
                    vt[nbr.getIdx()] = MKRTree(nbr, vt[aom.getIdx()])
                }
            }
        } while i != bv.endBit()
        
        if next.isEmpty() {
            break
        }
        
        curr = next
        level += 1
    } while level < MK_RTREE_CUTOFF

}
