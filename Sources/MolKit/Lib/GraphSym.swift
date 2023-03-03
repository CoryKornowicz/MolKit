//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/24/23.
//

import Foundation


class MKGraphSymPrivate { 

    var _frag_atoms: MKBitVec = MKBitVec()
    var _pmol: MKMol
    var _conanLabels: [Int]?
    var _stereoUnits: MKStereoUnitSet = MKStereoUnitSet()

    init(_ mol: MKMol) {
        _pmol = mol
    }

    /**
    * Like OBAtom::GetHvyDegree(): Counts the number non-hydrogen
    * neighbors, but doesn't count atoms not in the fragment.
    */
    func getHvyDegree(_ atom: MKAtom) -> Int {
        var count = 0 
        guard let bonds = atom.getBondIterator() else { return count }
        for bond in bonds {
                let nbr = bond.getNbrAtom(atom)
                if _frag_atoms.bitIsSet(nbr.getIndex()) && nbr.getAtomicNum() != MKElements.Hydrogen.atomicNum {
                    count += 1
                }
        }
        return count
    }


    /**
    * Sums the bond order over the bonds from this atom to other atoms
    * in the fragment.  Single = 1, double = 2, triple = 3, aromatic = 1.6,
    * but sum is rounded to nearest integer.
    *
    * This is used for fragment symmetry perception instead of the "implicit
    * valence" used by the standard OpenBabel symmetry perception.  It
    * has the same effect, but we don't have to worry about hydrogen counts,
    * EXCEPT for aromatic N, where the difference between n and [nH] is
    * critical.
    */
    func getHyvBondSum(_ atom: MKAtom) -> Int {
        var count: Float = 0.0
        guard let bonds = atom.getBondIterator() else { return Int(count) }
        for bond in bonds {
            let nbr = bond.getNbrAtom(atom)
            if _frag_atoms.bitIsSet(nbr.getIndex()) && nbr.getAtomicNum() != MKElements.Hydrogen.atomicNum {
                if bond.isAromatic() {
                    count += 1.6
                } else {
                    count += Float(bond.getBondOrder())
                }
            }
        }
        if atom.getAtomicNum() == 7 && atom.isAromatic() && atom.getTotalDegree() == 3 { // [nH] - add another bond
            count += 1.0
        }
        return Int(count + 0.5)  // round to nearest integer
    }

    /**
    * Finds all atoms that are part of a ring in the current fragment.
    * We start with the whole molecule's rings, and eliminate any that
    * have atoms not in the subset.  For the rings that are left, mark
    * each atom of the ring as a ring atom.
    *
    * @return A bit vector where TRUE means it's a ring atom.
    */
    func findRingAtoms(_ ring_atoms: inout MKBitVec) {
        ring_atoms.resize(UInt32(_pmol.numAtoms()))
        ring_atoms.clear()

        let sssRings: [MKRing] = _pmol.getSSSR()

        for ri in sssRings {
            let btmp: MKBitVec = _frag_atoms & ri._pathSet // intersection: fragment and ring
            if btmp == ri._pathSet { //                       all ring atoms are in the fragment?
                ring_atoms |= ri._pathSet //                  yes - add this ring's atoms
            }
        }
    }

    /**
    * Creates a new vector of symmetry classes based on an existing
    * vector.  (Helper routine to GetGIDVector.)  On return, vp2 will
    * have newly-extended connectivity sums, but the numbers (the class
    * IDs) are very large.
    *
    * (Comments by CJ) This appears to compute the "extended connectivity
    * sums" similar to those described by Weininger, Morgan, etc. It uses
    * vp1 as its starting point (the current connectivity sums), and puts
    * the new sums in vp2.  Note that vp1 is modified along the way.
    *
    * Note that, per Weininger's warning, this assumes the initial class
    * ID's are less than 100, which is a BAD assumption, e.g. OCC...CCN
    * would have more than 100 symmetry classes if the chain is more than
    * 98 carbons long.  Should change this to use Weininger's product of
    * corresponding primes.
    */
    func createNewClassVector(_ vp1: [Pair<MKAtom, Int>], _ vp2: inout [Pair<MKAtom, Int>]){
        // There may be fewer atoms than in the whole molecule, so we can't
        // index the vp1 array by atom->GetIdx().  Instead, create a quick
        // mapping vector of idx-to-index for vp1.

        var idx2index: [Int] = [Int](repeating: -1, count: _pmol.numAtoms()+1)
        var index = 0
        for vp_iter in vp1 {
            let idx = vp_iter.0.getIdx()
            idx2index[idx] = index
            index += 1
        }
        // vp2 will hold the newly-extended symmetry classes

        vp2.reserveCapacity(vp1.count)
        vp2.removeAll()

        // Loop over original atoms.
        // Create a new extended varient for each atom.  Get its neighbors' class ID's,
        // sort them into ascending order, and create a sum of (c0 + c1*10^2 + c2*10^4 + ...)
        // which becomes the new class ID (where c0 is the current classID).

        for vp_iter in vp1 {
            let atom = vp_iter.0
            var id = vp_iter.1

            var vtmp: [Int] = [Int]()
            guard let nbrs = atom.getNbrAtomIterator() else { continue }
            for nbr in nbrs {
                let idx = nbr.getIdx()
                if _frag_atoms.bitIsSet(idx) {
                    vtmp.append(vp1[idx2index[idx]].1)
                }
            }
            vtmp.sort(by: {$0 < $1})
            var m = 100 
            for k in vtmp {
                id += k * m
                m *= 100
            }
            vp2.append((atom, id))
        }
    }

    static func createNewClassVector(_ mol: MKMol, _ vp1: [Pair<MKAtom, Int>], _ vp2: inout [Pair<MKAtom, Int>]){
        // There may be fewer atoms than in the whole molecule, so we can't
        // index the vp1 array by atom->GetIdx().  Instead, create a quick
        // mapping vector of idx-to-index for vp1.

        var idx2index: [Int] = [Int](repeating: -1, count: mol.numAtoms()+1)
        var index = 0
        for vp_iter in vp1 {
            let idx = vp_iter.0.getIdx()
            idx2index[idx] = index
            index += 1
        }
        // vp2 will hold the newly-extended symmetry classes

        vp2.reserveCapacity(vp1.count)
        vp2.removeAll()

        // Loop over original atoms.
        // Create a new extended varient for each atom.  Get its neighbors' class ID's,
        // sort them into ascending order, and create a sum of (c0 + c1*10^2 + c2*10^4 + ...)
        // which becomes the new class ID (where c0 is the current classID).

        for vp_iter in vp1 {
            let atom = vp_iter.0
            var id = vp_iter.1

            var vtmp: [Int] = [Int]()
            guard let nbrs = atom.getNbrAtomIterator() else { continue }
            for nbr in nbrs {
                let idx = nbr.getIdx()
                vtmp.append(vp1[idx2index[idx]].1)
            }
            vtmp.sort(by: {$0 < $1})
            var m = 100 
            for k in vtmp {
                id += k * m
                m *= 100
            }
            vp2.append((atom, id))
        }
    }

    /**
    * Calculates a set of graph invariant indexes using the graph theoretical
    * distance, number of connected heavy atoms, aromatic boolean, ring
    * boolean, atomic number, and summation of bond orders connected to the
    * atom.
    *
    * We have to recalculate which atoms are in rings by taking the fragment's
    * atoms into account when we generate the graph invarients.
    *
    * Vector is indexed from zero (not one, like atom->GetIdx()).
    *
    * NOTE: This may need to be extended to include the bond-invariant properties,
    * particularly the size of all rings the bond is in (from a SSSR).
    */
    func getGIVector(_ vid: inout [Int]) { 
        // Prepare the vector...
        vid.removeAll()
        vid.reserveCapacity(_pmol.numAtoms())

        // The "graph theoretical distance" for each atom (see comments in the function)
        var v: [Int] = [Int]()
        getGTDVector(&v)

        // Compute the ring atoms for this particular fragment (set of atoms)
        var ring_atoms: MKBitVec = MKBitVec()
        findRingAtoms(&ring_atoms)
        var i = 0 
        for atom in _pmol.getAtomIterator() {
            vid[i] = Int(MKGraphSym.NoSymmetryClass)
            if _frag_atoms.bitIsSet(atom.getIdx()) {
                vid[i] = 
                v[i]                                                   // 10 bits: graph-theoretical distance
                | (getHvyDegree(atom) << 10)                           //  4 bits: heavy valence
                | ((atom.isAromatic() ? 1 : 0) << 14)                  //  1 bit:  aromaticity
                | ((ring_atoms.bitIsSet(atom.getIdx()) ? 1 : 0) << 15) //  1 bit:  ring atom
                | (atom.getAtomicNum() << 16)                          //  7 bits: atomic number
                | (getHyvBondSum(atom) << 23)                          //  4 bits: heavy bond sum 
                | (7 + atom.getFormalCharge() << 27)                   //  4 bits: formal charge
            }
            i += 1
        }

    }

    /**
    * Calculates the graph theoretical distance of each atom.
    * Vector is indexed from zero.
    *
    * NOTE: "Indexed from zero" means it's one off from the atom->GetIdx()
    * that's used to index atoms inside the molecule!
    *
    * NOTE: This function is hard to decipher, and seems to be misnamed.
    * A "distance" should be be between two atoms, but there's more here
    * than that.  It seems to be doing a breadth-first search to find the
    * most-distant atom from each atom, and reporting the number of steps
    * (which happens to be the graph-theoretical distance) to that atom.
    * The name "Graph Theoretical Distance" is thus misleading.
    */
    @discardableResult
    func getGTDVector(_ gtd: inout [Int]) -> Bool {
        gtd.removeAll()
        gtd.reserveCapacity(_pmol.numAtoms())

        let next: MKBitVec = MKBitVec()
        var curr: MKBitVec = MKBitVec()
        var used: MKBitVec = MKBitVec()
        var gtdcount: Int = 0
        var natom = 0

        next.clear()

        for atom in _pmol.getAtomIterator() {
            let idx = atom.getIdx()
            if !_frag_atoms.bitIsSet(idx) {
                gtd[idx-1] = Int(MKGraphSym.NoSymmetryClass)
                continue
            }
            gtdcount = 0 
            used.clear()
            curr.clear()
            used.setBitOn(UInt32(idx))
            curr.setBitOn(UInt32(idx))

            while !curr.isEmpty() {
                next.clear()
                natom = curr.nextBit(-1)
                while natom != -1 {
                    guard let atom1 = _pmol.getAtom(natom) else { return false }
                    if !_frag_atoms.bitIsSet(atom1.getIdx()) {
                        continue
                    }
                    guard let bonds = atom1.getBondIterator() else { continue }
                    for bond in bonds {
                        let nbr_idx = bond.getNbrAtomIdx(atom1)
                        if _frag_atoms.bitIsSet(nbr_idx) && 
                            !used.bitIsSet(nbr_idx) && 
                            !curr.bitIsSet(nbr_idx) &&
                            bond.getNbrAtom(atom1).getAtomicNum() != MKElements.Hydrogen.atomicNum {
                            next.setBitOn(UInt32(nbr_idx))
                        }
                    }
                    natom = curr.nextBit(natom)
                }
                used |= next 
                curr = next
                gtdcount += 1
            }
            gtd[idx-1] = gtdcount
        }
        return true
    }

    /**
    * Counts the number of unique symmetry classes in a list.
    *
    * (NOTE: CJ -- It also appears to MODIFY the list.  It sorts it in order
    * of class ID, then renumbers the ID's zero through N-1.  See the comments
    * in CreateNewClassVector() about how it returns very large numbers for the
    * class IDs it creates.  These are replaced by lower, sequential numbers here.)
    */
    static func countAndRenumberClasses(_ vp: inout [Pair<MKAtom, Int>], _ count: inout Int) { 
        count = 1
        vp.sort(by: {$0.1 < $1.1})
        var kiter = vp.makeIterator()
        var k = kiter.next()
        while k != nil {
            var id = k!.1
            if id != 0 || id != MKGraphSym.NoSymmetryClass {
                k!.1 = 1
                for var k in kiter {
                    if k.1 != id {
                        id = k.1
                        count += 1
                        k.1 = count
                    } else {
                        k.1 = count
                    }
                }
            }
            k = kiter.next()
        }
    }
    
    /**
    * This is the core of symmetry analysis.  Starting with a set of
    * classes on each atom, it "spreads" them using a sum-of-invariants
    * of each atom's class and its neighbors' classes.  This iterates
    * until a stable solution is found (further spreading doesn't
    * change the answer).
    *
    * @return The number of distinct symmetry classes found.
    */
    func extendInvariants(_ symmetry_classes: inout [Pair<MKAtom, Int>]) -> Int { 
        var nclasses1: Int = 0
        var nclasses2: Int = 0
        var tmp_classes: [Pair<MKAtom, Int>] = [Pair<MKAtom, Int>]()

        // How many classes are we starting with?  (The "renumber" part isn't relevant.)
        MKGraphSymPrivate.countAndRenumberClasses(&symmetry_classes, &nclasses1)

        let nfragatoms = _frag_atoms.countBits()

        // LOOP: Do extended sum-of-invarients until no further changes are
        // noted.  (Note: This is inefficient, as it re-computes extended sums
        // and re-sorts the entire list each time.  You can save a lot of time by
        // only recomputing and resorting within regions where there is a tie
        // initially.  But it's a lot more code.)
        if nclasses1 < nfragatoms {
            // TODO: why is this a for loop?
            for _ in 0..<100 { //sanity check - shouldn't ever hit this number
                createNewClassVector(symmetry_classes, &tmp_classes)
                MKGraphSymPrivate.countAndRenumberClasses(&tmp_classes, &nclasses2)
                symmetry_classes = tmp_classes
                if nclasses1 == nclasses2 {
                    break
                }
                nclasses1 = nclasses2
            }
        }

        createNewClassVector(symmetry_classes, &tmp_classes)
        MKGraphSymPrivate.countAndRenumberClasses(&tmp_classes, &nclasses2)

        if nclasses1 != nclasses2 {
            symmetry_classes = tmp_classes
            return extendInvariants(&symmetry_classes)
        }

        return nclasses1
    }   

    /**
    * Calculates a set of canonical symmetry identifiers for a molecule.
    * Atoms with the same symmetry ID are symmetrically equivalent.  By
    * "canonical", we mean it generates a repeatable labelling of the
    * atoms, i.e. the same fragment will get the same symmetry labels in
    * any molecule in which it occurs.
    *
    * Vector is indexed from zero, corresponding to (atom->GetIdx() - 1).
    *
    * The bit vector "_frag_atoms" specifies a fragment of the molecule,
    * where each bit represents the presence or absence of the atom in
    * the fragment.  Symmetry is computed as though the fragment is the
    * only part that exists.
    */
    func calculateSymmetry(_ atom_sym_classes: inout [UInt]) -> Int {

        // Get vector of graph invariants.  These are the starting "symmetry classes".

        var vgi: [Int] = []
        getGIVector(&vgi)
        var symmetry_classes: [Pair<MKAtom, Int>] = [Pair<MKAtom, Int>]()
        // Create a vector-of-pairs, associating each atom with its Class ID.
        for atom in _pmol.getAtomIterator() {
            let idx = atom.getIdx()
            if _frag_atoms.bitIsSet(idx) {
                symmetry_classes.append((atom, vgi[idx-1]))
            }
            // } else {
            //     symmetry_classes.append(MKGraphSym.NoSymmetryClass)
            // }
        }

        // The heart of the matter: Do extended sum-of-invariants until no further
        // changes are noted.

        let nclasses = extendInvariants(&symmetry_classes)

        // Convert to a vector indexed by Index
        // Atoms not in the fragment will have a value of OBGraphSym::NoSymmetryClass
        atom_sym_classes.removeAll()
        atom_sym_classes = [UInt](repeating: UInt(MKGraphSym.NoSymmetryClass), count: _pmol.numAtoms())
        for i in 0..<symmetry_classes.count {
            atom_sym_classes[symmetry_classes[i].0.getIdx()] = UInt(symmetry_classes[i].1)
        }

        // Store the symmetry classes in an OBPairData
        var temp: String = ""
        for sym_ter in atom_sym_classes {
            temp += "\(sym_ter) "
        }

        let symData = MKPairData<String>()
        symData.setAttribute("OpenBabel Symmetry Classes")
        symData.setOrigin(.local) // will not show as sdf or cml property
        symData.setValue(temp)
        _pmol.setData(symData)
        return nclasses
    }

    func iterate(_ symClasses: inout [Int]) -> Int {
        // Create a vector-of-pairs, associating each atom with its Class ID.
        var symmetry_classes = [Pair<MKAtom, Int>]()
        for atom in _pmol.getAtomIterator() {
            let idx = atom.getIdx()
            if _frag_atoms.bitIsSet(idx) {
                symmetry_classes.append(Pair(atom, symClasses[idx-1]))
            }
        }

        // The heart of the matter: Do extended sum-of-invariants until no further
        // changes are noted.
        let nclasses = extendInvariants(&symmetry_classes)

        // Convert to a vector indexed by Index
        // Atoms not in the fragment will have a value of OBGraphSym::NoSymmetryClass
        symClasses.removeAll()
        symClasses = [Int](repeating: Int(MKGraphSym.NoSymmetryClass), count: _pmol.numAtoms())
        for i in 0..<symmetry_classes.count {
            symClasses[symmetry_classes[i].0.getIndex()] = symmetry_classes[i].1
        }

        return nclasses
    }

    func canonicalLabels(_ symmetry_classes: [Int], _ canon_labels: [Int], _ maxSeconds: Int) { 

    }
}

class MKGraphSym {
    
    static let NoSymmetryClass: UInt = 0x7FFFFFFF
    
    private var d: MKGraphSymPrivate
    
    init(_ pmol: MKMol, _ frag_atoms: inout MKBitVec?) {
        self.d = MKGraphSymPrivate(pmol)
        d._pmol = pmol
        if frag_atoms != nil {
            d._frag_atoms = frag_atoms!
        } else {
            d._frag_atoms.resize(UInt32(d._pmol.numAtoms()))
            for atom in d._pmol.getAtomIterator() {
                d._frag_atoms.setBitOn(UInt32(atom.getIdx()))
            }
        }
    }
    
    /**
    * Calculate the symmetry classes for the molecule. The result will be
    * stored in @p symmetry_classes.
    *
    * The results in @p symmetry_classes will be ordered by symmetry
    * classes. Use the OBAtom* pointer in the std::pair to match the atoms
    * with the right symmetry classes.
    *
    * @return The number of symmetry classes.
    */
    @discardableResult
    func getSymmetry(_ symmetry_classes: inout [UInt]) -> Int {
        clearSymmetry() // For the moment just recalculate the symmetry classes
        // Check to see whether we have already calculated the symmetry classes
        let pd: MKPairData<String>? = d._pmol.getData("OpenBabel Symmetry Classes") as? MKPairData<String>
        
        var nclasses = 0 
        if pd == nil {
            nclasses = d.calculateSymmetry(&symmetry_classes)
        } else {
            let iss = pd!.getValue()!
            symmetry_classes.removeAll()
            //  parse numbers from iss into symmetry classes vector 
            // chuck the string on spaces,
            let parts = iss.components(separatedBy: .whitespaces)
            for part in parts {
                if let num = part.toInt() {
                    symmetry_classes.append(UInt(num))
                }
            }
            // Now find the number of unique elements
            symmetry_classes = symmetry_classes.unique()
            nclasses = symmetry_classes.count
            
        }
        return nclasses
    }
    
    /**
    * Clear the symmetry classes data stored in the molecule specified when
    * construting the OBGraphSym object.
    */
    func clearSymmetry() {
        d._pmol.deleteData("OpenBabel Symmetry Classes")
    }
    
}
//! \brief Handle and perceive graph symmtery for canonical numbering



