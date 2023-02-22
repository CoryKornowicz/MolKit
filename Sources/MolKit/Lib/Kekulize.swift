//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/21/23.
//

import Foundation



/**
   \brief Kekulize a molecule by assigning bond orders of 1 or 2 to aromatic bonds

   Some file formats describe bond orders as aromatic. Such bonds require kekulization
   before the molecule is returned by the reader. Normally, a user should never need
   to call this function themselves.

   This function takes an OBMol which has atoms and bonds marked as aromatic,
   aromatic bonds whose bond orders have been set to single, and aromaticity set as perceived.
   The function assumes that atoms joined by aromatic bonds have been marked as aromatic
   (if they are so intended).

   The purpose of the function is to set the bond orders of the aromatic bonds to either 1 or 2
   in such a way that the valencies of all of the aromatic atoms are satisfied. Failure to do
   this will result in one or more atoms having unsatisfied valences, indicated by a radical.
   Such a failure can only occur if an atom is incorrectly marked as aromatic, or is correctly
   marked as aromatic but has incorrect valence (e.g. 'n' instead of '[nH]' in SMILES).

   \return Whether kekulization was successful
   **/


func getMaxAtomIdx(_ mol: MKMol) -> Int {
    return mol.numAtoms() + 1
}

func getMaxBondIdx(_ mol: MKMol) -> Int {
    return mol.numBonds() + 1
}

class Kekulizer {
    
    var m_mol: MKMol
    var needs_dbl_bond: MKBitVec?
    var doubleBonds: MKBitVec?
    var kekule_system: MKBitVec?
    var atomArraySize: Int = 0
    var bondArraySize: Int = 0
    var m_path: [Int] = []
    
    init(_ mol: MKMol) {
        self.m_mol = mol
        self.needs_dbl_bond = nil
        self.doubleBonds = nil
        self.kekule_system = nil
        
        self.atomArraySize = getMaxAtomIdx(mol) + 1
        self.bondArraySize = getMaxBondIdx(mol) + 1
    }
    
    func greedyMatch() -> Bool {
        
        // What atoms need a double bond? The job of kekulization is
        // to give all of these atoms a single double bond.
        
        needs_dbl_bond = MKBitVec(UInt32(atomArraySize))// defaults to all False
        for atom in m_mol.getAtomIterator() {
            if needsDoubleBond(atom) {
                needs_dbl_bond!.setBitOn(UInt32(atom.getIdx()))
            }
        }
        
        // Make a copy of needs_dbl_bond, to restrict the traversal in BackTrack()
        kekule_system = MKBitVec(needs_dbl_bond!)
        
        // Create lookup of degrees
        var degrees: [Int] = []
        degrees.reserveCapacity(atomArraySize)
        var degreeOneAtoms: [MKAtom] = []
        
        for atom in m_mol.getAtomIterator() {
            let atom_idx = atom.getIdx()
            if !needs_dbl_bond!.bitIsSet(atom_idx) {
                degrees[atom_idx] = 0
                continue
            }
            var mdeg = 0
            guard let bonds = atom.getBondIterator() else { continue }
            for bond in bonds {
                if !bond.isAromatic() { continue }
                let nbr = bond.getNbrAtom(atom)
                if needs_dbl_bond!.bitIsSet(nbr.getIdx()) {
                    mdeg += 1
                }
            }
            degrees[atom_idx] = mdeg
            if mdeg == 1 {
                degreeOneAtoms.append(atom)
            }
        }
        
        // Location of assigned double bonds
        doubleBonds = MKBitVec(UInt32(bondArraySize)) // defaults all to false
        var finished = false
        repeat {
            // Complete all of the degree one nodes
            while !degreeOneAtoms.isEmpty {
                guard let atom = degreeOneAtoms.popLast() else { break }
                
                // some nodes may already have been handled
                if !needs_dbl_bond!.bitIsSet(atom.getIdx()) { continue }
                guard let bonds = atom.getBondIterator() else { continue }
                for bond in bonds {
                    if !bond.isAromatic() { continue }
                    let nbr = bond.getNbrAtom(atom)
                    if !needs_dbl_bond!.bitIsSet(nbr.getIdx()) { continue }
                    // create a double bond from atom -> nbr
                    doubleBonds!.setBitOn(UInt32(bond.getIdx()))
                    needs_dbl_bond!.setBitOff(UInt32(atom.getIdx()))
                    needs_dbl_bond!.setBitOff(UInt32(nbr.getIdx()))
                    // now update degree information for nbr's neighbors
                    guard let nbrbonds = nbr.getBondIterator() else { continue }
                    for nbrbond in nbrbonds {
                        if nbrbond == bond || !nbrbond.isAromatic() { continue }
                        let nbrnbr = nbrbond.getNbrAtom(nbr)
                        let nbrnbrIdx = nbrnbr.getIdx()
                        if !needs_dbl_bond!.bitIsSet(nbrnbrIdx) { continue }
                        degrees[nbrnbrIdx] -= 1
                        if degrees[nbrnbrIdx] == 1 {
                            degreeOneAtoms.append(nbrnbr)
                        }
                    }
                }
                // only a single double bond can be made to atom so we can break here
                break
            }
            
            if needs_dbl_bond!.isEmpty() {
                finished = true
                break
            }
            
            // Now handle any remaining degree 2 or 3 nodes
            // We handle deg 2 nodes first and then 3, and the iteration over these nodes
            // is abstracted away. Once a double-bond is added that generates more
            // degree one nodes, then the iterator is exited
            
            var iterator: NodeIterator = NodeIterator(m_degrees: degrees, m_atomArraySize: atomArraySize)
            var change = false
            var atomIdx = iterator.next()
            repeat {
                if !needs_dbl_bond!.bitIsSet(atomIdx) { continue }
            // The following is almost identical to the code above for deg 1 atoms
            // except for handling the variable 'change'
                guard let atom = m_mol.getAtom(atomIdx) else { break }
                guard let bonds = atom.getBondIterator() else { continue }
                
                for bond in bonds {
                    if !bond.isAromatic() { continue }
                    let nbr = bond.getNbrAtom(atom)
                    if !needs_dbl_bond!.bitIsSet(nbr.getIdx()) { continue }
                    // create a double bond from atom -> nbr
                    doubleBonds!.setBitOn(UInt32(bond.getIdx()))
                    needs_dbl_bond!.setBitOff(UInt32(atomIdx))
                    needs_dbl_bond!.setBitOff(UInt32(nbr.getIdx()))
                    // now update degree information for both atom's and nbr's neighbors
                    for N in 0..<2 {
                        let ref = N == 0 ? atom : nbr
                        guard let nbrbonds = ref.getBondIterator() else { continue }
                        for nbrbond in nbrbonds {
                            if nbrbond == bond || !nbrbond.isAromatic() { continue }
                            let nbrnbr = nbrbond.getNbrAtom(ref)
                            let nbrnbrIdx = nbrnbr.getIdx()
                            if !needs_dbl_bond!.bitIsSet(nbrnbrIdx) { continue }
                            degrees[nbrnbrIdx] -= 1
                            if degrees[nbrnbrIdx] == 1 {
                                degreeOneAtoms.append(nbrnbr)
                                change = true
                            }
                        }
                    }
                    // only a single double bond can be made to atom so we can break here
                    break
                }
                
                if change {
                    break // exit the iterator once we have actually set a double bond
                }
                
                atomIdx = iterator.next()
            } while atomIdx != 0
            // We exit if we are finished or if no degree 2/3 nodes can be set
            if !change {
                break
            }
        } while true // Main Loop
        
        return finished
    }
    
    func backTrack() -> Bool {
        // With an odd number of bits, it's never going to kekulize fully, but let's fill in as many as we can
        guard var needs_dbl_bond = needs_dbl_bond else { return false }
        guard var doubleBonds = doubleBonds else { return false }
        let cound = needs_dbl_bond.countBits()
        var total_handled = 0
        var idx = needs_dbl_bond.firstBit()
        repeat {
            total_handled += 1
            
            // If there is no additional bit available to match this bit, then terminate
            if total_handled == cound {
                return false
            }
            
            // Our goal is to find an alternating path to another atom
            // that needs a double bond
            needs_dbl_bond.setBitOff(UInt32(idx)) // to avoid the trivial null path being found
            var visited = MKBitVec(UInt32(atomArraySize))
            m_path.removeAll()
            let found_path = findPath(idx, false, &visited)
            if !found_path { // could only happen if not kekulizable
                needs_dbl_bond.setBitOn(UInt32(idx)) // reset
                continue
            }
            total_handled += 1
            m_path.append(idx)
            needs_dbl_bond.setBitOff(UInt32(m_path[0]))
            // Flip all of the bond orders on the path from double<-->single
            for i in 0..<m_path.count-1 {
                guard let bond = m_mol.getBond(m_path[i], m_path[i+1]) else { continue }
                if i % 2 == 0 {
                    doubleBonds.setBitOn(UInt32(bond.getIdx()))
                } else {
                    doubleBonds.setBitOff(UInt32(bond.getIdx()))
                }
            }
            idx = needs_dbl_bond.nextBit(idx)
        } while idx != -1
        
        return needs_dbl_bond.isEmpty()
    }
    
    func assignDoubleBonds() {
        guard var bit = doubleBonds?.firstBit() else { return }
        repeat {
            guard let bond = m_mol.getBond(bit) else { continue }
            bond.setBondOrder(2)
            bit = doubleBonds?.nextBit(bit) ?? -1
        } while bit != -1
    }
    
    // The isDoubleBond alternates between double and single, as we need to find
    // an alternating path
    private func findPath(_ atomidx: Int, _ isDoubleBond: Bool, _ visited: inout MKBitVec) -> Bool {
//        Nil checking
        if needs_dbl_bond == nil || doubleBonds == nil || kekule_system == nil {
            return false
        }
        
        if needs_dbl_bond!.bitIsSet(atomidx) {
            return true
        }
        
        visited.setBitOn(UInt32(atomidx))
        guard let atom = m_mol.getAtom(atomidx) else { return false }
        guard let bonds = atom.getBondIterator() else { return false }
        
        for bond in bonds {
            if !bond.isAromatic() { continue }
            let nbr = bond.getNbrAtom(atom)
            if !kekule_system!.bitIsSet(nbr.getIdx()) { continue }
            if doubleBonds!.bitIsSet(Int(bond.getIdx())) == isDoubleBond {
                if visited.bitIsSet(nbr.getIdx()) { continue }
                let found_path: Bool = findPath(nbr.getIdx(), !isDoubleBond, &visited)
                if found_path {
                    m_path.append(nbr.getIdx())
                    return true
                }
            }
        }
        visited.setBitOff(UInt32(atomidx))
        return false
    }
    
}

extension Kekulizer {
    
    func isSpecialCase(_ atom: MKAtom) -> Bool {
        switch atom.getAtomicNum() {
        case 7:
        // Any exo-cyclic double bond from a N
        // e.g. pyridine N-oxide as the double bond form
            if atom.getTotalDegree() == 3 && atom.getFormalCharge() == 0 {
                return true
            }
        case 16: // e.g. Cs1(=O)ccccn1 but not O=s1(=O)cccn1
            if atom.getTotalDegree() == 4 && atom.getFormalCharge() == 0 && atom.getTotalValence() < 6 {
                return true
            }
        default: break
        }
        
        return false
    }

    func needsDoubleBond(_ atom: MKAtom) -> Bool {
        if !atom.isAromatic() {
//            Maybe this is not always a good checkpoint?
            return false
        }
        
        // Does it already have an explicit double bond?
        guard let bonds = atom.getBondIterator() else { return false }
        for bond in bonds {
            if bond.isAromatic() { continue }
            let nbr = bond.getNbrAtom(atom)
            switch bond.getBondOrder() {
            case 0, 1:
                continue
            case 2:
                if isSpecialCase(atom) {
                    return true
                }
                return false
            default: // bond order > 2
                return false
            }
        }
        
        // Is it one of the cases where we know that it only has single bonds?
        
        let chg = atom.getFormalCharge()
        let deg = atom.getTotalDegree()
        
        switch atom.getAtomicNum() {
        case 6:
            if deg == 3 && (chg == 1 || chg == -1) { return false }
        case 5, 7, 15, 33, 51, 83:
            switch chg {
            case 0:  // e.g. a pyrrole-type nitrogen
                if deg == 3 || deg > 4 { return false }
            case -1:
                if deg == 2 { return false }
            case 1:
                if deg > 3 { return false }
            default: break
            }
        case 8, 16, 34, 52:
            switch chg {
            case 0:
                if deg == 2 || deg == 4 || deg > 5 { return false }
            case 1, -1:
                if deg == 3 || deg == 5 || deg > 6 { return false }
            default: break
            }
        default: break
        }
        
        return true // It needs a double bond
    }
    
    
}



class NodeIterator {
    
    var m_degrees: [Int]
    var m_atomArraySize: Int
    var m_counter: Int
    var finishedDegTwo: Bool
    
    init(m_degrees: [Int], m_atomArraySize: Int) {
        self.m_degrees = m_degrees
        self.m_atomArraySize = m_atomArraySize
        self.m_counter = 0
        self.finishedDegTwo = false
    }
        
    func next() -> Int {
        m_counter += 1
        if !finishedDegTwo {
            for m in m_counter..<m_atomArraySize {
                if m_degrees[m] == 2 {
                    return m
                }
            }
            finishedDegTwo = true
            m_counter = 1  // first atom has idx
        }
        
        // return nodes with degree > 2
        for m in m_counter..<m_atomArraySize {
            if m_degrees[m] > 2 {
                return m
            }
        }
        
        // Finished - return 0 signalling the end of iteration
        return 0
    }
    
}

// MKKekulize() implements a two-step kekulization:
//   Step one: try a greedy match
//   Step two: try an exhaustive backtracking (using the results of step one)
//
// The greedy match algorithm is outlined in the thesis of John May
// and indeed NeedsDoubleBond() is based on the implementation in Beam.
// The greedy algorithm almost always works. But when it doesn't, step two is needed.
//
// The goal of the exhaustive backtracking is find a path of alternating single/double
// bonds between two radicals and flip those bonds. For more information, read about
// augmenting paths in the context of perfect matching. John's thesis instead describes
// the use of Edmond's Blossom algorithm which scales better - this may or may not be
// faster in practice for typical chemical graphs.
//
// Potential speedups:
//   * Is OBBitVec performant? I don't know - it seems to do a lot of bounds checking.
//     You could try replacing all usages theoreof with
//     std::vector<char>, where the char could possibly handle several flags.
//   * Before trying the exhaustive search, try a BFS. I have a feeling that this would work
//     90% of the time.
//   * There's a lot of switching between atoms and atom indices (and similar for bonds).
//     Was this completely necessary?
//   * The iterator over degree 2 and 3 nodes may iterate twice - it would have been
//     faster if I just took the first degree 2 or 3 node I came across, but would
//     it have worked as well?
// I would like to thank Noel M. O'Boyle, who thanked John Mayfield, for their implementation.


func MKKekulize(_ mol: MKMol) -> Bool {
    let kekulizer = Kekulizer(mol)
    var success = kekulizer.greedyMatch()
    if !success {
        success = kekulizer.backTrack()
    }
    
    kekulizer.assignDoubleBonds()
    return success
}
