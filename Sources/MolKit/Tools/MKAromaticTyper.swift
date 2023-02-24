



import Foundation

enum ExocyclicAtom {
    case NO_EXOCYCLIC_ATOM
    case EXO_OXYGEN
    case EXO_NONOXYGEN
}


// Start of helper functions for AssignOBAromaticityModel

func findExocyclicAtom(_ atom: MKAtom) -> ExocyclicAtom {
    guard let bonds = atom.getBondIterator() else { return .NO_EXOCYCLIC_ATOM }
    for bond in bonds {
        if bond.getBondOrder() == 2 && !bond.isInRing() {
            var atomicnum = bond.getNbrAtom(atom).getAtomicNum()
            switch atomicnum {
                case MKElements.getAtomicNum("O"):
                    return .EXO_OXYGEN
                default: 
                    return .EXO_NONOXYGEN
            }
        }
    }
    return .NO_EXOCYCLIC_ATOM
}

func hasExocyclicBondToOxygenMinus(_ atom: MKAtom) -> Bool {
    guard let bonds = atom.getBondIterator() else { return false }
    for bond in bonds {
        if bond.getBondOrder() == 1 && !bond.isInRing() {
            var nbr = bond.getNbrAtom(atom)
            if nbr.getAtomicNum() == MKElements.getAtomicNum("O") && nbr.getFormalCharge() == -1 {
                return true
            }
        }
    }
    return false
}

func hasExocyclicDblBondToOxygen(_ atom: MKAtom) -> Bool {
    guard let bonds = atom.getBondIterator() else { return false }
    for bond in bonds {
        if bond.getBondOrder() == 2 && !bond.isInRing() && bond.getNbrAtom(atom).getAtomicNum() == MKElements.getAtomicNum("O") {
            return true
        }
    }
    return false
}

func hasExocyclicDblBondToHet(_ atom: MKAtom) -> Bool {
    guard let bonds = atom.getBondIterator() else { return false }
    for bond in bonds {
        if bond.getBondOrder() == 2 && !bond.isInRing() {
            var atomicnum = bond.getNbrAtom(atom).getAtomicNum()
            switch atomicnum {
            case MKElements.getAtomicNum("C"), MKElements.getAtomicNum("H"):
                break
            default: return true
            }
        }
    }
    return false
}

func assignMKAromaticityModel(_ atm: MKAtom, _ min: inout Int, _ max: inout Int) -> Bool {
    // The Open Babel aromaticity model
    //
    // Return minimum and maximum pi-electrons contributed to an aromatic system
    // The return value indicates a potentially aromatic atom (i.e. any patterns matched)
    //
    // The original code used SMARTS patterns organised in increasing order of
    // prioirity (i.e. later matches overrode earlier ones). These SMARTS patterns
    // are now implemented in the code below, but are included in the comments
    // for reference (Case 1->22).

    if !atm.isInRing() {
        min = 0
        max = 0
        return false
    }    

    let elem = atm.getAtomicNum()
    let chg = atm.getFormalCharge()
    let deg = atm.getExplicitDegree() + Int(atm.getImplicitHCount())
    let val = atm.getExplicitValence() + atm.getImplicitHCount()

    switch elem {
    case MKElements.getAtomicNum("C"): // carbon
        switch chg {
        case 0 :
            if val == 4 && deg  == 3 {
                if hasExocyclicDblBondToHet(atm) {
                    min = 0
                    max = 0
                    return true
                } else {
                    min = 1 
                    max = 1
                    return true
                }
            }
        case 1: 
            if val == 3 {
                switch deg {
                case 3: 
                    min = 0
                    max = 0
                    return true
                case 2:
                    min = 1
                    max = 1
                    return true
                default: break
                }
            }
        case -1:
            if val == 3 {
                switch deg {
                case 3: 
                    min = 2
                    max = 2
                    return true
                case 2:
                    min = 1
                    max = 1
                    return true
                default: break
                }
            }
        default: break
        }
    case MKElements.getAtomicNum("N"),
         MKElements.getAtomicNum("P"): // phosphorus, nitrogen
        switch chg {
        case 0:
            switch val {
            case 3:
                switch deg {
                    case 3:
                        min = 2
                        max = 2
                        return true
                    case 2:
                        min = 1
                        max = 1
                        return true
                    default: break
                }
            case 5:
                if deg == 3 {
                    var exoatom = findExocyclicAtom(atm)
                    switch exoatom {
                        case .EXO_OXYGEN:
                            min = 1
                            max = 1
                            return true
                        case .EXO_NONOXYGEN:
                            min = 2
                            max = 2
                            return true
                        default: break
                    }
                }
            default: break
            }
        case 1: 
            if val == 4 && deg == 3 {
                min = 1
                max = 1
                return true
            }
        case -1:
            if val == 2 && deg == 2 {
                min = 2
                max = 2
                return true
            }
        default: break
        }
    case MKElements.getAtomicNum("O"),
         MKElements.getAtomicNum("Se"): // oxygen, selenium
        switch chg {
        case 0:
            if val == 2 && deg == 2 {
                min = 2
                max = 2
                return true
            }
        case 1:
            if val == 3 && deg == 2 {
                min = 1
                max = 1
                return true
            } 
        default: break
        }
    case MKElements.getAtomicNum("S"): // sulfur
        switch chg {
        case 0:
            switch val {
            case 2: 
                if deg == 2 {
                    min = 2
                    max = 2
                    return true
                }
            case 4: 
                if deg == 3 && hasExocyclicDblBondToOxygen(atm) {
                    min = 2
                    max = 2
                    return true
                }
            default: break 
            }
        case 1:
            if val == 3 {
                switch deg {
                case 2: 
                    min = 1
                    max = 1
                    return true
                case 3: 
                    if hasExocyclicBondToOxygenMinus(atm) {
                        min = 2
                        max = 2
                        return true
                    }
                default: break
                }
            }
        default: break
        }
    case MKElements.getAtomicNum("B"): // boron
        if chg == 0 && val == 3 {
            switch deg {
            case 3: 
                min = 0
                max = 0
                return true
            case 2:
                min = 1
                max = 1
                return true
            default: break
            }
        }
    case MKElements.getAtomicNum("As"): // arsenic
        switch chg {
        case 0:
            if val == 3 {
                switch deg {
                case 3: 
                    min = 2
                    max = 2
                    return true
                case 2:
                    min = 1
                    max = 1
                    return true
                default: break
                }
            }
        case 1: 
            if val == 4 && deg == 3 {
                min = 1
                max = 1
                return true
            }
        default: break
        }
    case 0: //Asterisk
        if chg == 0 {
            switch val {
                case 2: 
                    switch deg {
                    case 2, 3: 
                        min = 0
                        max = 2
                        return true
                    default: break
                    }
                case 3: 
                    switch deg {
                    case 2, 3: 
                        min = 0
                        max = 1
                        return true
                    default: break
                    }
            default: break
            }
        }
    default: break
    }
    min = 0 
    max = 0
    return false // nothing matched
}


class MKAromaticTyperMolState {

    var mol: MKMol
    var _vpa: [Bool] = [] // //!< potentially aromatic atoms
    var _visit: [Bool] = []
    var _root: [Bool] = []
    var _velec: [Pair<Int, Int>] = [] //!< # electrons an atom contributes

    init(_ mol: MKMol) {
        self.mol = mol
        _vpa = [Bool](repeating: false, count: mol.numAtoms()+1)
        _visit = [Bool](repeating: false, count: mol.numAtoms()+1)
        _root = [Bool](repeating: false, count: mol.numAtoms()+1)
        _velec = [Pair<Int, Int>](repeating: Pair<Int, Int>(0, 0), count: mol.numAtoms()+1)
    }

    func assignAromaticFlags() {

        //unset all aromatic flags
        for atom in mol.getAtomIterator() {
            atom.setAromatic(false)
        }
        for bond in mol.getBondIterator() {
            bond.setAromatic(false)
        }

        // New code using lookups instead of SMARTS patterns
        for atom in mol.getAtomIterator() {
            let idx = atom.getIdx()
            _vpa[idx] = assignMKAromaticityModel(atom, &_velec[idx].0, &_velec[idx].1)
        }

        //propagate potentially aromatic atoms
        for atom in mol.getAtomIterator() {
            let idx = atom.getIdx()
            if _vpa[idx] {
                propagatePotentialAromatic(atom)
            }
        }

        //select root atoms
        selectRootAtoms()

        excludeSmallRing() // remove 3 membered rings from consideration

        //loop over root atoms and look for aromatic rings
        for atom in mol.getAtomIterator() {
            if _root[atom.getIdx()] {
                checkAromaticity(atom, 14)
            }
        }
    }

    //! "Anti-alias" potentially aromatic flags around a molecule
    //! (aromatic atoms need to have >= 2 neighboring ring atoms)
    func propagatePotentialAromatic(_ atom: MKAtom) {
        var count = 0
        if let bonds = atom.getBondIterator() {
            for bond in bonds {
                if bond.isInRing() && _vpa[bond.getNbrAtom(atom).getIdx()] {
                    count += 1
                }
            }
        }
        
        if count < 2 {
            _vpa[atom.getIdx()] = false
            if count == 1 {
                if let bonds = atom.getBondIterator() {
                    for bond in bonds {
                        if bond.isInRing() && _vpa[bond.getNbrAtom(atom).getIdx()] {
                            propagatePotentialAromatic(bond.getNbrAtom(atom))
                        }
                    }
                }
            }
        }
    }
    
    /**
    * \brief Select the root atoms for traversing atoms in rings.
    *
    * Picking only the begin atom of a closure bond can cause
    * difficulties when the selected atom is an inner atom
    * with three neighbour ring atoms. Why ? Because this atom
    * can get trapped by the other atoms when determining aromaticity,
    * because a simple visited flag is used in the
    * OBAromaticTyper::TraverseCycle() method.
    *
    * Ported from JOELib, copyright Joerg Wegner, 2003 under the GPL version 2
    * Improved by Fabian (fab5) in 2009 -- PR#2889708
    *
    * @param mol the molecule
    * @param avoidInnerRingAtoms inner closure ring atoms with more than 2 neighbours will be avoided
    *
    */
    func selectRootAtoms(_ avoidInnerRingAtoms: Bool = true) {
        
        var sssRings = mol.getSSSR()

        var cbonds: [MKBond] = []
        var tmpRootAtoms: [Int] = []
        var tmp: [Int] = []
        var hvyNbrs: Int = 0
        var ringNbrs: Int = 0
        var newRoot: Int = -1
        var ringAtoms: [[MKRing]] = [] // store ring pointers on an atom basis

        //generate list of closure bonds
        for bond in mol.getBondIterator() {
            if bond.isClosure() {
                cbonds.append(bond)
                if avoidInnerRingAtoms {
                    tmpRootAtoms.append(bond.getBeginAtomIdx())
                }
            }
        }

        if avoidInnerRingAtoms {
            //for every atom fill vector with ring pointer it's associated with
            for k in sssRings {
                var tmp = k._path
                for j in tmp {
                    ringAtoms[j].append(k)
                }
            }
        }

        //loop over closure bonds

        for bond in cbonds {

            // BASIC APPROACH
            // pick beginning atom at closure bond
            // this is really ready, isn't it ! ;-)

            var rootAtom = bond.getBeginAtomIdx()
            _root[rootAtom] = true

            // EXTENDED APPROACH
            if avoidInnerRingAtoms {
                // count the number of neighbor ring atoms
                guard let atom = mol.getAtom(rootAtom) else { break }
                ringNbrs = 0
                hvyNbrs = 0
                
                guard let nbrs = atom.getNbrAtomIterator() else { break }
                
                for nbr in nbrs {
                    // we can get this from atom->GetHvyDegree()
                    // but we need to find neighbors in rings too
                    // so let's save some time
                    if nbr.getAtomicNum() != MKElements.Hydrogen.atomicNum {
                        hvyNbrs += 1
                        if nbr.isInRing() {
                            ringNbrs += 1
                        }
                    }
                    
                    // if this atom has more than 2 neighbor ring atoms
                    // we could get trapped later when traversing cycles
                    // which can cause aromaticity false detection
                    
                    newRoot = -1
                    
                    if ringNbrs > 2 {
                        // try to find another root atom
                        // only loop over rings which contain rootAtom
                        for ring in ringAtoms[rootAtom] {
                            tmp = ring._path
                            
                            var checkThisRing = false
                            var rootAtomNumber = 0
                            var idx = 0
                            // avoiding two root atoms in one ring!
                            for j in 0..<tmpRootAtoms.count {
                                idx = tmpRootAtoms[j]
                                if ring.isInRing(idx) {
                                    rootAtomNumber += 1
                                    if rootAtomNumber >= 2 {
                                        break
                                    }
                                }
                            }
                            if rootAtomNumber < 2 {
                                for j in 0..<tmp.count {
                                    // find critical ring
                                    if tmp[j] == rootAtom {
                                        checkThisRing = true
                                    } else {
                                        // second root atom in this ring ?
                                        if _root[tmp[j]] == true {
                                            // when there is a second root
                                            // atom this ring can not be
                                            // used for getting an other
                                            // root atom
                                            checkThisRing = false
                                            break
                                        }
                                    }
                                }
                            }
                            // check ring for getting another
                            // root atom to avoid aromaticity typer problems
                            if checkThisRing {
                                // check if we can find another root atom
                                for m in 0..<tmp.count {
                                    ringNbrs = 0
                                    hvyNbrs = 0
                                    guard let atom = mol.getAtom(tmp[m]) else { break }
                                    guard let nbrs2 = atom.getNbrAtomIterator() else { break }
                                    for nbr2 in nbrs2 {
                                        if nbr2.getAtomicNum() != MKElements.Hydrogen.atomicNum {
                                            hvyNbrs += 1
                                            if nbr2.isInRing() {
                                                ringNbrs += 1
                                            }
                                        }
                                    }
                                    // if the number of neighboured heavy atoms is also
                                    // the number of neighboured ring atoms, the aromaticity
                                    // typer could be stuck in a local traversing trap
                                    if ringNbrs <= 2 && ring.isInRing(atom.getIdx()) {
                                        newRoot = tmp[m]
                                    }
                                }
                            }
                        }
                        
                        if (newRoot != -1) && (rootAtom != newRoot) {
                            // unset root atom
                            _root[rootAtom] = false
                            //pick new root atom
                            _root[newRoot] = true
                        }
                    } // if (ringNbrs > 2)
                }// end for
            } // if (avoid)
        }// end for(closure bonds)
    }
    
    //! Remove 3-member rings from consideration
    func excludeSmallRing() {
//        This needs to be refactor ASAP
        for atom in mol.getAtomIterator() {
            if _root[atom.getIdx()] {
                if let bonds = atom.getBondIterator() {
                    for bondj in bonds {
                        if bondj.isInRing() && _vpa[bondj.getNbrAtom(atom).getIdx()] {
                            if let nbrbonds = bondj.getNbrAtom(atom).getBondIterator() {
                                for bondk in nbrbonds {
                                    if bondk.getNbrAtom(bondj.getNbrAtom(atom)) != atom &&
                                        bondk.isInRing() && _vpa[bondk.getNbrAtom(bondj.getNbrAtom(atom)).getIdx()] {
                                        if atom.isConnected(bondk.getNbrAtom(bondj.getNbrAtom(atom))) {
                                            _root[atom.getIdx()] = false
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    //! Check aromaticity starting from the root atom, up to a specified depth
    func checkAromaticity(_ atom: MKAtom, _ depth: Int) {
        var erange: Pair<Int, Int>
        if let bonds = atom.getBondIterator() {
            for bond in bonds {
                if bond.isInRing() { // check all rings, regardless of assumed aromaticity
                    erange = _velec[atom.getIdx()]
                    if traverseCycle(atom, bond.getNbrAtom(atom), bond, &erange, depth - 1) {
                        atom.setAromatic()
                        bond.setAromatic()
                    }
                }
            }
        }
    }

    /** \brief Traverse a potentially aromatic cycle starting at @p root.
      \return  True if the cycle is likely aromatic
      \param root  The initial, "root" atom in traversing this ring
      \param atom  The current atom to visit and check
      \param prev  The bond traversed in moving to this @p atom
      \param er    The min and max number of pi electrons for this ring
      \param depth The maximum number of atoms to visit in a ring (e.g., 6)

      This method traverses a potentially aromatic ring, adding up the possible
      pi electrons for each atom. At the end (e.g., when @p atom == @p root)
      the Huekel 4n+2 rule is checked to see if there is a possible electronic
      configuration which corresponds to aromaticity.
    **/
    func traverseCycle(_ root: MKAtom, _ atom : MKAtom, _ prev: MKBond, _ er: inout Pair<Int, Int>, _ depth: Int) -> Bool {
        
        if atom == root {
            for i in er.0...er.1 {
                if i % 4 == 2 && i > 2 {
                    return true
                }
            }
            return false
        }
        
        if depth == 0 || !_vpa[atom.getIdx()] || _visit[atom.getIdx()] {
            return false
        }
        
        var result = false
        
        var depth = depth - 1
        
        er.0 += _velec[atom.getIdx()].0
        er.1 += _velec[atom.getIdx()].1
        
        _visit[atom.getIdx()] = true
        if let bonds = atom.getBondIterator() {
            for bond in bonds {
                let nbr = bond.getNbrAtom(atom)
                if bond != prev && bond.isInRing() && _vpa[nbr.getIdx()] {
                    if traverseCycle(root, nbr, bond, &er, depth) {
                        result = true
                        bond.setAromatic()
                    }
                }
            }
        }
        
        _visit[atom.getIdx()] = false
        if result {
            atom.setAromatic()
        }
        
        er.0 -= _velec[atom.getIdx()].0
        er.1 -= _velec[atom.getIdx()].1

        return result
      }

}


class MKAromaticTyper {

    init() { }
   
    func assignAromaticFlags(_ mol: MKMol) {
        if mol.hasAromaticPerceived() { return }
        mol.setAromaticPerceived()
        let molstate = MKAromaticTyperMolState(mol)
        molstate.assignAromaticFlags()
    }

}
