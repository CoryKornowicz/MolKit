



/**
* Calculate the canonical labels for the molecule. Stereochemistry is
* included in the algorithm and the canonical labels. The result will be
* stored in @p canonical_labels.
*
* @param mol The molecule.
* @param symmetry_classes The symmetry_classes for the molecule. These can
* be obtained using the OBGraphSym class.
* @param canonical_labels Reference to the object to store the results in.
* @param mask The fragment to label. When the bit for an atom is set, it is
* included in the fragment. If no bits are set, all atoms will be included.
* Atoms are indexed from 1 (i.e. OBAtom::GetIdx()).
* @param maxSeconds Timeout in seconds.
* @param onlyOne If true, the first found labels are returned. These are
* canonical labels without considering stereochemistry and other attributes
* not included in the symmetry classes.
*
* @return The canonical labels for the molecule in @p canonical_labels.
*
* @see @ref canonical_code_algorithm
* @since 2.3
*/

import Foundation
import Algorithms
import Collections

let MAX_IDENTITY_NODES: Int = 50

func totalHydrogenCount(_ atom: MKAtom) -> UInt {
    return atom.explicitHydrogenCount() + atom.getImplicitHCount()
}

/**
* Helper function for getFragment below.
*/
func addNbrs(_ fragment: MKBitVec, _ atom: MKAtom, _ mask: MKBitVec, _ metalloceneBonds: [MKBond]) {
    for nbr in atom.getNbrAtomIterator()! {
        if !mask.bitIsSet(nbr.getIdx()) { // interestingly, you mask __only__ the atoms you want to check, not the atoms to ignore
            continue
        }
        // skip visited atoms
        if fragment.bitIsSet(nbr.getIdx()) {
            continue
        }
        // skip mettalocene bonds
        if metalloceneBonds.contains(where: { $0 === atom.getParent()!.getBond(atom, nbr) }) {
            continue
        }
        // add the neighbor atom to the fragment
        fragment.setBitOn(UInt32(nbr.getIdx()))
        // recurse...
        addNbrs(fragment, nbr, mask, metalloceneBonds)
    }
}

/**
* Create an OBBitVec objects with bets set for the fragment consisting of all
* atoms for which there is a path to atom without going through skip. These
* fragment bitvecs are indexed by atom idx (i.e. OBAtom::GetIdx()).
 */
func getFragment(_ atom: MKAtom, _ mask: MKBitVec, _ metalloceneBonds: [MKBond] = [MKBond]()) -> MKBitVec {
    let fragment: MKBitVec = MKBitVec()
    fragment.setBitOn(UInt32(atom.getIdx()))
    // start the recursion
    addNbrs(fragment, atom, mask, metalloceneBonds)
    return fragment
}


struct getFragmentImpl {
    // OBBitVec &fragment, OBAtom *atom, OBAtom *skip, const OBBitVec &mask
    static func addNbrs(_ fragment: MKBitVec, _ atom: MKAtom, _ skip: MKAtom, _ mask: MKBitVec) {
        for nbr in atom.getNbrAtomIterator()! {
            // don't pass through skip
            if nbr.getIdx() == skip.getIdx() {
                continue
            }
            // skip visited atoms
            if fragment.bitIsSet(nbr.getIdx()) {
                continue
            }
            if !mask.bitIsSet(nbr.getIdx()) {
                continue
            }
            // add the neighbor atom to the fragment
            fragment.setBitOn(UInt32(nbr.getIdx()))
            // recurse...
            addNbrs(fragment, nbr, skip, mask)
        }
    }
}

func getFragment(_ atom: MKAtom, _ skip: MKAtom, _ mask: MKBitVec) -> MKBitVec {
    let fragment: MKBitVec = MKBitVec()
    fragment.setBitOn(UInt32(atom.getIdx()))
    // start the recursion
    getFragmentImpl.addNbrs(fragment, atom, skip, mask)
    return fragment
}

func isFerroceneBond(_ bond: MKBond) -> Bool {
    if bond.getBondOrder() != 1 {
        return false
    }

    var Fe: MKAtom? = nil
    var C: MKAtom? = nil

    let begin: MKAtom = bond.getBeginAtom()
    if begin.getAtomicNum() == 26 {
        Fe = begin
    }
    if begin.getAtomicNum() == 6 {
        C = begin
    }

    let end: MKAtom = bond.getEndAtom()
    if end.getAtomicNum() == 26 {
        Fe = end
    }
    if end.getAtomicNum() == 6 {
        C = end
    }

    if Fe == nil || C == nil {
        return false
    }

    if Fe!.getExplicitDegree() < 10 {
        return false
    }

    return C!.hasDoubleBond() && C!.isInRing()
}

func findMetalloceneBonds(_ bonds: inout [MKBond], _ mol: MKMol, _ symmetry_classes: [UInt]) {
    for atom in mol.getAtomIterator() {
        if !atom.isInRingSize(3) {
            continue
        }
        var nbrSymClasses: [UInt] = [UInt]()
        for nbr in atom.getNbrAtomIterator()! {
            if nbr.isInRingSize(3) {
                nbrSymClasses.append(symmetry_classes[nbr.getIdx()])
            }
        }
        if nbrSymClasses.count < 8 {
            continue
        }
        nbrSymClasses.sort()
        let numUnique: UInt = UInt(nbrSymClasses.unique().count)
        if numUnique > 1 {
            continue
        }
        for nbr in atom.getNbrAtomIterator()! {
            bonds.append(mol.getBond(atom, nbr)!)
        }
    }
}


// MARK: - Canonicalization Note* 

struct CanonicalLabelsImpl {

    typealias Orbit = [MKAtom]
    typealias Orbits = [Orbit]

    static func print_orbits(_ orbits: Orbits) {
        for orbit in orbits {
            for atom in orbit {
                print("(\(atom.getIdx()))")
            }
        }
    }

    static func print_orbits(_ label: String, _ orbits: Orbits) {
        print("\(label): ")
        print_orbits(orbits)
    }
//     TODO: these all neede to be connected with inout because structs are value types
//     ** technically with copy on write, the data should not change, but it will not hold its initial reference...which might not be an issue, but would cause a lot of extraneous writing.**
    /**
     * Structure to represent a completed labeling and it's associated canonical
     * candidate code.
     */
    struct FullCode {
        /**
         * The atom labels starting from 1. All excluded atoms have label 0.
         */
        var labels: [UInt]
        /**
         * The canonical candidate code resulting from the @p labels.
         */
        var code: [UInt]
        /**
         * Default constructor, used to create initial empty bestCode object.
         */
        init() {
            labels = [UInt]()
            code = [UInt]()
        }
        
        /**
         * Constructor specifying both the @p labels and @p from (the FROM part
         * of the canonical code). The other parts of the canonical code are
         * added using the CompleteCode function.
         */
        init(_ labels: [UInt], _ from: [UInt]) {
            self.labels = labels
            code = from
        }
        /**
         * Compare this object's canonical candidate code to the other object's
         * code and return true if this code is greater. This method is used to
         * select the final canonical code. In other words, the canonical labels
         * for a molecule are the labels with the greatest canonical code.
         */
        static func > (lhs: FullCode, rhs: FullCode) -> Bool {
            return compareArraysGreaterThan(lhs.code, rhs.code)
        }
    }
    
    /**
     * Structure to represent a partial labeling and the associated FROM part
     * of the canonical code.
     */
    struct PartialCode {
        
        /**
         * The atoms in canonical order. An atom is added when a label is assigned
         * to it. When all atoms are labeled, this vector will contain all atoms
         * in the fragment.
         */
        var atoms: [MKAtom] = []
        /**
         * The bonds in canonical order. Since the atoms are labeled from another
         * atom ("from"), these are the bonds connecting the "from" atom with the
         * next labeled atom. The FROM list forms a spanning tree that does not
         * include the ring closure bonds. These ring closure bonds are added at
         * the end of this vector, in the correct order, in the CompleteCode
         * function.
         */
        var bonds: [MKBond] = []
        /**
         * The FROM part of the canonical code. This is a vector containing the
         * label of the "from" atom for the atoms in canonical order. The initial
         * atom of a bonded fragment is not included in this list since there is
         * no "from" atom.
         */
        var from: [UInt] = []
        /**
         * The atom labels starting from 1. All atoms start with label 0.
         */
        var labels: [UInt] = []
        
        /**
         * Constructor specifying the number of atoms in the molecule. The @p
         * numAtoms are used to initialize the labels vector to contain label
         * 0 for all atoms.
         */
        init(_ numAtoms: Int) {
            self.labels = [UInt].init(repeating: 0, count: numAtoms)
        }
        
        /**
         * Add an atom to atoms. Used to add the initial atom.
         */
        mutating func add(_ atom: MKAtom) { atoms.append(atom) }
        /**
         * Add a bond to bonds. Used to add the ring closure bonds.
         */
        mutating func add(_ bond: MKBond) { bonds.append(bond) }
        
        /**
         * Add an atom that is labeled "from" another atom. This function adds
         * the @p atom to atoms, adds bond between @p fromAtom and @p atom to
         * bonds and updates the @p from vector.
         */
        mutating func add(_ fromAtom: MKAtom, _ atom: MKAtom) {
            from.append(labels[fromAtom.getIndex()])
            atoms.append(atom)
            bonds.append(atom.getParent()!.getBond(fromAtom, atom)!)
        }
        
        /**
         * Compare the @p from vector with a FullCode's code and return true if
         * this code is lower. This function is used to avoid completing a
         * labeling that will never result in a greatest canonical code.
         */ 
        
        static func < (lhs: PartialCode, rhs: FullCode) -> Bool {
            
            let numFrom = min(lhs.from.count, rhs.code.count)
            for i in 0..<numFrom {
                if lhs.from[i] > rhs.code[i] {
                    return false
                }
                if lhs.from[i] < rhs.code[i] {
                    return true
                }
            }
            return false
        }
    }
    
    /**
     * Structure to store and compute stereo center parities. Currently for
     * both tetrahedral and cistrans stereochemistry. Precomputing these
     * StereoCenter objects is done for performance reasons since parities
     * need to be computed for each canonical candidate labeling.
     */
    
    struct StereoCenter {
        
        /**
         * Indexes for the stereo center neighbor atoms. Tetrahedral centers have
         * all neighbor atoms in nbrIndexes1. CisTrans stereo centers store the
         * neighbor atoms for each double bond atom separately.
         */
        var nbrIndexes1: [UInt] = []
        var nbrIndexes2: [UInt] = []
        /**
         * Index(es) for the stereo center. This is used to sort the stereo
         * centers by label. Tetrahedral stereo centers have one index in this
         * vector and cistrans stereocenters have both double bond atom indexes
         * in this vector.
         */
        var indexes: [UInt] = []
        
        /**
         * Use the stored nbrIndexes to access the @p labels and compute a parity
         * or descriptor. This function returns 0 or 1 for both tetrahedral and
         * cistrans stereo centers.
         */
        func getDescriptor(_ symmetry_classes: [UInt], _ labels: [UInt]) -> UInt {
            // Unspecified stereo centers have their own descriptor.
            if nbrIndexes1.isEmpty {
                return 2
            }
            var refs1: [UInt] = []
            var refs2: [UInt] = []
            for i in 0..<nbrIndexes1.count {
                if nbrIndexes1[i] < labels.count {
                    refs1.append(labels[Int(nbrIndexes1[i])])
                } else {
                    refs1.append(nbrIndexes1[i])
                }
            }
            for i in 0..<nbrIndexes2.count {
                if nbrIndexes2[i] < labels.count {
                    refs2.append(labels[Int(nbrIndexes2[i])])
                } else {
                    refs2.append(nbrIndexes2[i])
                }
            }
            if indexes.count == 2 {
                let symOrder = symmetry_classes[Int(indexes[0])] < symmetry_classes[Int(indexes[1])]
                let canOrder = labels[Int(indexes[0])] < labels[Int(indexes[1])]
                if symOrder != canOrder && symmetry_classes[Int(indexes[0])] != symmetry_classes[Int(indexes[1])] {
                    refs1.swapAt(0, 1)
                }
            }
            return UInt(((MKStereo.numInversions(refs1) % 2 + MKStereo.numInversions(refs2) % 2) % 2))
        } 
    }
    
    /**
     * Sort StereoCenter objects by their label. Used in combination with
     * std::sort to create STEREO code.
     */
    struct SortStereoCenters: SortComparator {
        
        typealias Compared = StereoCenter
        var order: SortOrder
        var labels: [UInt]
        
        init(order: SortOrder, labels: [UInt]) {
            self.order = order
            self.labels = labels
        }
        
        static func == (lhs: CanonicalLabelsImpl.SortStereoCenters, rhs: CanonicalLabelsImpl.SortStereoCenters) -> Bool {
            return lhs.order == rhs.order && lhs.labels == rhs.labels
        }
        
        private func getLabel(_ c: CanonicalLabelsImpl.StereoCenter) -> UInt {
            switch c.indexes.count {
            case 2:
                return min(labels[Int(c.indexes[0])], labels[Int(c.indexes[1])])
            default:
                return labels[Int(c.indexes[0])]
            }
        }
        
        func compare(_ lhs: CanonicalLabelsImpl.StereoCenter, _ rhs: CanonicalLabelsImpl.StereoCenter) -> ComparisonResult {
            let lh = getLabel(lhs)
            let rh = getLabel(rhs)
            if lh < rh {
                return .orderedAscending
            } else if lh > rh {
                return .orderedDescending
            } else {
                return .orderedSame
            }
        }
        
        func hash(into hasher: inout Hasher) {
            hasher.combine(order.hashValue)
            hasher.combine(labels)
        }
        
    }
    
    /**
     * Sort FullCode objects by their code. Used in combination with
     * std::sort order disconnected fragments.
     */
    static func sortCode(_ code1: FullCode, _ code2: FullCode) -> Bool {
//        return code1.code < code2.code // original implementationwas assumed to be using lexicographical comparison...trailing zeros could be an issue, easy fix if needed
        return  compareArraysLessThan(code1.code, code2.code)
    }
    /**
     * Sort FullCode objects with associated int in std::pair. Used to sort
     * ligands (i.e. equivalent fragments connected to a symmetric atom).
     */
    static func sortCode2(_ code1: (Int, FullCode), _ code2: (Int, FullCode)) -> Bool {
        return compareArraysGreaterThan(code1.1.code, code2.1.code)
    }
    /**
     * Sort atoms by ascending ranks (e.g. labels).
     */
    struct SortAtomsAscending: SortComparator {
        
        typealias Compared = MKAtom
        var order: SortOrder
        var _ranks: [UInt]
        
        init(ranks: [UInt]) {
            self.order = .forward
            self._ranks = ranks
        }
        func compare(_ lhs: MKAtom, _ rhs: MKAtom) -> ComparisonResult {
            let ra = _ranks[lhs.getIndex()]
            let rb = _ranks[rhs.getIndex()]
            if ra < rb {
                return .orderedAscending
            } else if ra > rb {
                return .orderedDescending
            } else {
                return .orderedSame
            }
        }
    }

    /**
     * Sort atoms by descending ranks (e.g. symmetry classes)
     */
    struct SortAtomsDescending: SortComparator {
        
        typealias Compared = MKAtom
        var order: SortOrder
        var _ranks: [UInt]
        
        init(ranks: [UInt]) {
            self.order = .reverse
            self._ranks = ranks
        }
        func compare(_ lhs: MKAtom, _ rhs: MKAtom) -> ComparisonResult {
            let ra = _ranks[lhs.getIndex()]
            let rb = _ranks[rhs.getIndex()]
            if ra < rb {
                return .orderedDescending
            } else if ra > rb {
                return .orderedAscending
            } else {
                return .orderedSame
            }
        }
    }
    
    /**
     * Structure used while labeling a single connected fragment. All input is
     * specified using the constructor and @p code is generated as a result of
     * assigning labels.
     */
    struct State {
        
        var symmetry_classes: [UInt]
        /**
         * The connected fragment. This is a subset of the mask.
         */
        var fragment: MKBitVec
        var onlyOne: Bool = false
        /**
         * The pre-computed stereo centers. Non-const since it needs to be
         * sorted.
         */
        var stereoCenters: [StereoCenter]
        /**
         * The partial code with labels and FROM code.
         */
        var code: PartialCode
        /**
         * Identity nodes of the search tree.
         */
        var identityCodes: [FullCode]
        var backtrackDepth: UInt = 0 
        var orbits: Orbits 
        var mcr: MKBitVec 
        
        init(_ _symmetry_classes: [UInt], _ _fragment: MKBitVec, _ _stereoCenters: [StereoCenter], 
             _ _identityCodes: [FullCode], _ _orbits: Orbits, _ _mcr: MKBitVec, _ _onlyOne: Bool = false) {
            symmetry_classes = _symmetry_classes
            fragment = _fragment
            onlyOne = _onlyOne
            stereoCenters = _stereoCenters
            code = PartialCode(_symmetry_classes.count)
            identityCodes = _identityCodes
            orbits = _orbits
            mcr = _mcr
            
            mcr.clear()
            if mcr.isEmpty() {
                for i in 0..<_symmetry_classes.count {
                    mcr.setBitOn(UInt32(i+1))
                }
            }
        }
    } 
    
    /**
     * Simple struct to manage timeout.
     */
    struct Timeout {
        var startTime: Date
        var maxTime: TimeInterval
        init(_ _maxTime: TimeInterval) {
            startTime = Date()
            maxTime = _maxTime
        }
    }
    
    /**
    * Given a complete labeling (@p state.code.labels) and FROM code (@p state.code.from),
    * construct a FullCode with complete code.
    */
    static func completeCode(_ mol: MKMol, _ fullcode: inout FullCode, _ state: inout State) {
        var code: PartialCode = state.code
        // initialize the FullCode object with the found labels and the from code
        fullcode = FullCode(code.labels, code.from)
        var numClosures: UInt = 0
        //
        // the RING-CLOSURE list
        //
        var current_label: UInt = 1 // the second atom is always labeled from 1
        // iterate over all atoms in the order they are labeled
        // (this ensures [1 3] < [3 1])
        for j in 0..<code.atoms.count {
            let atom = code.atoms[j]
            // still need to sort [1 3] and [1 4]
            var closures: [Pair<MKBond, UInt>] = []
            for bond in atom.getBondIterator()! {
                // skip atoms not in the fragment 
                if !state.fragment.bitIsSet(bond.getNbrAtom(atom).getIdx()) {
                    continue
                }
                // a closure bond is a bond not found while generating the FROM spanning tree.
                if !code.bonds.contains(bond) {
                    closures.append((bond, code.labels[bond.getNbrAtom(atom).getIndex()]))
                }
            }
            // do the sorting: [1 3] < [1 4]
            closures.sort { (a, b) -> Bool in
                return a.1 < b.1
            }
            for k in closures {
                // add the closure bond to the code
                fullcode.code.append(current_label)
                fullcode.code.append(k.1)
                // add the bond to the list (needed for BOND-TYPES below)
                code.add(k.0)
                numClosures += 1
            }
            current_label += 1
        } 

        // Isotopes are only considered if there are isotopes.
        var hasIsotope = false
      // Charges are only considered if there is a formal charge.
        var hasCharge = false
        _ = MKStereoFacade(mol)
        for i in 0..<code.atoms.count {
            let atom = code.atoms[i]
            if atom.getIsotope() != 0 {
                hasIsotope = true
            }
            if atom.getFormalCharge() != 0 {
                hasCharge = true
            }
            // Include all hydrogens
            let hydrogens_to_include = totalHydrogenCount(atom)
            let c = 10000 * atom.getSpinMultiplicity() + 1000 * Int(hydrogens_to_include) + atom.getAtomicNum()
            // add the atomic number to the code
            fullcode.code.append(UInt(c))
        }
        
        // 
        // the (optional) ISOTOPE list 
        //
        if hasIsotope {
            for j in 0..<code.atoms.count {
                fullcode.code.append(code.atoms[j].getIsotope())
            }
        }
        //
        // the (optional) CHARGES list
        // 
        if hasCharge {
            for j in 0..<code.atoms.count {
                fullcode.code.append(UInt(code.atoms[j].getFormalCharge()))
            }
        }
        //
        // the BOND-TYPES list
        //
        for j in 0..<code.bonds.count {
            let bond = code.bonds[j]
            if bond.isAromatic() {
                fullcode.code.append(5)
            } else {
                fullcode.code.append(bond.getBondOrder())
            }
        }
        // the STEREO flag
        if state.stereoCenters.count > 0 {
            // sort the stereo centers 
            state.stereoCenters.sort(using: SortStereoCenters(order: .forward, labels: code.labels))
            
            for i in 0..<state.stereoCenters.count {
                var isInFragment: Bool = false
                for j in 0..<state.stereoCenters[i].indexes.count {
                    if state.fragment.bitIsSet(Int(state.stereoCenters[i].indexes[j]) + 1) {
                        isInFragment = true
                        break
                    }
                }
                // ignore stereo centers not in this fragment
                if isInFragment {
                    fullcode.code.append(state.stereoCenters[i].getDescriptor(state.symmetry_classes, code.labels))
                }
            }
        }
        // backtrack
        for _ in 0..<numClosures {
            _ = code.bonds.popLast()
        }
    }
    /**
     * This function implements optimization 1 as described above.
     */
    static func labelFragments(_ current: MKAtom, _ nbrs: inout [MKAtom], _ label: UInt, _ timeout: inout Timeout, _ bestCode: inout FullCode, _ state: inout State) {
        guard let mol = current.getParent() else { 
            fatalError("Parent Mol could not be retrieved")
        }
        var code: PartialCode = state.code

        // sort the neighbor atoms according to symmetry class 
        nbrs.sort(using: SortAtomsDescending(ranks: state.symmetry_classes))
        // Count the number of unique neighbor symmetry classes.
        var numUnique: UInt = 1
        var lastSymClass = state.symmetry_classes[nbrs[0].getIndex()]
        for i in 0..<nbrs.count {
            let symClass = state.symmetry_classes[nbrs[i].getIndex()]
            if symClass != lastSymClass {
                numUnique += 1
            }
            lastSymClass = symClass
        }

        if numUnique < nbrs.count {
            // The canonical codes for the equivalent ligands
            var lcodes: [Pair<Int, CanonicalLabelsImpl.FullCode>] = []
            var ligandSizes: [UInt] = []
            for i in 0..<nbrs.count {
                let ligand: MKBitVec = getFragment(current, nbrs[i], state.fragment)
                ligandSizes.append(UInt(ligand.countBits()))
                var lbestCode = FullCode()
                if ligandSizes.last! == 1 {
                    // Avoid additional state creation if the neighbor is a single terminal atom.
                    lbestCode.code.append(UInt(nbrs[i].getAtomicNum()))
                    lbestCode.labels.resize(state.symmetry_classes.count, with: 0)
                    lbestCode.labels[nbrs[i].getIndex()] = 1
                } else if ligandSizes.last! == 2 {
                    // Avoid additional state creation if the neighbor is a fragment of 2 atoms.
                    lbestCode.labels.resize(state.symmetry_classes.count, with: 0)
                    lbestCode.code.append(1) // FROM
                    lbestCode.code.append(UInt(nbrs[i].getAtomicNum()))
                    lbestCode.labels[nbrs[i].getIndex()] = 1
                    for nbr in nbrs[i].getNbrAtomIterator()! {
                        if !state.fragment.bitIsSet(nbr.getIdx()) { continue }
                        if code.labels[nbr.getIndex()] == 0 { continue }
                        lbestCode.code.append(UInt(nbr.getAtomicNum())) // ATOM-TYPES 2
                        let bond = mol.getBond(nbrs[i], nbr)!
                        if bond.isAromatic() {
                            lbestCode.code.append(5)
                        } else {
                            lbestCode.code.append(bond.getBondOrder())
                        }
                        lbestCode.labels[nbr.getIndex()] = 2
                    }
                } else {
                    // Start labeling from the ligand atom.
                    let identityCode: [FullCode] = []
                    let orbits = Orbits()
                    let mcr = MKBitVec()
                    var lstate = State(state.symmetry_classes, ligand, state.stereoCenters, identityCode, orbits, mcr, state.onlyOne)
                    lstate.code.add(nbrs[i])
                    lstate.code.labels[nbrs[i].getIndex()] = 1
                    canonicalLabelRecursive(nbrs[i], 1, &timeout, &lbestCode, &lstate)
                }
                // Store the canonical code (and labels) for the ligand.
                lcodes.append((i, lbestCode))
            }

            // Sort the codes for the fragments. Each neighbor symmetry class is sorted separately.
            var firstIndex: UInt = 0 
            lastSymClass = state.symmetry_classes[nbrs[0].getIndex()]
            for i in 0..<nbrs.count {
                let symClass = state.symmetry_classes[nbrs[i].getIndex()]
                if symClass != lastSymClass {
                    lcodes[Int(firstIndex)..<i].sort(by: sortCode2)
                    firstIndex = UInt(i)
                }
                lastSymClass = symClass
            }
            // Make sure to sort the last set too.
            lcodes[Int(firstIndex)...].sort(by: sortCode2)
            // Label the neighbor atoms by updating code.
            var atoms: [MKAtom] = []
            var nextLbl: UInt = label + 1
            for l in 0..<lcodes.count {
                //print_vector("LIG CODE", lcodes[l].second.code);
                let atom = nbrs[lcodes[l].0]
                code.add(current, atom)
                code.labels[atom.getIndex()] = nextLbl
                atoms.append(atom)
                nextLbl += 1
            }
            // Convert the labels from the ligands to labels in the whole fragment.
            let ligandSize = ligandSizes.max()!
            for lbl in 1..<ligandSize {
                for l in 0..<lcodes.count {
                    if lbl >= ligandSizes[lcodes[l].0] { continue }
                    var atom: MKAtom? = nil
                    for i in 0..<mol.numAtoms() {
                        if lcodes[l].1.labels[i] == lbl {
                            atom = mol.getAtom(i+1)
                            break
                        }
                    }
                    if atom == nil { continue }
                    var atomNbrs: [MKAtom] = []
                    guard let nbrs = atom!.getNbrAtomIterator() else { continue }
                    for nbr in nbrs {
                        if lcodes[l].1.labels[nbr.getIndex()] > lbl {
                            if atoms.contains(nbr) { continue }
                            atomNbrs.append(nbr)
                        }
                    }

                    atomNbrs.sort(using: SortAtomsAscending(ranks: lcodes[l].1.labels))

                    for i in 0..<atomNbrs.count {
                        code.add(atom!, atomNbrs[i])
                        code.labels[atomNbrs[i].getIndex()] = nextLbl
                        atoms.append(atomNbrs[i])
                        nextLbl += 1
                    }
                }
            }
            
            // Recurse...
            canonicalLabelRecursive(current, nextLbl-1, &timeout, &bestCode, &state)

            // Backtrack...
            for j in 0..<atoms.count {
                _ = code.atoms.popLast()
                _ = code.bonds.popLast()
                _ = code.from.popLast()
                code.labels[atoms[j].getIndex()] = 0
            }
        } else {
            // There are no duplicated symmetry classes for the neighbor atoms.
            // Label each neighbor atom in sorted sequence.
            var lbl: UInt = label
            for i in 0..<nbrs.count {
                lbl += 1
                code.add(current, nbrs[i])
                code.labels[nbrs[i].getIndex()] = lbl
            }

            //  Recurse... 
            canonicalLabelRecursive(current, lbl, &timeout, &bestCode, &state)

            // Backtrack...
            for i in 0..<nbrs.count {
                _ = code.atoms.popLast()
                _ = code.bonds.popLast()
                _ = code.from.popLast()
                code.labels[nbrs[i].getIndex()] = 0
            }
        }
    }
    
    static func findOrbits(_ orbits: inout Orbits, _ mol: MKMol, _ labels1: [UInt], _ labels2: [UInt]) {

        var newOrbits: Orbits = []
        var visited: [Bool] = Array(repeating: false, count: labels1.count)

        // Construct the orbits for automorphic permutation labels1 -> labels2
        for i in 0..<labels1.count {
            if visited[i] { continue }
            let vi = labels1[i]
            if vi == labels2[i] { continue }
            var j = i 
            var vj = labels2[j]
            newOrbits.resize(newOrbits.count + 1, with: Orbit())
            if var last = newOrbits.last {
                last.append(mol.getAtom(j+1)!)
            } else {
                fatalError("Could not get last element??")
            }
            visited[i] = true
            while vi != vj {
                for k in i..<labels1.count {
                    if vj == labels1[k] {
                        j = k 
                        vj = labels2[j]
                        if var last = newOrbits.last {
                            last.append(mol.getAtom(j+1)!)
                        } else {
                            fatalError("Could not get last element??")
                        }
                        visited[j] = true
                        break
                    }
                }
            }
        }

        // Merge the orbits with previously found orbits.
        for j in 0..<newOrbits.count {
            var merge: Bool = false 
            for k in 0..<orbits.count {
                for l in 0..<newOrbits[j].count {
                    if orbits[k].contains(newOrbits[j][l]) {
                        merge = true
                    }
                }
                if merge {
                    for l in 0..<newOrbits[j].count {
                        if !orbits[k].contains(newOrbits[j][l]) {
                            orbits[k].append(newOrbits[j][l])
                        }
                    }
                    break
                }
            }
            if !merge {
                orbits.append(newOrbits[j])
            }
        }
        newOrbits.removeAll()
        // var vivisited: [Bool] = Array(repeating: false, count: orbits.count)
        // TODO: check if this array is really needed and it was a bug in the original code
        for i in 0..<orbits.count {
            if visited[i] { continue }
            visited[i] = true
            var newOrbit: Orbit = orbits[i]
            for j in i..<orbits.count {
                if visited[j] { continue }
                newOrbit.sort { atom1, atom2 in
                    atom1.getIndex() < atom2.getIndex()
                }
                orbits[j].sort { atom1, atom2 in
                    atom1.getIndex() < atom2.getIndex()
                }
                var result = Set(newOrbit).intersection(orbits[j])
                if result.isEmpty {
                    continue
                }
                visited[j] = true
                result.removeAll()
                result = Set(newOrbit).union(orbits[j])
                newOrbit = Array(result)
            }
            newOrbits.append(newOrbit)
        }
        orbits = newOrbits
    }
    
    /**
     * Update the minimum cell representations (mcr).
     */
    static func updateMCR(_ mcr: inout MKBitVec, _ orbits: inout Orbits, _ bestLabels: [UInt]) {
        for i in 0..<bestLabels.count {
            mcr.setBitOn(UInt32(i+1))
        }
        for j in 0..<orbits.count {
            orbits[j].sort(using: SortAtomsAscending(ranks: bestLabels))
            for k in 0..<orbits[j].count {
                if k != 0 {
                    mcr.setBitOn(UInt32(orbits[j][k].getIdx()))
                }
            }
        }
    }
    /**
     * This is the recursive function implementing the labeling algorithm
     * outlined above (steps 1-3). This function works for single connected
     * fragments and uses @p state to hold most information. The @p bestCode
     * is the current greatest code (if any). given the @p current atom and
     * the last @p label, this function labels the neighbor atoms with
     * label+1, label+2, ...
     */
    static func canonicalLabelRecursive(_ current: MKAtom, _ label: UInt, _ timeout: inout Timeout, _ bestCode: inout FullCode, _ state: inout State) {
        guard let mol = current.getParent() else {
            fatalError("Could not get parent")
        }
        var code = state.code
        if state.backtrackDepth > 0 {
            if code.atoms.count > state.backtrackDepth {
                // Log here??
                return
            }
            if code.atoms.count == state.backtrackDepth {
                state.backtrackDepth = 0
                print("BACKTRACK DONE")
                return
            } else if code.atoms.count < state.backtrackDepth {
                state.backtrackDepth = 0
            }
        }
        
        // Check if there is a full mapping.
        if label == state.fragment.countBits() {
            // Complete the canonical code
            var fullcode = FullCode()
            completeCode(mol, &fullcode, &state)
            
            // Check previously found codes to find redundant subtrees.
            for i in stride(from: state.identityCodes.count, to: 0, by: -1) {
                if fullcode.code == state.identityCodes[i-1].code {
                    //
                    // An explicit automorphism has been found.
                    //
                    var v1: [UInt] = Array(repeating: 0, count: fullcode.labels.count)
                    for j in 0..<fullcode.labels.count {
                        if fullcode.labels[j] != 0{
                            v1[Int(fullcode.labels[j]) - 1] = UInt(j + 1)
                        }
                    }

                    var v2: [UInt] = Array(repeating: 0, count: state.identityCodes[i-1].labels.count)
                    for j in 0..<state.identityCodes[i-1].labels.count {
                        if state.identityCodes[i-1].labels[j] != 0 {
                            v2[Int(state.identityCodes[i-1].labels[j]) - 1] = UInt(j + 1)
                        }
                    }

                    state.backtrackDepth = 0 
                    for j in 0..<v1.count {
                        if v1[j] != v2[j] {
                            //  Yield backtrackDepth 
                            return 
                        }
                        if v1[j] != 0 {
                            state.backtrackDepth += 1
                        }
                    }
                }
            }
            if fullcode.code == bestCode.code {
                updateMCR(&state.mcr, &state.orbits, bestCode.labels)
                findOrbits(&state.orbits, mol, fullcode.labels, bestCode.labels)
            } else if fullcode > bestCode {
                // if fullcode is greater than bestCode, we have found a new greatest code
                bestCode = fullcode
            }
            if state.identityCodes.count < MAX_IDENTITY_NODES {
                state.identityCodes.append(FullCode())
                if var last = state.identityCodes.last {
                    swap(&last.labels, &fullcode.labels)
                    swap(&last.code, &fullcode.code)
                } else {
                    fatalError("Could not grab last element we just added")
                }
            } else {
                swap(&state.identityCodes[MAX_IDENTITY_NODES - 1].labels, &fullcode.labels)
                swap(&state.identityCodes[MAX_IDENTITY_NODES - 1].code, &fullcode.code)
            }
            
            return
        }
        
//         avoid endless loops
        if (Date().timeIntervalSince(timeout.startTime) > timeout.maxTime) {
            return
        }

        // If there is a bestCode and only one labeling is required, return.
        if state.onlyOne && !bestCode.code.isEmpty {
            return
        }

        // Abort early if this will not lead to a greatest canonical code. The
        // code.from vector is compared with the elements in the bestCode.
        if code < bestCode {
            return
        }

        // Find the neighbors of the current atom to assign the next label(s).
        var nbrs: [MKAtom] = []
        var nbrSymClasses: [UInt] = []
        guard let currentNbrs = current.getNbrAtomIterator() else { 
            fatalError("Could not get neighbor iterator???")
        }

        for nbr in currentNbrs {
            // skip atoms not in the fragment 
            if !state.fragment.bitIsSet(nbr.getIdx()) {
                continue
            }
            // skip atoms already labeled
            if code.labels[nbr.getIdx()] != 0 {
                continue
            }
            guard let bond = mol.getBond(current, nbr) else {
                // TODO: maybe throw an error here
                break
            }
            // Ugly, but it helps...
            if !isFerroceneBond(bond) {
                nbrSymClasses.append(state.symmetry_classes[nbr.getIndex()])
                nbrs.append(nbr)
            }
        }

        if nbrs.isEmpty {
            // If there are no neighbor atoms to label, recurse with the next
            // current atom.
            let nextLabel = code.labels[current.getIndex()] + 1
            for i in 0..<code.labels.count {
                if code.labels[i] == nextLabel {
                    canonicalLabelRecursive(mol.getAtom(i+1)!, label, &timeout, &bestCode, &state)
                    return
                }
            }
            return
        }

        if !current.isInRing() {
            // Optimization to avoid generating the n! permutations. If the current
            // atom is not in a ring, the neighbor atoms with the same symmetry
            // are not connected and it is possible to label without considering
            // the rest of the molecule. The found ligand labels can be sorted and
            // the labels can be calculated using offsets.
            //
            // This modifies the final canonical code since it produces a different
            // canonical order but although different, the code is still canonical.
            labelFragments(current, &nbrs, label, &timeout, &bestCode, &state)
        } else {
            // Create all possible labelings for the neighbors.
            //
            // Example: The current atom has 4 neighbors to label.
            //
            // neighbor atom   : 1 2 3 4
            // symmetry classes: 3 2 2 1
            //
            // The first neighbor with symmetry class 3 is added to allOrderedNbrs.
            //
            // allOrderedNbrs = [ [1] ]
            //
            // The two neighbor atoms with symmetry class 2 are added next. All
            // permutations are added (2 in this case).
            //
            // allOrderedNbrs = [ [1,2,3], [1,3,2] ]
            //
            // The last atom is similar to the first one.
            //
            // allOrderedNbrs = [ [1,2,3,4], [1,3,2,4] ]
            //
            // Note: If there are atoms with a large number of neighbor atoms
            // with the same symmetry class (n), this can result in a large number
            // of permutations (n!).
            var allOrderedNbrs: [[MKAtom]] = [[MKAtom]].init(repeating: [], count: 1)
            while !nbrs.isEmpty {
                // Select the next nbr atoms with highest symmetry classes.
                guard let maxSymClass = nbrSymClasses.min() else { return }
                var finalNbrs: [MKAtom] = []
                for i in 0..<nbrs.count {
                    if nbrSymClasses[i] == maxSymClass {
                        finalNbrs.append(nbrs[i])
                    }
                }
                // Remove the selected atoms from nbrs and nbrSymClasses (this could be made more efficient)
                for i in 0..<finalNbrs.count {
                    if let idx = nbrs.firstIndex(of: finalNbrs[i]) {
                        nbrs.remove(at: idx)
                    }
                    if let idx = nbrSymClasses.firstIndex(of: maxSymClass) {
                        nbrSymClasses.remove(at: idx)
                    }
                }

                if finalNbrs.count == 1 {
                    // If there is only one atom with the same symmetry class, label it
                    // and select the next group of neighbor atoms with the same symmetry
                    // class.
                    for i in 0..<allOrderedNbrs.count {
                        allOrderedNbrs[i].append(finalNbrs[0])
                    }
                } else {
                    
                    // Sort the atoms lexicographically. TODO: make sure this works
                    finalNbrs.sort { (a, b) -> Bool in
                        return a.getIdx() < b.getIdx()
                    }

                    // Copy the current labelings for the neighbor atoms.
                    let allOrderedNbrsCopy = allOrderedNbrs


                    // Add the first permutation for the neighbor atoms.
                    for j in 0..<allOrderedNbrs.count {
                        for i in 0..<finalNbrs.count {
                            allOrderedNbrs[j].append(finalNbrs[i])
                        }
                    }

                    // Add the other permutations.
                    for perm in finalNbrs.permutations() {
                        if state.mcr.bitIsSet(perm[0].getIdx()) {
                            for j in 0..<allOrderedNbrsCopy.count {
                                allOrderedNbrs.append(allOrderedNbrsCopy[j])
                                for i in 0..<perm.count {
                                    if var last = allOrderedNbrs.last {
                                        last.append(perm[i])
                                    }
                                }
                            }
                        }
                    }

                }// finalNbrs.size() != 1
            }// while (!nbrs.empty())

            for i in 0..<allOrderedNbrs.count {
                // convert the order stored in allorderedNbrs to labels.
                var lbl = label
                for j in 0..<allOrderedNbrs[i].count {
                    lbl += 1
                    code.add(current, allOrderedNbrs[i][j])
                    code.labels[allOrderedNbrs[i][j].getIndex()] = lbl
                } 

                // Recurse...
                canonicalLabelRecursive(current, lbl, &timeout, &bestCode, &state)
                
                // backtrack...
                for j in 0..<allOrderedNbrs[i].count {
                    _ = code.atoms.popLast()
                    _ = code.bonds.popLast()
                    _ = code.from.popLast()
                    code.labels[allOrderedNbrs[i][j].getIndex()] = 0
                }

                //  Optimization 
                if state.backtrackDepth != 0 {
                    if code.atoms.count <= state.backtrackDepth {
                        print("BACKTRACK DONE 3")
                        state.backtrackDepth = 0
                    }
                }
            }
        }// current is in ring end 
    }
    
    /**
     * Select an initial atom from a fragment to assign the first label.
     */
    static func findStartAtoms(_ mol: MKMol, _ fragment: MKBitVec, _ symmetry_classes: [UInt]) -> [MKAtom] {
        // find the symmetry class in the fragment using criteria 
        var ranks: [UInt] = []
        for i in 0..<mol.numAtoms() {
            if !fragment.bitIsSet(i+1) {
                continue
            }
            guard let atom = mol.getAtom(i+1) else { 
                // Thow error here? TODO: 
                continue 
            }
            var rank = 10000 * Int(symmetry_classes[i])
            rank += 1000 * atom.getSpinMultiplicity()
            rank += 10 * Int((atom.getFormalCharge() + 7))
            rank += Int(totalHydrogenCount(atom))
            ranks.append(UInt(rank))
        }
        guard let lowestRank = ranks.min() else {
            fatalError("ranks are empty")
        }
        var result: [MKAtom] = []
        for i in 0..<mol.numAtoms() {
            if !fragment.bitIsSet(i+1) {
                continue
            }
            guard let atom = mol.getAtom(i+1) else { 
                // Thow error here? TODO: 
                continue 
            }
            var rank = 10000 * Int(symmetry_classes[i])
            rank += 1000 * atom.getSpinMultiplicity()
            rank += 10 * Int((atom.getFormalCharge() + 7))
            rank += Int(totalHydrogenCount(atom))
            if rank == lowestRank {
                result.append(atom)
            }
        }
        return result
    }
    /**
     * Compute the canonical labels for @p mol. This function handles a whole molecule with
     * disconnected fragments. It finds a canonical code for each fragment, sorts these codes
     * and computes labels for the molecule as a whole.
     *
     * This is the CanonicalLabelsImpl entry point.
     */
    
    static func calcCanonicalLabels(_ mol: MKMol, _ symmetry_classes: [UInt], _ canonical_labels: inout [UInt], _ stereoUnits: MKStereoUnitSet, _ mask: MKBitVec, _ stereoFacade: MKStereoFacade?, _ maxSeconds: Int, _ onlyOne: Bool = false) {

        // Handle some special cases 
        if mol.numAtoms() == 0 {
            return 
        }
        if mol.numAtoms() == 1 {
            canonical_labels.resize(1, with: 1)
            return
        }

        // Set all labels to 0.
        canonical_labels.removeAll()
        canonical_labels.resize(mol.numAtoms(), with: 0)

        var metalloceneBonds: [MKBond] = []
        findMetalloceneBonds(&metalloceneBonds, mol, symmetry_classes)

        // Find the (dis)connected fragments.
        var visited: MKBitVec = MKBitVec()
        var fragments: [MKBitVec] = []
        for i in 0..<mol.numAtoms() {
            if (!mask.bitIsSet(i+1) || visited.bitIsSet(i+1)) {
                continue
            }
            fragments.append(getFragment(mol.getAtom(i+1)!, mask, metalloceneBonds))
            visited |= fragments.last!
        }
        // Pre-compute the stereo center information. (See StereoCenter)
        var stereoCenters: [StereoCenter] = []
        if stereoFacade != nil {
            for i in 0..<stereoUnits.count {
                let unit: MKStereoUnit = stereoUnits[i]
                if unit.type == .Tetrahedral {
                    guard let atom = mol.getAtomById(unit.id) else {
                        continue
                    }
                    // Add the StereoCenter indexes.
                    stereoCenters.resize(stereoCenters.count + 1, with: StereoCenter())
                    if var last = stereoCenters.last {
                        last.indexes.append(UInt(atom.getIndex()))
                    } else {
                        fatalError("stereoCenters is empty")
                    }
                    if !stereoFacade!.hasTetrahedralStereo(unit.id.intValue) {
                        continue
                    }
                    guard let config: MKTetrahedralStereo.Config = stereoFacade!.getTetrahedralStereo(unit.id.intValue)?.getConfig() else {
                        fatalError("config is nil")
                    }
                    if !config.specified {
                        continue 
                    }
                    // Add the neighbor atom indexes.
                    let from = mol.getAtomById(config.from_or_towrds.refValue)
                    if from != nil && from!.getAtomicNum() != MKElements.Hydrogen.atomicNum {
                        if var last = stereoCenters.last {
                            last.nbrIndexes1.append(UInt(from!.getIndex()))
                        }
                    } else {
                        if var last = stereoCenters.last {
                            last.nbrIndexes1.append(UInt(Int.max))
                        }
                    }
                    for j in 0..<config.refs.count {
                        let ref = mol.getAtomById(config.refs[j])
                        if ref != nil && ref?.getAtomicNum() != MKElements.Hydrogen.atomicNum {
                            if var last = stereoCenters.last {
                                last.nbrIndexes1.append(UInt(ref!.getIndex()))
                            }
                        } else {
                            if var last = stereoCenters.last {
                                last.nbrIndexes1.append(UInt(Int.max))
                            }
                        }
                    } 
                } else if unit.type == .CisTrans {
                    let bond = mol.getBondById(unit.id)
                    if bond != nil || bond!.isAromatic() {
                        continue
                    }
                    let begin = bond?.getBeginAtom()
                    let end = bond?.getEndAtom()
                    if begin == nil || end == nil {
                        continue
                    }
                    // Add the StereoCenter indexes.
                    stereoCenters.resize(stereoCenters.count + 1, with: StereoCenter())
                    if var last = stereoCenters.last {
                        last.indexes.append(UInt(begin!.getIndex()))
                        last.indexes.append(UInt(end!.getIndex()))
                    } else {
                        fatalError("stereoCenters is empty")
                    }

                    if !stereoFacade!.hasCisTransStereo(unit.id.intValue) {
                        continue
                    }
                    guard let config: MKCisTransStereo.Config = stereoFacade!.getCisTransStereo(unit.id.intValue)?.getConfig() else {
                        fatalError("config is nil")
                    }
                    if !config.specified {
                        continue 
                    }

                    // Add the neighbor atom indexes.
                    for j in 0..<config.refs.count {
                        let ref = mol.getAtomById(config.refs[j])
                        let r = (ref !=  nil && ref?.getAtomicNum() != MKElements.Hydrogen.atomicNum) ? UInt(ref!.getIndex()) : UInt(Int.max)
                        if stereoCenters.last!.nbrIndexes1.count < 2 {
                            if var last = stereoCenters.last {
                                last.nbrIndexes1.append(r)
                            }
                        } else {
                            if var last = stereoCenters.last {
                                last.nbrIndexes2.append(r)
                            }
                        }
                    }
                }
            }
        }

        // Find the canonical code for each fragment.
        var fcodes: [CanonicalLabelsImpl.FullCode] = []
        for f in 0..<fragments.count {
            var fragment = fragments[f]

            // Select the first atom.
            var startAtoms: [MKAtom] = findStartAtoms(mol, fragment, symmetry_classes)

            var timeout = Timeout(TimeInterval(maxSeconds))
            var bestCode: CanonicalLabelsImpl.FullCode = CanonicalLabelsImpl.FullCode()
            var identityCodes: [CanonicalLabelsImpl.FullCode] = []
            var orbits: Orbits = Orbits()
            var mcr: MKBitVec = MKBitVec()

            for i in 0..<startAtoms.count {
                let atom = startAtoms[i]
                // Start labeling of the fragment.
                var state: State = State(symmetry_classes, fragment, stereoCenters, identityCodes, orbits, mcr, onlyOne)
                //if (!state.mcr.BitIsSet(atom->GetIdx()) && atom->IsInRing())
                //  continue; // Original impl, might not be needed
                state.code.add(atom)
                state.code.labels[atom.getIndex()] = 1
                canonicalLabelRecursive(atom, 1, &timeout, &bestCode, &state)
            }
            // Throw an error if the timeout is exceeded.
            if Date().timeIntervalSince(timeout.startTime) > timeout.maxTime {
                // throw MKError.timeout
                print("ERROR: maximum time exceeded...")
                return 
            }
            // Store the canonical code for the fragment.
            fcodes.append(bestCode)
        }
        // Sort the codes for the fragments.
        fcodes.sort(by: sortCode)
        // Construct the full labeling from the sorted fragment labels.
        var offset: UInt = 0
        for f in 0..<fcodes.count {
            if fcodes[f].labels.count == 0 {
                continue // defensive programming ?
            }
            var max_label: UInt = 0
            for i in 0..<mol.numAtoms() {
                if fcodes[f].labels[i] != 0 {
                    canonical_labels[i] = fcodes[f].labels[i] + offset
                    max_label = max(max_label, canonical_labels[i])
                }
            }
            offset = max_label
        }
    }
}

/*
   * See below for detailed description of parameters.
   *
   * The main purpose of this function is calling CanonicalLabelsImpl::CalcCanonicalLabels
   * with the correct parameters regarding stereochemistry.
   */
func canonicalLabels(_ mol: MKMol, _ symmetry_classes: inout [UInt], _ canonical_labels: inout [UInt], _ mask: MKBitVec = MKBitVec(), _ maxSeconds: Int = 5, _ onlyOne: Bool = false) {
    // make sure the mask is valid: no mask = all atoms
    var maskCopy = mask 
    if maskCopy.countBits() == 0 {
        for atom in mol.getAtomIterator() {
            maskCopy.setBitOn(UInt32(atom.getIdx()))
        }
    }
    if onlyOne {
        // Only one labeling requested. This results in canonical labels that do not
        // consider stereochemistry. Used for finding stereo centers with automorphisms.
        CanonicalLabelsImpl.calcCanonicalLabels(mol, symmetry_classes, &canonical_labels, MKStereoUnitSet(), maskCopy, nil, maxSeconds, true)
    } else {
        var metalloceneBonds: [MKBond] = []
        findMetalloceneBonds(&metalloceneBonds, mol, symmetry_classes)

        // Check if there are specified stereo centers 
        var hasAtLeastOneDefined: Bool = false
        let sf: MKStereoFacade = MKStereoFacade(mol, m_perceive: false)
        for atom in mol.getAtomIterator() {
            for bond in atom.getBondIterator()! {
                if metalloceneBonds.contains(bond) {
                    continue
                }
                if sf.hasTetrahedralStereo(atom.getId().intValue) {
                    if sf.getTetrahedralStereo(atom.getId().intValue)?.getConfig().specified ?? false {
                        hasAtLeastOneDefined = true
                        break
                    }
                }
            }
        }
        
        for bonds in mol.getBondIterator() {
            if sf.hasCisTransStereo(bonds.getId().intValue) {
                if sf.getCisTransStereo(bonds.getId().intValue)?.getConfig().specified ?? false {
                    hasAtLeastOneDefined = true
                    break
                }
            }
        }

        if !hasAtLeastOneDefined {
            // If there are no specified stereo centers, we don't need to find stereogenic units.
            CanonicalLabelsImpl.calcCanonicalLabels(mol, symmetry_classes, &canonical_labels, MKStereoUnitSet(), maskCopy, nil, maxSeconds)
            return
        }
        // Find the stereogenic units
        var stereoUnits: MKStereoUnitSet = findStereogenicUnits(mol, symClasses: &symmetry_classes)
        // Mark all invalid stereo data as unspecified
        if var stereoData: [MKGenericData] = mol.getAllData(.StereoData) {
            for data in stereoData {
                let type: MKStereo.TType = (data as? MKStereoBase)!.getType()
                if type == .Tetrahedral {
                    let ts: MKTetrahedralStereo = data as! MKTetrahedralStereo
                    var config = ts.getConfig()
                    var valid = true
                    if !ts.isValid() {
                        valid = false
                    }
                    if let center = mol.getAtomById(config.center) {
                        if !isTetrahedral(center, stereoUnits) {
                            valid = false
                        }
                    } else {
                        valid = false
                    }
                    if !valid {
                        config.specified = false 
                        ts.setConfig(config)
                    }
                }
                if type == .CisTrans {
                    let ct: MKCisTransStereo = data as! MKCisTransStereo
                    let config = ct.getConfig()
                    var valid = true
                    if !ct.isValid() {
                        valid = false
                    }
                    let beginAtom = mol.getAtomById(config.begin)
                    let endAtom = mol.getAtomById(config.end)
                    if beginAtom == nil || endAtom == nil {
                        valid = false
                    } else {
                        if let bond = mol.getBond(beginAtom!, endAtom!) {
                            if !isCisTrans(bond, stereoUnits) {
                                valid = false
                            }
                        } else {
                            valid = false
                        }
                    }
                    if !valid {
                        config.specified = false 
                        ct.setConfig(config)
                    }
                }
            }

            // Determine stereochemistry from coordinates if needed
            if !mol.hasChiralityPerceived() {
                switch mol.getDimension() {
                case 2: 
                    mol.deleteData(.StereoData)
                    tetrahedralFrom2D(mol, stereoUnits)
                    cisTransFrom2D(mol, stereoUnits)
                case 3:
                    mol.deleteData(.StereoData)
                    tetrahedralFrom3D(mol, stereoUnits)
                    cisTransFrom3D(mol, stereoUnits)
                default:
                    tetrahedralFrom0D(mol, stereoUnits)
                    cisTransFrom0D(mol, stereoUnits)
                }
            }

            // The mol->DeleteData() above deleted some of the data that the OBStereoFacade sf()
            // is using.  We have to create a new one.  (CJ - fixes a bug.) // TODO: figure out if we really need this
            let newsf: MKStereoFacade = MKStereoFacade(mol, m_perceive: false)
            // Start the labeling process
            CanonicalLabelsImpl.calcCanonicalLabels(mol, symmetry_classes, &canonical_labels, stereoUnits, maskCopy, newsf, maxSeconds)
        }
    }
    // if the labeling failed, just return the identity labels to avoid craches
    if canonical_labels.isEmpty {
        for i in 0..<symmetry_classes.count {
            canonical_labels.append(UInt(i+1))
        }
    }
}

/* -*-C++-*-
+======================================================================
|
|                       !!!DEPRECATED!!!
|
|
| AUTHOR: Craig A. James, eMolecules, Inc.
|
| DESCRIPTION: CANONICALIZATION OF SMILES
|
|       This is a specialized SMILES canonicalization algorithm.  Although
|       it can be applied in the standard fashion to a whole molecule,
|       its real job is to generate canonical SMILES for fragments, or
|       "subsets", of the atoms of a molecule.
|
|       For example, consider the first three atoms of Oc1ccccc1.  With
|       a "normal" SMILES canonicalizer, you couldn't generate a SMILES
|       for Occ, because it's not a valid molecule.  However, this system
|       can do exactly that, by taking both the whole molecule (which
|       retains the aromaticity), and a "subset" bitmap that specifies
|       which atoms are to be included in the SMILES.
|
|       Canonicalization is carried out per Weininger et al (J. Chem.
|       Inf. Comput. Sci., Vol. 29, No. 2, 1989, pp 97-101), with some
|       modifications to handle bond symmetries not foreseen by Weininger
|       in that paper.
|
|       WARNING - KNOWN BUG: These functions make use of a bitmap vector
|       to represent a "fragment" -- a subset of the atoms in a molecule.
|       But this means the bonds of the fragment are implicit, not explicit,
|       which is incorrect.  For example, if you want to break one bond of
|       cyclehexane (C1CCCCC1), all six atoms will still be there, so the
|       "fragment" will be cyclic.  This is relevant when generating fragment
|       SMILES for ring systems where breaking a bond can reduce the number
|       of ring without removing any atoms.  We need to add a pair of bit
|       vectors, the atoms AND the bonds, to represent a fragment.  (Note
|       that this is also an ambiguity in OpenBabel itself, which represents
|       a ring as a set of atoms. This is only valid if the ring is a member
|       of a SSSR.)
+======================================================================
*/


// MARK: - Canonicalization Note*
/**
* @page canonical_code_algorithm Canonical Coding Algorithm
*
* @section canonical_introduction Introduction
* The aim of a canonical coding algorithm is to assign unique labels to
* atoms regardless of their input order. In practise, these labels or orders
* are atom indexes resulting from the order in which the atoms are defined
* in an input file. Although most chemical file formats could be used to
* store canonical ordered molecules, canonical single line notations are used
* more often since they allow two canonical molecules to be compared using a
* string comparison. Two well known examples are canonical smiles and inchi.
* While the canonical smiles for the same molecule should always be the same
* when using a specific implementation (i.e. toolkits), these canonical
* smiles are usually not transferable between implementations. There is only
* one inchi implementation and the canonical code, although not formally
* specified is always the same. A typical use case for canonical codes is to
* determine if a molecule is already in a database.
*
* The rest of this page documents the OpenBabel canonical coding algorithm.
*
* @section Topological Symmetry
* The starting point of the canonical coding algorithm is the topological
* symmetry or graph symmetry. This symmetry can be expressed by assigning
* symmetry classes to atoms. Symmetric atoms have the same symmetry class.
* The term color is often used in graph theory as a synonym for symmetry
* classes.
*
* Symmetry classes are used as graph invariants to order atoms. The symmetry
* classes can be used in two ways:
*
* - lowest symmetry class = lowest label (i.e. 1)
* - highest symmetry class = lowest label
*
* Both methods are valid. OpenBabel uses the second method. In the following
* sections there are more of these decisions that have to be made. The first,
* choice is used as example to illustrate that it doesn't make any real
* difference as long as the same choice is used every time. When trying to
* reproduce OpenBabel canonical smiles the same choices or conventions have
* to be used.
*
* The used initial graph invariant is specified below and the iteration
* algorithm is similar to the original Morgan algorithm. The OBGraphSym class
* should be consulted for a detailed description.
    @verbatim
        0                   1                   2                   3
        1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2
    GI = D D D D D D D D D D D V V V V A R N N N N N N N B B B B C C C C

    GI: Graph invariant (32 bits)
    D: Graph theoretical distance (10 bits)
    V: Heavy valence (4 bits)
    A: Aromaticity (1 bit)
    R: Ring atom (1 bit)
    N: Atomic number (7 bits)
    B: Heavy bond sum (4 bits)
    C: Formal charge + 7 (4 bits)
    @endverbatim
*
* (FIXME: ugly smiles)
*
* The symmetry classes encode various attributes for each atom. Consult the
* OBGraphSym documentation for details. (note: when trying to reproduce
* OpenBabel canonical smiles, make sure these graph invariants are the same)
*
* @section canonical_labeling Labeling
* The labeling starts by selecting the initial atom using a set of rules to
* ensure it starts at terminal atoms. (FIXME) Label 1 is assigned to this atom
* and the labeling continues using this atom as the current atom:
*
* -# Find the neighbor atoms of the current atom to label next. Already labeled
*   and atoms not in the fragment are ignored. If there are no such atoms, go
*   to 3.
* -# Sort the neighbors by their symmetry classes. Each group of atoms with the
*    same symmetry class is now handled in sequence. If there is only one atom
*    in the group, this atom is labeled and the process continues with the next
*    group. If there are multiple (n) atoms in a group, all n! permutations of
*    the group are tried. (note: The Optimization section below mentions some
*    exceptions which affect the final result. These optimizations are part of
*    the algorithm)
* -# When all neighbor atoms of the current atom are labeled, the next atom
*    (i.e. the atom with current_label+1 label) is considered the current atom
*    and step 1 through 3 are repeated. If there is no next atom, the labeling
*    is done.
*
* The algorithm just described results in a set of labellings that are candidates
* for the canonical labeling. For each labeling, a code is generated to select
* the canonical code and the associated labels.
*
* @section canonical_code Canonical Code
* A canonical code can be anything encoding the molecular graph with all the
* properties that should be canonicalized. The code should be invariant to the
* original atom order. At this point, it may seem that this is similar to the
* graph invariants from above (i.e. symmetry classes). However, generating these
* codes allows for additional attributes to be included. A simple example are
* isotopes. More importantly, stereochemistry is included in these codes which
* is an attribute that does not depend an isolated atom.
*
* The actual code is a list of numbers (unsigned int) constructed by joining
* smaller codes. These smaller codes usually encode specific attributes.
*
* @subsection canonical_from FROM code
* When labeling the molecule as described earlier, all but the initial atom in
* a connected fragment have a "from" atom (i.e. current atom in the algorithm).
* The FROM code contains the label of the "from" atom, for each atom (ignoring
* the initial atom) ordered by their assigned label. The FROM code is a
* spanning tree containing all bonds that are not ring closures.
*
* @subsection canonical_closure CLOSURE code
* The ring closure bonds are easily identified while creating the FROM code.
* The CLOSURE code contains two labels for each closure bond. The labels are
* sorted in the following way:
*
* - [1 6] < [6 1] (i.e. 1 6 is used, individual bond sorting)
* - [1 3][1 4] < [1 4][1 3] (i.e. 1 3 1 4 is used, sorting bonds)
* - [1 3][5 7] < [5 7][1 3]
*
* The square brackets are only to indicate individual bonds, the code is just
* a single list.
*
* @subsection canonical_atomtype ATOM-TYPES & ISOTOPES code
* The ATOM-TYPES code is simply the atomic number for each atom ordered by
* labels. If the molecule contains isotopes, an additional ISOTOPES code is
* used to correctly canonicalize these molecules.
*
* @subsection canonical_bondtype BOND-TYPES code
* The BOND-TYPES code contains, for each bond, the bond order or 5 if the
* bond is aromatic. The bonds are ordered using the order obtained during
* the FROM code construction. The closure bonds are last, ordered as
* described in the CLOSURE code section.
*
* @subsection canonical_stereo STEREO code
* The STEREO code contains a descriptor for each stereo center. The stereo
* centers are ordered using the central atom(s) (i.e. center for tetrahedral
* stereochemistry, lowest label for cistrans).
*
* @subsection canonical_canonical Canonical Code
* The canonical code is the greatest code (elements are compared, the first
* difference determines order). The complete canonical code can also be seen
* as layers. The first layer is the FROM code together with the CLOSURE code.
* Together these codes describe the graph without node (atom) or edge (bond)
* properties. The ATOM-TYPES code is the layer adding node properties. The
* optional ISOTOPES code also encodes node properties. Edge properties are
* specified in the BOND-TYPES code. The previous codes together are enough
* to specify the constitution of a molecule. Additional layers such as the
* STEREO code could be ignored for certain use cases. For example, if two
* molecules have the same canonical code up to the STEREO part, they are
* stereoisomers.
*
* @section canonical_fragments Disconnected Fragments
* Disconnected fragments are handled individually and the resulting canonical
* codes are sorted. The final canonical code for the whole molecule, are the
* sorted canonical codes joined together.
*
* @section canonical_optimization Optimizations
* This section documents optimizations to reduce the number of states which
* affect the final result. The aim of these optimizations is to avoid the n!
* permutations in step 2 of the labeling algorithm. Some of the optimizations
* affect the final result and should be implemented when trying to reproduce
* the same canonical code.
*
* @subsection canonical_opt1 Optimization 1
* If the current atom is not in a ring, it follows that the bonds to all it's
* neighbors are not ring bonds. It also follows that there is no path between
* two neighbor atoms with the same symmetry class without going through the
* current atom. Therefore, these ligands (i.e. the fragments that would be
* created by deleting the current atom) can be labeled individually. For
* example, for an atom with 4 neighbors with the same symmetry class, the
* number of permutations to label is reduced from 4! * 4 = 96 to 4. The times
* 4 is because after generating the 24 permutations, the ligands are not
* labeled yet.
*
* @note This optimization <b>affects the final result</b>.
*
* @subsection canonical_opt2 Optimization 2 (not implemented)
* If the bond(s) to the neighbor atom(s) with the same symmetry classes is not
* a ring bond, optimization 1 could be used. This is the more general version
* of optimization 1.
*
* @note This optimization <b>affects the final result</b> but is <b>not
* implemented</b>.
*
* @subsection canonical_opt3 Optimization 3 (not implemented)
* If the current atom is a stereo center, only permutations with the highest
* descriptor produce a greatest canonical code. This does not change the final
* result and can be implemented without changing canonical smiles etc.
*
* @subsection canonical_opt4 Other optimizations
* In essence, the canonical coding algorithm is very similar to the algorithm
* from the well-known nauty package [1, 2]. However, the optimizations mentioned
* above are far easier to implement, document and result in acceptable
* performance for most structures. This enables anyone to reproduce OpenBabel
* canonical smiles without the absolute requirement to implement the complex
* optimizations from nauty. Since the optimizations from the nauty paper do
* not have any affect on the final result, these can be implemented as required.
*
* The OpenBabel algorithm uses the molecule graph as the main object while
* searching for a canonical code. Nauty starts from a partition of the set
* of nodes which is refined. This is not the essential part though. The
* similarity between the algorithms is that both algorithms share the
* concept of a search tree. A completed canonical candidate code is a leaf
* node in this search tree. If a second labeling is found with the same
* canonical code, there exists a permutation mapping the first labels to
* the second. This permutation is an automorphism that takes all attributes
* of the canonical code into account. The found automorphism is not only
* an automorphism of the molecule but also an isomorphism of the search tree.
* These accidentally discovered automorphisms can be used to reduce the
* number of states to consider in various ways.
*
* @subsubsection canonical_opt5 Optimization 4
*
*
*
*
*
    @verbatim
    [1] Brendan D. McKay, Backtrack programming and isomorphism rejection on
        ordered subsets, ARS Combinatoria, Vol. 5, 1979, pp. 65-99
        http://cs.anu.edu.au/~bdm/papers/backtrack1978.pdf
    [2] Brendan D. McKay, Practical graph isomorphism, Congressus Numerantium,
        Vol. 30, 1981, pp. 45-87
        http://cs.anu.edu.au/~bdm/papers/pgi.pdf
    @endverbatim
*/
