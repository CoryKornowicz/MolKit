//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/18/23.
//

import Foundation

/**
* Convert any reference to @p atomId in a stereo object to an OBStereo::ImplicitRef.
* This function is called from OBMol::DeleteHydrogens()
* (via OBMol::DeleteHydrogen()) to remove any explicit references to a
* hydrogen atom that has been deleted. However, the code is not specific
* to hydrogen atoms and could be used for other atoms.
*
* @param mol The molecule
* @param atomId The Id of the atom to be converted to an OBStereo::ImplicitRef
* @since version 2.3
*/
func stereoRefToImplicit(_ mol: MKMol, _ atomId: Ref) {
    guard let vdata: [MKGenericData] = mol.getAllData(.StereoData) else {
        print("ERROR: no stereo data found in stereoRefToImplicit")
        return
    }
    for data in vdata {
        guard let datatype: MKStereo.TType = (data as? MKStereoBase)?.getType() else { continue }
        // maybe we should just unset the stereochemistry instead 
        if datatype != .CisTrans && datatype != .Tetrahedral {
            print("This function should be updated to handle additional stereo types.\nSome stereochemistry objects may contain explicit refs to hydrogens which have been removed.")
            continue
        }
        // replace and references to atomId with ImplicitRef
        if datatype == .CisTrans {
            let ct = data as! MKCisTransStereo
            let ct_config = ct.getConfig()
            ct_config.refs.replace(atomId, with: .ImplicitRef)
            ct.setConfig(ct_config)
        } else if datatype == .Tetrahedral {
            let ts = data as! MKTetrahedralStereo
            var ts_config = ts.getConfig()
            if ts_config.from_or_towrds == atomId {
                ts_config.from_or_towrds = .from(.ImplicitRef)
            }
            ts_config.refs.replace(atomId, with: .ImplicitRef)
            ts.setConfig(ts_config)
        }
    }
}
/**
* Convert any reference to an OBStereo::ImplicitRef attached to @p centerId
* in a stereo object to an explicit reference to @p newId.
* This function is called from OBMol::AddHydrogens() and
* OBMol::AddHydrogen() to convert any implicit references to a
* hydrogen atom that has just been added. However, the code is not specific
* to hydrogen atoms and could be used for other atoms.
*
* @param mol The molecule
* @param centerId The Id of the atom to which the new explicit atom is attached
* @param newId The Id of the atom which was previously an OBStereo::ImplicitRef
* @since version 2.4
*/
func implicitRefToStereo(_ mol: MKMol, _ centerId: Ref, _ newId: Ref) {
    guard let vdata: [MKGenericData] = mol.getAllData(.StereoData) else {
        print("ERROR: no stereo data found in stereoRefToImplicit")
        return
    }
    for data in vdata {
        guard let datatype: MKStereo.TType = (data as? MKStereoBase)?.getType() else { continue }
        // maybe we should just unset the stereochemistry instead
        if datatype != .CisTrans && datatype != .Tetrahedral {
            print("This function should be updated to handle additional stereo types.\nSome stereochemistry objects may contain explicit refs to hydrogens which have been removed.")
            continue
        }
        // Replace any references to ImplicitRef (attached to centerId) with newId
        if datatype == .CisTrans {
            let ct = data as! MKCisTransStereo
            let ct_config = ct.getConfig()
            if ct_config.begin == centerId || ct_config.end == centerId {
                // Assumption: the first two refs are on the begin atom, the last two on the end atom
                if ct_config.begin == centerId {
                    ct_config.refs.replaceInRange(ct_config.refs.startIndex...ct_config.refs.startIndex.advanced(by: 2), .ImplicitRef, with: newId)
                }
                if ct_config.end == centerId {
                    ct_config.refs.replaceInRange(ct_config.refs.startIndex.advanced(by: 2)...ct_config.refs.endIndex, .ImplicitRef, with: newId)
                }
                ct.setConfig(ct_config)
            }
        } else if datatype == .Tetrahedral {
            let ts = data as! MKTetrahedralStereo
            var ts_config = ts.getConfig()
            if ts_config.center == centerId {
                if case let .from(ref) = ts_config.from_or_towrds, ref == .ImplicitRef {
                    ts_config.from_or_towrds = .from(newId)
                }
                ts_config.refs.replace(.ImplicitRef, with: newId)
            }
            ts.setConfig(ts_config)
        }
    }
}

/**
* Perform a quick check for tetrahedral stereo centers. Used by
* FindStereogenicUnits to return quickly if there are no stereogenic units.
*/
func mayHaveTetrahedralCenter(_ mol: MKMol) -> Bool {
    for atom in mol.getAtomIterator() {
        if (atom.getHyb() == 3 && atom.getHeavyDegree() >= 3) {
            return true
        }
    }
    return false
}

/**
* Perform a quick check for stereogenic bonds. Used by FindStereogenicUnits
* to return quickly if there are no stereogenic units.
*/
func mayHaveCisTransBond(_ mol: MKMol) -> Bool {
    for bond in mol.getBondIterator() {
        if (bond.getBondOrder() == 2) {
            return true
        }
    }
    return false
}

/**
   * Check if the specified atom is a potential stereogenic atom.
   *
   * Criteria:
   * - sp3 hybridization (or P and sp3d hybridization)
   * - not connected to more than 4 atoms
   * - at least 3 "heavy" neighbors
   *
   * Nitrogen (neutral) is treated as a special case since the barrier of inversion is
   * low in many cases making the atom non-stereogenic. Only bridge-head
   * nitrogen atoms (i.e. nitrogen has 3 neighbors in rings) will be
   * considered stereogenic.
   */
func isPotentialTetrahedral(_ atom: MKAtom) -> Bool {
    // consider only potential steroecenters
    if (atom.getHyb() != 3 && !(atom.getHyb() == 5 && atom.getAtomicNum() == MKElements.Phosphorus.atomicNum)) || 
        (atom.getTotalDegree() > 4 || atom.getHeavyDegree() < 3 || atom.getHeavyDegree() > 4) {
        return false
    }
    // skip non-chiral N
    if (atom.getAtomicNum() == 7 && atom.getFormalCharge() == 0) {
        var nbrRingAtomCount = 0
        guard let nbors = atom.getNbrAtomIterator() else { return false }
        for nbr in nbors {
            if (nbr.isInRing()) {
                nbrRingAtomCount += 1
            }
        }
        if (nbrRingAtomCount < 3) {
            return false
        }
    }
    if (atom.getAtomicNum() == 6) {
        if (atom.getFormalCharge() != 0) {
            return false
        }
        guard let nbors = atom.getNbrAtomIterator() else { return false }
        for nbr in nbors {
            if (nbr.getAtomicNum() == 26 && nbr.getExplicitDegree() > 7) {
                return false
            }
        }
    }
    return true
}
  /**
   * Check if the specified bond is a potential stereogenic bond.
   *
   * Criteria:
   * - must be a double bond
   * - must not be in a ring
   * - both begin and end atom should have at least one single bond
   */
func isPotentialCisTrans(_ bond: MKBond) -> Bool {
    if (bond.getBondOrder() != 2) {
        return false
    }
    if (bond.isInRing()) {
        return false
    }
    if (!bond.getBeginAtom().hasSingleBond() || !bond.getEndAtom().hasSingleBond()) {
        return false
    }
    if (bond.getBeginAtom().getHeavyDegree() == 1 || bond.getEndAtom().getHeavyDegree() == 1) {
        return false
    }
    if (bond.getBeginAtom().getHeavyDegree() > 3 || bond.getEndAtom().getHeavyDegree() > 3) {
        return false
    }
    return true
}

/**
* Check if the specified stereogenic unit is in a fragment.
*/
func isUnitInFragment(_ mol: MKMol, _ unit: MKStereoUnit, _ fragment: MKBitVec) -> Bool {
    if (unit.type == .Tetrahedral) {
        if (fragment.bitIsSet(unit.id.intValue)) {
            return true
        }
    } else if (unit.type == .CisTrans) {
        guard let bond = mol.getBondById(unit.id) else { return false }
        let begin = bond.getBeginAtom()
        let end = bond.getEndAtom()
        if (fragment.bitIsSet(begin.getId().intValue) || fragment.bitIsSet(end.getId().intValue)) {
            return true
        }
    }
    return false
}

/**
   * Check if the specified atom is a tetrahedral center (i.e. there is a Tetrahedral
   * OBStereoUnit in units with the same id)
   */
func isTetrahedral(_ atom: MKAtom, _ units: [MKStereoUnit]) -> Bool {
    for unit in units {
        if (unit.type == .Tetrahedral && unit.id.intValue == atom.getId().intValue) {
            return true
        }
    }
    return false
}

//   /**
//    * Check if the specified bond is a double bond stereocenter (i.e. there is a CisTrans
//    * OBStereoUnit in units with the same id)
//    */
func isCisTrans(_ bond: MKBond, _ units: [MKStereoUnit]) -> Bool {
    for unit in units {
        if (unit.type == .CisTrans && unit.id.intValue == bond.getId().intValue) {
            return true
        }
    }
    return false
}

/**
* Classify the tetrahedral atom using the NeighborSymmetryClasses types.
*/
func classifyTetrahedralNbrSymClasses(_ symClasses: [UInt], _ atom: MKAtom) -> NeighborSymmetryClasses {
    var nbrClasses = [UInt]()
    var nbrClassesCopy = [UInt]()
    var uniqueClasses = [UInt]()
    for nbr in atom.getNbrAtomIterator()! {
        nbrClasses.append(symClasses[nbr.getId().intValue])
    }
    //     // add an implicit ref if there are only 3 explicit
    if nbrClasses.count == 3 {
        nbrClasses.append(UInt(RefValue.ImplicitRef.intValue))
    }
//     // use some STL to work out the number of unique classes
    nbrClassesCopy = nbrClasses // keep copy for count below
    nbrClasses.sort()
    uniqueClasses = nbrClasses.unique()
    switch uniqueClasses.count {
    case 4:
        return .T1234 // e.g. 1 2 3 4
    case 3:
        return .T1123 // e.g. 1 1 2 3
    case 2:
        if nbrClassesCopy.filter({ $0 == uniqueClasses[0] }).count == 2 {
            return .T1122 // e.g. 1 1 2 2
        } else {
            return .T1112 // e.g. 1 1 1 2
        }
    case 1:
        return .T1111 // e.g. 1 1 1 1
    default:
        return .T1111 // e.g. 1 1 1 1
    }
}


  /**
   * Classify the cis/trans bond using the NeighborSymmetryClasses types.
   */
func classifyCisTransNbrSymClasses(_ symClasses: [UInt], _ doubleBond: MKBond, _ atom: MKAtom) -> NeighborSymmetryClasses {
    var nbrClasses = [UInt]()
    for nbr in atom.getNbrAtomIterator()! {
        if (nbr.getId() != doubleBond.getNbrAtom(atom).getId()) {
            nbrClasses.append(symClasses[nbr.getId().intValue])
        }
    }
    if nbrClasses.count == 1 {
        nbrClasses.append(UInt(RefValue.ImplicitRef.intValue))
    }
    if nbrClasses[0] == nbrClasses[1] {
        return .C11 // e.g. 1 1
    } else {
        return .C12 // e.g. 1 2
    }
}

/**
* Merge the rings in a molecule and return the result as OBBitVec objects.
* Rings are merged if they share at least one atom (e.g. bridged, spiro,
* adjacent, ...).
*/
func mergeRings(_ mol: MKMol, _ symClasses: [UInt]) -> [MKBitVec] {
    var result = [MKBitVec]()
    let rings = mol.getSSSR()
    for ring in rings {
        // check if ring shares atom with previously found ring
        var found = false
        for (_, r) in result.enumerated() {
            // foreach ring atom
            var shared = [UInt]()
            for path in ring._path {
            // check if the ring atom is in the current result bitvec
                if (r.bitIsSet(path)) {
                    shared.append(UInt(path))
                }
            }
            if (shared.count > 1) {
                found = true
            } else if (shared.count == 1) {
                let classification = classifyTetrahedralNbrSymClasses(symClasses, mol.getAtom(Int(shared[0]))!)
                if (classification == .T1122 || classification == .T1111) {
                    found = true
                }
            }
            if (found) {
                // add bits for the atoms in the ring
                for path in ring._path {
                    r.setBitOn(UInt32(path))
                }
                break
            }
        }
        // add the ring as a new bitvec if it shares no atom with a previous ring
        if (!found) {
            let r = MKBitVec()
            for path in ring._path {
                r.setBitOn(UInt32(path))
            }
            result.append(r)
        }
    }
    return result
}


// TODO: I have c++ overloading...can we shrink these into two separate functions? 

/**
* Helper function for getFragment below.
*/
func addNbrs(_ fragment: inout MKBitVec, _ atom: MKAtom, _ skip: MKAtom) {
    guard let nbors = atom.getNbrAtomIterator() else {
        // TODO: throw error
        return
    }
    for nbr in nbors {
        // don't pass through skip
        if (nbr.getId() == skip.getId()) {
            continue
        }
        // skip visited atoms
        if (fragment.bitIsSet(nbr.getId().intValue)) {
            continue
        }
        // add the neighbor atom to the fragment
        fragment.setBitOn(UInt32(nbr.getId().intValue))
        // recurse...
        addNbrs(&fragment, nbr, skip)
    }
}

/**
* Create an OBBitVec objects with bets set for the fragment consisting of all
* atoms for which there is a path to atom without going through skip. These
* fragment bitvecs are indexed by unique id (i.e. OBAtom::GetId()).
*/
func getFragment(_ atom: MKAtom, _ skip: MKAtom) -> MKBitVec {
    var fragment: MKBitVec = MKBitVec()
    fragment.setBitOn(UInt32(atom.getId().intValue))
    // start the recursion
    addNbrs(&fragment, atom, skip)
    return fragment
}

/**
 * Helper function for getFragment below.
 */
func addNbrs(_ fragment: inout MKBitVec, _ atom: MKAtom, _ mask: inout MKBitVec) {
    guard let nbors = atom.getNbrAtomIterator() else {
        print("ERROR: no neighbors")
        return
    }
    for nbr in nbors {
        if (!mask.bitIsSet(nbr.getId().intValue)) {
            continue
        }
        // skip visited atoms
        if (fragment.bitIsSet(nbr.getId().intValue)) {
            continue
        }
        // add the neighbor atom to the fragment
        fragment.setBitOn(UInt32(nbr.getId().intValue))
        // recurse...
        addNbrs(&fragment, nbr, &mask)
    }
}

/**
 * Create an OBBitVec objects with bets set for the fragment consisting of all
 * atoms for which there is a path to atom without going through skip. These
 * fragment bitvecs are indexed by atom idx (i.e. OBAtom::GetIdx()).
 */
func getFragment(_ atom: MKAtom, _ mask: inout MKBitVec) -> MKBitVec
{
    var fragment: MKBitVec = MKBitVec()
    fragment.setBitOn(UInt32(atom.getIdx()))
    // start the recursion
    addNbrs(&fragment, atom, &mask)
    return fragment
}


private protocol StereoRingType {
    var outsideNbrs: [MKAtom] { get set }
    var insideNbrs: [MKAtom] { get set }
    var outIdx: Ref { get set }
    var inIdx: Ref { get set }
    
    func isInRing(_ ring: StereoRing) -> Bool
}
extension StereoRingType {
     
}

struct StereoRing {

    struct ParaAtom: StereoRingType {
        
        // typealias of CenterType was declared here but not mentioned elsewhere in the code
        var id: Ref
        // union {
        //     unsigned int inIdx, outIdx;
        // };
        // replaced with one variable that is generic to each
        var Idx: Ref
        var inIdx: Ref {
            get {
                return Idx
            }
            set {
                Idx = newValue
            }
        }
        var outIdx: Ref {
            get {
                return Idx
            }
            set {
                Idx = newValue
            }
        }
        
        var insideNbrs: [MKAtom] = []
        var outsideNbrs: [MKAtom] = []

        init(id: Ref, idx: Ref) {
            self.id = id
            self.Idx = idx
        }

        func getCenter(_ mol: MKMol) -> MKAtom? {
            return mol.getAtomById(self.id)
        }
        
        func isInRing(_ ring: StereoRing) -> Bool {
            for atom in ring.paraAtoms {
                if (atom.Idx == self.Idx) {
                    return true
                }
            }
            return false
        }
    }

    struct ParaBond: StereoRingType {
        var id: Ref
        var inIdx: Ref
        var outIdx: Ref
        var insideNbrs: [MKAtom] = []
        var outsideNbrs: [MKAtom] = []

        init(id: Ref, inIdx: Ref, outIdx: Ref) {
            self.id = id
            self.inIdx = inIdx
            self.outIdx = outIdx
        }

        func getCenter(_ mol: MKMol) -> MKBond? {
            return mol.getBondById(self.id)
        }

        func isInRing(_ ring: StereoRing) -> Bool {
            for bond in ring.paraBonds {
                if (bond.inIdx == self.inIdx) {
                    return true
                }
            }
            return false
        }
    }

    var paraAtoms: [ParaAtom] = []
    var paraBonds: [ParaBond] = []
    var trueCount: UInt = 0

    init() {
        self.trueCount = 0
    }
}

private func checkLigands<T: StereoRingType>(_ currentPara: T, _ units: MKStereoUnitSet) -> Bool {
    if currentPara.outsideNbrs.count == 1 {
        print("OK: ligands")
        return true
    }
    guard let mol = currentPara.insideNbrs[0].getParent() else {
        print("ERROR: no parent")
        return false
    }
    guard let outAtom = mol.getAtomById(currentPara.outIdx) else {
        print("ERROR: no atom")
        return false
    }
    let ligand = getFragment(currentPara.outsideNbrs[0], outAtom)
    for unit in units {
        if isUnitInFragment(mol, unit, ligand) {
            print("OK: ligands")
            return true
        }
    }
    print("NOT OK: ligands")
    return false
}

// bool ApplyRule1(const Type &currentPara, const std::vector<unsigned int> &symmetry_classes,
//       const std::vector<StereoRing> &rings, std::vector<bool> &visitedRings, const OBStereoUnitSet &units,
//       std::vector<unsigned int> stereoAtoms)
private func applyRule1<T: StereoRingType>(_ currentPara: T, _ symmetry_classes: inout [UInt], _ rings: inout [StereoRing],
                                           _ visitedRings: inout [Bool], _ units: MKStereoUnitSet, _ stereoAtoms: inout [UInt]) -> Bool {
    var foundRing = false
    let idx = currentPara.inIdx
    for i in 0..<rings.count {
        // skip visited rings
        if visitedRings[i] {
            continue
        }
        // Check if currentPara is in this ring
        if !currentPara.isInRing(rings[i]) {
            continue
        }
        //
        // A new ring containing currentPara is found
        //
        foundRing = true
        // if there are one or more true stereo centers, currentPara is a stereo center
        if rings[i].trueCount > 0 {
            return true
        }
        // check if there is at least one other potential atom
        for j in 0..<rings[i].paraAtoms.count {
            let paraAtom = rings[i].paraAtoms[j]
            // skip idx
            if paraAtom.inIdx == idx {
                continue
            }
            // there is another atom already identified as stereo atom
            if stereoAtoms.contains(UInt(paraAtom.inIdx.intValue)) {
                return true
            }
            if paraAtom.outsideNbrs.count == 1 {
                // only 1 ring substituent, the other is implicit H -> topologically different
                return true
            } else {
                if paraAtom.outsideNbrs.count != 2 {
                    return false
                }
                // two ring substituents, need to check for topological difference
                if symmetry_classes[paraAtom.outsideNbrs[0].getIndex()] != symmetry_classes[paraAtom.outsideNbrs[1].getIndex()] {
                    // they are different
                    return true
                } else {
                    // they are the same and they might also be in a ring -> apply rule 1 recursive
                    visitedRings[i] = true
                    if applyRule1(paraAtom, &symmetry_classes, &rings, &visitedRings, units, &stereoAtoms) {
                        return true
                    }
                }
            }
        }
        // check if there is at least one other potential bond
        for j in 0..<rings[i].paraBonds.count {
            let paraBond = rings[i].paraBonds[j]
            // skip idx
            if paraBond.inIdx == idx {
                continue
            }
            // there is another atom already identified as stereo atom
            if stereoAtoms.contains(UInt(paraBond.inIdx.intValue)) {
                return true
            }
            if paraBond.outsideNbrs.count == 1 {
                // only 1 ring substituent, the other is implicit H -> topologically different
                return true
            } else {
                if paraBond.outsideNbrs.count != 2 {
                    continue
                }
                // two ring substituents, need to check for topological difference
                if symmetry_classes[paraBond.outsideNbrs[0].getIndex()] != symmetry_classes[paraBond.outsideNbrs[1].getIndex()] {
                    // they are different
                    return true
                } else {
                    // they are the same and they might also be in a ring -> apply rule 1 recursive
                    visitedRings[i] = true
                    if applyRule1(paraBond, &symmetry_classes, &rings, &visitedRings, units, &stereoAtoms) {
                        return true
                    }
                }
            }
        }
    }
    // if a non-visited ring was found and true was not returned -> it does not
    // contain any stereocenters other than idx
    if foundRing {
        print("NOT OK: no ring found")
        return false
    }
    print("NOT OK: no ring found")
    return false
}


// void StartRule1(const std::vector<unsigned int> &symmetry_classes, const std::vector<StereoRing> &rings,
//       OBStereoUnitSet &units, std::vector<unsigned int> &stereoAtoms)
func startRule1(_ symmetry_classes: inout [UInt], _ rings: inout [StereoRing], _ units: inout MKStereoUnitSet, _ stereoAtoms: inout [UInt]) {
    for i in 0..<rings.count {
        print("Checking ring: \(i)")
        // tetrahedral atoms
        for j in 0..<rings[i].paraAtoms.count {
            let paraAtom = rings[i].paraAtoms[j]
            // skip the atom if it is already in stereoAtoms
            if stereoAtoms.contains(UInt(paraAtom.inIdx.intValue)) {
                continue
            }
            var visitedRings = [Bool](repeating: false, count: rings.count)
            //visitedRings[i] = true;
            if applyRule1(paraAtom, &symmetry_classes, &rings, &visitedRings, units, &stereoAtoms) {
                var isStereoUnit = false
                if paraAtom.outsideNbrs.count == 1 {
                    isStereoUnit = true
                }
                if paraAtom.outsideNbrs.count == 2 {
                    if symmetry_classes[paraAtom.outsideNbrs[0].getIndex()] == symmetry_classes[paraAtom.outsideNbrs[1].getIndex()] {
                        // check for spiro atom
                        var isSpiro = false
                        for k in 0..<rings[i].paraAtoms.count {
                            let paraAtom2 = rings[i].paraAtoms[k]
                            if paraAtom.inIdx == paraAtom2.outIdx && paraAtom.insideNbrs == paraAtom2.outsideNbrs {
                                isSpiro = true
                                if applyRule1(paraAtom2, &symmetry_classes, &rings, &visitedRings, units, &stereoAtoms) {
                                    isStereoUnit = true
                                }
                            }
                        }
                        if !isSpiro {
                            isStereoUnit = checkLigands(paraAtom, units)
                        }
                        //cout << "isStereoUnit = " << isStereoUnit << endl;
                    } else {
                        isStereoUnit = true
                    }
                }
                if isStereoUnit {
                    stereoAtoms.append(UInt(paraAtom.inIdx.intValue))
                    let atom = paraAtom.insideNbrs[0].getParent()!.getAtomById(paraAtom.id)
                    units.append(MKStereoUnit(.Tetrahedral, atom!.getId(), true))
                }
            }
        }
        // cistrans bonds
        for j in 0..<rings[i].paraBonds.count {
            let paraBond = rings[i].paraBonds[j]
            // skip the atom if it is already in stereoAtoms
            if stereoAtoms.contains(UInt(paraBond.inIdx.intValue)) {
                continue
            }
            var visitedRings = [Bool](repeating: false, count: rings.count)
            //visitedRings[i] = true;
            if applyRule1(paraBond, &symmetry_classes, &rings, &visitedRings, units, &stereoAtoms) {
                var isStereoUnit = false
                if paraBond.outsideNbrs.count == 1 {
                    isStereoUnit = true
                }
                if paraBond.outsideNbrs.count == 2 {
                    if symmetry_classes[paraBond.outsideNbrs[0].getIndex()] == symmetry_classes[paraBond.outsideNbrs[1].getIndex()] {
                        // check for spiro bond
                        var isSpiro = false
                        for k in 0..<rings[i].paraBonds.count {
                            let paraBond2 = rings[i].paraBonds[k]
                            if paraBond.inIdx == paraBond2.outIdx && paraBond.insideNbrs == paraBond2.outsideNbrs {
                                isSpiro = true
                                if applyRule1(paraBond2, &symmetry_classes, &rings, &visitedRings, units, &stereoAtoms) {
                                    isStereoUnit = true
                                }
                            }
                        }
                        if !isSpiro {
                            isStereoUnit = checkLigands(paraBond, units)
                        }
                        //cout << "isStereoUnit = " << isStereoUnit << endl;
                    } else {
                        isStereoUnit = true
                    }
                }
                if isStereoUnit {
                    stereoAtoms.append(UInt(paraBond.inIdx.intValue))
                    stereoAtoms.append(UInt(paraBond.outIdx.intValue))
                    let bond = paraBond.insideNbrs[0].getParent()!.getBondById(paraBond.id)
                    units.append(MKStereoUnit(.CisTrans, bond!.getId(), true))
                }
            }
        }
    }
}

///@name Stereogenic unit identification
///@{
/**
* Find the stereogenic units in a molecule using a set of rules.<sup>1</sup>
*
* The potential stereocenters are identified first. A potential tetrahedral
* stereogenic atom is any atom meeting the following criteria:
*
* - sp3 hybridization
* - at least 3 "heavy" neighbors
*
* Nitrogen is treated as a special case since the barrier of inversion is
* low in many cases making the atom non-stereogenic. Only bridge-head
* nitrogen atoms (i.e. nitrogen has 3 neighbors in rings) will be
* considered stereogenic.
*
* Potential stereogenic double bonds are identified using another set of
* simple criteria:
*
* - must be a double bond
* - must not be in a ring
* - both begin and end atom should have at least one single bond
*
* True stereocenters (i.e. stereocenters with topologically different
* ligands) are identified first. For tetrahedral stereocenters, true
* stereocenters will have 4 different neighbor atom symmetry classes
* and this can be expressed using T1234 to classify these stereocenters.
* For stereogenic bonds, a similar classification C12 can be used but
* both begin and end atom have their own classification and the bond
* is only a true stereocenter if both atoms are C12.
*
* Para stereocenters are all stereocenters where there are at least two
* equivalent neighbor atom symmetry classes. These are T1123, T1112, T1111
* and T1122 for tetrahedral stereocenters and C11 for double bonds. To
* determine which of the remaining potential stereocenters really are
* stereocenters, a set of rules is used.<sup>1</sup>
*
* Rule 1 is applied recusively:
*
* All rings are merged "mergedRings". A merged ring is simply a fragment consisting
* of all atoms of a ring system (bridged, spiro, adjacent, ...). If two rings in the
* SSSR set share an atom, they are merged.
*
* Each merged must at least have two para-stereocenters (or 1 true + 1 para) in order
* for the para-stereocenter to be valid. This is repeated until no new stereocenters
* are identified.
*
* rule 1a for double bonds:
* - bond atom in ring has two identical symmetry classes for it's neighbor atoms (-> para)
* - other bond atom:
*   - has two different symmetry classes for it's neighbours -> new stereocenter
*   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
*
* rule 1b for tetracoord atoms:
* - at least two neighbour symmetry classes are the same (-> para)
* - other pair:
*   - has two different symmetry classes for it's neighbours -> new stereocenter
*   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
*
* Rules 2 and 3 are applied sequential (i.e. only once).
*
* Rule 2a for tetracoordinate carbon:
* - 1 or 2 pair identical ligands
* - each ligand contains at least 1 true-stereocenter or 2 para-stereocenters (from rule 1)
*
* Rule 2b for tetracoordinate carbon:
* - 3 or 4 identical ligands with at least
*   - 2 true-stereocenters
*   - 2 separate assemblies of para-stereocenters (from rule 1)
*
* Rule 3 for double bonds:
* - 1 or 2 pair identical ligands (on begin and end atom)
* - each pair contains at least 1 true-stereocenter or 2 para-stereocenters (from rule 1)
*
*
* @verbatim
    Reference:
    [1] M. Razinger, K. Balasubramanian, M. Perdih, M. E. Munk, Stereoisomer
    Generation in Computer-Enhanced Structure Elucidation, J. Chem. Inf.
    Comput. Sci. 1993, 33, 812-825
    @endverbatim
*/
func findStereogenicUnits(_ mol: MKMol, symClasses: inout [UInt]) -> MKStereoUnitSet {
    
    var units: MKStereoUnitSet = MKStereoUnitSet()

    // do quick test to see if there are any possible stereogenic units
    if !mayHaveTetrahedralCenter(mol) && !mayHaveCisTransBond(mol) {
        return units
    }

    // make sure we have symmetry classes for all atoms
    if symClasses.count != mol.numAtoms() {
        return units
    }

    // para-stereocenters candidates
    var stereoAtoms: [UInt] = [UInt]() // Tetrahedral = idx, CisTrans = begin & end idx
    var paraAtoms: [UInt] = [UInt]()
    var paraBonds: [UInt] = [UInt]()

    /**
     * true Tetrahedral stereocenters:
     * - have four different symmetry classes for the ligands to the central atom
     */
    var isChiral: Bool = false
    for atom in mol.getAtomIterator() {
        if !isPotentialTetrahedral(atom) {
            continue
        }

        // list containing neighbor symmetry classes
        var tlist: [UInt] = [UInt]()
        isChiral = true

        // check neighbors to see if this atom is stereogenic
        guard let nbors = atom.getNbrAtomIterator() else {
//           TODO: throw error gracefully and return
            return units
        }
        for nbr in nbors {
            // check if we already have a neighbor with this symmetry class
            for k in 0..<tlist.count {
                if symClasses[nbr.getIndex()] == tlist[k] {
                    isChiral = false
                    // if so, might still be a para-stereocenter
                    paraAtoms.append(UInt(atom.getIdx()))
                }
            }

            if isChiral {
                // keep track of all neighbors, so we can detect duplicates
                tlist.append(symClasses[nbr.getIndex()])
            } else {
                break
            }
        }

        if isChiral {
            // true-stereocenter found
            stereoAtoms.append(UInt(atom.getIndex()))
            units.append(MKStereoUnit(.Tetrahedral, atom.getId()))
        }
    }

    /**
     * true CisTrans stereocenters:
     * - each terminal has two different symmetry classes for it's ligands
     */
    var isCisTransBond: Bool = false
    for bond in mol.getBondIterator() {
        if bond.isInRing() && bond.isAromatic() {
            continue
        }
        if bond.getBondOrder() == 2 {
            let begin = bond.getBeginAtom()
            let end = bond.getEndAtom()
            if begin.getTotalDegree() > 3 || end.getTotalDegree() > 3 {
                continue // e.g. C=Ru where the Ru has four substituents
            } 
            // Needs to have at least one explicit single bond at either end
            // FIXME: timvdm: what about C=C=C=C
            if !begin.hasSingleBond() || !end.hasSingleBond() {
                continue
            }
            isCisTransBond = true

            if begin.getExplicitDegree() == 2 {
                // Begin atom has two explicit neighbors. One is the end atom. The other should
                // be a heavy atom - this is what we test here.
                // (There is a third, implicit, neighbor which is either a hydrogen
                // or a lone pair.)
                if begin.explicitHydrogenCount() == 1 {
                    isCisTransBond = false
                }
            } else if begin.getExplicitDegree() == 3 {
                var tlist: [UInt] = [UInt]()
                for nbr in begin.getNbrAtomIterator()! {
                    // skip end atom
                    if nbr.getId() == end.getId() {
                        continue
                    }
                    // do we already have an atom with this symmetry class?
                    if tlist.count > 0 {
                        // compare second with first
                        if symClasses[nbr.getIndex()] == tlist[0] {
                            isCisTransBond = false
                            // if same, might still be a para-stereocenter
                            paraBonds.append(UInt(bond.getIdx()))
                        }
                        break
                    }
                    // save first symmetry class
                    tlist.append(symClasses[nbr.getIndex()])
                }
            } else {
                // Valence is not 2 or 3, for example SR3=NR
                isCisTransBond = false
            }
            if !isCisTransBond {
                continue
            }

            if end.getExplicitDegree() == 2 {
                // see comment above for begin atom
                if end.explicitHydrogenCount() == 1 {
                    isCisTransBond = false
                }
            } else if end.getExplicitDegree() == 3 {
                var tlist: [UInt] = [UInt]()
                for nbr in end.getNbrAtomIterator()! {
                    // skip end atom
                    if nbr.getId() == begin.getId() {
                        continue
                    }
                    // do we already have an atom with this symmetry class?
                    if tlist.count > 0 {
                        // compare second with first
                        if symClasses[nbr.getIndex()] == tlist[0] {
                            // if same, might still be a para-stereocenter
                            paraBonds.append(UInt(bond.getIdx()))
                            isCisTransBond = false
                        }
                        break
                    }
                    // save first symmetry class
                    tlist.append(symClasses[nbr.getIndex()])
                }
            } else {
                // Valence is not 2 or 3, for example SR3=NR
                isCisTransBond = false
            }
            if isCisTransBond {
            //   true-stereocenter found
                units.append(MKStereoUnit(.CisTrans, bond.getId()))
            }
        }
    }
    
    /**
     * Apply rule 1 from the Razinger paper recusively:
     *
     * All rings are merged "mergedRings". A merged ring is simply a fragment consisting
     * of all atoms of a ring system (bridged, spiro, adjacent, ...). If two rings in the
     * SSSR set share an atom, they are merged.
     *
     * Each merged must at least have two para-stereocenters (or 1 true + 1 para) in order
     * for the para-stereocenter to be valid. This is repeated until no new stereocenters
     * are identified.
     *
     * rule 1a for double bonds:
     * - bond atom in ring has two identical symmetry classes for it's neighbor atoms (-> para)
     * - other bond atom:
     *   - has two different symmetry classes for it's neighbours -> new stereocenter
     *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
     *
     * rule 1b for tetracoord atoms:
     * - at least two neighbour symmetry classes are the same (-> para)
     * - other pair:
     *   - has two different symmetry classes for it's neighbours -> new stereocenter
     *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
     *
     * NOTE: there must always be at least 2 new stereocenters (or one existing + 1 newly found) in order for them to be valid
     */
    let lssr = mol.getLSSR()
    var rings = [StereoRing]()
    for i in 0..<lssr.count {
        rings.append(StereoRing())
        var ring = rings.last!
        for j in 0..<stereoAtoms.count {
            if lssr[i]._pathSet.bitIsSet(Int(stereoAtoms[j])) {
                ring.trueCount += 1
            }
        }
        for j in 0..<paraAtoms.count {
            if lssr[i]._pathSet.bitIsSet(Int(paraAtoms[j])) {
                guard let atom = mol.getAtom(Int(paraAtoms[j])) else {
                    fatalError("ERROR: atom is not found!!")
                }
                ring.paraAtoms.append(StereoRing.ParaAtom(id: atom.getId(), idx: .Ref(Int(paraAtoms[j]))))
                for nbr in atom.getNbrAtomIterator()! {
                    if lssr[i]._pathSet.bitIsSet(nbr.getIndex()) {
                        guard var lastAtom = ring.paraAtoms.last else { break }
                        lastAtom.insideNbrs.append(nbr)
                    } else {
                        guard var lastAtom = ring.paraAtoms.last else { break }
                        lastAtom.outsideNbrs.append(nbr)
                    }
                }
                if ring.paraAtoms.last!.insideNbrs.count != 2 {
                    ring.paraAtoms.removeLast()
                }
            }
        }
        for j in 0..<paraBonds.count {
            guard let bond = mol.getBond(Int(paraBonds[j])) else {
                fatalError("ERROR: Bond is not found")
            }
            let beginIdx = bond.getBeginAtomIdx()
            let endIdx = bond.getEndAtomIdx()
            if lssr[i]._pathSet.bitIsSet(Int(beginIdx)) {
                ring.paraBonds.append(StereoRing.ParaBond(id: bond.getId(), inIdx: .Ref(beginIdx), outIdx: .Ref(endIdx)))
                for nbr in bond.getBeginAtom().getNbrAtomIterator()! {
                    if nbr.getIndex() == endIdx {
                        continue
                    }
                    guard var lastBond = ring.paraBonds.last else { break }
                    lastBond.insideNbrs.append(nbr)
                }
                for nbr in bond.getEndAtom().getNbrAtomIterator()! {
                    if nbr.getIndex() == beginIdx {
                        continue
                    }
                    guard var lastBond = ring.paraBonds.last else { break }
                    lastBond.outsideNbrs.append(nbr)
                }
                if ring.paraBonds.last!.insideNbrs.count != 2 {
                    ring.paraBonds.removeLast()
                }
            }
            if lssr[i]._pathSet.bitIsSet(Int(endIdx)) {
                ring.paraBonds.append(StereoRing.ParaBond(id: bond.getId(), inIdx: .Ref(endIdx), outIdx: .Ref(beginIdx)))
                for nbr in bond.getEndAtom().getNbrAtomIterator()! {
                    if nbr.getIndex() == beginIdx {
                        continue
                    }
                    guard var lastBond = ring.paraBonds.last else { break }
                    lastBond.insideNbrs.append(nbr)
                }
                for nbr in bond.getBeginAtom().getNbrAtomIterator()! {
                    if nbr.getIndex() == endIdx {
                        continue
                    }
                    guard var lastBond = ring.paraBonds.last else { break }
                    lastBond.outsideNbrs.append(nbr)
                }
                if ring.paraBonds.last!.insideNbrs.count != 2 {
                    ring.paraBonds.removeLast()
                }
            }
            if ring.paraAtoms.count + ring.paraBonds.count == 1 {
                ring.paraAtoms.removeAll()
                ring.paraBonds.removeAll()
            }
        }

        


    }
        
    var numStereoUnits: Int

    repeat {
        numStereoUnits = units.count
        startRule1(&symClasses, &rings, &units, &stereoAtoms)
    } while units.count > numStereoUnits
    
    
    let mergedRings = mergeRings(mol, symClasses)
    
    /**
     * Apply rule 2a for tetracoordinate carbon:
     * - 1 or 2 pair identical ligands
     * - each pair contains at least 1 true-stereocenter or 2 para-stereocenters
     *
     * Apply rule 2b for tetracoordinate carbon:
     * - 3 or 4 identical ligands with at least
     *   - 2 true-stereocenters
     *   - 2 separate assemblies of para-stereocenters
     */
    
    for idx in paraAtoms {
        guard let atom = mol.getAtom(Int(idx)) else {
            fatalError("ERROR: atom is not discovered in paraAtoms")
        }
        // make sure we didn't add this atom already from rule 1
        var alreadyAdded: Bool = false 
        for u2 in units {
            if u2.type == .Tetrahedral {
                if atom.getId() == u2.id {
                    alreadyAdded = true
                    break
                }
            }
        }
        if alreadyAdded {
            continue
        }

        let classification = classifyTetrahedralNbrSymClasses(symClasses, atom)
        switch classification {
            case .T1123:
            // rule 2a with 1 pair
            let duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses: symClasses)
            guard let ligandAtom = findAtomWithSymmetryClass(atom, symClass: duplicatedSymClass, symClasses: symClasses) else { break }
            if containsAtLeast_1true_2para(ligandAtom, atom: atom, units: units) {
                units.append(MKStereoUnit(.Tetrahedral, atom.getId(), true))
            }
        case .T1122:
            // rule 2a with 2 pairs
            var duplicatedSymClass1: UInt = 0
            var duplicatedSymClass2: UInt = 0
            findDuplicatedSymmetryClasses(atom, symClasses: symClasses, duplicated1: &duplicatedSymClass1, duplicated2: &duplicatedSymClass2)
            guard let ligandAtom1 = findAtomWithSymmetryClass(atom, symClass: duplicatedSymClass1, symClasses: symClasses) else { break }
            guard let ligandAtom2 = findAtomWithSymmetryClass(atom, symClass: duplicatedSymClass2, symClasses: symClasses) else { break }
            if (containsAtLeast_1true_2para(ligandAtom1, atom: atom, units: units) &&
                containsAtLeast_1true_2para(ligandAtom2, atom: atom, units: units)) {
                units.append(MKStereoUnit(.Tetrahedral, atom.getId(), true))
            }
        case .T1112:
            // rule 2b with 3 identical
            let duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses: symClasses)
            guard let ligandAtom = findAtomWithSymmetryClass(atom, symClass: duplicatedSymClass, symClasses: symClasses) else { break }
            if containsAtLeast_2true_2paraAssemblies(ligandAtom, atom: atom, units: units, mergedRings: mergedRings) {
                units.append(MKStereoUnit(.Tetrahedral, atom.getId(), true))
            }
        case .T1111:
            // rule 2b with 4 identical
            let duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses: symClasses)
            guard let ligandAtom = findAtomWithSymmetryClass(atom, symClass: duplicatedSymClass, symClasses: symClasses) else { break }
            if containsAtLeast_2true_2paraAssemblies(ligandAtom, atom: atom, units: units, mergedRings: mergedRings) {
                units.append(MKStereoUnit(.Tetrahedral, atom.getId(), true))
            }
        default: break
        }
    }
    
    /**
     * Apply rule 3 for double bonds.
     * - 1 or 2 pair identical ligands (on begin and end atom)
     * - each pair contains at least 1 true-stereocenter or 2 para-stereocenters (from rule1)
     */
    for idx in paraBonds {
        guard let bond = mol.getBond(Int(idx)) else { continue }
        var alreadyAdded: Bool = false
        for u2 in units {
            if u2.type == .CisTrans {
                if bond.getId() == u2.id {
                    alreadyAdded = true
                    break
                }
            }
        }
        
        if alreadyAdded {
            continue
        }
        
        let begin = bond.getBeginAtom()
        let end = bond.getEndAtom()
        //         begin classification
        let beginClassification = classifyCisTransNbrSymClasses(symClasses, bond, begin)
        var beginValid = false
        switch beginClassification {
        case .C12:
            beginValid = true
        case .C11:
            // find the ligand 
            var ligandAtom: MKAtom? = nil 
            for nbr in begin.getNbrAtomIterator()! {
                if (nbr.getIdx() != bond.getBeginAtomIdx()) && (nbr.getIdx() != bond.getEndAtomIdx()) {
                    ligandAtom = nbr
                    break
                }
            }
            guard ligandAtom != nil else {
                print("ERROR: ligand atom is nil??")
                break 
            } // TODO: Throw Error Here 
            let ligand: MKBitVec = getFragment(ligandAtom!, begin)
            for u2 in units {
                if u2.type == .Tetrahedral {
                    if ligand.bitIsSet(u2.id.intValue) {
                        beginValid = true
                        break
                    }
                } else if u2.type == .CisTrans {
                    guard let bond2 = mol.getBondById(u2.id) else { continue }
                    let begin2 = bond2.getBeginAtom()
                    let end2 = bond2.getEndAtom()
                    if ligand.bitIsSet(begin2.getId()) || ligand.bitIsSet(end2.getId()) {
                        beginValid = true
                        break
                    }
                }
            }
        default: break 
        }
        
        if !beginValid {
            continue
        }
        //         end classification
        let endClassification = classifyCisTransNbrSymClasses(symClasses, bond, end)
        var endValid = false
        switch endClassification {
        case .C12:
            endValid = true
        case .C11:
            // find the ligand
            var ligandAtom: MKAtom? = nil
            for nbr in end.getNbrAtomIterator()! {
                if (nbr.getIdx() != bond.getBeginAtomIdx()) && (nbr.getIdx() != bond.getEndAtomIdx()) {
                    ligandAtom = nbr
                    break
                }
            }
            guard ligandAtom != nil else {
                print("ERROR: ligand atom is nil??")
                break
            } // TODO: Throw Error Here
            let ligand: MKBitVec = getFragment(ligandAtom!, end)
            for u2 in units {
                if u2.type == .Tetrahedral {
                    if ligand.bitIsSet(u2.id.intValue) {
                        endValid = true
                        break
                    }
                } else if u2.type == .CisTrans {
                    guard let bond2 = mol.getBondById(u2.id) else { continue }
                    let begin2 = bond2.getBeginAtom()
                    let end2 = bond2.getEndAtom()
                    if ligand.bitIsSet(begin2.getId()) || ligand.bitIsSet(end2.getId()) {
                        endValid = true
                        break
                    }
                }
            }
        default: break
        }
        
        if endValid {
            units.append(MKStereoUnit(.CisTrans, bond.getId(), true))
        }
    }
    
//    MARK: FOR DEBUG ONLY
    
    for unit in units {
        if unit.type == .Tetrahedral {
            print("Tetrahedral(center = \(unit.id), para = \(unit.para))")
        }
        if unit.type == .CisTrans {
            print("CisTrans(center = \(unit.id), para = \(unit.para))")
        }
        if unit.type == .SquarePlanar {
            print("SquarePlanar(center = \(unit.id), para = \(unit.para))")
        }
    }

    return units
}


/**
* Helper function for FindStereogenicUnits using automorphisms.
*
* Find the duplicated symmetry class for neighbors of atom. This method only works if there is
* only one duplicated symmetry class (i.e. T1123, T1112, T1111).
*/
func findDuplicatedSymmetryClass(_ atom: MKAtom, symClasses: [UInt]) -> UInt {
    // find the duplicated symmetry class
    var duplicatedSymClass = MKGraphSym.NoSymmetryClass
    var nbrSymClasses = [UInt]()
    for nbr in atom.getNbrAtomIterator()! {
        nbrSymClasses.append(symClasses[nbr.getIndex()])
    }
    for i in 0..<nbrSymClasses.count {
        if nbrSymClasses.filter({$0 == nbrSymClasses[i]}).count >= 2 {
            duplicatedSymClass = nbrSymClasses[i]
            break
        }
    }
    return duplicatedSymClass
}

/**
* Helper function for FindStereogenicUnits using automorphisms.
*
* Find the duplicated symmetry classes for neighbors of atom. This method only works for the
* T1122 case.
*/
func findDuplicatedSymmetryClasses(_ atom: MKAtom, symClasses: [UInt], duplicated1: inout UInt, duplicated2: inout UInt) {
    var nbrSymClasses = [UInt]()
    for nbr in atom.getNbrAtomIterator()! {
        nbrSymClasses.append(symClasses[nbr.getIndex()])
    }
    nbrSymClasses.sort()
    duplicated1 = nbrSymClasses[0]
    duplicated2 = nbrSymClasses[2]
}

/**
* Helper function for FindStereogenicUnits using automorphisms.
*
* Find the duplicated symmetry classes for neighbors of atoms. This method works for all
* cases (i.e. T1234, T1123, T1112, T1111 and T1122).
*/
func findDuplicatedSymmetryClasses(_ atom: MKAtom, symClasses: [UInt]) -> [UInt] {
    var nbrSymClasses = [UInt]()
    for nbr in atom.getNbrAtomIterator()! {
        nbrSymClasses.append(symClasses[nbr.getIndex()])
    }
    nbrSymClasses.sort()
    var result = [UInt]()
    for i in 0..<nbrSymClasses.count {
        if nbrSymClasses.filter({$0 == nbrSymClasses[i]}).count > 1 {
            if !result.contains(nbrSymClasses[i]) {
                result.append(nbrSymClasses[i])
            }
        }
    }
    return result
}

/**
* Helper functions for FindStereogenicUnits (using automorphisms).
*
* These functions determine if an automorphism permutation invert the
* configuration of stereocenters by exchanging equivalent neighbor atoms
* (i.e. neighbor atoms with the same topological symmetry class).
*
* @note: The molecule should be ordered by topological canonical labels.
*/


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
//
//  FindStereogenicUnits using automorphisms
//
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

/**
* Find an atom with the specified symmetry class. The first found atom is returned
* or 0 when there is no such atom. This function is intended to be used in cases
* where any atom with the specified symmetry class can be used. For example, when
* checking a fragments for stereocenters, the result will be the same for any atom
* with a specified (duplicated) symmetry class.
*/
func findAtomWithSymmetryClass(_ atom: MKAtom, symClass: UInt, symClasses: [UInt]) -> MKAtom? {
    for nbr in atom.getNbrAtomIterator()! {
        if symClasses[nbr.getIndex()] == symClass {
            return nbr
        }
    }
    return nil
}

/**
* Helper function to determine if a stereogenic center with duplicated symmetry classes
* really is a stereogenic center.
*
* Check if the ligandAtom's fragment (see getFragment()) contains at least one
* true- or 1 para-stereocenter. This is rule 1 (a & b) in the Razinger paper on
* stereoisomer generation.
*/
func containsAtLeast_1true_1para(_ ligandAtom: MKAtom, skip: MKAtom, units: [MKStereoUnit]) -> Bool {
    guard let mol = skip.getParent() else { return false }
    // create the fragment bitvec
    let ligand = getFragment(ligandAtom, skip)
    for u2 in units {
        if isUnitInFragment(mol, u2, ligand) {
            return true
        }
    }
    return false
}

/**
* Helper function to determine if a stereogenic center with duplicated symmetry classes
* really is a stereogenic center.
*
* Check if the ligandAtom's fragment (see getFragment()) contains at least one
* true- or 2 para-stereocenter. This is rule 2a and rule 3 in the Razinger
* paper on stereoisomer generation.
*/
func containsAtLeast_1true_2para(_ ligandAtom: MKAtom, atom: MKAtom, units: [MKStereoUnit]) -> Bool {
    guard let mol = atom.getParent() else { return false }
    // check if ligand contains at least:
    // - 1 true-stereocenter
    // - 2 para-stereocenters
    let ligand = getFragment(ligandAtom, atom)
    var foundTrueStereoCenter = false
    var paraStereoCenterCount = 0
    for u2: MKStereoUnit in units {
        if isUnitInFragment(mol, u2, ligand) {
            if u2.para {
                paraStereoCenterCount += 1
            } else {
                foundTrueStereoCenter = true
            }
        }
    }
    if foundTrueStereoCenter || paraStereoCenterCount >= 2 {
        return true
    }
    if ligandAtom.isInRing() && atom.isInRing() && paraStereoCenterCount > 0 {
        return true
    }
    return false
}

/**
* Helper function to determine if a stereogenic center with duplicated symmetry classes
* really is a stereogenic center.
*
* Check if the ligandAtom's fragment (see getFragment()) contains at least one
* true- or 2 separate assemblies of at least 2 para-stereocenter. This is rule
* 2b in the Razinger paper on stereoisomer generation.
*/
func containsAtLeast_2true_2paraAssemblies(_ ligandAtom: MKAtom, atom: MKAtom, units: [MKStereoUnit], mergedRings: [MKBitVec]) -> Bool {
    guard let mol = atom.getParent() else { return false }
    // check if ligand contains at least:
    // - 2 true-stereocenter
    // - 2 separate para-stereocenters assemblies
    let ligand = getFragment(ligandAtom, atom)
    var trueStereoCenterCount = 0
    var ringIndices = [UInt]()
    for u2 in units {
        if u2.type == .Tetrahedral {
            if ligand.bitIsSet(u2.id.intValue) {
                if u2.para {
                    if let paraAtom: MKAtom = mol.getAtomById(u2.id) {
                        for ringIdx in 0..<mergedRings.count {
                            if mergedRings[ringIdx].bitIsSet(paraAtom.getIdx()) {
                                if !ringIndices.contains(UInt(ringIdx)) {
                                    ringIndices.append(UInt(ringIdx))
                                }
                            }
                        }
                    }
                } else {
                    trueStereoCenterCount += 1
                }
            }
        } else if u2.type == .CisTrans {
            guard let bond = mol.getBondById(u2.id) else { return false }
            let begin = bond.getBeginAtom() 
            let end = bond.getEndAtom()
            if ligand.bitIsSet(begin.getIdx()) || ligand.bitIsSet(end.getIdx()) {
                if u2.para {
                    for ringIdx in 0..<mergedRings.count {
                        if mergedRings[ringIdx].bitIsSet(begin.getIdx()) || mergedRings[ringIdx].bitIsSet(end.getIdx()) {
                            if !ringIndices.contains(UInt(ringIdx)) {
                                ringIndices.append(UInt(ringIdx))
                            }
                        }
                    }
                } else {
                    trueStereoCenterCount += 1
                }
            }
        }
    }
    if (trueStereoCenterCount >= 2 || ringIndices.count >= 2) {
        return true
    }
    return false
}



/**
 * Helper functions for FindStereogenicUnits (using automorphisms).
 *
 * These functions determine if an automorphism permutation invert the
 * configuration of stereocenters by exchanging equivalent neighbor atoms
 * (i.e. neighbor atoms with the same topological symmetry class).
 *
 * @note: The molecule should be ordered by topological canonical labels.
 */

struct StereoInverted {
    struct Entry {
        var p: Automorphism?
        var invertedAtoms: [MKAtom] = []
        var invertedBonds: [MKBond] = []
    }
    
    /**
     * Check if the specified automorphism causes an inversion of configuration
     * for the specified tetrahedral stereogenic center.
     */
    static func permutationInvertsTetrahedralCenter(_ p: Automorphism, _ center: MKAtom, _ symmetry_classes: [UInt], _ canon_labels: [UInt]) -> Bool {
        // Find the duplicated ligand symmetry class(es)
        let duplicatedSymClasses = findDuplicatedSymmetryClasses(center, symClasses: symmetry_classes)
        var duplicatedAtoms: [[MKAtom]] = []
        var permutated: Int = 0
        for i in 0..<duplicatedSymClasses.count {
            let duplicatedSymClass: UInt = duplicatedSymClasses[i]
            
            // Store the ligand indexes for the atoms with the duplicated symmetry class
            var tlist1: [Pair<UInt, UInt>] = []
            for nbr in center.getNbrAtomIterator()! {
                if symmetry_classes[nbr.getIndex()] == duplicatedSymClass {
                    tlist1.append((UInt(nbr.getIndex()), canon_labels[nbr.getIndex()]))
                    if duplicatedAtoms.count > 1 {
                        guard var last = duplicatedAtoms.last else { return false }
//                        TODO: throw error here
                        last.append(nbr)
                    } else {
                        duplicatedAtoms.append([nbr])
                    }
                }
            }
//             sort the indices
            tlist1.sort { $0.1 < $1.1 }
            // Translate the sorted indexes using the automorphism
            var tlist2: [UInt] = []
            for j in 0..<tlist1.count {
                var t: UInt = 0
                if mapsTo(p, tlist1[j].1, &t) {
                    tlist2.append(canon_labels[Int(t)])
                }
            }
            // Permute the flag
//            quickly map tlist2 to refs
            let tlist2Ref = tlist2.map { RefValue.Ref(Int($0)) }
            if ((MKStereo.numInversions(tlist2Ref) % 2) != 0) {
                //permutated = !permutated;
                permutated += 1
            }
        }
        if (permutated == 2) {
            let lssr = center.getParent()!.getLSSR()
            assert( duplicatedAtoms.count == 2 )
            assert( duplicatedAtoms[0].count == 2 )
            assert( duplicatedAtoms[1].count == 2 )
            for i in 0..<lssr.count {
                if lssr[i]._pathSet.bitIsSet(duplicatedAtoms[0][0].getIdx()) &&
                    lssr[i]._pathSet.bitIsSet(duplicatedAtoms[0][1].getIdx()) {
                    return false
                }
                if lssr[i]._pathSet.bitIsSet(duplicatedAtoms[1][0].getIdx()) &&
                    lssr[i]._pathSet.bitIsSet(duplicatedAtoms[1][1].getIdx()) {
                    return false
                }
            }
            return true
        }
        return false
    }
    
    static func permutationInvertsCisTransBeginOrEndAtom(_ p: Automorphism, _ bond: MKBond, _ beginOrEnd: MKAtom, _ canon_labels: [UInt]) -> Bool {
        let otherAtom = bond.getNbrAtom(beginOrEnd)
        var tlist1: [Pair<UInt, UInt>] = []
        for nbr in beginOrEnd.getNbrAtomIterator()! {
            if nbr != otherAtom {
                tlist1.append((UInt(nbr.getIndex()), canon_labels[nbr.getIndex()]))
            }
        }
        // sort the indices
        tlist1.sort { $0.1 < $1.1 }
        // Translate the sorted indexes using the automorphism
        var tlist2: [UInt] = []
        for j in 0..<tlist1.count {
            var t: UInt = 0
            if mapsTo(p, tlist1[j].1, &t) {
                tlist2.append(canon_labels[Int(t)])
            }
        }
        let tlist2Ref = tlist2.map { RefValue.Ref(Int($0)) }
        return ((MKStereo.numInversions(tlist2Ref) % 2) != 0)
    }

    /**
     * Check if the specified automorphism causes an inversion of configuration
     * for the specfied stereogenic double bond.
     */
     static func permutationInvertsCisTransCenter(_ p: Automorphism, _ bond: MKBond, _ canon_labels: [UInt]) -> Bool {
        // begin atom
        let beginInverted = permutationInvertsCisTransBeginOrEndAtom(p, bond, bond.getBeginAtom(), canon_labels)
        // end atom 
        let endInverted = permutationInvertsCisTransBeginOrEndAtom(p, bond, bond.getEndAtom(), canon_labels)
        // combine result using xor operation
        if (!beginInverted && endInverted) || (beginInverted && !endInverted) {
            return true
        }
        return false
     }
    
    /**
     * Perform the computation.
     */
    static func compute(_ mol: MKMol, _ symClasses: inout [UInt], _ automorphisms: Automorphisms) -> [Entry] {
        // We need topological canonical labels for this
        var canon_labels: [UInt] = []
        canonicalLabels(mol, &symClasses, &canon_labels, MKBitVec(), 5, true)
        // the result
        var result: [Entry] = []
        // make a list of stereogenic centers inverted by the automorphism permutations
        for i in 0..<automorphisms.count {
            var entry = Entry()
            entry.p = automorphisms[i]
            // Check the atoms
            for atom in mol.getAtomIterator() {
                // consider only potential stereo centers
                if !isPotentialTetrahedral(atom) {
                    continue
                }
                // add the atom to the inverted list if the automorphism inverses it's configuration
                if permutationInvertsTetrahedralCenter(automorphisms[i], atom, symClasses, canon_labels) {
                    entry.invertedAtoms.append(atom)
                }
            }
            // Check the bonds
            for bond in mol.getBondIterator() {
                // consider only potential stereo centers
                if !isPotentialCisTrans(bond) {
                    continue
                }
                // add the bond to the inverted list if the automorphism inverses it's configuration
                if permutationInvertsCisTransCenter(entry.p!, bond, canon_labels) {
                    entry.invertedBonds.append(bond)
                }
            }
            result.append(entry)
        }
        return result
    }
}


/**
* @brief Find the stereogenic units in a molecule making use of the automorphisms.
*
* The potential stereocenters are identified first. A potential tetrahedral
* stereogenic atom is any atom meeting the following criteria:
*
* - sp3 hybridization
* - at least 3 "heavy" neighbors
*
* Nitrogen is treated as a special case since the barrier of inversion is
* low in many cases making the atom non-stereogenic. Only bridge-head
* nitrogen atoms (i.e. nitrogen has 3 neighbors in rings) will be
* considered stereogenic.
*
* Potential stereogenic double bonds are identified using another set of
* simple criteria:
*
* - must be a double bond
* - must not be in a ring
* - both begin and end atom should have at least one single bond
*
* Once the potential stereocenters are found, the automorphisms are the key
* to identifying real stereogenic units. Automorphisms can be seen as
* permutations that permutate a graph back to the same graph. Such a
* permutation can only exchange atoms with the same symmetry class and it
* follows that the use of automorphisms takes symmetry into account. The
* definitions below use a concept where the automorphisms cause inversions
* of configuration to potential stereocenters. Such an inversion occurs
* whenever an automorphism exchanges two equivalent (i.e. with the same
* symmetry class) neighbor atoms attached to the potential stereogenic unit.
*
* @par Definition for tetrahedral stereocenters:
* A potential stereocenter really is a stereocenter if there exists no automorphic
* permutation causing an inversion of the configuration of only the potential
* stereogenic unit under consideration.
* If there exists at least one automorphic permutation causing an inversion of
* the configuration, then the potential stereogenic center can be a stereogenic
* center if the number of topologically equivalent neighbors (ligands) of the
* potential stereogenic center is less than or equal to the number of different
* configurations of these ligands.<sup>1</sup>
*
* The actual number of configurations needed for the ligands depends on the
* classification (i.e. T1234, T1123, ...) of the stereo center. These classes
* reflect the symmetry classes of the neighbor atoms of the center.
*
* - T1123: 1 true stereocenter OR 2 para stereocenters
* - T1122: 1 true stereocenter OR 2 para stereocenters (for both)
* - T1112: 2 true stereocenters OR 2 para stereocenter assemblies
* - T1111: 2 true stereocenters OR 2 para stereocenter assemblies
*
* @par Definition for double bond stereocenters:
* A potential stereogenic double bond really is a stereogenic bond if there
* exists no automorphic permutation causing an inversion of the configuration
* of only the potential stereogenic unit under consideration. The bond can still
* be a stereogenic bond if there exists such an automorphism when the number of
* configurations of the pair of topologically equivalent geminal ligands, which
* are exchanged by the automorphism, is greater than or equal to two (i.e. the
* number of topologically equivalent geminal ligands.<sup>1</sup>
*
* For stereogenic bonds, there is only one case but both begin and end atom
* have to be checked.
*
* - C11: 1 true stereocenter OR 1 para stereocenter
*
* These criteria are analogous to the rules from the Razinger paper on
* stereoisomer generation. Since the existence of stereocenters can depend
* on the existence of other stereocenters (in the ligands), the stereocenters
* are found by iterating until no new stereocenters are found.
*
* @verbatim
    Reference:
    [1] M. Perdih, M. Razinger, Stereochemistry and Sequence Rules:
    A Proposal for Modification of Cahn-Ingold-Prelog System,
    Tetrahedron: Asymmetry, 1994, Vol. 5, No. 5, 835-861
    @endverbatim
*/
func findStereogenicUnits(_ mol: MKMol, _ symClasses: inout [UInt], _ automorphisms: Automorphisms) -> MKStereoUnitSet {
    var units: MKStereoUnitSet = MKStereoUnitSet()
    // do quick test to see if there are any possible stereogenic units
    // do quick test to see if there are any possible stereogenic units
    if !mayHaveTetrahedralCenter(mol) && !mayHaveCisTransBond(mol) {
        return units
    }
    // make sure we have symmetry classes for all atoms
    if symClasses.count != mol.numAtoms() {
        return units
    }
    // Compute which automorphisms cause inversion of configuration
    // for the stereogenic units
    let inverted = StereoInverted.compute(mol, &symClasses, automorphisms)
    // std::vector<OBBitVec> mergedRings = mergeRings(mol, symClasses);
    let mergedRings = mergeRings(mol, symClasses)
    // std::vector<unsigned long> doneAtoms, doneBonds;
    var doneAtoms: [UInt] = []
    var doneBonds: [UInt] = []
    // unsigned int lastSize = units.size();
    var lastSize = units.count
    // while (true) {
    while true {
        
        for atom in mol.getAtomIterator() {
            if doneAtoms.contains(UInt(atom.getId().intValue)) {
                continue
            }
            // consider only potential steroecenters
            if !isPotentialTetrahedral(atom) {
                continue
            }
            // A potential stereocenter is really a stereocenter if there exists no automorphic
            // permutation causing an inversion of the configuration of only the potential
            // stereogenic unit under consideration.
            var foundPermutation: Bool = false
            for i in 0..<inverted.count {
                let atoms = inverted[i].invertedAtoms
                if atoms.count != 1 {
                    continue
                }
                let bonds = inverted[i].invertedBonds
                if bonds.count != 0 {
                    continue
                }
                if atoms[0] == atom {
                    foundPermutation = true
                    break
                }
            }
            
            var classification = classifyTetrahedralNbrSymClasses(symClasses, atom)
            
            if !foundPermutation {
                // true-stereocenter found
                let isParaCenter = (classification == .T1234) ? false : true
                units.append(MKStereoUnit(.Tetrahedral, atom.getId(), isParaCenter))
                doneAtoms.append(UInt(atom.getId().intValue))
            } else {
                // count ligand configurations:
                // If there exists at least one automorphic permutation causing the inversion of the
                // configuration of only the stereogenic unit under consideration, then the potential
                // stereocenter can be a stereocenter if the number of topologically equivalent neighbors
                // (ligands) of potential stereogenic is less than or equal to the number of configurations
                // of these ligands.
                //
                // In practise:
                //    T1123 -> 1 true stereocenter OR 2 para stereocenters
                //    T1122 -> 1 true stereocenter OR 2 para stereocenters (for both)
                //    T1112 -> 2 true stereocenters OR 2 para stereocenter assemblies
                //    T1111 -> 2 true stereocenters OR 2 para stereocenter assemblies
                switch classification {
                case .T1123:
                    // rule 2a with 1 pair
                    let duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses: symClasses)
                    guard let ligandAtom = findAtomWithSymmetryClass(atom, symClass: duplicatedSymClass, symClasses: symClasses) else { break }
                    if containsAtLeast_1true_2para(ligandAtom, atom: atom, units: units) {
                        units.append(MKStereoUnit(.Tetrahedral, atom.getId(), true))
                        doneAtoms.append(UInt(atom.getId().intValue))
                    }
                case .T1122:
                    // rule 2a with 2 pairs
                    var duplicatedSymClass1: UInt = 0
                    var duplicatedSymClass2: UInt = 0
                    findDuplicatedSymmetryClasses(atom, symClasses: symClasses, duplicated1: &duplicatedSymClass1, duplicated2: &duplicatedSymClass2)
                    guard let ligandAtom1 = findAtomWithSymmetryClass(atom, symClass: duplicatedSymClass1, symClasses: symClasses) else { break }
                    guard let ligandAtom2 = findAtomWithSymmetryClass(atom, symClass: duplicatedSymClass2, symClasses: symClasses) else { break }
                    if (containsAtLeast_1true_2para(ligandAtom1, atom: atom, units: units) &&
                        containsAtLeast_1true_2para(ligandAtom2, atom: atom, units: units)) {
                        units.append(MKStereoUnit(.Tetrahedral, atom.getId(), true))
                        doneAtoms.append(UInt(atom.getId().intValue))
                    }
                case .T1112, .T1111:
                    // rule 2b with 4 identical
                    let duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses: symClasses)
                    guard let ligandAtom = findAtomWithSymmetryClass(atom, symClass: duplicatedSymClass, symClasses: symClasses) else { break }
                    if containsAtLeast_2true_2paraAssemblies(ligandAtom, atom: atom, units: units, mergedRings: mergedRings) {
                        units.append(MKStereoUnit(.Tetrahedral, atom.getId(), true))
                        doneAtoms.append(UInt(atom.getId().intValue))
                    }
                default: break
                }
            }
        }

        for bond in mol.getBondIterator() {
            if doneBonds.contains(UInt(bond.getId().intValue)) {
                continue
            }
            // consider only potential steroecenters
            if !isPotentialCisTrans(bond) {
                continue
            }
            // A double bond is a stereogenic bond if there exists no automorphic
            // permutation causing an inversion of the configuration of only the potential
            // stereogenic unit under consideration.
            var foundPermutation: Bool = false
            for i in 0..<inverted.count {
                let atoms = inverted[i].invertedAtoms
                if atoms.count != 0 {
                    continue
                }
                let bonds = inverted[i].invertedBonds
                if bonds.count != 1 {
                    continue
                }
                if bonds[0] == bond {
                    foundPermutation = true
                    break
                }
            }

            let beginClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond.getBeginAtom())
            let endClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond.getEndAtom())
            if !foundPermutation {
                //  true stereocenter found
                var isParaCenter = (beginClassification == .C12) && (endClassification == .C12) ? false : true
                units.append(MKStereoUnit(.CisTrans, bond.getId(), isParaCenter))
                doneBonds.append(UInt(bond.getId().intValue))
            } else {
                // count ligand configurations:
                var beginValid: Bool = false
                switch beginClassification {
                case .C12:
                    beginValid = true
                    break
                case .C11:
                    // find the ligand
                    var ligandAtom: MKAtom? = nil
                    for nbr in bond.getBeginAtom().getNbrAtomIterator()! {
                        if (nbr.getIdx() != bond.getBeginAtomIdx()) && (nbr.getIdx() != bond.getEndAtomIdx()) {
                            ligandAtom = nbr
                            break
                        }
                    }
                    if ligandAtom != nil {
                        beginValid = containsAtLeast_1true_1para(ligandAtom!, skip: bond.getBeginAtom(), units: units)
                    }
                default: break
                }
                if !beginValid {
                    continue
                }
                var endValid: Bool = false
                switch endClassification {
                case .C12:
                    endValid = true
                    break
                case .C11:
                    // find the ligand
                    var ligandAtom: MKAtom? = nil
                    for nbr in bond.getEndAtom().getNbrAtomIterator()! {
                        if (nbr.getIdx() != bond.getBeginAtomIdx()) && (nbr.getIdx() != bond.getEndAtomIdx()) {
                            ligandAtom = nbr
                            break
                        }
                    }
                    if ligandAtom != nil {
                        endValid = containsAtLeast_1true_1para(ligandAtom!, skip: bond.getEndAtom(), units: units)
                    }
                default: break
                }
                if endValid {
                    units.append(MKStereoUnit(.CisTrans, bond.getId(), true))
                    doneBonds.append(UInt(bond.getId().intValue))
                }
            }
        }
        if units.count == lastSize {
            break
        }
        lastSize = units.count
    }

//    MARK: FOR DEBUG ONLY
    
    for unit in units {
        if unit.type == .Tetrahedral {
            print("Tetrahedral(center = \(unit.id), para = \(unit.para))")
        }
        if unit.type == .CisTrans {
            print("CisTrans(center = \(unit.id), para = \(unit.para))")
        }
        if unit.type == .SquarePlanar {
            print("SquarePlanar(center = \(unit.id), para = \(unit.para))")
        }
    }

    return units
}


/**
* Perform symmetry analysis.
*
* @return vector containing symmetry classes index by OBAtom::GetIndex().
*/
func findSymmetry(_ mol: MKMol) -> [Ref] {
    var symVec: MKBitVec? = nil
    let symmetry = MKGraphSym(mol, &symVec)
    var symClasses = [Ref](repeating: .NoRef, count: mol.numAtoms())
    symmetry.getSymmetry(&symClasses)
    return symClasses
}

