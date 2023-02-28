

import Foundation


class MKChemTsfm {
    
    var _bgn: MKSmartsPattern = MKSmartsPattern()
    var _end: MKSmartsPattern = MKSmartsPattern()
    var _vadel: [Int] = []
    var _vele: [Pair<Int, Int>] = []
    var _vchrg: [Pair<Int, Int>] = []
    var _vbdel: [Pair<Int, Int>] = []
    var _vbond: [Pair<Pair<Int, Int>, Int>] = []
    
    init () { }
    
    //! Initialize this transformation with the supplied SMARTS patterns
    func initialize(_ bgn: String, _ end: String) -> Bool {
        if !_bgn.initialize(bgn) { return false }
        if !end.isEmpty { if !_end.initialize(end) { return false } }
        
        //find atoms to be deleted
        var vb: Int? = 0
        var found: Bool = false
        
        for i in 0..<_bgn.numAtoms() {
            if let vb = _bgn.getVectorBinding(i) {
                found = false
                for j in 0..<_end.numAtoms() {
                    if vb == _end.getVectorBinding(j) {
                        found = true
                        break
                    }
                }
                if !found {
                    _vadel.append(i)
                }
            }
        }
        
        //find elements to be changed
        for i in 0..<_bgn.numAtoms() {
            // Allow single-atom transformations without vector bindings
            vb = _bgn.getVectorBinding(i)
            if vb != nil || _bgn.numAtoms() == 1 {
                guard let ele = _bgn.getAtomicNum(i) else { break }
                for j in 0..<_end.numAtoms() {
                    if vb == _end.getVectorBinding(j) {
                        guard let elej = _end.getAtomicNum(j) else { break }
                        if ele != elej {
                            _vele.append((i, elej))
                            break
                        }
                    }
                }
            }
        }
        
        //find charges to modify
        for i in 0..<_bgn.numAtoms() {
            // Allow single-atom transformations without vector bindings
            // PR#2802980.
            vb = _bgn.getVectorBinding(i)
            if vb != nil || _bgn.numAtoms() == 1 {
                guard let chrg = _bgn.getCharge(i) else { break }
                for j in 0..<_end.numAtoms() {
                    if vb == _end.getVectorBinding(j) {
                        guard let chrgj = _end.getCharge(j) else { break }
                        if chrg != chrgj {
                            _vchrg.append((i, chrgj))
                            break
                        }
                    }
                }
            }
        }
        
        //find bonds to be modified
        
        var bsrc = 0
        var bdst = 0
        var bord = 0
        var esrc = 0
        var edst = 0
        var eord = 0
        
        var bvb1, bvb2: Int?
        var evb1, evb2: Int?
        
        for i in 0..<_bgn.numBonds() {
            _bgn.getBond(&bsrc, &bdst, &bord, i)
            bvb1 = _bgn.getVectorBinding(bsrc)
            bvb2 = _bgn.getVectorBinding(bdst)
            if bvb1 == nil || bvb2 == nil { continue }
            
            for j in 0..<_end.numBonds() {
                _end.getBond(&esrc, &edst, &eord, j)
                evb1 = _bgn.getVectorBinding(esrc)
                evb2 = _bgn.getVectorBinding(edst)
                if ((bvb1 == evb1 && bvb2 == evb2) || (bvb1 == evb2 && bvb2 == evb1)) {
                    if bord == eord {
                        break //nothing to modify if bond orders identical
                    }
                    _vbond.append(((bsrc, bdst), eord))
                    break
                }
            }
        }
        
        //make sure there is some kind of transform to do here
        if _vadel.isEmpty && _vchrg.isEmpty && _vbond.isEmpty && _vele.isEmpty {
            return false
        }
        
        return true
    }
    
    func apply(_ mol: MKMol) -> Bool {
        
        if !_bgn.match(mol) { return false }
        
        mol.beginModify()
        
        let mlist = _bgn.getUMapList()
        
        if !_vchrg.isEmpty {  //modify charges
            for i in mlist {
                for j in _vchrg {
                    if j.0 < i.count { //goof proofing
//                        MARK: maybe throw error here instead
                        guard let atom = mol.getAtom(i[j.0]) else { continue }
                        let old_charge = atom.getFormalCharge()
                        if j.1 != old_charge {
                            atom.setFormalCharge(j.1)
                            MKAtomAssignTypicalImplicitHydrogens(atom) //update with new charge info
                        }
                    }
                }
            }
        }
        
        if !_vbond.isEmpty { //modify bond orders
            for i in mlist {
                for j in _vbond {
//                    MARK: need to print and error check here better
                    guard let bond = mol.getBond(i[j.0.0], i[j.0.1]) else { continue }
                    let old_bond_order = bond.getBondOrder()
                    bond.setBondOrder(UInt(j.1))
                    for k in 0..<2 {
                        let atom = k == 0 ? bond.getBeginAtom() : bond.getEndAtom()
                        let new_hcount = atom.getImplicitHCount() - (UInt(j.1) - old_bond_order)
//                        MARK: using UInt here could throw and underflow error!!
                        atom.setImplicitHCount(UInt(new_hcount))
                    }
                }
            }
        }
        
        if !_vele.isEmpty || !_vadel.isEmpty { // delete atoms and change elements
            
            if !_vele.isEmpty {
                for i in mlist {
                    for k in _vele {
                        guard let atom = mol.getAtom(i[k.0]) else { continue }
                        atom.setAtomicNum(k.1)
                    }
                }
            }
            
            //make sure same atom isn't deleted twice
            var vda: [Bool] = [Bool].init(repeating: false, count: mol.numAtoms()+1)
            var vdel: [MKAtom] = []
            
            for i in mlist {
                for j in _vadel {
                    if !vda[i[j]] {
                        vda[i[j]] = true
                        guard let atom = mol.getAtom(i[j]) else { break }
                        vdel.append(atom)
                    }
                }
            }
            
            for k in vdel {
                mol.deleteAtom(k)
            }
        }
        
        mol.endModify()
        return true
    }
    
    /*! Is this transformation an acid dissociation?
     *  \code
     *      Ka
     *  HA ----> A(-)         (the H(+) will be deleted)
     *  \endcode
     *
     *  IsAcid() will check the charge in the end SMARTS pattern.
     *  \return true if the charge is less than 0 (-1).
     */
    func isAcid() -> Bool {
        
        if _bgn.numAtoms() > _end.numAtoms() { return true } // O=CO[#1:1] >> O=CO
        
        for i in 0..<_end.numAtoms() {
            if let chrg = _end.getCharge(i) {
                if chrg < 0 {
                    return true
                }
            }
        }
        return false
    }
    /*! Is this a transformation to the conjugated acid from a base?
     *  \code
     *      Ka
     *  HA ----> A(-)         (the H(+) will be deleted)
     *  \endcode
     *
     *  IsBase() will check the charge in the end SMARTS pattern.
     *  \return true if the charge is higher than 0 (+1).
     */
    func isBase() -> Bool {
        for i in 0..<_end.numAtoms() {
            if let chrg = _end.getCharge(i) {
                if chrg > 0 {
                    return true
                }
            }
        }
        return false
    }
    
}


/*! \brief Corrections for pH used by OBMol::CorrectForPH()
 *
 *  The data/phmodel.txt file contains transformations which are applied
 *  to correct the charges for a given pH. This function uses the
 *  Henderson-Hasselbalch equation to calculate which species (protonated/
 *  unprotonated) is present in the highest concentration at the given pH.
 *
 *  For acids an entry would look like:
 *  \code
 *  # carboxylic acid
 *  O=C[OD1:1] >> O=C[O-:1]    4.0
 *  \endcode
 *
 *  The 4.0 is the pKa for the dissociation [HA] -> [H+] + [A-]. To
 *  calculate [HA]/[A-] we use:
 *  \code
 *  [HA] / [A-] = 10^(pKa - pH)
 *
 *  [HA]/[A-] > 1  :  [HA] > [A-]
 *  [HA]/[A-] < 1  :  [A-] > [HA]
 *  \endcode
 *
 *  For a base, an entry would look be:
 *  \code
 *  # methyl amine
 *  C[N:1] >> C[N+:1]    10.7
 *  \endcode
 *
 *  Here, the 10.7 is the pKa for the dissociation [BH+] -> [H+] + [B:]. To
 *  calculate [BH+]/[B:] we use:
 *  \code
 *  [BH+] / [B:] = 10^(pKa - pH)
 *
 *  [BH+]/[B:] > 1  :  [BH+] > [B:]
 *  [BH+]/[B:] < 1  :  [B:] > [BH+]
 *  \endcode
 *
 *  The transformations are all applied (if needed at the specified pH value) in
 *  the same order they are found in data/phmodel.txt.
 */

class MKPhModel: MKGlobalDataBase {
    
    var _vtsfm: [MKChemTsfm] = []
    var _vpka: [Double] = []
    var _vschrg: [Pair<MKSmartsPattern, [Double]>] = []

    init() {
        super.init(fileName: "phmodel", subDir: "Data")
        self.readFile()
    }
    
    override func readFile() {
        guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }
        
        filePath.foreachRow { rowContent, lineNum in
            if rowContent.starts(with: "TRANSFORM") {
                var vs = rowContent.components(separatedBy: .whitespaces)
                vs.removeAll(where: { $0 == "" })
                if vs.count < 5 {
                    //                  TODO: Handle error, log
                    print("Could not parse line in phmodel table from phmodel.txt")
                }
                let tsfm = MKChemTsfm()
                if !tsfm.initialize(vs[1], vs[3]) {
                    print("Could not parse line in phmodel table from phmodel.txt")
                }
                _vtsfm.append(tsfm)
                _vpka.append(String(vs[4]).toDouble()!)
            } else if rowContent.starts(with: "SEEDCHARGE") {
                var vs = rowContent.components(separatedBy: .whitespaces)
                vs.removeAll(where: { $0 == "" })
                if vs.count < 2 {
                    print("Could not parse line in phmodel table from phmodel.txt")
                }
                let sp = MKSmartsPattern()
                if !sp.initialize(vs[1]) {
                    print("Could not parse line in phmodel table from phmodel.txt")
                }
                if (vs.count-2) != sp.numAtoms() {
                    print("Could not parse line in phmodel table from phmodel.txt")
                }
                var vf : [Double] = []
                for i in vs[2...] {
                    vf.append(i.toDouble()!)
                }
                _vschrg.append((sp, vf))
            }
        }
    }
    
    override func getSize() -> Int {
        return _vschrg.count
    }
    
    
    func assignSeedPartialCharge(_ mol: MKMol) {
        mol.setPartialChargesPerceived()
        
        if !mol.automaticPartialCharge() { return }
        
        for i in _vschrg {
            var mlist: [[Int]] = []
            if i.0.match(mol, &mlist, .AllUnique) {
                for j in mlist {
                    for k in 0..<j.count{
                        guard let atom = mol.getAtom(j[k]) else { break }
                        atom.setPartialCharge(i.1[k])
                    }
                }
            }
        }
    }
    
    func correctForPH(_ mol: MKMol, _ pH: Double = 7.4) {
        
        if mol.isCorrectedForPH() { return }
        if mol.getDimension() > 0 && !mol.automaticPartialCharge() { return }
        
        mol.setCorrectedForPH()
        
        mol.deleteHydrogens()
        
        for i in 0..<_vtsfm.count {
            if _vpka[i] > 1e+9 {
                // always apply when pKa is > 1e+9
                _vtsfm[i].apply(mol)
            } else {
                // 10^(pKa - pH) = [HA] / [A-]
                //
                // > 1 : [HA] > [A-]
                // < 1 : [HA] < [A-]
                if _vtsfm[i].isAcid() {
                    if pow(10, (_vpka[i] - pH)) < 1.0 {
                        _vtsfm[i].apply(mol)
                    }
                }
                // 10^(pKa - pH) = [BH+] / [B:]
                //
                // > 1 : [BH+] > [B:]
                // < 1 : [BH+] < [B:]
                if _vtsfm[i].isBase() {
                    if pow(10, (_vpka[i] - pH)) > 1.0 {
                        _vtsfm[i].apply(mol)
                    }
                }
            }
        }
    }
}
