//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/20/23.
//

import Foundation

class MKResidueData: MKGlobalDataBase {
    
    var _resnum: Int = 0
    var _resname: [String] = []
    var _resatoms: [[String]] = []
    var _resbonds: [[Pair<String, Int>]] = []
    
    //variables used only temporarily for parsing resdata.txt
    var _vatmtmp: [String] = []
    var _vtmp: [Pair<String, Int>] = []
        
    init() {
        super.init(fileName: "resdata", subDir: "Data")
        self.readFile()
    }
    
    override func getSize() -> Int {
        return self._resname.count
    }
    
    override func readFile() {
        //  Try to load contents of file
        guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }

        filePath.foreachRow { rowContent, lineNum in
            if !rowContent.starts(with: "#") {
                var lineContents = rowContent.components(separatedBy: .whitespaces)
                lineContents.removeAll(where: { $0 == "" })
                if !lineContents.isEmpty {
                    
                    if lineContents[0] == "BOND" {
                        let b1 = lineContents[1]
                        let b2 = lineContents[2]
                        let s = b1 < b2 ? "\(b1) \(b2)" : "\(b2) \(b1)"
                        guard var bo = lineContents[3].toInt() else {
                            print("Could not parse bond order \(lineContents[3]) in resdata.txt, line \(lineNum)")
                            return
                        }
                        _vtmp.append((s, bo))
                    }
                    
                    if lineContents[0] == "ATOM" && lineContents.count == 4 {
                        _vatmtmp.append(lineContents[1])
                        _vatmtmp.append(lineContents[2])
                        _vatmtmp.append(lineContents[3])
                    }
                    
                    if lineContents[0] == "RES" {
                        _resname.append(lineContents[1])
                    }
                    
                    if lineContents[0] == "END" {
                        _resatoms.append(_vatmtmp)
                        _resbonds.append(_vtmp)
                        _vatmtmp = []
                        _vtmp = []
                    }
                    
                }
            }
        }
    }
    
    //! Sets the table to access the residue information for a specified
    //!  residue name
    //! \return whether this residue name is in the table
    func setResName(_ name: String) -> Bool {
        for i in 0..<self._resname.count {
            if self._resname[i] == name {
                self._resnum = i
                return true
            }
        }
        _resnum = -1 
        return false
    }
    
    //! \return the bond order for the bond specified in the current residue
    //! \deprecated Easier to use the two-argument form
    func lookupBO(_ name: String) -> Int {
        if _resnum == -1 { return 0 }
        for i in _resbonds[_resnum] {
            if i.0 == name {
                return i.1
            }
        }
        return 0
    }
    //! \return the bond order for the bond specified between the two specified
    //! atom labels
    func lookupBO(_ name: String, _ name1: String) -> Int {
        if _resnum == -1 { return 0 }
        let s = name < name1 ? "\(name) \(name1)" : "\(name1) \(name)"
        for i in _resbonds[_resnum] {
            if i.0 == s {
                return i.1
            }
        }
        return 0
    }
    
    //! Look up the atom type and hybridization for the atom label specified
    //! in the first argument for the current residue
    //! \return whether the atom label specified is found in the current residue
    func lookupType(_ atmid: String, _ type: inout String, _ hyb: inout Int) -> Bool {
        if _resnum == -1 { return false }

        for i in _resatoms[_resnum] {
            if i[0] == atmid {
                type = i[1]
                hyb = i[2].toInt()!
                return true
            }
        }
        return false
    }
    
    //! Assign bond orders, atom types and residues for the supplied OBMol
    //! based on the residue information assigned to atoms
    func assignBonds(_ mol: inout MKMol) -> Bool {
        
        var bo: Int = 0
        var skipres: String = "" // Residue Number to skip
        var rname: String = "" 
        //assign residue bonds
        for i in 0..<mol.numAtoms() {
            guard let a1 = mol.getAtom(i) else {
                print("Could not get atom \(i)")
                return false
            }
            guard let r1 = a1.getResidue() else { 
                // print("Could not get residue for atom \(a1)") // atoms may not have residues
                return false
            }
            if skipres.length > 0 && skipres == r1.getNumString() {
                continue
            }
            if r1.getName() != rname {
                skipres = setResName(r1.getName()) ? "" : r1.getNumString()
                rname = r1.getName()
            }
            //assign bonds for each atom
            for j in i..<mol.numAtoms() {
                guard let a2 = mol.getAtom(j) else {
                    print("Could not get atom \(j)")
                    return false
                }
                guard let r2 = a2.getResidue() else { 
                    // print("Could not get residue for atom \(a2)") // atoms may not have residues
                    return false
                }
                if r1.getName() != r2.getName() || r1.getNumString() != r2.getNumString() || r1.getChain() != r2.getChain() {
                                                                                        // Fixes PR#2889763 - Fabian
                    break
                }
                bo = lookupBO(r1.getAtomID(a1), r2.getAtomID(a2))
                if bo != 0 {
                    if bo > 0 {
                        // Suggested by Liu Zhiguo 2007-08-13
                        // for predefined residues, don't perceive connection
                        // by distance
                        //                v = a1->GetVector() - a2->GetVector();
                        //                if (v.length_2() < 3.5) //check by distance
                        mol.addBond(a1.getIdx(), a2.getIdx(), bo)
                    }
                }
            }
        }

        var hyb: Int = 0
        var type: String = ""
        //types and hybridization
        skipres = ""
        rname = ""

        for a1 in mol.getAtomIterator() {
            if a1.getAtomicNum() == MKElements.Oxygen.atomicNum && a1.getExplicitDegree() == 0 {
                a1.setType("O3")
                continue
            }
            if a1.getAtomicNum() == MKElements.Hydrogen.atomicNum {
                a1.setType("H")
                continue
            }
            // valence rule for O-
            if a1.getAtomicNum() == MKElements.Oxygen.atomicNum && a1.getExplicitDegree() == 1 {
                
                guard let bond = a1.getBondIterator()!.next() else {
                    print("Could not get bond for atom \(a1)")
                    return false
                }
                if bond.getBondOrder() == 2 {
                    a1.setType("O2")
                    a1.setHyb(2)
                } else if bond.getBondOrder() == 3 {
                    // Leave the protonation/deprotonation to phmodel
                    a1.setType("O3")
                    a1.setHyb(3)
                    // PR#3203039 -- Fix from Magnus Lundborg
                    //                a1->SetFormalCharge(0);
                }
                continue
            }

            guard let r1 = a1.getResidue() else { 
                // print("Could not get residue for atom \(a1)") // atoms may not have residues
                return false
            }
            if skipres.length > 0 && skipres == r1.getNumString() {
                continue
            }
            if r1.getName() != rname {
                // if SetResName fails, skip this residue
                skipres = setResName(r1.getName()) ? "" : r1.getNumString()
                rname = r1.getName()
                // Should there be a continue here??
            }
            if lookupType(r1.getAtomID(a1), &type, &hyb) {
                a1.setType(type)
                a1.setHyb(hyb)
            } else {} // try to figure it out by bond order ???
        }
        return true
    }
}
