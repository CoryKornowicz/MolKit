//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/21/23.
//

import Foundation

class MKBondTyper: MKGlobalDataBase {
    
    var _fgbonds: [Pair<MKSmartsPattern, [Int]>] = []
    
    init() {
        super.init(fileName: "bondtyp", subDir: "Data")
        self.readFile()
    }
    
    override func readFile() {
        guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }
        
        filePath.foreachRow { rowContent, lineNum in
            if !rowContent.starts(with: "#") {
                var bovector: [Int] = []
                var vs = rowContent.components(separatedBy: .whitespaces)
                vs.removeAll(where: { $0 == "" })
                // Make sure we actually have a SMARTS pattern plus at least one triple
                // and make sure we have the correct number of integers
                if vs.count < 4 { return } // just ignore empty (or short lines)
                else if vs.count >= 4 && (vs.count % 3 != 1) {
//                    TODO: probably need to check error handling here
                    print("Error in OBBondTyper. Pattern is incorrect, found \(vs.count) tokens.")
                    print("Buffer is: \(rowContent)")
                    return
                }
                
                let sp = MKSmartsPattern()
                if sp.initialize(vs[0]) {
                    for i in 1..<vs.count {
                        guard let bo = String(vs[i]).toInt() else { return }
                        bovector.append(bo)
                    }
                    _fgbonds.append((sp, bovector))
                }
            }
        }
    }
    
    override func getSize() -> Int {
        return self._fgbonds.count
    }
    
    //! \name Bond Perception Routines
    //@{
    //! Assign bonds to functional groups based on the bond typer database
    func assignFunctionalGroupBonds(_ mol: MKMol) {
     
        var currentPattern: MKSmartsPattern
        var assignments: [Int]
        var mlist: [[Int]]
        // Loop through for all the functional groups and assign bond orders
        for i in _fgbonds {
            
            currentPattern = i.0
            assignments = i.1
            
            if currentPattern.match(mol) {
                mlist = currentPattern.getUMapList()
                for matchs in mlist {
                    // Now loop through the bonds to assign from _fgbonds chucked into three elements
                    // at a time
                    for j in assignments.chunks(ofCount: 3) {
                        // along the assignments vector: atomID1 atomID2 bondOrder
                        guard let a1 = mol.getAtom(matchs[assignments[j[0]]]) else { continue }
                        guard let a2 = mol.getAtom(matchs[assignments[j[1]]]) else { continue }
                        
                        guard let b1 = a1.getBond(a2) else { continue }
                        b1.setBondOrder(UInt(j[2]))
                    } // bond order assignments
                }  // each match
            } // current pattern matches
        } // for(functional groups)
        
        // FG with distance and/or bond criteria
        // Carbonyl oxygen C=O (O must be neutral)
        let carbo: MKSmartsPattern = MKSmartsPattern()
        carbo.initialize("[#8D1;!-][#6](*)(*)")
        
        if carbo.match(mol) {
            mlist = carbo.getUMapList()
            for l in mlist {
                guard let a1 = mol.getAtom(l[0]) else { continue }
                guard let a2 = mol.getAtom(l[1]) else { continue }
                
                let angle = a2.averageBondAngle()
                let dist1 = a1.getDistance(a2)
                
                // carbonyl geometries ?/
                if angle > 115 && angle < 150 && dist1 < 1.28 {
                    if !a1.hasDoubleBond() { // no double bond already assigned
                        guard let b1 = a1.getBond(a2) else { continue }
                        b1.setBondOrder(2)
                    }
                }
            }
        }// Carbonyl oxygen
        
        // thione C=S
        let thione: MKSmartsPattern = MKSmartsPattern()
        thione.initialize("[#16D1][#6](*)(*)")
        
        if thione.match(mol) {
            mlist = thione.getUMapList()
            for l in mlist {
                guard let a1 = mol.getAtom(l[0]) else { continue }
                guard let a2 = mol.getAtom(l[1]) else { continue }
                
                let angle = a2.averageBondAngle()
                let dist1 = a1.getDistance(a2)
                
                // thione geometries ?/
                if angle > 115 && angle < 150 && dist1 < 1.72 {
                    if !a1.hasDoubleBond() { // no double bond already assigned
                        guard let b1 = a1.getBond(a2) else { continue }
                        b1.setBondOrder(2)
                    }
                }
            }
        } // thione
        
        // Isocyanate N=C=O or Isothiocyanate
        var dist1OK: Bool = false
        let isocyanate: MKSmartsPattern = MKSmartsPattern()
        isocyanate.initialize("[#8,#16;D1][#6D2][#7D2]")
        
        if isocyanate.match(mol) {
            mlist = isocyanate.getUMapList()
            for l in mlist {
                guard let a1 = mol.getAtom(l[0]) else { continue }
                guard let a2 = mol.getAtom(l[1]) else { continue }
                guard let a3 = mol.getAtom(l[2]) else { continue }
                
                let angle = a2.averageBondAngle()
                let dist1 = a1.getDistance(a2)
                let dist2 = a2.getDistance(a3)
                
                // isocyanate geometry or Isothiocyanate geometry ?
                if a1.getAtomicNum() == MKElements.getAtomicNum("O") {
                    dist1OK = dist1 < 1.28
                } else {
                    dist1OK = dist1 < 1.72
                }
                
                if angle > 150 && dist1OK && dist2 < 1.34 {
                    guard let b1 = a1.getBond(a2) else { continue }
                    guard let b2 = a2.getBond(a3) else { continue }
                    
                    b1.setBondOrder(2)
                    b2.setBondOrder(2)
                }
            }
        } // isocyanate
        
        // oxime C=S
        let oxime: MKSmartsPattern = MKSmartsPattern()
        oxime.initialize("[#6D3][#7D2][#8D2]")

        if oxime.match(mol) {
            mlist = oxime.getUMapList()
            for l in mlist {
                guard let a1 = mol.getAtom(l[0]) else { continue }
                guard let a2 = mol.getAtom(l[1]) else { continue }
                
                let angle = a2.averageBondAngle()
                let dist1 = a1.getDistance(a2)
                
                // oxime geometries ?/
                if angle > 110 && angle < 150 && dist1 < 1.4 {
                    if !a1.hasDoubleBond() { // no double bond already assigned
                        guard let b1 = a1.getBond(a2) else { continue }
                        b1.setBondOrder(2)
                    }
                }
            }
        } // oxime
        
        // oxido-n+ (e.g., pyridine-N-oxide)
        let oxidopyr: MKSmartsPattern = MKSmartsPattern()
        oxidopyr.initialize("[#8D1][#7D3r6]")

        if oxidopyr.match(mol) {
            mlist = oxidopyr.getUMapList()
            for l in mlist {
                guard let a1 = mol.getAtom(l[0]) else { continue }
                guard let a2 = mol.getAtom(l[1]) else { continue }
                
                let angle = a2.averageBondAngle()
                let dist1 = a1.getDistance(a2)
                
                if angle > 110 && angle < 150 && dist1 < 1.35 {
                    //  TODO: Assert these are being added to the correct atoms
                    a1.setFormalCharge(-1) // oxygen
                    a2.setFormalCharge(+1) // nitrogen
                }
            }
        } // oxido-n+
    }
}
