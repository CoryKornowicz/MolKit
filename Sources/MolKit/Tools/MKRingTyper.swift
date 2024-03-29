//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/21/23.
//

import Foundation

class MKRingTyper: MKGlobalDataBase {
    
    var _ringtyp: [Pair<MKSmartsPattern, String>] = [] //!< ring type rules
    
    public init() {
        super.init(fileName: "ringtyp", subDir: "Data")
        self.readFile()
    }
    
    //! \return the number of SMARTS patterns
    override func getSize() -> Int {
        return self._ringtyp.count
    }
    
    override func readFile() {
        guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }
        
        filePath.foreachRow { rowContent, lineNum in
            if rowContent.starts(with: "RINGTYP") {
                var vs = rowContent.components(separatedBy: .whitespaces)
                vs.removeAll(where: { $0 == "" })
                if vs.count < 3 {
//                  TODO: Handle error, log
                    print("Could not parse RING line in ring type table from ringtyp.txt")
                }
                let sp = MKSmartsPattern()
                if sp.initialize(vs[2]) {
                    _ringtyp.append((sp, vs[1]))
                } else {
//                  TODO: Handle error, log
                    print("Could not parse RING line in ring type table from ringtyp.txt")
                }
            }
        }
    }
    
    //! Assign external atomic types (ringtyp.txt)
    func assignTypes(_ mol: MKMol) {
        mol.setRingTypesPerceived()

        let rlist = mol.getSSSR()
        var member_count = 0 

        for i2 in _ringtyp { // for each ring type
            var mlist: [[Int]] = []
            if i2.0.match(mol, &mlist) {
                for j2 in mlist { // for each found match
                    for i in rlist { // for each ring 
                        member_count = 0
                        for j in j2 { // for each atom in the match
                            guard let atom = mol.getAtom(j) else { break }
                            if i.isMember(atom) {
                                member_count += 1
                            }
                        }
                        if i.size() == member_count {
                            i.setType(i2.1)
                        }
                    }
                }
            }
        }
    }
}
