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
    
    public static let sharedInstance = MKResidueData()
    
    private init() {
        super.init(fileName: "resdata", subDir: "Data")
    }
    
    override func getSize() -> Int {
        return self._resname.count
    }
    
    override func readFile() {
        //        Try to load contents of file
                guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }
        filePath.foreachRow { rowContent, lineNum in
            
        }
    }
    
    //! Sets the table to access the residue information for a specified
    //!  residue name
    //! \return whether this residue name is in the table
    func setResName(_ name: String) -> Bool {
        return false
    }
    
    //! \return the bond order for the bond specified in the current residue
    //! \deprecated Easier to use the two-argument form
    func lookupBO(_ name: String) -> Int {
        return 0
    }
    //! \return the bond order for the bond specified between the two specified
    //! atom labels
    func lookupBO(_ name: String, _ name1: String) -> Int {
        return 0
    }
    
    //! Look up the atom type and hybridization for the atom label specified
    //! in the first argument for the current residue
    //! \return whether the atom label specified is found in the current residue
    func lookupType(_ name: String, _ atom: String, _ i: Int) -> Bool {
        return false
    }
    
    //! Assign bond orders, atom types and residues for the supplied OBMol
    //! based on the residue information assigned to atoms
    func assignBonds(_ mol: inout MKMol) -> Bool {
        return false
    }
    
}
