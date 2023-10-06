//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation

  /*\brief Assigns atom types, hybridization, and formal charges
    The MKAtomTyper class is designed to read in a list of atom typing
    rules and apply them to molecules. The code that performs atom
    typing is not usually used directly as atom typing, hybridization
    assignment, and charge are all done
    automatically when their corresponding values are requested of
    atoms:
    \code
    atom->GetType();
    atom->GetFormalCharge();
    atom->GetHyb();
    \endcode
  */

public class MKAtomTyper: MKGlobalDataBase {

    private var _vinthyb: [Pair<MKSmartsPattern, Int>] = []    //!< internal hybridization rules
    private var _vexttyp: [Pair<MKSmartsPattern, String>] = [] //!< external atom type rules
    
    public init() {
        super.init(fileName: "atomtyp", subDir: "Data")
        self.readFile()
    }
    
    public func get_vinthyb() -> [Pair<MKSmartsPattern, Int>] {
        return _vinthyb
    }
    
    public func get_vettyp() -> [Pair<MKSmartsPattern, String>] {
        return _vexttyp
    }
    
    //! \return the number of internal hybridization rules
    override func getSize() -> Int {
        return _vinthyb.count
    }
    
    override func readFile() {
        guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }
        
        filePath.foreachRow { rowContent, lineNum in
//            INTHYB Line
            if rowContent.starts(with: "INTHYB") {
                var vs = rowContent.components(separatedBy: .whitespaces)
                vs.removeAll(where: { $0 == "" })
                if vs.count < 3 {
                    print("Could not parse INTHYB line in atom type table from atomtyp.txt")
                }
                let sp = MKSmartsPattern()
                if sp.initialize(vs[1]) {
                    _vinthyb.append((sp, String(vs[2]).toInt()!))
                } else {
                    print("Could not parse INTHYB line in atom type table from atomtyp.txt")
                }
            }
//            EXTTYP Line
            else if rowContent.starts(with: "EXTTYP") {
                var vs = rowContent.components(separatedBy: .whitespaces)
                vs.removeAll(where: { $0 == "" })
                if vs.count < 3 {
                    print("Could not parse EXTTYP line in atom type table from atomtyp.txt")
                }
                let sp = MKSmartsPattern()
                if sp.initialize(vs[1]) {
                    _vexttyp.append((sp, String(vs[2])))
                } else {
                    print("Could not parse EXTTYP line in atom type table from atomtyp.txt")
                }
            } 
        }
    }
    
    func assignTypes(_ mol: MKMol) {
        mol.setAtomTypesPerceived()

        for i in _vexttyp {
            var mlist: [[Int]] = []
            if i.0.match(mol, &mlist) {
                for j in mlist {
                    guard let atom = mol.getAtom(j[0]) else { break }
                    atom.setType(i.1)
                }
            }
        }

        // Special cases
        for atom in mol.getAtomIterator() {
            // guanidinium. Fixes PR#1800964
            if atom.getType() == "C2" {
                var guanidineN = 0 
                guard let nbrs = atom.getNbrAtomIterator() else { break }
                for nbr in nbrs {
                    if nbr.getType() == "Ng+" || nbr.getType() == "Npl" || nbr.getType() == "N2" {
                        guanidineN += 1
                    }
                }
                if guanidineN == 3 {
                    atom.setType("C+")
                }
            }// end C2 carbon for guanidinium
        }// end special cases
    }
    
    func assignHyb(_ mol: MKMol) {
        
        MolKit._AromTyper.assignAromaticFlags(mol)

        mol.setHybridizationPerceived()

        mol.getAtomIterator().forEach({ $0.setHyb(0) })

        for i in _vinthyb {
            var mlist: [[Int]] = []
            if i.0.match(mol, &mlist) {
                for j in mlist {
                    guard let atom = mol.getAtom(j[0]) else { break }
                    atom.setHyb(i.1)
                }
            }
        }

        // check all atoms to make sure *some* hybridization is assigned 
        for atom in mol.getAtomIterator() {
            if atom.getHyb() == 0 {
                switch atom.getExplicitDegree() {
                    case 0,1,2:
                        atom.setHyb(1)
                    case 3:
                        atom.setHyb(2)
                    case 4:
                        atom.setHyb(3)
                    default:
                        atom.setHyb(atom.getExplicitDegree())
                }
            }
        }
    }
}

