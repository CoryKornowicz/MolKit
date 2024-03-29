//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/20/23.
//

import Foundation

/** \class OBTypeTable data.h <openbabel/data.h>
      \brief Atom Type Translation Table

      Molecular file formats frequently store information about atoms in an
      atom type field. Some formats store only the element for each atom,
      while others include hybridization and local environments, such as the
      Sybyl mol2 atom type field. The OBTypeTable class acts as a translation
      table to convert atom types between a number of different molecular
      file formats. The constructor for OBTypeTable automatically reads the
      text file types.txt. An instance of
      OBTypeTable (ttab) is declared external in data.cpp and is referenced as
      extern OBTypeTable ttab in mol.h.  The following code demonstrates how
      to use the OBTypeTable class to translate the internal representation
      of atom types in an OBMol Internal to Sybyl Mol2 atom types.

      \code
      ttab.SetFromType("INT");
      ttab.SetToType("SYB");
      OBAtom *atom;
      vector<OBAtom*>::iterator i;
      string src,dst;
      for (atom = mol.BeginAtom(i);atom;atom = mol.EndAtom(i))
      {
         src = atom->GetType();
         ttab.Translate(dst,src);
         cout << "atom number " << atom->GetIdx() << "has mol2 type " << dst << endl;
      }
      \endcode

      Current atom types include (defined in the top line of the data file types.txt):
      - INT (Open Babel internal codes)
      - ATN (atomic numbers)
      - HYB (hybridization)
      - MMD (MacroModel)
      - MM2 (MM2 force field)
      - XYZ (element symbols from XYZ file format)
      - ALC (Alchemy)
      - HAD (H added)
      - MCML (MacMolecule)
      - C3D (Chem3D)
      - SYB (Sybyl mol2)
      - MOL (Sybyl mol)
      - MAP (Gasteiger partial charge types)
      - DRE (Dreiding)
      - XED (XED format)
      - DOK (Dock)
      - M3D (Molecular Arts M3D)
      - SBN (Sybyl descriptor types for MPD files)
      - PCM (PC Model)
  */
class MKTypeTable: MKGlobalDataBase {
    
    var _linecount: Int = 0
    var _ncols: Int = 0
    var _nrows: Int = 0
    var _from: Int = -1
    var _to: Int = -1
    var _colnames: [String] = []
    var _table: [[String]] = []
    
    init() {
        super.init(fileName: "types", subDir: "Data")
        self.readFile()
    }
    
    override func getSize() -> Int {
        return self._table.count
    }
    
    override func readFile() {
        //        Try to load contents of file
        guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }
        
        filePath.foreachRow { rowContent, lineNum in
            
            if !rowContent.starts(with: "#") {
                if self._linecount == 0 {
                    self._colnames = rowContent.components(separatedBy: .whitespaces)
                    self._ncols = self._colnames.count
                } else {
                    var newRow = rowContent.components(separatedBy: .whitespaces)
                    newRow.removeAll(where: { $0 == "" })
                    if newRow.count == self._ncols {
                        self._table.append(newRow)
                    } else {
                        print("Could not parse line in type translation table types.txt -- incorect number of columns")
                        print("found \(newRow.count) expected \(self._ncols)")
                        print(rowContent)
//                        TODO: Fix this reasonably, patch
                        self._table.append(Array(newRow[0..<self._ncols - 1]))
                    }
                }
                self._linecount += 1
            }
        }
    }
    
    func setFromType(_ from: String) -> Bool {
        var tmp: String = from 
        for i in 0..<_colnames.count {
            if _colnames[i] == tmp {
                _from = i
                return true
            }
        }
        print("Could not find type \(from) in type translation table types.txt")
        return false
    }
    
    func setToType(_ to: String) -> Bool {
        var tmp: String = to
        for i in 0..<_colnames.count {
            if _colnames[i] == tmp {
                _to = i
                return true
            }
        }
        print("Could not find type \(to) in type translation table types.txt")
        return false
    }
    
    func translate(_ from: String, _ to: inout String) -> Bool {
        if from == "" {
            return false
        }
        if _from >= 0 && _to >= 0 && _from < _table.count && _to < _table.count {
            for i in _table {
                if i.count > _from && i[_from] == from {
                    to = i[_to]
                    return true
                }
            }
        }
        //  throw werror, copy the string to to and return false
        print("Cannot perform atom type translation: table cannot find requested types.")
        to = from
        return false
    }
    
    func translate(_ from: String) -> String {
         if from == "" {
            return ""
        }
        if _from >= 0 && _to >= 0 && _from < _table.count && _to < _table.count {
            for i in _table {
                if i.count > _from && i[_from] == from {
                    return i[_to]
                }
            }
        }
        //  throw werror, copy the string to to and return false
        print("Cannot perform atom type translation: table cannot find requested types.")
        return ""
    }
    
    func getFromType() -> String {
        if _from > 0 && _from < _table.count {
            return _colnames[_from]
        } else {
            return _colnames[0]
        }
    }
    
    func getToType() -> String {
        if _to > 0 && _to < _table.count {
            return _colnames[_to]
        } else {
            return _colnames[0]
        }
    }
    
}
