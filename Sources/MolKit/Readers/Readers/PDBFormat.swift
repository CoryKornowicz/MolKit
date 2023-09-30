//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 9/26/23.
//

import Foundation

class PDBFormat: MKMoleculeFormat {
    
    override init() {
        super.init()
        MKConversion.registerFormat("pdb", self, "chemical/x-pdb")
        MKConversion.registerFormat("ent", self, "chemical/x-pdb")
        
        MKConversion.registerOptionParam("s", self, 0, .INOPTIONS)
        MKConversion.registerOptionParam("b", self, 0, .INOPTIONS)
        MKConversion.registerOptionParam("c", self, 0, .INOPTIONS)
        
        MKConversion.registerOptionParam("o", self, 0, .OUTOPTIONS)
        MKConversion.registerOptionParam("n", self, 0, .OUTOPTIONS)
    }
    
    required init(_ id: String, _ isDefault: Bool) {
        fatalError("init(_:_:) has not been implemented")
    }
    
    override func description() -> String? {
        return "Protein Data Bank format\n" +
                "Read Options e.g. -as\n" +
                "  s  Output single bonds only\n" +
                "  b  Disable bonding entirely\n" +
                "  c  Ignore CONECT records\n\n" +
                "Write Options, e.g. -xo\n" +
                "  n  Do not write duplicate CONECT records to indicate bond order\n" +
                "  o  Write origin in space group label (CRYST1 section)\n\n"
    }
    
    override func specificationURL() -> String {
        return "http://www.wwpdb.org/docs.html"
    }
    
    override func getMIMEType() -> String {
        return "chemical/x-pdb"
    }
    
    override func skipObjects(_ n: inout Int, _ pConv: MKConversion) -> Int {
        if n == 0 {
            ++n
        }
        if let ifs = pConv.getInStream() {
            var buffer: String = String()
            while n > 0 && ifs.getLine(&buffer) {
                if buffer.starts(with: "ENDMDL") {
                    --n
                }
            }
            return ifs.isEOF ? -1 : 1
        } else {
            return -1
        }
    }
    
    override func readMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        fatalError()
    }
    
    override func writeMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        fatalError()
    }
    
}
