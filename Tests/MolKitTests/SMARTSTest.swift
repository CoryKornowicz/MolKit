//
//  SMARTSTest.swift
//  
//
//  Created by Cory Kornowicz on 10/5/23.
//

import XCTest
@testable import MolKit

final class SMARTSTest: XCTestCase {

    enum SMARTSErrors: Error {
        case noPatterns
        case noResults
        case noSmiles
        case emptyMol
    }
    
    
    override class func setUp() {
        let smi: SMIFormat = SMIFormat() // initialize format so that MKConversion finds its input options
    }

    let subDir = "Data/testdata"
    
    let smarts_file = "smartstest"       // txt
    let smarts_results = "smartsresults" // txt
    let test_smiles = "attype.00"        // smi
    
    func testSMARTSPattern() throws {
        
        guard let smarts_test_patterns_url = Bundle.module.url(forResource: smarts_file, withExtension: "txt", subdirectory: subDir) else { throw SMARTSErrors.noPatterns }
        guard let smarts_results = Bundle.module.url(forResource: smarts_results, withExtension: "txt", subdirectory: subDir) else { throw SMARTSErrors.noResults }
        guard let smiles_url = Bundle.module.url(forResource: test_smiles, withExtension: "smi", subdirectory: subDir) else { throw SMARTSErrors.noSmiles }
        
        var numPatterns: Int = 0
        // get firstLine of smarts_results
        smarts_results.foreachRow { line, rowNum in
            if rowNum == 0 {  
                // split and grab first contents
                let contents = line.components(separatedBy: .whitespaces)
                if let npats = Int(contents[0]) {
                    numPatterns = npats
                }
            }
        }
        
        var vsp: [MKSmartsPattern] = []
        
        // create map of known patterns
        smarts_test_patterns_url.foreachRow { line, rowNum in
            if !line.starts(with: "#") {
                let sp = MKSmartsPattern()
                if sp.initialize(line) {
                    vsp.append(sp)
                }
            }
        }
        
        XCTAssert(numPatterns == vsp.count)
        
        let inputSmilesStream = try InputFileHandler(path: smiles_url, mode: "r")
        let outputStream = OutputStringStream()
        let conv: MKConversion = MKConversion(inputSmilesStream, outputStream)
        
        XCTAssert(conv.setInAndOutFormats("smi", "smi"))
        
        var mol = MKMol()
        
        smiles_url.foreachRow { smile, rowNum in
            if !smile.starts(with: "#") || !smile.isEmpty  {
                mol.clear()
                conv.read(&mol)
//                print(mol.numAtoms())
                if mol.isEmpty() {
                    print(smile)
                }
//                XCTAssert(!mol.isEmpty())
            }
        }
        
    }

//    func testBitVec() throws {
//        
//        var bv1 = MKBitVec(501)
//        var bits: [Int] = [Int].init(repeating: 0, count: bv1.getSize())
//        
//        for j in [0,1,22,33,4,500,6,7] {
//            bv1 |= j
//            bv1.toVecInt(&bits)
//            print(bits)
//            bits.removeAll()
//        }
//        
////        bv1.fromVecInt([0,1,2,3,4,500,6,7])
//        
////        bv1.toVecInt(&bits)
////        print(bits)
////        bits.removeAll()
//                
//        var i = 0
//        i = bv1.nextBit(0)
//        while i != bv1.endBit() {
//            print(i)
//            i = bv1.nextBit(i)
//        }
//        
//        
//        
//    }
    

}
