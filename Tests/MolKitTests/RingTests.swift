//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 10/13/23.
//

import XCTest
@testable import MolKit

final class RingTest: XCTestCase {
    
    let subDir = "Data/testdata"
    
    let ring_results = "ringresults"     // txt
    let test_smiles = "attype.00"        // smi
    
    enum RingErrors: Error {
        case noSmiles
        case noResults
        case noMoreLines
        case integerConversion
    }
    
    override class func setUp() {
        let smi: SMIFormat = SMIFormat() // adds it to MKPlugins
    }
    
    func testRings() throws {
        
        guard let ring_results = Bundle.module.url(forResource: ring_results, withExtension: "txt", subdirectory: subDir) else { throw RingErrors.noResults }
        guard let smiles_url = Bundle.module.url(forResource: test_smiles, withExtension: "smi", subdirectory: subDir) else { throw RingErrors.noSmiles }
                                
        let inputSmilesStream = try InputFileHandler(path: smiles_url, mode: "r")
        let outputStream = OutputStringStream()
        let conv: MKConversion = MKConversion(inputSmilesStream, outputStream)
        
        let data = try String(contentsOfFile: ring_results.relativePath, encoding: .utf8)
        let myStrings = data.components(separatedBy: .newlines)
        
        var myStringsOffset: Int = 0
                
        XCTAssert(conv.setInAndOutFormats("smi", "smi"))
        
        var currentTest: Int = 0
        var buffer: [String] = []
        
        try smiles_url.foreachRow { smile, rowNum in
            
            var mol = MKMol()
//            var molPassed: Bool = false
            var vb: [Bool] = []
//            var vi: [Int] = []
            var vr: [MKRing] = []

            conv.read(&mol)
            
            if myStringsOffset >= myStrings.count { throw RingErrors.noMoreLines }
            buffer = myStrings[myStringsOffset++].components(separatedBy: .whitespaces).filter({ !$0.isEmpty })
            
            vb.resize(mol.numBonds(), with: false)
            
            //MARK: Check ring bond data if buffer is not empty
            
            if !buffer.isEmpty {
                for i in buffer {
                    //convert i to integer and store at
                    guard let index = Int(i) else { throw RingErrors.integerConversion }
                    vb[index] = true
                }
            }
            
            for bond in mol.getBondIterator() {
                if vb[Int(bond.getIdx())] != bond.isInRing() {
                    print("no ok \(++currentTest)")
                    print("# ring bond data different than reference")
                    print("$ Molecule: \(mol.getTitle())")
                } else {
                    print("ok \(++currentTest)")
                }
            }
            
            //MARK: Check SSSR size
            
            vr = mol.getSSSR()
            
            if myStringsOffset >= myStrings.count { throw RingErrors.noMoreLines }
            buffer = myStrings[myStringsOffset++].components(separatedBy: .whitespaces).filter({ !$0.isEmpty })
            
            guard let vr_size = Int(buffer[0]) else { throw RingErrors.integerConversion }
            
            if vr.count != vr_size {
                print("not ok \(++currentTest)")
                print("# SSSR size different than reference")
                print("# Molecule: \(mol.getTitle())")
            } else {
                print("ok \(++currentTest)")
            }
            
            //MARK: Check atom ids match
            
            if myStringsOffset >= myStrings.count { throw RingErrors.noMoreLines }
            buffer = myStrings[myStringsOffset++].components(separatedBy: .whitespaces).filter({ !$0.isEmpty })
                        
            var vsIterator = MKIterator<String>(buffer)
            var i = vsIterator.next()
            
            for atom in mol.getAtomIterator() {
                
                if i == nil {
                    print("no ok \(++currentTest)")
                    print("# error in SSSR count")
                    print("# Molecule: \(mol.getTitle())")
                } else {
                    print("ok \(++currentTest)")
                }
                
                var count = 0
                for m in vr {
                    if m._pathSet[atom.getIdx()] {
                        count += 1
                    }
                }
                
                guard i != nil, let count_size = Int(i!) else {
                    continue
                }
                
                if count_size != count {
                    print("no ok \(++currentTest) # ring memebership test failed")
                    print("# Molecule: \(mol.getTitle())")
                } else {
                    print("ok \(++currentTest)")
                }
                
                i = vsIterator.next()
            }
        }
        print("executed \(currentTest) tests")
    }
}
