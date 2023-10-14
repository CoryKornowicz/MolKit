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
        let data = try String(contentsOfFile: smarts_results.relativePath, encoding: .utf8)
        let myStrings = data.components(separatedBy: .newlines)
        
        var myStringsOffset: Int = 0
        
        let contents = myStrings[myStringsOffset++].components(separatedBy: .whitespaces)
        
        if let npats = Int(contents[0]) {
            numPatterns = npats
        }
                
        XCTAssert(numPatterns != 0)
        
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
        
        var vs: [String] = []
        let inputSmilesStream = try InputFileHandler(path: smiles_url, mode: "r")
        let outputStream = OutputStringStream()
        let conv: MKConversion = MKConversion(inputSmilesStream, outputStream)
        
        XCTAssert(conv.setInAndOutFormats("smi", "smi"))
        
        var mol = MKMol()
        var currentMol: Int = 0
        var molPassed: Bool = false
        
        smiles_url.foreachRow { smile, rowNum in
            
            mol.clear()
            conv.read(&mol)
                        
            currentMol += 1
            molPassed = true
            
            for i in vsp {
                if myStringsOffset >= myStrings.count { continue }
                vs = myStrings[myStringsOffset++].components(separatedBy: .whitespaces).filter({ !$0.isEmpty })
                if vs.isEmpty { continue }
                
                i.match(mol)
                let mlist = i.getUMapList()
                
                if mlist.count != vs.count {
                    print("not ok \(currentMol)")
                    print("number of matches different than reference")
                    print("Expected \(vs.count) matches, found \(mlist.count)")
                    print("Error on molecule \(mol.getTitle())")
                    print("SMILES: \(smile)")
                    print("on pattern \(i.getSMARTS())")
                    if !mlist.isEmpty {
                        print("First match: atom# \(mlist[0][0])")
                    } else {
                        print("No matches found")
                    }
                    molPassed = false
                    continue
                }
                    
                if !mlist.isEmpty {
                    for (n, k) in vs.enumerated() {
                        if Int(k) != mlist[n][0] {
                            print("not ok \(currentMol)")
                            print("matching atom numbers different than reference")
                            print("expected \(k) but found \(mlist[n][0])")
                            print("# Molecule: \(mol.getTitle())")
                            print("# Pattern: \(i.getSMARTS())")
                        }
                    }
                }
            }
            
            if molPassed {
                print("ok \(currentMol) molecule passed tests")
            }
            
        }
        
        print("1..\(currentMol) ok")
    }

    

}
