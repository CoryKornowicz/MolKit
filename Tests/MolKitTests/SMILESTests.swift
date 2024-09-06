//
//  SmilesTest.swift
//  
//
//  Created by Cory Kornowicz on 9/29/23.
//

import XCTest
import MolKit

final class SmilesTest: XCTestCase {
    
    enum Errors: Error {
        case noStereoData
        case incorrectTetrahedralStereo
    }
    
    override class func setUp() {
        let smi: SMIFormat = SMIFormat() // initialize format so that MKConversion finds its input options
    }
    
    func testTetrahedralStereo() throws {
        
        var mol: MKMol = MKMol()
        let conv: MKConversion = MKConversion()
        XCTAssert(conv.setInFormat("smi"))
        print("Parsing smiles: C[C@H](O)N")
        XCTAssert(conv.readString(&mol, "C[C@H](O)N"))
        
        // Get the stereo data
        XCTAssert(mol.hasData(.StereoData))
        guard let stereoData = mol.getAllData(.StereoData) else {
            throw Errors.noStereoData
        }
        XCTAssert(stereoData.count == 1)
        
        // convert to tetrahedral data
        guard let ts = (stereoData[0] as? MKTetrahedralStereo) else {
            throw Errors.incorrectTetrahedralStereo
        }
        XCTAssert((ts as MKStereoBase).getType() == MKStereo.TType.Tetrahedral)
        
        print(ts.getConfig())
        
        // construct a valid configuration here
        //
        // C[C@H](O)N
        // 0 1 2  3 4  <- ids
        //
        
        let cfg: MKTetrahedralStereo.Config = MKTetrahedralStereo.Config(center: Ref(integerLiteral: 1), 
                                                                         from_or_towrds: .from(Ref(integerLiteral: 0)),
                                                                         winding: .Clockwise,
                                                                         view: .ViewFrom, specified: true, refs: MKStereo.makeRefs(4, 3, 2))
        print(cfg)
        // compare stereochemistry
        XCTAssert(ts.getConfig() == cfg)
    }
    
    func testComplexSMILES() throws {
        
        var mol: MKMol = MKMol()
        let conv: MKConversion = MKConversion()
        XCTAssert(conv.setInFormat("smi"))
//        let smi = "CC(=C)[C@@H]1CCC2(CC[C@]3(C)C(CCC4[C@@]5(C)CCC(=O)C(C)(C)C5CC[C@]43C)C12)N=C=O"
//        let smi = "CO[C@@H]1C[C@@H](C[C@H]2CC[C@H](C)[C@@H](O2)[C@@H](C)C(=O)O)O[C@]3(O[C@@](C)(C[C@H]3C)[C@@H]4CC[C@@](C)(O4)[C@H]5O[C@@H](C[C@H]5C)[C@@H]6O[C@](O)(CO)[C@@H](C)C[C@H]6C)[C@@H]1C"
        let smi = "O=C1C(=C[C@@H]2[C@@H](O1)CCCC2)C1CCCCC1"
        print("Parsing smiles: \(smi)")
        
        XCTAssert(conv.readString(&mol, smi))
//        mol.percieveBondOrders()
        // Get the stereo data
        let data = mol.getAllData(.StereoData)
        
        for dat in data! {
            print((dat as! MKTetrahedralStereo).getConfig())
        }
        
        for atom in mol.getAtomIterator() {
            print("Atom \(atom.getAtomicNum()) \(atom.getType()) ")
        }
        
        // PATTERN
        
        let pattern = MKSmartsPattern()
        pattern.initialize("[C@@H]")
        
        pattern.match(mol) 
        
        for i in pattern.getUMapList() {
            print(i)
        }
                
    }
    

}
