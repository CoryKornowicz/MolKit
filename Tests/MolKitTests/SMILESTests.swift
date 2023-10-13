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
        
        let cfg: MKTetrahedralStereo.Config = MKTetrahedralStereo.Config(center: Ref(integerLiteral: 1), from_or_towrds: .from(Ref(integerLiteral: 0)), winding: .Clockwise, view: .ViewFrom, specified: true, refs: MKStereo.makeRefs(4, 3, 2))
        print(cfg)
        // compare stereochemistry
        XCTAssert(ts.getConfig() == cfg)
    }
    
    func testComplexSMILES() throws {
        
        var mol: MKMol = MKMol()
        let conv: MKConversion = MKConversion()
        XCTAssert(conv.setInFormat("smi"))
        let smi =  "OP(=O)(=O)c1[se]c2ccccc2c1Br" //"OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1"
        print("Parsing smiles: \(smi)")
        
        XCTAssert(conv.readString(&mol, smi))
        
        // Get the stereo data
        /*XCTAssert*/(mol.hasData(.StereoData))
        
        for atom in mol.getAtomIterator() {
            print("Atom \(atom.getAtomicNum()) \(atom.getType()) ")
        }
        
        
    }
    

}
