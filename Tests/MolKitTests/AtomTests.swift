//
//  AtomTests.swift
//  
//
//  Created by Cory Kornowicz on 9/19/23.
//

import XCTest
import MolKit

final class AtomTests: XCTestCase {

    func testSetIdx() throws {
        // OBAtom isolation tests (no connection to residue, bond, molecule...)
        // Beat up on SetIdx
        let testAtom1 = MKAtom()
        
        testAtom1.setIdx(0)
        print(testAtom1.getIdx())
        XCTAssert(testAtom1.getIdx() == 0)
        testAtom1.setIdx(-1)
        print(testAtom1.getIdx())
        XCTAssert(testAtom1.getIdx() == -1)
        testAtom1.setIdx(1)
        print(testAtom1.getIdx())
        XCTAssert(testAtom1.getIdx() == 1)
    }
    
    func testSetAtomicNum() throws {
        // OBAtom isolation tests (no connection to residue, bond, molecule...)
        // Beat up on SetIdx
        let testAtom1 = MKAtom()
        
        testAtom1.setAtomicNum(0)
        print(testAtom1.getAtomicNum())
        XCTAssert(testAtom1.getAtomicNum() == 0)
        testAtom1.setAtomicNum(-1)
        print(testAtom1.getAtomicNum())
        XCTAssert(testAtom1.getAtomicNum() == -1)
        testAtom1.setAtomicNum(200)
        print(testAtom1.getAtomicNum())
        XCTAssert(testAtom1.getAtomicNum() == 200)
        testAtom1.setAtomicNum(300)
        print(testAtom1.getAtomicNum())
        XCTAssert(testAtom1.getAtomicNum() == 300)
        testAtom1.setAtomicNum(1)
        print(testAtom1.getAtomicNum())
        XCTAssert(testAtom1.getAtomicNum() == 1)
        
        testAtom1.getCoordinate()
    }

//    func testPerformanceExample() throws {
//        // This is an example of a performance test case.
//        self.measure {
//            // Put the code you want to measure the time of here.
//        }
//    }

}
