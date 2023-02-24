import XCTest
@testable import MolKit

final class MolKitTests: XCTestCase {
    
//    func testExample() throws {
//        // This is an example of a functional test case.
//        // Use XCTAssert and related functions to verify your tests produce the correct
//        // results.
//        XCTAssert(MolKit != nil)
//    }

    // Test MKAtom Class 

    func testMKAtom() throws {
        let atom = MKAtom()
        XCTAssertNotNil(atom)
    }

    func testMKAtomGetAtomicNum() throws {
        let atom = MKAtom()
        XCTAssertEqual(atom.getAtomicNum(), 0)
    }
    
    
    func testMKElementsGetAtomicNum() throws {
        XCTAssertEqual(MKElements.getAtomicNum("S"), 16)
    }
    
    func testMKElementsGetExactMass() throws {
        XCTAssertEqual(MKElements.getExactMass(5, 10), 10.012937000)
    }
    
    // Test MKGlobalDataBase subclasses
    
    
    func testMKSpaceGroups() throws {
        let spacegroups = MolKit._SpaceGroups
        print(spacegroups.getSize())
        XCTAssert(spacegroups.getSize() > 0)
    }
    
    func testMKAtomicHeatOfFormationTable() throws {
        let atomicHOFTable = MolKit._AtomicHeatOfFormationTable
        print(atomicHOFTable.getSize())
        XCTAssert(atomicHOFTable.getSize() > 0)
    }
    
    func testMKTypeTable() throws {
        let typeTable = MolKit._TypeTable
        print(typeTable.getSize())
        XCTAssert(typeTable.getSize() > 0)
    }
    
    func testMKRingTyper() throws {
        let ringTyper = MolKit._RingTyper
        print(ringTyper.getSize())
        XCTAssert(ringTyper.getSize() == 40)
    }
    
}
