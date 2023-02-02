import XCTest
@testable import MolKit

final class MolKitTests: XCTestCase {
    
    func testExample() throws {
        // This is an example of a functional test case.
        // Use XCTAssert and related functions to verify your tests produce the correct
        // results.
        XCTAssertEqual(MolKit().text, "Hello, World!")
    }

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
    
}
