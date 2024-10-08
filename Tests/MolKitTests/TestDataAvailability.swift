//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/4/23.
//

import XCTest
@testable import MolKit

final class DataAvailabilityTests: XCTestCase {
    
    public func testFindSpaceGroup() throws {
        
        // Define bundle path to Data folder
        let path = Bundle.module.url(forResource: "space-groups", withExtension: "txt", subdirectory: "Data")
        let data = try String(contentsOf: path!)
                
        XCTAssert(!data.isEmpty)
    }
    
    func testFindAtomHOF() throws {
        
        // Define bundle path to Data folder
        let path = Bundle.module.url(forResource: "atomization-energies", withExtension: "txt", subdirectory: "Data")
        let data = try String(contentsOf: path!)
        
        XCTAssert(!data.isEmpty)
    }
    
    
    func testFindTypeTable() throws {
        
        // Define bundle path to Data folder
        let path = Bundle.module.url(forResource: "types", withExtension: "txt", subdirectory: "Data")
        let data = try String(contentsOf: path!)
        
        XCTAssert(!data.isEmpty)
    }
    
    func testFindAllDataFiles() throws {
        
        // Define bundle path to Data folder
        let path = Bundle.module.paths(forResourcesOfType: "txt", inDirectory: "Data")
        XCTAssert(path.count > 0)
    }
    
}
