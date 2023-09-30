//
//  SpeedBenchmarks.swift
//  
//
//  Created by Cory Kornowicz on 9/26/23.
//

import XCTest
import MolKit

final class SpeedBenchmarks: XCTestCase {

    // MKMol 1: create/destroy Mol with 10000 atoms
    func testPerformanceMKMol1() throws {
        self.measure {
            let mol = MKMol()
            for _ in 0..<1000 {
                mol.newAtom()
            }
        }
    }
    
    

}
