//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation

class MKMol: MKBase {
    
    var atoms: [MKAtom]
    
    init(atoms: [MKAtom]) {
        self.atoms = atoms
    }
    
    func getAtom(_ id: Int) -> MKAtom {
        return self.atoms[id]
    }

    func newAtom() -> MKAtom {
        return MKAtom()
    }

    func addBond(_ start_idx: Int, _ end_idx: Int, _ type: Int) {}

    func numAtoms() -> Int {
        return self.atoms.count
    }

    func automaticPartialCharge() -> Bool {
        return false
    }

    func beginModify() {}
    func endModify() {}
    
    func findChildren(_ idxA: Int, _ idx2: Int) -> [Int] {
        return []
    }
    
    func hasPartialChargesPerceived() -> Bool {
        return false
    }
    
    func hasAtomTypesPerceived() -> Bool {
        return false 
    }
    
    func hasAromaticPerceived() -> Bool {
        return false
    }

    func hasChainsPerceived() -> Bool {
        return false
    }
    
    func hasHybridizationPerceived() -> Bool {
        return false
    }

    func hasRingAtomsAndBondsPerceived() -> Bool {
        return false
    }

    func hasSSSRPerceived() -> Bool {
        return false
    }

    func findRingAtomsAndBonds() {
        
    }

    func findSSSR() {
        
    }

    func isPeriodic() -> Bool {
        return false
    }

    func getTorsion(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom, _ d: MKAtom) -> Double {
        return 0.0
    }

    func getSSSR() -> [MKRing] {
        return []
    }

}
