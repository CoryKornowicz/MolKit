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
    
    func automaticPartialCharge() -> Bool {
        return false
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

    func findRingAtomsAndBonds() {
        
    }

    func getTorsion(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom, _ d: MKAtom) -> Double {
        return 0.0
    }

}
