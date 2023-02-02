//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation
import simd

public let OB_SSSR_MOL: UInt = 1<<1
public let OB_RINGFLAGS_MOL: UInt = 1<<2
public let OB_AROMATIC_MOL: UInt = 1<<3
public let OB_ATOMTYPES_MOL: UInt = 1<<4
public let OB_CHIRALITY_MOL: UInt = 1<<5
public let OB_PCHARGE_MOL: UInt = 1<<6
public let OB_HYBRID_MOL: UInt = 1<<8
public let OB_CLOSURE_MOL: UInt = 1<<11
public let OB_H_ADDED_MOL: UInt = 1<<12
public let OB_PH_CORRECTED_MOL: UInt = 1<<13
public let OB_CHAINS_MOL: UInt = 1<<15
public let OB_TCHARGE_MOL: UInt = 1<<16
public let OB_TSPIN_MOL: UInt = 1<<17
public let OB_RINGTYPES_MOL: UInt = 1<<18
public let OB_PATTERN_STRUCTURE: UInt = 1<<19
public let OB_LSSR_MOL: UInt = 1<<20
public let OB_ATOMSPIN_MOL: UInt = 1<<21
public let OB_REACTION_MOL: UInt = 1<<22
public let OB_PERIODIC_MOL: UInt = 1<<23

// Flags 24-32 are unspecified 

public let OB_CURRENT_CONFORMER	= -1

enum HydrogenType {
    case AllHydrogen, PolarHydrogen, NonPolarHydrogen 
}

class MKMol: MKBase {
    
    private var _autoPartialCharge: Bool = false
    private var _autoFormalCharge: Bool = false

    private var _title: String = ""
    private var _vatom: [MKAtom]? = nil
    
    init(atoms: [MKAtom]) {
        self._vatom = atoms
    }
    
    func getAtom(_ id: Int) -> MKAtom {
        guard let atoms = self._vatom else { return MKAtom() }
        return atoms[id]
    }
    
    func getAtoms() -> [MKAtom] {
        return self._vatom!
    }

    func newAtom() -> MKAtom {
        return MKAtom()
    }

    func addBond(_ start_idx: Int, _ end_idx: Int, _ type: Int) {}

    func numAtoms() -> Int {
        return (self._vatom != nil) ? self._vatom!.count : 0
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

    func hasClosureBondsPerceived() -> Bool {
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
