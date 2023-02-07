//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation
import Surge
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

enum HydrogenType: Int {
    case AllHydrogen = 0
    case PolarHydrogen = 1
    case NonPolarHydrogen = 2
}

class MKMol: MKBase {
    
    private var _autoPartialCharge: Bool = false
    private var _autoFormalCharge: Bool = false

    private var _title: String = ""
    private var _vatom: [MKAtom]? = []
    private var _vatomIds: [MKAtom]? = []
    private var _vbond: [MKBond]? = []
    private var _vbondIds: [MKBond]? = []
    private var _dimension: UInt16 = 0
    private var _totalCharge: Int = 0
    private var _totalSpin: UInt = 0
    private var _c: Array<SIMD3<Double>> = []
    private var _vconf: Array<SIMD3<Double>> = []
    private var _energy: Double = 0.0
    private var _natoms: UInt = 0
    private var _nbonds: UInt = 0
    private var _residue: [MKResidue]? = []
    private var _internals: [MKInternalCoord]? = []
    private var _mod: UInt16 = 0 

    
    init(atoms: [MKAtom]) {
        self._vatom = atoms
    }
    
    func getTitle() -> String {
        return self._title
    }

    func setTitle(_ title: String) {
        self._title = title
    }

    func getAtom(_ id: Int) -> MKAtom {
        guard let atoms = self._vatom else { return MKAtom() }
        return atoms[id]
    }
    
    func getAtoms() -> [MKAtom] {
        return self._vatom!
    }

    func getAtomIterator() -> MKIterator<MKAtom>? {
        guard let atoms = self._vatom else { return nil }
        return MKIterator<MKAtom>(atoms)
    }

    func newAtom() -> MKAtom {
        return MKAtom()
    }
    
    func deleteAtom(_ atom: MKAtom) {
        
    }

    func numAtoms() -> Int {
        return (self._vatom != nil) ? self._vatom!.count : 0
    }

    func addBond(_ start_idx: Int, _ end_idx: Int, _ type: Int) {}


    func automaticPartialCharge() -> Bool {
        return false
    }

    func beginModify() {}
    func endModify() {}
    
    func findChildren(_ idxA: Int, _ idx2: Int) -> [Int] {
        return []
    }
    
    func has2D(_ not3D: Bool = false) -> Bool {
        return false
    }
    
    func has3D() -> Bool {
        return false
    }
    
    func hasNonZeroCoords() -> Bool {
        return false
    }
    
    func hasAromaticPerceived() -> Bool {
        return self.hasFlag(OB_AROMATIC_MOL)
    }
    
    func hasSSSRPerceived() -> Bool {
        return self.hasFlag(OB_SSSR_MOL)
    }
    
    func hasLSSRPerceived() -> Bool {
        return self.hasFlag(OB_LSSR_MOL)
    }
    
    func hasRingAtomsAndBondsPerceived() -> Bool {
        return self.hasFlag(OB_RINGFLAGS_MOL)
    }
    
    func hasAtomTypesPerceived() -> Bool {
        return self.hasFlag(OB_ATOMTYPES_MOL)
    }
    
    func hasRingTypesPerceived() -> Bool {
        return self.hasFlag(OB_RINGTYPES_MOL)
    }
    
    func hasChiralityPerceived() -> Bool {
        return self.hasFlag(OB_CHIRALITY_MOL)
    }
    
    func hasPartialChargesPerceived() -> Bool {
        return self.hasFlag(OB_PCHARGE_MOL)
    }
    
    func hasHybridizationPerceived() -> Bool {
        return self.hasFlag(OB_HYBRID_MOL)
    }
    
    func hasClosureBondsPerceived() -> Bool {
        return self.hasFlag(OB_CLOSURE_MOL)
    }

    func hasChainsPerceived() -> Bool {
        return self.hasFlag(OB_CHAINS_MOL)
    }
    
    func hasHydrogensAdded() -> Bool {
        return self.hasFlag(OB_H_ADDED_MOL)
    }
    
    func isReaction() -> Bool {
        return self.hasFlag(OB_REACTION_MOL)
    }
    
    func isEmpty() -> Bool {
        return self._natoms == 0
    }
    
    
    

    func findRingAtomsAndBonds() {
        
    }

    func findSSSR() {
        
    }

    func isPeriodic() -> Bool {
        return false
    }

    func getAngle(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom) -> Double {
        return a.getAngle(b, c)
    }

    func getTorsion(_ a: Int, _ b: Int, _ c: Int, _ d: Int) -> Double {
        guard let i = self._vatom?[a-1], let j = self._vatom?[b-1], let k = self._vatom?[c-1], let m = self._vatom?[d-1] else { return 0.0 }
        return self.getTorsion(i,j,k,m)
    }
    
    func getTorsion(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom, _ d: MKAtom) -> Double {
        
        if !self.isPeriodic() {
            return calculateTorsionAngle(a.getVector(), b.getVector(), c.getVector(), d.getVector())
        }else {
            return 0.0
        }
        
    }

    func getSSSR() -> [MKRing]? {
        return nil
    }

    func setPeriodicMol(_ value: Bool = true) {
        if value {
            self.setFlag(OB_PERIODIC_MOL)
        } else {
            self.unsetFlag(OB_PERIODIC_MOL)
        }
    }
    
    
}
