//
//  MKGastChrg.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation

public let OB_GASTEIGER_DENOM: Double = 20.02
public let OB_GASTEIGER_DAMP: Double =  0.5
public let OB_GASTEIGER_ITERS: Double = 6

struct MKGasteigerState {
    var a: Double, b: Double, c: Double, denom: Double, chi: Double, q: Double = 0.0
    
    init() {
        self.a = 0.0
        self.b = 0.0
        self.c = 0.0
        self.denom = 0.0
        self.chi = 0.0
        self.q = 0.0
    }

    init(a: Double, b: Double, c: Double, denom: Double, chi: Double, q: Double) {
        self.a = a
        self.b = b
        self.c = c
        self.denom = denom
        self.chi = chi
        self.q = q
    }
}

class MKGastChrg {
    
    var _gsv: Array<MKGasteigerState> = []
    
    init() {}

    private func initialPartialCharges(_ mol: MKMol) { }

    private func gasteigerSigmaChi(_ atom: MKAtom, _ a: Double ...) {
        
    }

    func gsvResize(_ int: Int) { 
        self._gsv = Array<MKGasteigerState>(repeating: MKGasteigerState(), count: int)
    }

    func assignPartialCharges(_ mol: MKMol) {
        
    }

    deinit {
        self._gsv.removeAll()
    }
    
}
