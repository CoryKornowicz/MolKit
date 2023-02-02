//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 1/31/23.
//

import Foundation

class MKAtom: MKBase {
    
    var _parent: MKBase?
    var _residue: MKResidue?
    
    
    public override init() {
        super.init()
        self._parent = nil
        self.clear()
    }
    

    private func clear() {
    
        
    }
    
    deinit {
        if let res = self._residue {
            _ = res.removeAtom(self)
        }
    }
    
}
