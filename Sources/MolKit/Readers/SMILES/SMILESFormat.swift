//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/24/23.
//

import Foundation



// This code uses the old OpenEye SMILES parser
// but replaces the SMILES export with Craig James canonical smiles
// (For regular SMILES, the canonical order is not computed and ignored)

let IMPLICIT_CIS_RING_SIZE = 8

// some constant variables
let BondUpChar: Character = "\\"
let BondDownChar: Character = "/"

// This function return true for sulfur and nitrogen
// (I'm not sure that is the right approach, longterm)
// TODO: Oxygen should be able to have lone pair as well, no?
func canHaveLonePair(_ elem: Int) -> Bool {
    switch elem {
    case MKElements.Nitrogen.atomicNum, MKElements.Sulfur.atomicNum:
        return true
    default: return false
    }
}


