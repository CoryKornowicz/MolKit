//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/24/23.
//

import Foundation

public func isOutsideOrganicSubset(_ elem: Int) -> Bool {
  switch (elem) {
  case  0, // *
      5, // B
      6, // C
      7, // N
      15, // P
      8, // O
      16, // S
      9, // F
      17, // Cl
      35, // Br
      53: // I
    return false
  default:
    return true
  }
}


/* Return the implicit Smiles valence for element "elem" with neutral charge
 * and "val" explicit nbrs. When writing SMILES, a return value of 0 indicates a hypervalent structure
 */

func smilesValence(_ elem: Int, _ val: Int, _ reading: Bool = true) -> Int {
    switch elem {
    case 5: // B
        if val <= 3 { return 3 }
    case 6: // C
        if val <= 4 { return 4 }
    case 7, 15: // N, P
        switch val {
        case 0, 1, 2, 3:
            return 3
        case 4,5:
            return 5
        default: break
        }
    case 8:
        if val <= 2 { return 2 }
    case 16:
        switch val {
        case 0, 1, 2:
            return 2
        case 3, 4:
            return 4
        case 5, 6:
            return 6
        default: break
        }
    case 9,17,35,53:
        if val <= 1 { return 1 }
    default: break
    }

    return reading ? val : 0
}
