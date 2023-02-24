//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/7/23.
//

import Foundation
import Collections


extension MutableCollection {
  mutating func updateEach(_ update: (inout Element) -> Void) {
    for i in indices {
      update(&self[i])
    }
  }
}

enum SubscriptError: Error {
    case outOfBounds
    case greaterThanZero
    case lessThanLastIndex
}

extension Collection where Indices.Iterator.Element == Index {
    
    public subscript(safelyAccess ind: Index) -> Iterator.Element {
        get { return (ind as! Int) < 0 ? self[index(endIndex, offsetBy: (ind as! Int) - 1)] : self[ind]}
        
        
    }

}
