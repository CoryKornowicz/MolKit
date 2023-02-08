//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/7/23.
//

import Foundation


extension MutableCollection {
  mutating func updateEach(_ update: (inout Element) -> Void) {
    for i in indices {
      update(&self[i])
    }
  }
}
