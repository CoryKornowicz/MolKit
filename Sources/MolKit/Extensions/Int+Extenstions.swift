//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 3/19/23.
//

import Foundation

extension Int {
    static func & (_ lhs: Int, _ rhs: UInt16) -> Int {
        return lhs & Int(rhs)
    }
}


extension UInt16 {
    static func & (_ lhs: UInt16, _ rhs: Int) -> Int {
        return Int(lhs) & rhs
    }
}
