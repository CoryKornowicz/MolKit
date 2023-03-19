//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 3/17/23.
//

import Foundation


extension ExpressibleByIntegerLiteral {
    init(_ booleanLiteral: BooleanLiteralType) {
        self = booleanLiteral ? 1 : 0
    }
}

extension BooleanLiteralType {
    init(_ intergerLiteral: any ExpressibleByIntegerLiteral) {
        self = intergerLiteral as! Int > 0
    }
}

extension Bool {
    
    static func ^ (_ lhs: Bool, _ rhs: Bool) -> Bool {
        return Bool(Int(lhs) ^ Int(rhs))
    }
    
    static func | (_ lhs: Bool, _ rhs: Bool) -> Bool {
        return Bool(Int(lhs) | Int(rhs))
    }
    
}
