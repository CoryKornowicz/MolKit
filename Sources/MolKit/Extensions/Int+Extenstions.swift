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

// Add support for incrementing underlying variable after value return
// e.g.
// var j: Int = 0
// print(j++) // prints 0, THEN increments j by 1
// print(j)   // prints 1

postfix operator ++
postfix operator --

@discardableResult 
@inline(__always)
postfix func ++ <A: FixedWidthInteger>(val: inout A) -> A {
    defer {
        val += 1
    }
    return val
}

@discardableResult
@inline(__always)
postfix func -- <A: FixedWidthInteger>(val: inout A) -> A {
    defer {
        val -= 1
    }
    return val
}


// Add support for incrementing underlying variable before value return
// e.g.
// var j: Int = 0
// print(--j) // decrements j by -1, THEN returns it
// print(j)   // prints -1

prefix operator ++
prefix operator --

@discardableResult
@inline(__always)
prefix func ++ <A: FixedWidthInteger>(val: inout A) -> A {
    val += 1
    return val
}

@discardableResult
@inline(__always)
prefix func -- <A: FixedWidthInteger>(val: inout A) -> A {
    val -= 1
    return val
}

