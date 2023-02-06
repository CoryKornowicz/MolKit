//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/6/23.
//

import Foundation

// Make Typing SOOOOO Much Easier and use this Vectorable Type instead

protocol Vectorable: FloatingPoint & ExpressibleByFloatLiteral {}

extension Float: Vectorable {}
extension Double: Vectorable {}

// Harder than anticipated -- Defaulting to using Double everywhere 
