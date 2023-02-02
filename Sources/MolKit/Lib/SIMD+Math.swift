//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation
import simd
    
extension SIMD3 {
    
    static func vector_angle(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Double {
        
        var dp: Double = simd_dot(simd_normalize(a), simd_normalize(b))

        if (dp < -0.999999) { dp = -0.9999999 }

        if (dp > 0.9999999) { dp = 0.9999999 }
        // Convert from radians to degrees
        return acos(dp) * 180.0 / Double.pi
    }
    
}
    
func isNearZero<U: FloatingPoint>(_ x: U, _ epsilon: U) -> Bool {
    return abs(x) < epsilon
}
