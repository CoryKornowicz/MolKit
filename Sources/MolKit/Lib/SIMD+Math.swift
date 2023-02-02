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
    
    // Generate a random unit vector uniformly distributed over the unit sphere
    static func randomUnitVector() -> SIMD3<Double> {
        let x: Double = Double.random(in: -1.0...1.0)
        let y: Double = Double.random(in: -1.0...1.0)
        let z: Double = Double.random(in: -1.0...1.0)
        let v: SIMD3<Double> = SIMD3<Double>(x, y, z)
        return simd_normalize(v)
    }
}
    
func isNearZero<U: FloatingPoint>(_ x: U, _ epsilon: U) -> Bool {
    return abs(x) < epsilon
}
