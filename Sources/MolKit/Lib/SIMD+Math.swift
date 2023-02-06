//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation
import Surge
import simd
import Accelerate


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

// MARK: - Cross Product


public func cross3x3(_ lhs: Vector<Float>, _ rhs: Vector<Float>) -> Vector<Float> {
    var res: Vector<Float> = Vector<Float>.init(dimensions: 3, repeatedValue: 0.0)
    res[0] = (lhs[1] * rhs[2]) - (lhs[2] * rhs[1])
    res[1] = (lhs[0] * rhs[2]) - (lhs[2] * rhs[0])
    res[2] = (lhs[0] * rhs[1]) - (lhs[1] * rhs[0])
    return res
}


public func calculateTorsionAngle(_ a: Vector<Float>,
                           _ b: Vector<Float>,
                           _ c: Vector<Float>,
                           _ d: Vector<Float>) -> Float {
    
    let b1: Vector<Float> = a - b
    let b2: Vector<Float> = b - c
    let b3: Vector<Float> = c - d
    
    let rb2: Float = sqrt(dot(b2, b2))
    
    let b2xb3: Vector<Float> = cross3x3(b2, b3)
    let b1xb2: Vector<Float> = cross3x3(b1, b2)
    let torsion = -atan2(dot(rb2 * b1, b2xb3), dot(b1xb2, b2xb3))
    return torsion * (180.0/Float.pi)
}

    
public func cross3x3(_ lhs: Vector<Double>, _ rhs: Vector<Double>) -> Vector<Double> {
    var res: Vector<Double> = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
    res[0] = (lhs[1] * rhs[2]) - (lhs[2] * rhs[1])
    res[1] = (lhs[0] * rhs[2]) - (lhs[2] * rhs[0])
    res[2] = (lhs[0] * rhs[1]) - (lhs[1] * rhs[0])
    return res
}


public func calculateTorsionAngle(_ a: Vector<Double>,
                           _ b: Vector<Double>,
                           _ c: Vector<Double>,
                           _ d: Vector<Double>) -> Double {
    
    let b1: Vector<Double> = a - b
    let b2: Vector<Double> = b - c
    let b3: Vector<Double> = c - d
    
    let rb2: Double = sqrt(dot(b2, b2))
    
    let b2xb3: Vector<Double> = cross3x3(b2, b3)
    let b1xb2: Vector<Double> = cross3x3(b1, b2)
    let torsion = -atan2(dot(rb2 * b1, b2xb3), dot(b1xb2, b2xb3))
    return torsion * (180.0/Double.pi)
}
    
@inline(__always)
func withVector<V: Vectorable>(from vector: Vector<V>, _ closure: (inout Vector<V>) -> ()) -> Vector<V> {
    var copy = vector
    closure(&copy)
    return copy
}

public func normalize(_ lhs: Vector<Float>) -> Vector<Float> {
    let mag = sqrt(sumsq(lhs.scalars))
    return withVector(from: lhs) { $0 /= mag }
}

public func normalize(_ lhs: Vector<Double>) -> Vector<Double> {
    let mag = sqrt(sumsq(lhs.scalars))
    return lhs / mag
}

public func length(_ lhs: Vector<Float>) -> Float {
    return sqrt(sumsq(lhs.scalars))
}

public func length(_ lhs: Vector<Double>) -> Double {
    return sqrt(sumsq(lhs.scalars))
}

public func vector_angle(_ lhs: Vector<Double>, _ rhs: Vector<Double>) -> Double {
    
    var dp: Double = dot(normalize(lhs), normalize(rhs))

    if (dp < -0.999999) { dp = -0.9999999 }

    if (dp > 0.9999999) { dp = 0.9999999 }
    // Convert from radians to degrees
    return acos(dp) * 180.0 / Double.pi
}

public func vector_angle(_ lhs: Vector<Float>, _ rhs: Vector<Float>) -> Float {
    
    var dp: Float = dot(normalize(lhs), normalize(rhs))

    if (dp < -0.999999) { dp = -0.9999999 }

    if (dp > 0.9999999) { dp = 0.9999999 }
    // Convert from radians to degrees
    return acos(dp) * 180.0 / Float.pi
}
