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

public let VZero : Vector<Double> = [0.0, 0.0, 0.0]

// MARK: - Float

extension Float {

    // Degree to Radian
    public var degreesToRadians: Float {
        return self * .pi / 180
    }

    // Radian to Degree

    public var radiansToDegrees: Float {
        return self * 180 / .pi
    }

    // Square function 
    public var square: Float {
        return self * self
    }

}

// MARK: - Double 

extension Double {

    // Degree to Radian
    public var degreesToRadians: Double {
        return self * .pi / 180
    }

    // Radian to Degree

    public var radiansToDegrees: Double {
        return self * 180 / .pi
    }
    
    // Square Function
    public var square: Double {
        return self * self
    }

}

// MARK: Array Extenstions

//extension Array where Element == any Sequence {
//    
////    static func - (_ lhs: Array<Double>, _ rhs: Array<Double>) -> Array<Double> {
////        var arr = Array(repeating: 0.0, count: lhs.count)
////        arr = zip(lhs, rhs).compactMap({ (l, r) in
////            l - r
////        })
////        return arr
////    }
//    
//}

// MARK: - SIMD3

extension SIMD3 {
    
    static func vector_angle(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Double {
        
        var dp: Double = simd_dot(simd_normalize(a), simd_normalize(b))

        if (dp < -0.999999) { dp = -0.9999999 }

        if (dp > 0.9999999) { dp = 0.9999999 }
        // Convert from radians to degrees
        return acos(dp).radiansToDegrees
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


func isApprox(_ val: Vector<Float>, to: Vector<Float>, epsilon: Float) -> Bool {
    return (distSq(val, to) <= ((epsilon * epsilon) * min(length(val), length(to))))
}

func isApprox(_ val: Vector<Double>, to: Vector<Double>, epsilon: Double) -> Bool {
    return (distSq(val, to) <= ((epsilon * epsilon) * min(length(val), length(to))))
}



func isApprox(_ val: Float, to: Float, epsilon: Float) -> Bool {
    return( abs(val - to) <= epsilon * min(abs(val), abs(to)))
}

func isApprox(_ val: Double, to: Double, epsilon: Double) -> Bool {
    return( abs(val - to) <= epsilon * min(abs(val), abs(to)))
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
    return torsion.radiansToDegrees
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
    return torsion.radiansToDegrees
}
    

// MARK: Inlineable function blocks for Vector and Matrix that mutate the reference instead of makeing a copy

@inline(__always)
func withVector<V: Vectorable>(from vector: Vector<V>, _ closure: (inout Vector<V>) -> ()) -> Vector<V> {
    var copy = vector
    closure(&copy)
    return copy
}


@inline(__always)
func withMatrix<Scalar>(from matrix: Matrix<Scalar>, _ closure: (inout Matrix<Scalar>) -> ()) -> Matrix<Scalar> {
    var copy = matrix
    closure(&copy)
    return copy
}

// MARK: Vector Extensions

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
    return acos(dp).radiansToDegrees
}

public func vector_angle(_ lhs: Vector<Float>, _ rhs: Vector<Float>) -> Float {
    
    var dp: Float = dot(normalize(lhs), normalize(rhs))

    if (dp < -0.999999) { dp = -0.9999999 }

    if (dp > 0.9999999) { dp = 0.9999999 }
    // Convert from radians to degrees
    return acos(dp).radiansToDegrees
}

// MARK: Matrix Extensions

extension Matrix where Scalar == Double {
    
    mutating func fillOrtho(alpha: Double, beta: Double, gamma: Double, A: Double, B: Double, C: Double) {
        
        let alpha = alpha.degreesToRadians
        let beta = beta.degreesToRadians
        let gamma = gamma.degreesToRadians
        let cosAlpha = cos(alpha)
        let cosBeta = cos(beta)
        let cosGamma = cos(gamma)
        let sinGamma = sin(gamma)
        
        // from the PDB specification:
        // http://www.rcsb.org/pdb/docs/format/pdbguide2.2/part_75.html
        // since we'll ultimately divide by (a * b), we've factored those out here
        let V = C * sqrt(1 - cosAlpha.square - cosBeta.square - cosGamma.square + 2 * cosAlpha * cosBeta * cosGamma)
        
        self[0, 0] = A
        self[0, 1] = B * cosGamma
        self[0, 2] = C * cosBeta

        self[1, 0] = 0
        self[1, 1] = B * sinGamma
        self[1, 2] = C * (cosAlpha - cosBeta * cosGamma) / sinGamma

        self[2, 0] = 0
        self[2, 1] = 0
        self[2, 2] = V / sinGamma
    }
    
}

extension Matrix where Scalar == Float {
    
    mutating func fillOrtho(alpha: Float, beta: Float, gamma: Float, A: Float, B: Float, C: Float) {
        
        let alpha = alpha.degreesToRadians
        let beta = beta.degreesToRadians
        let gamma = gamma.degreesToRadians
        let cosAlpha = cos(alpha)
        let cosBeta = cos(beta)
        let cosGamma = cos(gamma)
        let sinGamma = sin(gamma)
        
        // from the PDB specification:
        // http://www.rcsb.org/pdb/docs/format/pdbguide2.2/part_75.html
        // since we'll ultimately divide by (a * b), we've factored those out here
        let V = C * sqrt(1 - cosAlpha.square - cosBeta.square - cosGamma.square + 2 * cosAlpha * cosBeta * cosGamma)
        
        self[0, 0] = A
        self[0, 1] = B * cosGamma
        self[0, 2] = C * cosBeta

        self[1, 0] = 0
        self[1, 1] = B * sinGamma
        self[1, 2] = C * (cosAlpha - cosBeta * cosGamma) / sinGamma

        self[2, 0] = 0
        self[2, 1] = 0
        self[2, 2] = V / sinGamma
    }
    
}

// MARK: Extension Vector

extension Vector {
    
    var x: Element {
        scalars[0]
    }
    
    var y: Element {
        scalars[1]
    }
    
    var z: Element {
        scalars[2]
    }
    
}
