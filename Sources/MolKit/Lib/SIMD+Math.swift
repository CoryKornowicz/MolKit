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
public let VX    : Vector<Double> = [1.0, 0.0, 0.0]
public let VY    : Vector<Double> = [0.0, 1.0, 0.0]
public let VZ    : Vector<Double> = [0.0, 0.0, 1.0]


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

public func square(_ x: Double) -> Double {
    return x * x
}

public func square(_ x: Float) -> Float {
    return x * x
}

public func length_2(_ vec: Vector<Float>) -> Float {
    return vec.map { val in square(val) }.reduce(0, +)
}

public func length_2(_ vec: Vector<Double>) -> Double {
    return vec.map { val in square(val) }.reduce(0, +)
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

func isNearZero(_ x: Float) -> Bool {
    isNearZero(x, 1.0e-6)
}

func isNearZero(_ x: Double) -> Bool {
    isNearZero(x, 1.0e-6)
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


func isNegligible(_ val: Double, _ to: Double, precision: Double = 1e-11) -> Bool {
    return (abs(val) <= precision * abs(to))
}

func isNegligible(_ val: Float, _ to: Float, precision: Float = 1e-11) -> Bool {
    return (abs(val) <= precision * abs(to))
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
    
    mutating func rotAboutAxisByAngle(_ v: Vector<Double>, _ angle: Double) {
        let theta = angle.degreesToRadians
        let s = sin(theta)
        let c = cos(theta)
        let t = 1 - c
        
        let vtmp = normalize(v)
        
        let x = vtmp.x
        let y = vtmp.y
        let z = vtmp.z
        
        self[0, 0] = t * x * x + c
        self[0, 1] = t * x * y + s * z
        self[0, 2] = t * x * z - s * y

        self[1, 0] = t * x * y - s * z
        self[1, 1] = t * y * y + c
        self[1, 2] = t * y * z + s * x

        self[2, 0] = t * x * z + s * y
        self[2, 1] = t * y * z - s * x
        self[2, 2] = t * z * z + c
    }
    
    mutating func setupRotMat(_ phi: Double, _ theta: Double, _ psi: Double) {
        let p = phi.degreesToRadians
        let h = theta.degreesToRadians
        let b = psi.degreesToRadians
        
        let cx = cos(p)
        let sx = sin(p)
        let cy = cos(h)
        let sy = sin(h)
        let cz = cos(b)
        let sz = sin(b)

        self[0, 0] = cy * cz
        self[0, 1] = cy * sz
        self[0, 2] = -sy

        self[1, 0] = sx * sy * cz - cx * sz
        self[1, 1] = sx * sy * sz + cx * cz
        self[1, 2] = sx * cy

        self[2, 0] = cx * sy * cz + sx * sz
        self[2, 1] = cx * sy * sz - sx * cz
        self[2, 2] = cx * cy
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

    mutating func rotAboutAxisByAngle(_ v: Vector<Float>, _ angle: Float) {
        let theta = angle.degreesToRadians
        let s = sin(theta)
        let c = cos(theta)
        let t = 1 - c
        
        let vtmp = normalize(v)
        
        let x = vtmp.x
        let y = vtmp.y
        let z = vtmp.z
        
        self[0, 0] = t * x * x + c
        self[0, 1] = t * x * y + s * z
        self[0, 2] = t * x * z - s * y

        self[1, 0] = t * x * y - s * z
        self[1, 1] = t * y * y + c
        self[1, 2] = t * y * z + s * x

        self[2, 0] = t * x * z + s * y
        self[2, 1] = t * y * z - s * x
        self[2, 2] = t * z * z + c
    }
    
    mutating func setupRotMat(_ phi: Float, _ theta: Float, _ psi: Float) {
        let p = phi.degreesToRadians
        let h = theta.degreesToRadians
        let b = psi.degreesToRadians
        
        let cx = cos(p)
        let sx = sin(p)
        let cy = cos(h)
        let sy = sin(h)
        let cz = cos(b)
        let sz = sin(b)

        self[0, 0] = cy * cz
        self[0, 1] = cy * sz
        self[0, 2] = -sy

        self[1, 0] = sx * sy * cz - cx * sz
        self[1, 1] = sx * sy * sz + cx * cz
        self[1, 2] = sx * cy

        self[2, 0] = cx * sy * cz + sx * sz
        self[2, 1] = cx * sy * sz - sx * cz
        self[2, 2] = cx * cy
    }
    
}

extension Matrix {
    
    init(vecs: [Vector<Scalar>]) {
        let rows = vecs.count
        var cols: Int = 0
        vecs.forEach { vec in
            if vec.count > cols {
                cols = vec.count
            }
        }
        self.init(rows: rows, columns: cols) { row, col in
            vecs[row][col]
        }
    }
    
    
    var diagonal: Vector<Scalar> {
        let minStride = Swift.min(self.rows, self.columns)
        let res = Vector.init(dimensions: minStride) { idx in
            self[idx, idx]
        }
        return res
    }
    
}

public let MAX_SWEEPS = 50

public func mk_make_rmat<Scalar: ExpressibleByFloatLiteral & FloatingPoint>(_ a: inout Matrix<Scalar>, _ rmat: inout Matrix<Scalar>) {
    /*
  	onorm, dnorm - hold the sum of diagonals and off diagonals to check Jacobi completion
  	d[3] - holds the diagonals of the input vector, which transofrm to become the Eigenvalues
  	r1, r2 - hold 1st two Eigenvectors
  	v1,v2,v3 - hold orthogonal unit vectors derived from Eigenvectors
  	
  	The junction performs a Jacobi Eigenvalue/vector determination 
  	(https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm) on the supplied
  	Inertial Tensor in a, and returns a unit transform matrix rmat as a row matrix.
  	To work, a must be diagonally symmetric (i.e a[i][j] = a[j][i])
  	v starts out holding the unit matrix (i.e. no transform in co-ordinate frame), 
  	and undergoes the same rotations as applied to the Inertial Tensor in the Jacobi
  	process to arrive at the new co-ordinate frame.
  	Finally, the eigenvalues are sorted in order that the largest principal moment aligns to the 
  	new x-axis
  	*/

    var onorm, dnorm: Scalar
    var b, dma, q, t, c, s: Scalar
    var d: Vector<Scalar> = a.diagonal
    var atemp, vtemp, dtemp: Scalar
    var v: Matrix<Scalar> = Matrix.eye(rows: 3, columns: 3)
    var r1: Vector<Scalar> = Vector.init([0.0,0.0,0.0])
    var r2: Vector<Scalar> = Vector.init([0.0,0.0,0.0])
    var v1: Vector<Scalar> = Vector.init([0.0,0.0,0.0])
    var v2: Vector<Scalar> = Vector.init([0.0,0.0,0.0])
    var v3: Vector<Scalar> = Vector.init([0.0,0.0,0.0])
    var k: Int
    
sweep_loop: for _ in 1...MAX_SWEEPS {
        dnorm = 0.0
        onorm = 0.0
        for j in 0..<3 {
            dnorm += abs(d[j])
            // MARK: Potential error here if the comparator is backwards
            for i in 0...j-1 {
                onorm += abs(a[i, j])
            }
        }
//        Test for convergence
        if (onorm/dnorm <= 1.0e-12) {
            break sweep_loop
        }
        
        for j in 1..<3 {
            for i in 0...j-1 {
                b = a[i, j]
                if abs(b) > 0.0 {
                    dma = d[j] - d[i]
                    if ((abs(dma) + abs(b)) <= abs(dma)) {
                        t = b / dma
                    } else {
                        q = 0.5 * dma / b
                        t = 1.0 / abs(q) + sqrt(1.0+q*q)
                        if q < 0.0 {
                            t = -t
                        }
                    }
                    c = 1.0/sqrt(t * t + 1.0)
                    s = t * c
                    a[i, j] = 0.0
                    /* Perform a Jacobi rotation on the supplied matrix*/
                    for k in 0...i-1 {
                        atemp = c * a[k, i] - s * a[k, j]
                        a[k, j] = s * a[k, i] + c * a[k, j]
                        a[k, i] = atemp 
                    }
                    for k in i+1...j-1 {
                        atemp = c * a[i, k] - s * a[k, j]
                        a[k, j] = s * a[i, k] + c * a[k, j]
                        a[i, k] = atemp
                    }
                    for k in j+1..<3 {
                        atemp = c * a[i, k] - s * a[j, k]
                        a[j, k] = s * a[i, k] + c * a[j, k]
                        a[i, k] = atemp
                    }
                    /* Rotate the reference frame */
                    for k in 0..<3 {
                        vtemp = c * v[k, i] - s * v[k, j]
                        v[k, j] = s * v[k, i] + c * v[k, j]
                        v[k, i] = vtemp
                    }
                    let ccdi = (c * c * d[i])
                    let ccdj = (c * c * d[j])
                    let ssdi = (s * s * d[i])
                    let ssdj = (s * s * d[j])
                    dtemp = ccdi + ssdj - (2.0 * c * s * b)
                    d[j] = ssdi + ccdj + (2.0 * c * s * b)
                    d[i] = dtemp
                }// end if
            } // end for i
        }// end for j
    }// end for l
    
    /* Now sort the eigenvalues and eigenvectors*/

    for j in 0..<3-1 {
        k = j
        dtemp = d[k]
        for i in j+1..<3 {
            if d[i] < dtemp {
                k = i
                dtemp = d[k]
            }
        }

        if k > j {
            d[k] = d[j]
            d[j] = dtemp
            for i in 0..<3 {
                dtemp = v[i, k]
                v[i, k] = v[i, j]
                v[i, j] = dtemp
            }
        }
    }

    /* Transfer the 1st two eigenvectors into r1 and r2*/
    r1[0] = v[0, 0]
    r1[1] = v[1, 0]
    r1[2] = v[2, 0]
    r2[0] = v[0, 1]
    r2[1] = v[1, 1]
    r2[2] = v[2, 1]
    
    /* Generate the 3rd unit vector for the new coordinate frame by cross product of r1 and r2*/
    v3[0] =  (r1[1] * r2[2])
    v3[0] -= (r1[2] * r2[1])
    v3[1] = (-r1[0] * r2[2])
    v3[1] += (r1[2] * r2[0])
    v3[2] =  (r1[0] * r2[1])
    v3[2] -= (r1[1] * r2[0])
    
    /* Ensure it is normalised |v3|=1 */
    let v302 = v3[0]*v3[0]
    let v312 = v3[1]*v3[1]
    let v322 = v3[2]*v3[2]
    s = sqrt( v302 + v312 + v322 )
    v3[0] /= s
    v3[1] /= s
    v3[2] /= s
    
    /* Generate the 2nd unit vector for the new co-ordinate frame by cross product of v3 and r1*/
    v2[0] =  v3[1]*r1[2] - v3[2]*r1[1]
    v2[1] = -v3[0]*r1[2] + v3[2]*r1[0]
    v2[2] =  v3[0]*r1[1] - v3[1]*r1[0]

    /* Ensure it is normalised |v2|=1 */
    let v202 = v2[0]*v2[0]
    let v212 = v2[1]*v2[1]
    let v222 = v2[2]*v2[2]
    s = sqrt( v202 + v212 + v222 )
    v2[0] /= s
    v2[1] /= s
    v2[2] /= s

    /* Generate the 1st unit vector for the new co-ordinate frame by cross product of v2 and v3*/
    v1[0] =  v2[1]*v3[2] - v2[2]*v3[1]
    v1[1] = -v2[0]*v3[2] + v2[2]*v3[0]
    v1[2] =  v2[0]*v3[1] - v2[1]*v3[0]

    /* Ensure it is normalised |v1|=1 */
    let v102 = v1[0]*v1[0]
    let v112 = v1[1]*v1[1]
    let v122 = v1[2]*v1[2]
    s = sqrt( v102 + v112 + v122 )
    v1[0] /= s
    v1[1] /= s
    v1[2] /= s
    /* Transfer to the row matrix form for the result*/ 
    rmat[0, 0] = v1[0]
    rmat[0, 1] = v1[1]
    rmat[0, 2] = v1[2]
    rmat[1, 0] = v2[0]
    rmat[1, 1] = v2[1]
    rmat[1, 2] = v2[2]
    rmat[2, 0] = v3[0]
    rmat[2, 1] = v3[1]
    rmat[2, 2] = v3[2]

}

// MARK: Extension Vector

extension Vector {
    
    var x: Element {
        get {
            scalars[0]
        }
        set {
            scalars[0] = newValue
        }
    }
    
    var y: Element {
        get {
            scalars[1]
        }
        set {
            scalars[1] = newValue
        }
    }
    
    var z: Element {
        get {
            scalars[2]
        }
        set {
            scalars[2] = newValue
        }
    }

}

extension Vector where Scalar == Float {
    func canBeNormalized() -> Bool {
        if x == 0 && y == 0 && z == 0 {
            return false
        }
        return (canBeSquared(x) && canBeSquared(y) && canBeSquared(z))
    }
    /*! \brief Construct a unit vector orthogonal to *this.

    It requires that *this is normalizable; otherwise it just
    returns false. See CanBeNormalized()

    @param res reference by which to pass the result.

    @returns always true. (Return value kept for compatibility,
                as old versions of OpenBabel used to check for
                normalizability).
    */
    func createOrthoVector(_ res: inout Vector<Scalar>) -> Bool {
        /* Let us compute the crossed product of *this with a vector
        that is not too close to being colinear to *this.
        */
        /* unless the _vx and _vy coords are both close to zero, we can
        * simply take ( -_vy, _vx, 0 ) and normalize it.
        */
        if !isNegligible(x, z) || !isNegligible(y, z) {
            let norm = sqrt(x.square + y.square)
            res.x = -y / norm
            res.y = x / norm
            res.z = 0
        } else {
            /* if both _vx and _vy are close to zero, then the vector is close
            * to the z-axis, so it's far from colinear to the x-axis for instance.
            * So we take the crossed product with (1,0,0) and normalize it.
            */
            let norm = sqrt(y.square + z.square)
            res.x = 0
            res.y = -z / norm
            res.z = y / norm
        }
        return true
    }

    mutating func normalize() {
        let norm = sqrt(x.square + y.square + z.square)
        x /= norm
        y /= norm
        z /= norm
    }
    
    /*! Replaces *this with a random unit vector, which is (supposed
        to be) uniformly distributed over the unit sphere. Uses the
        system number generator with a time seed.
    */
    mutating func randomUnitVector() {
        // obtain a random vector with 0.001 <= length^2 <= 1.0, normalize
        // the vector to obtain a random vector of length 1.0.
        var l: Scalar = 0
        repeat {
            x = 2 * Float.random(in: 0...1) - 1
            y = 2 * Float.random(in: 0...1) - 1
            z = 2 * Float.random(in: 0...1) - 1
            l = x.square + y.square + z.square
        } while l > 1.0 || l < 1e-4
        normalize()
    }


}

extension Vector where Scalar == Double {
    
    func canBeNormalized() -> Bool {
        if x == 0 && y == 0 && z == 0 {
            return false
        }
        return (canBeSquared(x) && canBeSquared(y) && canBeSquared(z))
    }

    mutating func normalize() {
        let norm = sqrt(x.square + y.square + z.square)
        x /= norm
        y /= norm
        z /= norm
    }

    /*! Replaces *this with a random unit vector, which is (supposed
        to be) uniformly distributed over the unit sphere. Uses the
        system number generator with a time seed.
    */
    mutating func randomUnitVector() {
        // obtain a random vector with 0.001 <= length^2 <= 1.0, normalize
        // the vector to obtain a random vector of length 1.0.
        var l: Scalar = 0
        repeat {
            x = 2 * Double.random(in: 0...1) - 1
            y = 2 * Double.random(in: 0...1) - 1
            z = 2 * Double.random(in: 0...1) - 1
            l = x.square + y.square + z.square
        } while l > 1.0 || l < 1e-4
        normalize()
    }

    /*! \brief Construct a unit vector orthogonal to *this.

    It requires that *this is normalizable; otherwise it just
    returns false. See CanBeNormalized()

    @param res reference by which to pass the result.

    @returns always true. (Return value kept for compatibility,
                as old versions of OpenBabel used to check for
                normalizability).
    */
    func createOrthoVector(_ res: inout Vector<Scalar>) -> Bool {
        /* Let us compute the crossed product of *this with a vector
        that is not too close to being colinear to *this.
        */
        /* unless the _vx and _vy coords are both close to zero, we can
        * simply take ( -_vy, _vx, 0 ) and normalize it.
        */
        if !isNegligible(x, z) || !isNegligible(y, z) {
            let norm = sqrt(x.square + y.square)
            res.x = -y / norm
            res.y = x / norm
            res.z = 0
        } else {
            /* if both _vx and _vy are close to zero, then the vector is close
            * to the z-axis, so it's far from colinear to the x-axis for instance.
            * So we take the crossed product with (1,0,0) and normalize it.
            */
            let norm = sqrt(y.square + z.square)
            res.x = 0
            res.y = -z / norm
            res.z = y / norm
        }
        return true
    }
    
    static prefix func - (_ vec: Vector) -> Vector {
        for var scalar in vec.scalars {
            scalar *= -1
        }
        return vec
    }
    
}

// MARK: MATH Lib

//! Tests whether its argument can be squared without triggering an overflow or
//! underflow.
func canBeSquared(_ a: Double) -> Bool {
    if a == 0 { return true }
    let max_squarable_double = 1e150
    let min_squarable_double = 1e-150
    let abs_a = abs(a)
    return abs_a < max_squarable_double && abs_a > min_squarable_double
}

func canBeSquared(_ a: Float) -> Bool {
    if a == 0 { return true }
    let max_squarable_float: Float = 1e37
    let min_squarable_float: Float = 1e-37
    let abs_a: Float = abs(a)
    return abs_a < max_squarable_float && abs_a > min_squarable_float
}

/* Calculate the signed distance of point a to the plane determined by b,c,d */
func point2PlaneSigned(_ a: Vector<Double>, b: Vector<Double>, c: Vector<Double>, d: Vector<Double>) -> Double {
    let v_ba = a - b
    let v_normal = cross3x3(c - b, d - b)
    return Surge.dot(v_normal, v_ba) / length(v_normal)
}

/* Calculate the distance of point a to the plane determined by b,c,d */
func point2Plane(_ a: Vector<Double>, b: Vector<Double>, c: Vector<Double>, d: Vector<Double>) -> Double {
    return abs(point2PlaneSigned(a, b: b, c: c, d: d))
}

/* Calculate the angle between point a and the plane determined by b,c,d */
func point2PlaneAngle(_ a: Vector<Double>, b: Vector<Double>, c: Vector<Double>, d: Vector<Double>) -> Double {
    let v_ac = a - c
    let v_bc = b - c
    let v_cd = c - d
    let v_normal = cross3x3(v_bc, v_cd)
    return 90.0 - vector_angle(v_normal, v_ac)
}

func point2Line(_ a: Vector<Double>, b: Vector<Double>, c: Vector<Double>) -> Double {
    //http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    let v_cb = c - b
    let v_normal = cross3x3(a - b, a - c)
    return abs(length(v_normal) / length(v_cb))
}
