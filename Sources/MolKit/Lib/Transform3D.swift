

import Foundation 
import Surge
import Accelerate

class MKTransform3D {
    
    var _m: Matrix<Double> = Matrix.init(rows: 3, columns: 3, repeatedValue: 0.0)
    var _v: Vector<Double> = Vector.init(dimensions: 3, repeatedValue: 0.0)

    public init() {
        self._m = Matrix<Double>.init(rows: 3, columns: 3, repeatedValue: 0.0)
        self._v = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
    }

    public init(m: Matrix<Double>, v: Vector<Double>) {
        assert({ m.rows == 3 && m.columns == 3 }(), "Matrix is not 3x3")
        assert({ v.dimensions == 3 }(), "Vector does not 3 dimensions")
        
        self._m = m
        self._v = v
        self.normalize()
    }

    public func normalize() { 
        self._v[0] -= floor(self._v[0] + 0.01)
        self._v[1] -= floor(self._v[1] + 0.01)
        self._v[2] -= floor(self._v[2] + 0.01)
    }
    
    public func describeAsString() -> String {
        let v = self._v
        let m = self._m
        var outString: String = ""
        var first: Bool = false
        var neg: Bool = false
        var j: Int = 0
        var n: Int = 0
        
        for i in 0..<3 {
            
            if i != 0 { outString += "," }
            
            n = Int(floor(v[i] * 12.0 + 0.1))
            j = 0
            while (m[i,j] == 0) {
                j+=1
            }
            neg = m[i, j] < 0.0
            switch (n) {
            case 2:
                outString += ((neg) ? "1/6" : "1/6+")
            case 3:
                outString += ((neg) ? "1/4" : "1/4+")
            case 4:
                outString += ((neg) ? "1/3" : "1/3+")
            case 6:
                outString += ((neg) ? "1/2" : "1/2+")
            case 8:
                outString += ((neg) ? "2/3" : "2/3+")
            case 9:
                outString += ((neg) ? "3/4" : "3/4+")
            case 10:
                outString += ((neg) ? "5/6" : "5/6+")
            default:
                continue
            }
            first = true
            while (j < 3) {
                if (m[i, j] != 0.0) {
                    neg = m[i, j] < 0.0
                    switch (j) {
                    case 0:
                        outString += ((neg) ? "-x": (first ? "x": "+x"))
                    case 1:
                        outString += ((neg) ? "-y": (first ? "y": "+y"))
                    case 2:
                        outString += ((neg) ? "-z": (first ? "z": "+z"))
                    default:
                        continue
                    }
                    first = false
                }
                j+=1
            }
            
        }
        return outString
        
    }
    
    public func describeAsValue() -> String {
        let v = self._v
        let m = self._m
        var outString: String = ""
        
        outString += "\(m[0,0]) \(m[0,0]) \(m[0,0]) \(v[0]) "
        outString += "\(m[1,0]) \(m[1,0]) \(m[1,0]) \(v[1]) "
        outString += "\(m[2,0]) \(m[2,0]) \(m[2,0]) \(v[2])" 
        
        return outString
    }

    static func * (_ m: MKTransform3D, _ v: Vector<Double>) -> Vector<Double> {
        return m._m * v + m._v
    }
    
    static func * (_ m: MKTransform3D, _ v: MKTransform3D) -> MKTransform3D {
        return MKTransform3D(m: m._m * v._m, v: v._v)
    }
}

