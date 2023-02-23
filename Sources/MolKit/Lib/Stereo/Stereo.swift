
import Foundation
import Collections
/**
 * Special case Ref values.
 */
public enum RefValue: Equatable {
    case NoRef            //!< No Ref set (invalid Ref)
    case ImplicitRef      //!< Implicit Ref (i.e. hydrogen, N lone pair, ...).
    case Ref(_ value: Int)
    
    init?(rawValue: Int) {
        self = .Ref(rawValue)
    }
    
    public static func == (_ lhs: RefValue, _ rhs: RefValue) -> Bool {
        switch lhs {
        case .NoRef:
            return false
        case .ImplicitRef:
            return false
        case .Ref(let value):
            switch rhs {
            case .NoRef:
                return false
            case .ImplicitRef:
                return false
            case .Ref(let rvalue):
                return value == rvalue
            }
        }
    }
    
    static func < (_ lhs: RefValue, _ rhs: RefValue) -> Bool {
        switch lhs {
        case .NoRef:
            return false
        case .ImplicitRef:
            return false
        case .Ref(let value):
            switch rhs {
            case .NoRef:
                return false
            case .ImplicitRef:
                return false
            case .Ref(let rvalue):
                return value < rvalue
            }
        }
    }
}

typealias Ref = RefValue
typealias Refs = [Ref]

public struct MKStereo {
    
    /**
     * The various types of stereochemistry
     */
    enum TType: Int {
        case CisTrans = 0                //!< cis/trans double bond
        case ExtendedCisTrans = 2        //!< allene, biphenyl, ...
        case SquarePlanar = 4            //!< Square-planar stereochemistry
        case Tetrahedral = 8             //!< tetrahedral
        case ExtendedTetrachedral = 16   //!< extended tetrahedral
        case TrigonalBipyramidal = 32    //!< Trigonal-bipyramidal stereochemistry
        case Octahedral = 64             //!< Octahedral stereochemistry
    }
    
    /**
     * Bond directions used by StereoFrom0D() to translate to
     * internal CisTransStereo representation.
     */
    enum BondDirection: Int { // Values taken from MDL format
        case NotStereo =   0
        case UpBond =      1
        case DownBond =    6
        case UnknownDir =  4
    }
    
    /**
     * Shapes used by MKTetraPlanarStereo subclasses for
     * setting/getting reference ids.
     *
     * @image html SPshapes.png
     * @sa MKTetraPlanarStereo
     */
    enum Shape: Int {
        case ShapeU = 1
        case ShapeZ = 2
        case Shape4 = 3
    }
    
    /**
     * Views used by MKTetraNonPlanarStereo subclasses for
     * setting/getting reference ids.
     * @sa MKTetraNonPlanarStereo
     */
    enum View: Int
    {
        case ViewFrom = 1 //!< view from the atom (id parameter) towards the center atom
        case ViewTowards = 2 //!< view from center atom towards the atom (id parameter)
    }
    
    /**
     * Windings used by MKTetraNonPlanar subclasses for
     * setting/getting reference ids.
     * @sa MKTetraNonPlanar
     */
    enum Winding: Int {
        case Clockwise = 1      //!< Clockwise winding
        case AntiClockwise = 2  //!< AntiClockwise winding (or CounterClockwise)
        case UnknownWinding = 3 //!< The configuration is specified as unknown (squiggly line in depiction)
    }
    
    
    static func makeRefs(_ ref1: RefValue, _ ref2: RefValue, _ ref3: RefValue, _ ref4: RefValue = .NoRef) -> Refs {
        if ref4 == .NoRef {
            return [ref1, ref2, ref3].map { $0 }
        } else {
            return [ref1, ref2, ref3, ref4].map { $0 }
        }
    }
    
    static func containsSameRef(_ refs1: Refs, _ refs2: Refs) -> Bool {
        if refs1.count != refs2.count {
            return false
        }
        let count = zip(refs1, refs2).compactMap { (r1, r2) in r1 == r2 ? 1 : 0 }.reduce(0, +)
        return refs1.count == count
    }
    
    static func containsRef(_ refs: Refs, _ ref: Ref) -> Bool {
        for rref in refs {
            if rref == ref {
                return true
            }
        }
        return false
    }
    
    ///@name Low-level functions used by implementation.
    //@{
    /**
     * Compute the inversion vector for @p refs and return the sum of it's
     * elements. The ith element in the inversion vector is the number of
     * element to the right of element i with a lower value.
     *
     * The number of inversions is the same as the number of interchanges
     * of consecutive elements.
     *
     * When working with 3 refs from a tetrahedral configuration:
     * @code
     * permutation   inversion vector    sum
     * -------------------------------------
     * 123           0 0 0               0 (even) -> clockwise
     * 132           0 1 0               1 (odd)  -> anti-clockwise
     * 213           1 0 0               1 (odd)  -> anti-clockwise
     * 231           1 1 0               2 (even) -> clockwise
     * 312           2 0 0               2 (even) -> clockwise
     * 321           2 1 0               3 (odd)  -> anti-clockwise
     * @endcode
     */
    static func numInversions(_ refs: Refs) -> Int {
        var invVec: Refs = []
        for i in 0..<refs.count {
            var e = 0  // ith element
            // loop over elements to the right
            for j in i..<refs.count {
                if refs[j] < refs[i] { e += 1 }
            }
            invVec.append(RefValue.Ref(e))
        }
//        return sum of invVec
        return invVec.reduce(into: 0) { partialResult, ref in
            switch ref {
            case .Ref(let value):
                partialResult += value
            default:
                break
            }
        }
    }
    
    /**
     * Permutate element @p i with @p j in @p refs.
     *
     * @param refs The sequence with N elements to permutate.
     * @param i Element i (0...N-1) will be mutated to j and vice versa.
     * @param j Element j (0...N-1) will be mutated to i and vice versa.
     *
     * @note This method does nothing if i equals j.
     */
    static func permutate(_ refs: inout Refs, _ i: Int, _ j: Int) {
        if i >= refs.count { return }
        if j >= refs.count { return }
        let id = refs[i]
        refs[i] = refs[j]
        refs[j] = id
    }
    
    /**
     * Get @p refs with element @p i and @p j permutated.
     *
     * @param refs The sequence with N elements to permutate.
     * @param i Element @p i (0...N-1) will be mutated to @p j and vice versa.
     * @param j Element @p j (0...N-1) will be mutated to @p i and vice versa.
     *
     * @return @p refs with elements @p i and @p j permutated.
     *
     * @note This method does nothing if @p i equals @p j.
     */
    static func permuted(_ refs: Refs, _ i: Int, _ j: Int) -> Refs {
        if i >= refs.count { return refs }
        if j >= refs.count { return refs }
        var result = refs
        result[i] = refs[j]
        result[j] = refs[i]
        return result
    }
    
}

public struct MKStereoUnit {
    
    var type: MKStereo.TType
    var id: RefValue
    var para: Bool
    
    init() {
        self.type = .CisTrans
        self.id = .NoRef
        self.para = false
    }
    
    init(_ type: MKStereo.TType, _ id: Ref, _ para: Bool) {
        self.type = type
        self.id = id
        self.para = para
    }
}

/**
 * @brief A single set of MKStereoUnit objects.
 *
 * This type can be used to represent all stereogenic units in a molecule and
 * is used as return type of FinStereogenicUnits(). This set is also the input
 * for many functions requiring this information (e.g. StereoFrom2D, ...).
 */
typealias MKStereoUnitSet = [MKStereoUnit]
/**
 * @brief A set of sets of MKStereoUnit objects.
 *
 * This type is used for cases where there is some relationship between
 * individual MKStereoUnit objects.
 */
typealias MKStereoUnitSets = [MKStereoUnitSet]

class MKStereoBase: MKGenericData {

    var m_mol: MKMol 
    var m_specified: Bool

    init(_ mol: MKMol) {
        self.m_mol = mol
        self.m_specified = true
        super.init("StereoData", .StereoData, .perceived)
    }

    func getMolecule() -> MKMol {
        return self.m_mol
    }

    func setSpecified(_ specified: Bool) {
        self.m_specified = specified
    }

    func isSpecified() -> Bool {
        return self.m_specified
    }

    func getType() -> MKStereo.TType {
        return .CisTrans
    }
}

class MKStereoFacade {
    
    var m_mol: MKMol
    var m_init: Bool = false
    var m_perceive: Bool
    
    var m_tetrahedralMap: OrderedDictionary<Int, MKTetrahedralStereo> = [:]
    var m_cistransMap: OrderedDictionary<Int, MKCisTransStereo> = [:]
    var m_squarePlanarMap: OrderedDictionary<Int, MKSquarePlanarStereo> = [:]
    
    
    init(_ m_mol: MKMol, m_perceive: Bool = true) {
        self.m_mol = m_mol
        self.m_perceive = m_perceive
    }
    
    func numTetrahedralStereo() -> UInt {
        return UInt(0)
    }
    
    func getAllTetrahedralStereo() -> [MKTetrahedralStereo] {
        return []
    }
    
    func hasTetrahedralStereo(_ atomId: Int) -> Bool {
        return false
    }
        
    func getTetrahedralStereo(_ atomId: Int) -> MKTetrahedralStereo? {
        return nil
    }
    
    
}

