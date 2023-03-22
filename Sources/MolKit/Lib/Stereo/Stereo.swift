
import Foundation
import Collections
import Algorithms
import Surge

fileprivate let DELTA_ANGLE_FOR_OVERLAPPING_BONDS = 4.0

/**
 * Special case Ref values.
 */
public enum RefValue: Equatable, Hashable {
    case NoRef            //!< No Ref set (invalid Ref)
    case ImplicitRef      //!< Implicit Ref (i.e. hydrogen, N lone pair, ...).
    case Ref(_ value: Int)
    
    init?(rawValue: Int) {
        self = .Ref(rawValue)
    }
    
    var intValue: Int? {
        switch self {
        case .NoRef:
            return nil
        case .ImplicitRef:
            return -9999 // SMARTS ImplicitRef number
        case .Ref(var value):
            return value
        }
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
    
    static func == (_ lhs: RefValue, _ rhs: from_or_towrds) -> Bool {
        return lhs.intValue == rhs.value
    }
    
    static func -= (_ lhs: RefValue, _ rhs: Int) {
        if var refVal = lhs.intValue {
            refVal -= rhs
        }
    }
    
    static func - (_ lhs: RefValue, _ rhs: Int) -> RefValue {
        if case .Ref(var val) = lhs {
            val -= rhs
            return lhs
        } else if case .ImplicitRef = lhs {
            return .ImplicitRef
        }
        return lhs
    }
    
    public func hash(into hasher: inout Hasher) {
        hasher.combine(self.intValue ?? 0)
    }
}

extension RefValue: ExpressibleByIntegerLiteral {
    public typealias IntegerLiteralType = Int
    
    public init(integerLiteral value: Int) {
        self = .Ref(value)
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
    
    static func containsSameRefs(_ refs1: Refs, _ refs2: Refs) -> Bool {
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
    
    static func numInversions(_ refs: [UInt]) -> Int {
//        convert refs to Refs value and call original
        let refsIn = refs.map { Ref(rawValue: Int($0))! }
        return MKStereo.numInversions(refsIn)
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
    
    init(_ type: MKStereo.TType, _ id: Ref, _ para: Bool = false) {
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
        fatalError("not implemented in base class")
    }
}

/**
* @class OBStereoFacade stereo.h <openbabel/stereo/stereo.h>
* @brief Facade to simplify retrieval of OBStereoBase derived objects.
*
* The OBStereoFacade helps with retrieving OBStereoBase derived objects
* (i.e. OBTetrahedralStereo, OBCisTransStereo, ...) from an OBMol. This
* is done by iterating over all OBGenericData objects with data type
* OBGenericDataType::StereoData and checking the OBStereo::Type using
* OBStereoBase::GetType.
*
* @sa OBStereo OBStereoBase
* @since version 2.3
*/

class MKStereoFacade {
    
    var m_mol: MKMol
    var m_init: Bool = false
    var m_perceive: Bool = false
    
    var m_tetrahedralMap: OrderedDictionary<Int, MKTetrahedralStereo> = [:]
    var m_cistransMap: OrderedDictionary<Int, MKCisTransStereo> = [:]
    var m_squarePlanarMap: OrderedDictionary<Int, MKSquarePlanarStereo> = [:]
    
    
    init(_ m_mol: MKMol, m_perceive: Bool = true) {
        self.m_mol = m_mol
        self.m_perceive = m_perceive
    }

    private func ensureInit() {
        if !m_init { intialize() }
    }
    
    func intialize() {
        
        if (m_perceive && !m_mol.hasChiralityPerceived()) {
            perceiveStereo(&m_mol)
        }
        
        guard let stereoData = m_mol.getDataVector(.StereoData) else {
            fatalError("FAILED to retrieve stereo data for mol")
        }

        for data in stereoData {
            let type: MKStereo.TType = (data as! MKStereoBase).getType()
            if (type == .Tetrahedral) {
                let ts: MKTetrahedralStereo = (data as! MKTetrahedralStereo)
                let config: MKTetrahedralStereo.Config = ts.getConfig()
                if (config.center == .NoRef) {
                    continue
                }
                m_tetrahedralMap[config.center.intValue!] = ts
            } else {
                if (type == .SquarePlanar) {
                    let sp: MKSquarePlanarStereo = (data as! MKSquarePlanarStereo)
                    let config: MKSquarePlanarStereo.Config = sp.getConfig()
                    if (config.center == .NoRef) {
                        continue
                    }
                    m_squarePlanarMap[config.center.intValue!] = sp
                } else {
                    if (type == .CisTrans) {
                        let ct: MKCisTransStereo = (data as! MKCisTransStereo)
                        let config: MKCisTransStereo.Config = ct.getConfig()
                        // find the bond id from begin & end atom ids
                        var id: RefValue = .NoRef
                        let a: MKAtom? = m_mol.getAtomById(config.begin)
                        if (a == nil) {
                            continue
                        }
                        guard let bonds = a?.getBondIterator() else {
                            print("ERROR: No bonds for atom in stereo data??")
                            continue
                        }
                        for bond in bonds {
                            let beginId: Ref = bond.getBeginAtom().getId().ref
                            let endId: Ref = bond.getEndAtom().getId().ref
                            if ((beginId == config.begin && endId == config.end) ||
                                (beginId == config.end && endId == config.begin)) {
                                id = bond.getId()
                                break
                            }
                        }
                        if (id == .NoRef) {
                            continue
                        }
                        m_cistransMap[id.intValue!] = ct
                    }
                }
            }
        }
        m_init = true
    }
    
    ///@name Tetrahedral stereochemistry
      ///@{
      /**
       * Get the number of tetrahedral stereocenters.
       */
    func numTetrahedralStereo() -> UInt {
        ensureInit()
        return UInt(m_tetrahedralMap.count)
    }
    
    /**
       * Get all the OBTetrahedralStereo objects.
       */
    func getAllTetrahedralStereo() -> [MKTetrahedralStereo] {
        ensureInit()
        var result: [MKTetrahedralStereo] = []
        for (_, value) in m_tetrahedralMap {
            result.append(value)
        }
        return result
    }
    
    /**
       * Check if atom with @p id is a tetrahedral center.
       * @return True if the atom with @p id has tetrahedral stereochemistry.
       */
    func hasTetrahedralStereo(_ atomId: Int) -> Bool {
        ensureInit()
        return m_tetrahedralMap.contains { (key: Int, value: MKTetrahedralStereo) in
            key == atomId
        }
    }
    
    /**
       * Get the OBTetrahedralStereo object with @p atomId as center. This
       * function returns 0 if there is no OBTetrahedralStereo object found
       * with the specified center.
       */
    func getTetrahedralStereo(_ atomId: Int) -> MKTetrahedralStereo? {
        if !hasTetrahedralStereo(atomId) {
            return nil
        }
        return m_tetrahedralMap[atomId]
    }
    
    ///@name Cis/Trans stereochemistry
      ///@{
      /**
       * Get the number of cis/trans stereocenters.
       */
    func numCisTransStereo() -> UInt {
        ensureInit()
        return UInt(m_cistransMap.count)
    }

    /**
       * Get all the OBCisTransStereo objects.
       */
    func getAllCisTransStereo() -> [MKCisTransStereo] {
        ensureInit()
        var result: [MKCisTransStereo] = []
        for (_, value) in m_cistransMap {
            result.append(value)
        }
        return result
    }

    /**
       * Check if bond with @p id is a stereogenic cis/trans double bond.
       * @return True if the bond with @p id has cis/trans stereochemistry.
       */
    func hasCisTransStereo(_ bondId: Int) -> Bool {
        ensureInit()
        return m_cistransMap.contains { (key: Int, value: MKCisTransStereo) in
            key == bondId
        }
    }

    /**
       * Get the OBTetrahedralStereo object with @p bondId as double bond.
       * This function returns 0 if there is no OBCisTransStereo object found
       * with the specified bond.
       */
    func getCisTransStereo(_ bondId: Int) -> MKCisTransStereo? {
        if !hasCisTransStereo(bondId) {
            return nil
        }
        return m_cistransMap[bondId]
    }

    ///@name SquarePlanar stereochemistry
      ///@{
      /**
       * Get the number of square-planar stereocenters.
       */
    func numSquarePlanarStereo() -> UInt {
        ensureInit()
        return UInt(m_squarePlanarMap.count)
    }

    /**
    * Get all the OBSquarePlanarStereo objects.
    */
    func getAllSquarePlanarStereo() -> [MKSquarePlanarStereo] {
        ensureInit()
        var result: [MKSquarePlanarStereo] = []
        for (_, value) in m_squarePlanarMap {
            result.append(value)
        }
        return result
    }

    /**
       * Check if atom with @p id is a stereogenic square-planar atom.
       * @return True if the atom with @p id has square-planar stereochemistry.
       */
    func hasSquarePlanarStereo(_ atomId: Int) -> Bool {
        ensureInit()
        return m_squarePlanarMap.contains { (key: Int, value: MKSquarePlanarStereo) in
            key == atomId
        }
    }

    /**
       * Get the OBSquarePlanarStereo object with @p atomId as center. This
       * function returns 0 if there is no OBSquarePlanarStereo object found
       * with the specified center.
       */
    func getSquarePlanarStereo(_ atomId: Int) -> MKSquarePlanarStereo? {
        if !hasSquarePlanarStereo(atomId) {
            return nil
        }
        return m_squarePlanarMap[atomId]
    }

    
    func hasStereo(_ T: MKStereo.TType, _ id: Int) -> Bool {
        switch T {
        case .Tetrahedral:
            return self.hasTetrahedralStereo(id)
        case .CisTrans:
            return self.hasCisTransStereo(id)
        case .SquarePlanar:
            return self.hasSquarePlanarStereo(id)
        default:
            return false
        }
    }

    func getStereo<U: MKStereoBase>(_ T: MKStereo.TType, _ id: Int) -> U? {
        switch T {
        case .Tetrahedral:
            return self.getTetrahedralStereo(id) as? U
        case .CisTrans:
            return self.getCisTransStereo(id) as? U
        case .SquarePlanar:
            return self.getSquarePlanarStereo(id) as? U
        default:
            return nil
        }
    }
}

///@name High level functions
///@{
/**
* Convert 0D/2D/3D coordinates to OBStereo objects. The right function will
* be selected based on the molecule's dimensionality
* (i.e. OBMol::GetDimension()).
*
* @sa StereoFrom3D StereoFrom2D StereoFrom0D
* @since version 2.3
*/
func perceiveStereo(_ mol: inout MKMol,_ force: Bool = false) {

    switch mol.getDimension() {
    case 3: 
        stereoFrom3D(&mol, force)
    case 2: 
        stereoFrom2D(mol, nil, force)
    default: 
        stereoFrom0D(mol)
    }
    //  TODO: add logging here
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
//
//  From3D
//
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

/**
* Convert the 3D coordinates of molecule @p mol to OBStereo objects. This
* function makes use of the lower level functions TetrahedralFrom3D(),
* CisTransFrom3D(), SquarePlanarFrom3D(), ...
*
* Unless perception is forced, this function does nothing if stereochemistry
* has already been perceived (i.e. OBMol::HasChiralityPerceived()). Before
* doing the actual perception, any data of the OBGenericDataType::StereoData
* type will be deleted.
*
* @param mol The molecule containing 3D coordinates.
* @param force Force to run the perception even if the results are cached.
*
* @sa StereoFrom3D StereoFrom0D PerceiveStereo
* @since version 2.3
*/
func stereoFrom3D(_ mol: inout MKMol,_ force: Bool = false) {
    if mol.hasChiralityPerceived() && !force {
        return
    }
    var symmetryClasses = findSymmetry(mol).map { UInt($0.intValue!) }
    let stereogenicUnits = findStereogenicUnits(mol, symClasses: &symmetryClasses)
    mol.deleteData(.StereoData)
    tetrahedralFrom3D(mol, stereogenicUnits)
    cisTransFrom3D(mol, stereogenicUnits)
    mol.setChiralityPerceived()
}


//! Calculate the "sign of a volume" given by a set of 4 coordinates

 func volumeSign(_ a: Vector<Double>, _ b: Vector<Double>, _ c: Vector<Double>, _ d: Vector<Double>) -> Int {
//     let v = (b - a).cross(c - a).dot(d - a)
     let v = dot(cross3x3((b-a), (c-a)), (d-a))
     return v > 0 ? 1 : v < 0 ? -1 : 0
 }
// TODO: Make sure these are equivalent...
//func volumeSign(_ a: Vector<Double>, _ b: Vector<Double>, _ c: Vector<Double>, _ d: Vector<Double>) -> Double {
//
//    let A = b - a
//    let B = c - a
//    let C = d - a
//    let m: Matrix<Double> = Matrix<Double>.init(vecs: [A,B,C])
//    return det(m)!
//}

///@name Low level functions
///@{
/**
* Get a vector with all OBTetrahedralStereo objects for the molecule. This
* function is used by StereoFrom3D() with the @p addToMol parameter is set
* to true.
*
* The algorithm to convert the 3D coordinates to OBTetrahedralStereo object
* uses the sign of the volume described by the 4 center atom neighbors. Given
* 4 points \f$a\f$, \f$b\f$, \f$c\f$ and \f$d\f$, the signed volume \f$S_v\f$
* is defined as:
*
    \f[ S_v = \left| \begin{array}{ccc}
    x_b - x_a & y_b - y_a & z_b - z_a \\
    x_c - x_a & y_c - y_a & z_c - z_a \\
    x_d - x_a & y_d - y_a & z_d - z_a
    \end{array} \right| \f]
*
* The sign of \f$S_v\f$ changes when any of the points cross the plane defined
* by the other 3 points. To make this less abstract one could say that
* a change of sign is equal to inverting the tetrahedral stereochemistry.
*
* In case there are only 3 neighbor atoms for the tetrahedral center, the
* center atom itself is used as 4th point. This only changes the magnitude
* and not the sign of \f$S_v\f$ because the center atom is still on the same
* side of the plane.
*
* This function is also used for symmetry analysis to handle cases where
* there are two atoms in the same symmetry class that don't have the same
* stereochemistry. In this situation, the @p addToMol parameter is set to
* false and the returned objects will need to be deleted explicitly.
*
* @param mol The molecule.
* @param stereoUnits The stereogenic units.
* @param addToMol If true, the OBTetrahedralStereo objects will be added
* to the molecule using OBBase::SetData().
*
* @sa StereoFrom3D FindStereogenicUnits
* @since version 2.3
*/
@discardableResult
func tetrahedralFrom3D(_ mol: MKMol, _ stereoUnits: MKStereoUnitSet, _ addToMol: Bool = true) -> [MKTetrahedralStereo] {
    var configs: [MKTetrahedralStereo] = []
    let unc: MKUnitCell? = mol.getData(.UnitCell) as? MKUnitCell

    // find all tetrahedral centers
    let centers = stereoUnits.filter { $0.type == .Tetrahedral }.map { $0.id }

    for i in centers {
        guard let center = mol.getAtomById(i) else {
            continue
        }
        // make sure we have at least 3 heavy atom neighbors
        if center.getHeavyDegree() < 3 {
            print("Cannot calculate a signed volume for an atom with a heavy atom valence of \(center.getHeavyDegree())")
            continue 
        }

        var config: MKTetrahedralStereo.Config = MKTetrahedralStereo.Config()
        config.center = i
        for nbr in center.getNbrAtomIterator()! {
            if config.from_or_towrds.refValue == .NoRef {
                config.from_or_towrds = .from(nbr.getId().ref)
            } else {
                config.refs.append(nbr.getId().ref)
            }
        }
        let useCentralAtom = false
        // Create a vector with the coordinates of the neighbor atoms
        // and check for a bond that indicates unspecified stereochemistry
        var nbrCoords: [Vector<Double>] = []
        let from = mol.getAtomById(config.from_or_towrds.refValue)
        let bond = mol.getBond(from!, center)
        if bond!.isWedgeOrHash() && bond!.getBeginAtom() == center {
            config.specified = false
        }
        let centerCoord = center.getVector()
        if unc != nil {
            nbrCoords.append(unc!.unwrapCartesianNear(from!.getVector(), centerCoord))
        } else {
            nbrCoords.append(from!.getVector())
        }
        for id in config.refs {
            let nbr = mol.getAtomById(id)
            if unc != nil {
                nbrCoords.append(unc!.unwrapCartesianNear(nbr!.getVector(), centerCoord))
            } else {
                nbrCoords.append(nbr!.getVector())
            }
            let bond = mol.getBond(nbr!, center)
            if bond!.isWedgeOrHash() && bond!.getBeginAtom() == center {
                config.specified = false
            }
        }
        // Checks for a neighbour having 0 co-ords (added hydrogen etc)
        /* FIXME: needed? if the molecule has 3D coords, additional
         * hydrogens will get coords using OBAtom::GetNewBondVector
        for (std::vector<vector3>::iterator coord = nbrCoords.begin(); coord != nbrCoords.end(); ++coord) {
          // are the coordinates zero to 6 or more significant figures
          if (coord->IsApprox(VZero, 1.0e-6)) {
            if (!use_central_atom) {
              use_central_atom = true;
            } else {
              obErrorLog.ThrowError(__FUNCTION__,
                  "More than 2 neighbours have 0 co-ords when attempting 3D chiral calculation", obInfo);
            }
          }
        }
        */ // original code note

        // If we have three heavy atoms we can use the chiral center atom itself for the fourth
        // will always give same sign (for tetrahedron), magnitude will be smaller.
        if config.refs.count == 2 || useCentralAtom {
            nbrCoords.append(centerCoord)
            config.refs.append(.ImplicitRef) // need to add largest number on end to work
        }
        let sign = volumeSign(nbrCoords[0], nbrCoords[1], nbrCoords[2], nbrCoords[3])
        if sign < 0 {
            config.winding = .AntiClockwise
        }
        let th = MKTetrahedralStereo(mol)
        th.setConfig(config)
        configs.append(th)
        // add the data to the molecule if needed
        if addToMol {
            mol.setData(th)
        }
    }

    return configs 
}

/**
* Get a vector with all OBCisTransStereo objects for the molecule. This
* function is used by StereoFrom3D() with the @p addToMol parameter is set
* to true.
*
* The algorithm to convert the 3D coordinates to OBCisTransStereo objects
* considers the signed distance between the attached atoms and the plane
* through the double bond at right angles to the plane of the attached
* atoms. Bonds on the same side (cis) will share the same sign for the
* signed distance.
*
* Missing atom coordinates (OBStereo::ImplicitRef) and their bond
* vectors will be computed if needed.
*
@verbatim
        0      3     Get signed distance of 0 and 2 to the plane
        \    /      that goes through the double bond and is at
        C==C       right angles to the stereo bonds.
        /    \
        1      2     If the two signed distances have the same sign
                    then they are cis; if not, then trans.
@endverbatim
*
* This function is also used for symmetry analysis to handle cases where
* there are two atoms in the same symmetry class that don't have the same
* stereochemistry. In this situation, the @p addToMol parameter is set to
* false and the returned objects will need to be deleted explicitly.
*
* @param mol The molecule.
* @param stereoUnits The stereogenic units.
* @param addToMol If true, the OBCisTransStereo objects will be added
* to the molecule using OBBase::SetData().
*
* @sa StereoFrom3D FindStereogenicUnits
* @since version 2.3
*/
@discardableResult
func cisTransFrom3D(_ mol: MKMol, _ stereoUnits: MKStereoUnitSet, _ addToMol: Bool = true) -> [MKCisTransStereo] {
    var configs: [MKCisTransStereo] = []
    let uc = mol.getData(.UnitCell) as? MKUnitCell
    let bonds = stereoUnits.filter { $0.type == .CisTrans }.map { $0.id }

    for i in bonds {
        guard let bond = mol.getBondById(i) else {
            continue
        }
        let begin = bond.getBeginAtom()
        let end = bond.getEndAtom()

        // Create a vector with the coordinates of the neighbor atoms
        var bondVecs: [Vector<Double>] = []
        let config: MKCisTransStereo.Config = MKCisTransStereo.Config()
        //begin 
        config.begin = begin.getId().ref
        for nbr in begin.getNbrAtomIterator()! {
            if nbr.getId() == end.getId() {
                continue
            }
            config.refs.append(nbr.getId().ref)
            if let uc = uc {
                bondVecs.append(uc.minimumImageCartesian(nbr.getVector() - begin.getVector()))
            } else {
                bondVecs.append(nbr.getVector() - begin.getVector())
            }
        }
        if config.refs.count == 1 {
            config.refs.append(.ImplicitRef)
            var pos = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
            pos = begin.getNewBondVector(1.0)
            // WARNING: GetNewBondVector code has not yet been checked, since it's part of builder.cpp
            if let uc = uc {
                bondVecs.append(uc.minimumImageCartesian(pos - begin.getVector()))
            } else {
                bondVecs.append(pos - begin.getVector())
            }
        }
        //end 
        config.end = end.getId().ref
        var end_vec = end.getVector()
        if let uc = uc {
            end_vec = uc.unwrapCartesianNear(end_vec, begin.getVector())
        }
        for nbr in end.getNbrAtomIterator()! {
            if nbr.getId() == begin.getId() {
                continue
            }
            config.refs.append(nbr.getId().ref)
            if let uc = uc {
                bondVecs.append(uc.minimumImageCartesian(nbr.getVector() - end_vec))
            } else {
                bondVecs.append(nbr.getVector() - end_vec)
            }
        }
        if config.refs.count == 3 {
            config.refs.append(.ImplicitRef)
            var pos = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
            pos = end.getNewBondVector(1.0)
            // WARNING: GetNewBondVector code has not yet been checked, since it's part of builder.cpp
            if let uc = uc {
                bondVecs.append(uc.minimumImageCartesian(pos - end_vec))
            } else {
                bondVecs.append(pos - end_vec)
            }
        }
        var tor02: Double = 0.0
        var tor03: Double = 0.0
        var tor12: Double = 0.0
        var tor13: Double = 0.0
        if let uc = uc {
            let v0 = begin.getVector() + bondVecs[0]
            let v1 = begin.getVector() + bondVecs[1]
            let v2 = end.getVector() + bondVecs[2]
            let v3 = end.getVector() + bondVecs[3]

            var b = uc.unwrapCartesianNear(begin.getVector(), v0)
            var c = uc.unwrapCartesianNear(end.getVector(), b)
            var d = uc.unwrapCartesianNear(v2, c)
            tor02 = calculateTorsionAngle(v0, b, c, d)

            d = uc.unwrapCartesianNear(v3, c)
            tor03 = calculateTorsionAngle(v0, b, c, d)

            b = uc.unwrapCartesianNear(begin.getVector(), v1)
            c = uc.unwrapCartesianNear(end.getVector(), b)
            d = uc.unwrapCartesianNear(v2, c)
            tor12 = calculateTorsionAngle(v1, b, c, d)

            d = uc.unwrapCartesianNear(v3, c)
            tor13 = calculateTorsionAngle(v1, b, c, d)
        } else {
            tor02 = calculateTorsionAngle(begin.getVector() + bondVecs[0], begin.getVector(), end.getVector(), end.getVector() + bondVecs[2])
            tor03 = calculateTorsionAngle(begin.getVector() + bondVecs[0], begin.getVector(), end.getVector(), end.getVector() + bondVecs[3])
            tor12 = calculateTorsionAngle(begin.getVector() + bondVecs[1], begin.getVector(), end.getVector(), end.getVector() + bondVecs[2])
            tor13 = calculateTorsionAngle(begin.getVector() + bondVecs[1], begin.getVector(), end.getVector(), end.getVector() + bondVecs[3])
        }

        if abs(tor02) < 90.0 && abs(tor03) > 90.0 {
            // 0      2 //
            //  \    /  //
            //   C==C   //
            //  /    \  //
            // 1      3 //
            config.shape = .ShapeZ
            if abs(tor12) < 90.0 || abs(tor13) > 90.0 {
                print("Could not determine cis/trans from 3D coordinates, using unspecified")
                config.specified = false
            }
        } else if (abs(tor02) > 90.0 && abs(tor03) < 90.0) {
            // 0      3 //
            //  \    /  //
            //   C==C   //
            //  /    \  //
            // 1      2 //
            config.shape = .ShapeU
            if abs(tor12) > 90.0 || abs(tor13) < 90.0 {
                print("Could not determine cis/trans from 3D coordinates, using unspecified")
                config.specified = false
            }
        } else {
            print("Could not determine cis/trans from 3D coordinates, using unspecified")
            config.shape = .ShapeU
            config.specified = false
        }
        let ct = MKCisTransStereo(mol)
        ct.setConfig(config)
        configs.append(ct)
        // add the data to the molecule if needed
        if addToMol {
            mol.setData(ct)
        }
    }

    return configs
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
//  From2D
//
//  Reference:
//  [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the
//  Unambiguous Identification of the Stereochemical Characteristics of
//  Compounds During Their Registration in Databases. Molecules 2000, 6,
//  915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

/**
* Convert the 2D depiction of molecule @p mol to OBStereo objects.
* This function makes use of the lower level functions
* TetrahedralFrom2D(), CisTransFrom2D(), SquarePlanarFrom2D(), ...
*
* First, symmetry analysis taking stereochemistry into account is
* performed iteratively (see OBGraphSym). Next the 2D coordinates,
* OBBond::Wedge, OBBond::Hash, OBBond::WedgeOrHash and OBBond::CisOrTrans
* are used to construct OBStereoBase derived objects to store the
* stereochemistry. These objects will be added to @p mol.
*
* Unless perception is forced, this function does nothing if stereochemistry
* has already been perceived (i.e. OBMol::HasChiralityPerceived()). Before
* doing the actual perception, any data of the OBGenericDataType::StereoData
* type will be deleted.
*
    @verbatim
    Reference:
    [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the
    Unambiguous Identification of the Stereochemical Characteristics of
    Compounds During Their Registration in Databases. Molecules 2000, 6,
    915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
    @endverbatim
*
* @param mol The molecule containing 2D coordinates.
* @param updown A map of OBStereo::BondDirection for cis/trans bonds
* @param force Force to run the perception even if the results are cached.
*
* @sa StereoFrom3D StereoFrom0D PerceiveStereo
* @since version 2.3
*/
func stereoFrom2D(_ mol: MKMol, _ updown: [MKBond: MKStereo.BondDirection]? = nil,_ force: Bool = false) {
    if mol.hasChiralityPerceived() && !force {
        return
    }
    // TODO: Log here 
    var symmetry_classes = findSymmetry(mol).map { UInt($0.intValue!) }
    let stereogenicUnits = findStereogenicUnits(mol, symClasses: &symmetry_classes)
    mol.deleteData(.StereoData)
    tetrahedralFrom2D(mol, stereogenicUnits)
    cisTransFrom2D(mol, stereogenicUnits, updown)
    mol.setChiralityPerceived()
}

//! Calculate the "sign of a triangle" given by a set of 3 2D coordinates
func triangleSign(_ a: Vector<Double>, _ b: Vector<Double>, _ c: Vector<Double>) -> Double {
    return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x)
}


//! Calculate whether three vectors are arranged in order of increasing
//! angle anticlockwise (true) or clockwise (false) relative to a central point.

func angleOrder(_ a: Vector<Double>, _ b: Vector<Double>, _ c: Vector<Double>, _ center: Vector<Double>) -> Bool {

    var t = a - center 
    t = normalize(t)
    var u = b - center
    u = normalize(u)
    var v = c - center
    v = normalize(v)
    return triangleSign(t, u, v) > 0.0
}
//! Get the angle between three atoms (from -180 to +180)
//! Note: OBAtom.GetAngle just returns 0->180
func getAngle(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom) -> Double {
    
    var v1 = a.getVector() - b.getVector()
    var v2 = c.getVector() - b.getVector()
    if a.isPeriodic() { // Adapted from OBAtom.GetAngle
        guard let mol = a.getParent() else {
            fatalError("parent unretrievable")
        }
        let box = mol.getData(.UnitCell) as? MKUnitCell
        v1 = box!.minimumImageCartesian(v1)
        v2 = box!.minimumImageCartesian(v2)
    }
    if isNearZero(length(v1), 1.0e-3) || isNearZero(length(v2), 1.0e-3) {
        return 0.0
    }
    var angle = (atan2(v2.y, v2.x) - atan2(v1.y, v1.x)).radiansToDegrees
    while angle < -180.0 { angle += 360 }
    while angle > 180 { angle -= 360 }
    return angle
    
}

/**
* Get a vector with all OBTetrahedralStereo objects for the molecule. This
* function is used by StereoFrom2D() with the @p addToMol parameter is set
* to true.
*
* The algorithm to convert the 2D coordinates and bond properties
* (i.e. OBBond::Wedge, OBBond::Hash, OBBond::WedgeOrHash and OBBond::CisOrTrans)
* uses the sign of a triangle. Given 3 points \f$a\f$, \f$b\f$ and \f$c\f$, the
* sign of the trianle \f$S_t\f$ is defined as:
*
    \f[ S_t = (x_a - x_c) (y_b - y_c) - (y_a - y_c) (x_b - x_c) \f]
*
* This is equation 6 from on the referenced web page. The 3 points used
* to calculate the triangle sign always remain in the same plane (i.e. z = 0).
* The actual meaning of \f$S_t\f$ (i.e. assignment of OBStereo::Winding) depends
* on the 4th atom. When the atom is in front of the plane, the sign should be
* changed to have the same absolute meaning for an atom behind the plane and the
* same triangle. It is important to note that none of the z coordinates is ever
* changed, the molecule always stays 2D (unlike methods which set a pseudo-z
* coordinate).
*
* @todo document bond property interpretation!
*
* This function is also used for symmetry analysis to handle cases where
* there are two atoms in the same symmetry class that don't have the same
* stereochemistry. In this situation, the @p addToMol parameter is set to
* false and the returned objects will need to be deleted explicitly.
*
    @verbatim
    Reference:
    [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the
    Unambiguous Identification of the Stereochemical Characteristics of
    Compounds During Their Registration in Databases. Molecules 2000, 6,
    915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
    @endverbatim
*
* @param mol The molecule.
* @param stereoUnits The stereogenic units.
* @param addToMol If true, the OBTetrahedralStereo objects will be added
* to the molecule using OBBase::SetData().
*
* @sa StereoFrom2D FindStereogenicUnits
* @since version 2.3
*/
@discardableResult
func tetrahedralFrom2D(_ mol: MKMol, _ stereoUnits: MKStereoUnitSet, _ addToMol: Bool = true) -> [MKTetrahedralStereo] {
    var configs: [MKTetrahedralStereo] = []

    // find all tetrahedral centers
    let centers = stereoUnits.filter { $0.type == .Tetrahedral }.map { $0.id }
    for i in centers {
        guard let center = mol.getAtomById(i) else {
            fatalError("atom \(i) not found")
        }
        // make sure we have at least 3 heavy atom neighbors
        if center.getHeavyDegree() < 3 {
            print("Cannot calculate a signed volume for an atom with a heavy atom valence of \(center.getHeavyDegree())")
            continue
        }
        var config: MKTetrahedralStereo.Config = MKTetrahedralStereo.Config()
        config.center = i
        // We assume the 'tip-only' convention. That is, wedge or hash bonds only
        // determine the stereochemistry at their thin end (the BeginAtom)
        let tiponly: Bool = true
        // find the hash, wedge and 2 plane atoms
        var planeAtoms: [MKAtom] = []
        var wedgeAtoms: [MKAtom] = []
        var hashAtoms: [MKAtom] = []
        for bond in center.getBondIterator()! {
            let nbr = bond.getNbrAtom(center)
            // hash bonds
            if bond.isHash() {
                if bond.getBeginAtom().getId() == center.getId() {
                    // this is a 'real' hash bond going from the center to the neighbor
                    hashAtoms.append(nbr)
                } else {
                    // this is an 'inverted' hash bond going from the neighbor to the center
                    planeAtoms.append(nbr)
//                    if tiponly {
//                        planeAtoms.append(nbr)
//                    } else {
//                        wedgeAtoms.append(nbr)
//                    }
                } 
            } else if bond.isWedge() {
                // wedge bonds 
                if bond.getBeginAtom().getId() == center.getId() {
                    // this is a 'real' wedge bond going from the center to the neighbor
                    wedgeAtoms.append(nbr)
                } else {
                    // this is an 'inverted' wedge bond going from the neighbor to the center
                    planeAtoms.append(nbr)
//                    if tiponly {
//                        planeAtoms.append(nbr)
//                    } else {
//                        hashAtoms.append(nbr)
//                    }
                }
            } else if bond.isWedgeOrHash() {
                if !tiponly || (tiponly && bond.getBeginAtom().getId() == center.getId()) {
                    config.specified = true
                    config.winding = .UnknownWinding
                    break
                } else {
                    planeAtoms.append(nbr)
                }
            } else {
                // plane bonds
                planeAtoms.append(nbr)
            }
        }
        // Handle the case of a tet center with four plane atoms or
        //        3 plane atoms with the fourth bond implicit
        if planeAtoms.count == 4 || (planeAtoms.count == 3 && center.getExplicitDegree() == 3) {
            config.specified = false 
        }
        var success: Bool = true 
        if (!config.specified || (config.specified && config.winding == .UnknownWinding)) {
            // unspecified or specified as unknown 
            for nbr in center.getNbrAtomIterator()! {
                if config.from_or_towrds == .NoRef {
                    config.from_or_towrds = .from(nbr.getId().ref)
                } else {
                    config.refs.append(nbr.getId().ref)
                }
            }
            while config.refs.count < 3 {
                config.refs.append(.ImplicitRef)
            }
        } else {
            // config specified
            if hashAtoms.count == 4 || wedgeAtoms.count == 4 {
                success = false
            } else if (planeAtoms.count + hashAtoms.count + wedgeAtoms.count) == 4 {
                // Handle all explicit tetra with at least one stereobond
                var order: [MKAtom] = []
                // First of all, handle the case of three wedge (or three hash) and one other bond
                //          by converting it into a single hash (or single wedge) and three planes
                if wedgeAtoms.count == 3 || hashAtoms.count == 3 {
                    var pwedge: [MKAtom] = []
                    var phash: [MKAtom] = []
                    if wedgeAtoms.count == 3 {
                        pwedge = wedgeAtoms
                        phash = hashAtoms
                    } else {
                        phash = wedgeAtoms
                        pwedge = hashAtoms
                    }
                    if planeAtoms.count == 0 { // Already has the hash bond
                        planeAtoms.append(contentsOf: pwedge)
                        pwedge.removeAll()
                    } else { // Does not already have the hash bond
                        phash.append(planeAtoms[0])
                        planeAtoms.removeAll()
                        pwedge.removeAll()
                    }
                }
                var wedge: Bool = wedgeAtoms.count > 0
                order.append(wedge ? wedgeAtoms[0] : hashAtoms[0])
                var nbrs: [MKAtom] = []
                for nbr in center.getNbrAtomIterator()! {
                    if nbr.getId().ref != order[0].getId().ref {
                        nbrs.append(nbr)
                    }
                }
                // Add "nbrs" to "order" in order of anticlockwise stereo
                order.append(nbrs[0])
                if angleOrder(order[0].getVector(), order[1].getVector(), nbrs[1].getVector(), center.getVector()) {
                    order.append(nbrs[1])
                } else {
                    order.insert(nbrs[1], at: 1)
                }
                if angleOrder(order[0].getVector(), order[2].getVector(), nbrs[2].getVector(), center.getVector()) {
                    order.append(nbrs[2])
                } else {
                    if angleOrder(order[0].getVector(), order[1].getVector(), nbrs[2].getVector(), center.getVector()) {
                        order.insert(nbrs[2], at: 2)
                    } else {
                        order.insert(nbrs[2], at: 1)
                    }
                }
                // Handle the case of two planes with a wedge and hash bond opposite each other.
                // This is handled as in the InChI TechMan (Figure 9) by marking it ambiguous if
                // the (small) angle between the plane bonds is > 133, and basing the stereo on
                // the 'inner' bond otherwise. This is commonly used for stereo in rings.
                // See also Get2DTetrahedralAmbiguity() in ichister.c (part of InChI)
                if planeAtoms.count == 2 && wedgeAtoms.count == 1 {  // Two planes, 1 wedge, 1 hash
                    if order[2] == hashAtoms[0] { // The wedge and hash are opposite
                        let angle = getAngle(order[1], center, order[3]) // The anticlockwise angle between the plane atoms
                        if angle > -133 && angle < 133 { // This value is from the InChI TechMan Figure 9
                            if angle > 0 { // Change to three planes and the hash bond
                                order.rotate(toStartAt: 2) // Change the order so that it begins with the hash bond
                                //TODO: maybe replace this with a search index and then rotate to there to make sure it is really the hash bond
                                wedge = false
                                planeAtoms.append(wedgeAtoms[0])
                                wedgeAtoms.removeAll()
                            } else { // change to three planes and the wedge bond (note: order is already correct)
                                planeAtoms.append(hashAtoms[0])
                                hashAtoms.removeAll()
                            }
                        } // No need for "else" statement, as this will be picked up as ambiguous stereo below
                    }
                }

                config.from_or_towrds = .from(order[0].getId().ref)
                config.refs.reserveCapacity(3) // make sure this actually adds elements into the array 
                for i in 0..<3 {
                    config.refs[i] = order[i+1].getId().ref
                }
                if wedge {
                    config.winding = .AntiClockwise
                }
                // Check for ambiguous stereo based on the members of "order".
                // If the first is a wedge bond, then the next should be a plane/hash, then plane/wedge, then plane/hash
                // If not, then the stereo is considered ambiguous.

                var pwedge: [MKAtom] = []
                var phash: [MKAtom] = []
                if wedge {
                    pwedge = wedgeAtoms
                    phash = hashAtoms
                } else {
                    phash = wedgeAtoms
                    pwedge = hashAtoms
                }
                if pwedge.contains(order[1]) || phash.contains(order[2]) || pwedge.contains(order[3]) { // Ambiguous stereo
                    success = false
                }
            } else if (hashAtoms.count == 0 || wedgeAtoms.count == 0) { // 3 explicit bonds from here on
                // Composed of just wedge bonds and plane bonds, or just hash bonds and plane bonds
                // Pick a stereobond on which to base the stereochemistry:
                var order: [MKAtom] = []
                var wedge = wedgeAtoms.count > 0
                order.append(wedge ? wedgeAtoms[0] : hashAtoms[0])
                var nbrs: [MKAtom] = []
                for nbr in center.getNbrAtomIterator()! {
                    if nbr != order[0] {
                        nbrs.append(nbr)
                    }
                }
                // Add "nbrs" to "order" in order of anticlockwise stereo
                order.append(nbrs[0])
                if angleOrder(order[0].getVector(), order[1].getVector(), nbrs[1].getVector(), center.getVector()) {
                    order.append(nbrs[1])
                } else {
                    order.insert(nbrs[1], at: 1)
                }
                 // Handle the case of two planes with a wedge/hash in the small angle between them.
                // This is handled similar to the InChI TechMan (Figure 10) by treating the stereo bond
                // as being in the large angle. This is consistent with Symyx Draw.
                if planeAtoms.count == 2 { // two planes, one stereo
                    let angle = getAngle(order[1], center, order[2]) // The anticlockwise angle between the plane atoms
                    if angle < 0 { // INvert the stereo of the stereobond 
                        wedge = !wedge
                    }
                }
                config.from_or_towrds = .from(.ImplicitRef)
                config.refs.reserveCapacity(3) // make sure this actually adds elements into the array
                for i in 0..<3 {
                    // interestingly they only used the first 3 elements not the last three
                    config.refs[i] = order[i].getId().ref
                }
                if !wedge {
                    config.winding = .AntiClockwise
                }
            } else { // 3 explicit bonds with at least one hash and at least one wedge
                success = false
            }
        }  // end of config.specified

        if !success {
            print("Symmetry analysis found atom with id \(center.getId()) to be a tetrahedral atom but the wedge/hash bonds can't be interpreted.")
            print(" # in-plane bonds = \(planeAtoms.count)")
            print(" # wedge bonds = \(wedgeAtoms.count)")
            print(" # hash bonds = \(hashAtoms.count)")
            continue
        }
        let th = MKTetrahedralStereo(mol)
        th.setConfig(config)
        configs.append(th)
        // add the data to the molecule if needed
        if addToMol {
            mol.setData(th)
        }
    }

    return configs
}

/**
* Get a vector with all OBCisTransStereo objects for the molecule. This
* function is used by StereoFrom2D() with the @p addToMol parameter is set
* to true.
*
* This function is also used for symmetry analysis to handle cases where
* there are two atoms in the same symmetry class that don't have the same
* stereochemistry. In this situation, the @p addToMol parameter is set to
* false and the returned objects will need to be deleted explicitly.
*
* The algorithm for converting the 2D coordinates uses the same triangle
* sign as TetrahedralFrom2D(). Depending on sign of 2 triangles, the right
* OBStereo::Shape is selected.
@verbatim
   0      3
    \    /        2 triangles: 0-1-b & 2-3-a
     a==b    -->  same sign: U
    /    \        opposite sign: Z
   1      2
@endverbatim
*
* @param mol The molecule.
* @param stereoUnits The stereogenic units.
* @param updown A map of OBStereo::BondDirection for cis/trans bonds
* @param addToMol If true, the OBCisTransStereo objects will be added
* to the molecule using OBBase::SetData().
*
* @sa StereoFrom2D FindStereogenicUnits
* @since version 2.3
*/
@discardableResult
func cisTransFrom2D(_ mol: MKMol, _ stereoUnits: MKStereoUnitSet, _ updown: [MKBond: MKStereo.BondDirection]? = nil, _ addToMol: Bool = true) -> [MKCisTransStereo] {
    var configs: [MKCisTransStereo] = []

    // find all cis/trans bonds
    let bonds = stereoUnits.filter { $0.type == .CisTrans }.map { $0.id }
    for i in bonds {
        guard let bond = mol.getBondById(i) else {
            continue
        }
        let begin = bond.getBeginAtom()
        let end = bond.getEndAtom()

        // Create a vector with the coordinates of the neighbor atoms
        var bondVecs: [Vector<Double>] = []
        let config: MKCisTransStereo.Config = MKCisTransStereo.Config()
        config.specified = true

        // begin 
        config.begin = begin.getId().ref
        for nbr in begin.getNbrAtomIterator()! {
            if nbr.getId() == end.getId() {
                continue
            }
            config.refs.append(nbr.getId().ref)
            bondVecs.append(nbr.getVector())

            // Check whether a single bond with unknown dir starts at the dbl bond (tip-only convention)
            if let b = mol.getBond(begin, nbr) {
                if updown != nil {
                    let ud_cit = updown![b]
                    if ud_cit != nil && ud_cit == .UnknownDir && b.getBeginAtom() == begin {
                        config.specified = false
                    }
                }
            }
        }
        if config.refs.count == 1 {
            config.refs.append(.ImplicitRef)
            var pos = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
            pos = begin.getNewBondVector(1.0)
            bondVecs.append(pos)
        }
        // end
        config.end = end.getId().ref
        for nbr in end.getNbrAtomIterator()! {
            if nbr.getId() == begin.getId() {
                continue
            }
            config.refs.append(nbr.getId().ref)
            bondVecs.append(nbr.getVector())

            // Check whether a single bond with unknown dir starts at the dbl bond (tip-only convention)
            if let b = mol.getBond(end, nbr) {
                if updown != nil {
                    let ud_cit = updown![b]
                    if ud_cit != nil && ud_cit == .UnknownDir && b.getBeginAtom() == end {
                        config.specified = false
                    }
                }
            }
        }
        if config.refs.count == 3 {
            config.refs.append(.ImplicitRef)
            var pos = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
            pos = end.getNewBondVector(1.0)
            bondVecs.append(pos)
        }

        // hande the case where the dbl bond is marked as unknown stereo
        if updown != nil {
            let ud_cit = updown![bond]
            if ud_cit != nil && ud_cit == .UnknownDir {
                config.specified = false
            }
        }
        if config.specified == true { // Work out the stereochemistry
            // 0      3
            //  \    /        2 triangles: 0-1-b & 2-3-a
            //   a==b    -->  same sign: U
            //  /    \        opposite sign: Z
            // 1      2
            /*
            double sign1 = TriangleSign(begin->GetVector(), end->GetVector(), bondVecs[0]);
            double sign2 = TriangleSign(begin->GetVector(), end->GetVector(), bondVecs[2]);
            */
            let sign1 = triangleSign(bondVecs[0], bondVecs[1], end.getVector())
            let sign2 = triangleSign(bondVecs[2], bondVecs[3], begin.getVector())
            let sign = sign1 * sign2

            if sign < 0.0 { // opposite sign
                config.shape = .ShapeZ
            }
        }
        let ct = MKCisTransStereo(mol)
        ct.setConfig(config)
        configs.append(ct)
        // add to molecule if requested
        if addToMol {
            mol.setData(ct)
        }
    }
    return configs
}
/**
* Convert a molecule's OBTetrahedralStereo objects to a series of hash or
* wedge bonds. Note that the molecule itself is not modified; the result
* is returned in the maps @p updown and @p from, which indicate
* the origin and direction of each hash or wedge bond.
*
* When converting, the following guidelines are followed when trying to
* find the best candidate bond to set up/down for each OBTetrahedralStereo
* object:
* -# Should not already be set
* -# Should not be connected to a 2nd tet center
*    (this is acceptable in theory as the wedge is only at one end, but
*     in practice it may cause confusion and thus we avoid it)
* -# Preferably is not in a cycle
* -# Preferably is a terminal H
*
* If no bond can be found that matches rules 1 and 2 (and in theory this is possible)
* then an error message is logged and the function returns false. (If you find an
* example where this occurs, please file a bug.)
*
* @param mol The molecule.
* @param updown A map of OBStereo::BondDirection for each hash/wedge bond
* @param from A map of OBStereo::Ref indicating the origin of each hash/wedge bond
* @return True or False depending on whether the conversion was successful
* @since version 2.3
*/

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
//
//  From0D
//
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

/**
* Add missing OBStereo objects. Unlike StereoFrom3D() and StereoFrom2D(), this
* method only adds objects for previously unidentified objects since we
* don't want to loose any information. The Config::specified flag for the
* newly added structs is always set to false.
*
* For example, a smiles is read which has two tetrahedral centers. Only one has
* stereochemisrty specified using a '@' character. StereoFrom0D() will detect the
* second tetrahedral atom and add an OBTetrahedralStereo object to the molecule.
*
* @param mol The molecule.
*
* @sa StereoFrom3D StereoFrom2D PerceiveStereo
* @since version 2.3
*/
func stereoFrom0D(_ mol: MKMol) {
    if mol.hasChiralityPerceived() {
        return
    }
    var symmetry_classes = findSymmetry(mol).map { UInt($0.intValue!) }
    let stereogenicUnits = findStereogenicUnits(mol, symClasses: &symmetry_classes)
    tetrahedralFrom0D(mol, stereogenicUnits)
    cisTransFrom0D(mol, stereogenicUnits)
    mol.setChiralityPerceived()
}

/**
* Get a vector with all OBTetrahedralStereo objects for the molecule. This
* function is used by StereoFrom0D() with the @p addToMol parameter is set
* to true. There is no algorithm used here, all specified flags will be
* set to false.
*
* This function is also used for symmetry analysis to handle cases where
* there are two atoms in the same symmetry class that don't have the same
* stereochemistry. In this situation, the @p addToMol parameter is set to
* false and the returned objects will need to be deleted explicitly.
*
* @param mol The molecule.
* @param stereoUnits The stereogenic units.
* @param addToMol If true, the OBTetrahedralStereo objects will be added
* to the molecule using OBBase::SetData().
*
* @sa StereoFrom0D FindStereogenicUnits
* @since version 2.3
*/
@discardableResult
func tetrahedralFrom0D(_ mol: MKMol, _ stereoUnits: MKStereoUnitSet, _ addToMol: Bool = true) -> [MKTetrahedralStereo] {
    var configs: [MKTetrahedralStereo] = []
    
    // Delete any existing stereo objects that are not a member of 'centers'
    // and make a map of the remaining ones
    var existingMap: OrderedDictionary<Int, MKTetrahedralStereo> = [:]    
    let stereoData: [MKGenericData] = mol.getDataVector(.StereoData)!
    for data in stereoData {
        if (data as! MKStereoBase).getType() == .Tetrahedral {
            let ts = data as! MKTetrahedralStereo
            let center = ts.getConfig().center
            // check if the center is really stereogenic
            var isStereogenic = false
            for u in stereoUnits {
                if u.type == .Tetrahedral {
                    if u.id == center {
                        isStereogenic = true
                    }
                }
            }
            if isStereogenic {
                existingMap[center.intValue!] = ts
                configs.append(ts)
            } else {
                // According to OpenBabel, this is not a tetrahedral stereo
                print("Removed spurious TetrahedralStereo object")
                mol.deleteData(ts)
            }
        }
    }
    for u in stereoUnits {
        if u.type != .Tetrahedral {
            continue
        }
        if existingMap[u.id.intValue!] != nil {
            continue
        }
        let center = mol.getAtomById(u.id)
        var config = MKTetrahedralStereo.Config()
        config.specified = false
        config.center = u.id
        for nbr in center!.getNbrAtomIterator()! {
            if config.from_or_towrds.refValue == .NoRef {
                config.from_or_towrds = .from(nbr.getId().ref)
            } else {
                config.refs.append(nbr.getId().ref)
            }
        }
        if config.refs.count == 2 {
            config.refs.append(.ImplicitRef) // need to add largest number on end to work
            // TODO: Make sure this doesn't break logic down the line by having an implicit ref involved
        }
        let th = MKTetrahedralStereo(mol)
        th.setConfig(config)
        configs.append(th)
        // add the data to the molecule if needed
        if addToMol {
            mol.setData(th)
        }
    }
    
    return configs
}

/**
* Get a vector with all OBCisTransStereo objects for the molecule. This
* function is used by StereoFrom0D() with the @p addToMol parameter is set
* to true. There is no algorithm used here, all specified flags will be
* set to false.
*
* This function is also used for symmetry analysis to handle cases where
* there are two atoms in the same symmetry class that don't have the same
* stereochemistry. In this situation, the @p addToMol parameter is set to
* false and the returned objects will need to be deleted explicitly.
*
* @param mol The molecule.
* @param stereoUnits The stereogenic units.
* @param addToMol If true, the OBCisTransStereo objects will be added
* to the molecule using OBBase::SetData().
*
* @sa StereoFrom0D FindStereogenicUnits
* @since version 2.3
*/
@discardableResult
func cisTransFrom0D(_ mol: MKMol, _ stereoUnits: MKStereoUnitSet, _ addToMol: Bool = true) -> [MKCisTransStereo] {
    var configs: [MKCisTransStereo] = []
    
    let bonds = stereoUnits.filter { $0.type == .CisTrans }.map { $0.id }
    // Delete any existing stereo objects that are not a member of 'bonds'
    // and make a map of the remaining ones
    var existingMap: OrderedDictionary<Int, MKCisTransStereo> = [:]    
    let stereoData: [MKGenericData] = mol.getDataVector(.StereoData)!
    for data in stereoData {
        if (data as! MKStereoBase).getType() == .CisTrans {
            let ct = data as! MKCisTransStereo
            let config = ct.getConfig()
            // find the bond iff from begin & end atom ids 
            var id: Ref = .NoRef
            let a = mol.getAtomById(config.begin)
            if a == nil {
                continue 
            }
            for bond in a!.getBondIterator()! {
                let beginID = bond.getBeginAtom().getId().ref
                let endID = bond.getEndAtom().getId().ref
                if beginID == config.begin && endID == config.end ||
                    beginID == config.end && endID == config.begin {
                    id = bond.getId()
                    break
                }
            }

            if !bonds.contains(id) {
                // According to OpenBabel, this is not a cis trans stereo
                print("Removed spurious CisTransStereo object")
                mol.deleteData(ct)
            } else {
                existingMap[id.intValue!] = ct
                configs.append(ct)
            }
        }
    }

    for i in bonds {
        // If there already exists a OBCisTransStereo object for this
        // bond, leave it alone unless it's in a ring of small size
        let alreadyExists = existingMap.contains { (key: Int, value: MKCisTransStereo) in
            key == i.intValue!
        }
        guard let bond = mol.getBondById(i) else {
            continue
        }
        var ct: MKCisTransStereo
        var config: MKCisTransStereo.Config = MKCisTransStereo.Config()
        if alreadyExists {
            ct = existingMap[i.intValue!]!
            config = ct.getConfig()
        } else {
            let begin = bond.getBeginAtom()
            let end = bond.getEndAtom()
            config.specified = false
            config.begin = begin.getId().ref
            for nbr in begin.getNbrAtomIterator()! {
                if nbr.getId() == end.getId() {
                    continue
                }
                config.refs.append(nbr.getId().ref)
            }
            if config.refs.count == 1 {
                config.refs.append(.ImplicitRef)
            }
            config.end = end.getId().ref
            for nbr in end.getNbrAtomIterator()! {
                if nbr.getId() == begin.getId() {
                    continue
                }
                config.refs.append(nbr.getId().ref)
            }
            if config.refs.count == 3 {
                config.refs.append(.ImplicitRef)
            }
            ct = MKCisTransStereo(mol)
            ct.setConfig(config)
        }

        guard let ring = bond.findSmallestRing() else {
            continue
        }
        if ring.size() <= IMPLICIT_CIS_RING_SIZE {
            var ringRefs: [Ref] = []
            for i in 0..<2 {
                if config.refs[i*2] != .ImplicitRef && ring.isMember(mol.getAtomById(config.refs[i*2])!) {
                    ringRefs.append(config.refs[i*2])
                } else {
                    ringRefs.append(config.refs[i*2 + 1])
                }
            }
            if !ct.isCis(ringRefs[0], ringRefs[1]) {
                config.shape = .ShapeZ
            }
            config.specified = true
            ct.setConfig(config)
        }
        configs.append(ct)
        // add the data to the molecule if needed
        if addToMol && !alreadyExists {
            mol.setData(ct)
        }
    }
    return configs
}


func tetStereoToWedgeHash(_ mol: MKMol, _ updown: inout [MKBond: MKStereo.BondDirection], _ from: inout [MKBond: Ref]) -> Bool {
    // Store the tetcenters for the second loop (below)
    var tetcenters: Set<Ref> = []
    guard let vdata: [MKGenericData] = mol.getDataVector(.StereoData) else {
        fatalError("No stereo data to work with???")
    }
    for data in vdata {
        if (data as! MKStereoBase).getType() == .Tetrahedral {
            let ts = data as! MKTetrahedralStereo
            let cfg = ts.getConfig()
            tetcenters.insert(cfg.center)
        }
    }
    // This loop sets one bond of each tet stereo to up or to down (2D only)
    var alreadyset: Set<MKBond> = []
    let uc = mol.getData(.UnitCell) as? MKUnitCell
    for data in vdata {
        if (data as! MKStereoBase).getType() == .Tetrahedral {
            let ts = data as! MKTetrahedralStereo
            var cfg: MKTetrahedralStereo.Config = ts.getConfig()
            if cfg.specified {
                var chosen: MKBond? = nil
                guard let center: MKAtom = mol.getAtomById(cfg.center) else { continue }
                var center_coord: Vector<Double> = center.getVector()
                // Find the two bonds closest in angle and remember them if
                // they are closer than DELTA_ANGLE_FOR_OVERLAPPING_BONDS
                var nbrs: [MKAtom] = []
                for a in center.getNbrAtomIterator()! {
                    nbrs.append(a)
                }
                var min_angle: Double = 359.0
                var close_bond_a: MKBond? = nil
                var close_bond_b: MKBond? = nil
                for i in 0..<nbrs.count - 1 {
                    for j in i+1..<nbrs.count {
                        var angle = abs(nbrs[i].getAngle(center, nbrs[j]))
                        if angle < min_angle {
                            min_angle = angle
                            close_bond_a = mol.getBond(center, nbrs[i])!
                            close_bond_b = mol.getBond(center, nbrs[j])!
                        }
                    }
                }
                if min_angle > DELTA_ANGLE_FOR_OVERLAPPING_BONDS {
                    close_bond_a = nil
                    close_bond_b = nil
                }
                
                // Find the best candidate bond to set to up/down
                // 1. **Should not already be set**
                // 2. Should not be connected to a 2nd tet center
                //    (this is acceptable, as the wedge is only at one end, but will only confuse things)
                // 3. Preferably is not in a cycle
                // 4. Prefer neighbor with fewer bonds over neighbor with more bonds
                // 5. Preferably is a terminal H, C, or heteroatom (in that order)
                // 6. If two bonds are overlapping, choose one of these
                //    (otherwise the InChI code will mark it as ambiguous)
                var max_bond_score: Int = 0 // The test below (score > max_bond_score)
                // gave incorrect results when score < 0 and max_bond_score was an unsigned int
                // see https://stackoverflow.com/questions/5416414/signed-unsigned-comparisons#5416498
                for b in center.getBondIterator()! {
                    if alreadyset.contains(b) {
                        continue
                    }
                    let nbr = b.getNbrAtom(center)
                    var nbr_bonds: Int = nbr.getExplicitDegree()
                    var score: Int = 0
                    if !b.isInRing() {
                        if !nbr.isInRing() {
                            score += 8 // non-ring bond to a non-ring atom is good
                        } else {
                            score += 2 // non-ring bond to ring atom is bad
                        }
                    }
                    if !tetcenters.contains(nbr.getId().ref) {
                        score += 4
                    }
                    if nbr_bonds == 1 { // terminal atom...
                        score += 8 // strongly prefer terminal atoms
                        // TODO: oh how might this come back to haunt us
                    } else {
                        score -= nbr_bonds - 2 // bond to atom with many bonds is penalized
                    }
                    if nbr.getAtomicNum() == MKElements.Hydrogen.atomicNum {
                        score += 2 // prefer H
                    } else if nbr.getAtomicNum() == MKElements.Carbon.atomicNum {
                        score += 1 // then C
                    }
                    if (b == close_bond_a || b == close_bond_b) {
                        score += 16
                    }
                    if (score > max_bond_score) {
                        max_bond_score = score
                        chosen = b
                    }
                }
                guard chosen != nil else { // There is a remote possibility of this but let's worry about 99.9% of cases first
                    print("Failed to set stereochemistry as unable to find an available bond")
                    return false
                }
                alreadyset.insert(chosen!)
                // Determine whether this bond should be set hash or wedge (or indeed unknown)
                // (Code inspired by perception.cpp, TetrahedralFrom2D: plane1 + plane2 + plane3, wedge)
                var bonddir: MKStereo.BondDirection = .UnknownDir
                if cfg.winding != .UnknownWinding {
                    var test_cfg: MKTetrahedralStereo.Config = cfg
                    // If there is an implicit ref; let's make that the 'from' atom
                    // otherwise use the atom on the chosen bond
                    var implicit: Bool = false
                    if case let .from(ref) = test_cfg.from_or_towrds, ref != .ImplicitRef {
                        if test_cfg.refs.contains(where: { refVal in
                            refVal == .ImplicitRef
                        }) {
                            test_cfg = MKTetrahedralStereo.toConfig(test_cfg, .from(.ImplicitRef))
                            implicit = true
                        }
                    } else {
                        implicit = true
                    }
                    var anticlockwise_order: Bool = false
                    var useup: Bool = false
                    if implicit {
                        // Put the ref for the stereo bond second
                        while test_cfg.refs[1] != chosen?.getNbrAtom(center).getId().ref {
                            test_cfg.refs.rotate(toStartAt: 2)
                        }
                        if uc != nil {
                            anticlockwise_order = angleOrder(
                                uc!.unwrapCartesianNear(mol.getAtomById(test_cfg.refs[0])!.getVector(), center_coord),
                                uc!.unwrapCartesianNear(mol.getAtomById(test_cfg.refs[1])!.getVector(), center_coord),
                                uc!.unwrapCartesianNear(mol.getAtomById(test_cfg.refs[2])!.getVector(), center_coord),
                                center_coord)
                        } else {
                            anticlockwise_order = angleOrder(
                                mol.getAtomById(test_cfg.refs[0])!.getVector(),
                                mol.getAtomById(test_cfg.refs[1])!.getVector(),
                                mol.getAtomById(test_cfg.refs[2])!.getVector(),
                                center.getVector())
                        }
                        // Get the angle between the plane bonds
                        let angle = getAngle(mol.getAtomById(test_cfg.refs[0])!, center, mol.getAtomById(test_cfg.refs[2])!)
                        if ((angle < 0 && anticlockwise_order) || (angle > 0 && !anticlockwise_order)) { // Is the stereobond in the bigger angle?
                            // If the bonds are in anticlockwise order, a clockwise angle (<180) between plane bonds
                            // implies that the stereo bond is in the bigger angle. Otherwise it has the opposite meaning.
                            useup = anticlockwise_order
                        } else {
                            useup = !anticlockwise_order
                        }
                    } else {
                        test_cfg = MKTetrahedralStereo.toConfig(test_cfg, .from(chosen!.getNbrAtom(center).getId().ref))
                        if uc != nil {
                            anticlockwise_order = angleOrder(
                                uc!.unwrapCartesianNear(mol.getAtomById(test_cfg.refs[0])!.getVector(), center_coord),
                                uc!.unwrapCartesianNear(mol.getAtomById(test_cfg.refs[1])!.getVector(), center_coord),
                                uc!.unwrapCartesianNear(mol.getAtomById(test_cfg.refs[2])!.getVector(), center_coord),
                                center_coord)
                        } else {
                            anticlockwise_order = angleOrder(
                                mol.getAtomById(test_cfg.refs[0])!.getVector(),
                                mol.getAtomById(test_cfg.refs[1])!.getVector(),
                                mol.getAtomById(test_cfg.refs[2])!.getVector(),
                                center.getVector())
                        }
                        if anticlockwise_order {
                            useup = false 
                        } else {
                            useup = true
                        }
                    }
                    // Set to UpBond (filled wedge from cfg.center to chosen_nbr) or DownBond
                    bonddir = useup ? .UpBond : .DownBond
                }
                updown[chosen!] = bonddir
                from[chosen!] = cfg.center
            }
        }
    }
    
    return true
}
/**
* Return a set of double bonds corresponding to the OBCisTransStereo objects
* for which the stereochemistry is undefined.
*
* Note that this functions just iterates over the existing OBCisTransStereo
* objects - it does not try to identify new ones.
*
* @param mol The molecule
* @return A set of bonds with unspecified cis/trans stereochemistry
* @since version 2.3
*/
func getUnspecifiedCisTrans(_ mol: MKMol) -> Set<MKBond> {
    // get double bonds with unspecified CisTrans stereochemistry
    var unspec_ctstereo: Set<MKBond> = []
    guard let vdata: [MKGenericData] = mol.getDataVector(.StereoData) else {
        fatalError("No stereo data to work with???")
        // TODO: maybe this should just return
    }
    for data in vdata {
        if (data as! MKStereoBase).getType() == .CisTrans {
            let ct = data as! MKCisTransStereo
            let config = ct.getConfig()
            if !config.specified {
                guard let bond = mol.getBond(mol.getAtomById(config.begin)!, mol.getAtomById(config.end)!) else {
                    //  TODO: throw error? 
                    print("ERROR: bond is missing from stereo data")
                    continue 
                }
                unspec_ctstereo.insert(bond)
            }
        }
    }
    return unspec_ctstereo
}



/**
* @page Stereochemistry
* @section overview Overview of classes
*
* There are many molecules which contain stereogenic elements. However,
* certain cases (i.e. tetrahedral, cis/trans) are more common than others
* (i.e. allene, biphenyl, octrahedral, ...). For the common stereogenic
* units, classes are provided. The inheritance of these classes resembles
* the way they are split into groups.
*
* - OBStereoBase
*   - OBTetraNonPlanarStereo
*     - OBTetrahedralStereo
*     - OBExtendedTetrahedralStereo
*   - OBTetraPlanarStereo
*     - OBCisTransStereo
*     - OBExtendedCisTransStereo
*     - OBSquarePlanarStereo
*   - OBAxialStereo
*     - OBTrigonalBipyrimidalStereo
*     - OBOctahedralStereo
*
* @image html tetranonplanar.png
* @image html tetraplanar.png
*
* All specific classes (i.e. OBTetrahedralStereo, ...) have embedded Config
* structs which define the actual stereochemistry. All these Config structs
* use OBStereo::Ref values to reference or uniquely identify atoms. Make sure
* to read about OBStereo::Ref and the related functions (in OBStereo). OBStereo
* is also a placeholder for various enums with predefined values for parameters
* etc. These enums are used throughout the different stereo classes but having
* these enums in a single location makes it easier to remember. When working
* with stereo classes, you normally don't need to use any of the parent classes
* directly. Only OBStereo and the specific class are needed.
*
* @section usage Basic usage
*
* The OBStereoFacade hides the complexity of working with stereochemistry. When
* using openbabel as a library, this is by far the easiest way to access
* stereochemistry information.
* The header for the specific OBStereo::Type type is all you need to include.
* These are:
* - @em openbabel/stereo/tetrahedral.h
* - @em openbabel/stereo/cistrans.h
* - @em openbabel/stereo/squareplanar.h
*
* All these headers also include @em openbabel/stereo/stereo.h providing
* declarations for OBStereo & OBStereoFacade.
*
    @code
    #include <iostream>
    #include <openbabel/mol.h>
    #include <openbabel/obconversion.h>

    #include <openbabel/stereo/tetrahedral.h>

    using namespace OpenBabel;

    int main()
    {
    OBMol mol;
    OBConversion conv;
    conv.SetInFormat("smi");
    conv.ReadString(&mol, "C[C@H](Cl)Br");

    OBStereoFacade facade(&mol);

    FOR_ATOMS_OF_MOL(atom, mol) {
        if (facade.HasTetrahedralStereo(atom->GetId()))
        std::cout << facade.GetTetrahedralStereo(atom->GetId()) << std::endl;
    }
    }
    @endcode
*
* All specific stereo classes and their embedded Config struct have an
* operator<< function which allows them to be used with std::ostream objects
* (e.g. std::cout, std::err, ...). These functions are often useful when
* debugging code.
*
* @section details Details on implementation
*
* The detection of stereogenic units start with symmetry analysis. However, a
* complete symmetry analysis also needs to take stereochemistry into account.
* In practice, this means stereochemistry will be found iteratively. At each
* iteration, the current atom symmetry classes are used to identify stereogenic
* units. The details about how the symmetry classes are used depends on the type
* (OBStereo::Type) of stereogenic unit. For tetrahedral centers, having 3 heavy
* atom neighbors with different symmetry classes or 4 neighbors with different
* symmetry classes means the atom is chiral. See FindStereogenicUnits() for
* details.
*
* After identifying the stereogenic units, Config structs with all the
* information on the spacial arrangement of the groups still have to be
* created. This involves interpreting various ways to represent
* stereochemisrty:
*
* - 3D coordinates: StereoFrom3D()
* - 2D coordinates: StereoFrom2D()
* - 0D coordinates: StereoFrom0D()
*
* Both StereoFrom3D() and StereoFrom2D() delete all existing stereochemistry objects
* before adding new ones. For molecules with 3D coordinates, it is evident that
* all information is specified by the coordinates itself. However, if a file format
* uses stereo parity flags, Config structs must be constructed using lower level
* functions and StereoFrom3D() should not be called. In these cases information
* could be lost by calling StereoFrom3D() after reading the file (the stereo flag might have
* indicated the stereochemistry was unspecified or the flag might not match the
* coordinates). In the case of 2D molecules, the coordinates together with bond
* properties (OBBond::Hash, OBBond::Wedge, OBBond::WedgeOrHash and
* OBBond::CisOrTrans) define the stereochemistry. Again, lower level functions
* can be used when stereo flags need to be used.
*
* StereoFrom0D() works slightly different than 3D/2D. Here, deleting the
* stereochemistry would always result in lost information. Instead StereoFrom0D()
* only adds new objects for stereogenic units which were previously not found.
* For example, a smiles is read which has two tetrahedral centers. Only one has
* stereochemistry specified using a '@' character. StereoFrom0D() will detect the
* second tetrahedral atom and add an OBTetrahedralStereo object to the molecule.
* The Config::specified flag for the newly added structs is always set to false.
*
* Assuming the format code has correctly set the molecule dimensions (OBMol::GetDimesions),
* PerceiveStereo() will automatically select the correct function to call.
* When StereoFrom3D(), StereoFrom2D() or StereoFrom0D() are not used, make sure to always
* set OBMol::HasChiralityPerceived() before returning from the format's ReadMolecule().
*
*
* @section formats Guidelines for formats
* @subsection input Reading files
*
* - Read the section above
* - The MDL format (mdlformat.cpp) is a good example for 2D/3D formats with or
*   without parity flags.
* - The SMILES format (smilesformat.cpp) is a good example for 0D formats.
*
* @subsection output Writing files
*
* For many file formats no additional code is needed. For example, if a 3D format
* doesn't require stereo parity flags, writing the coordinates is enough. For 2D
* file formats it will often suffice to write the coordinates and bond properties.
* If parity flags are needed, the OBStereoFacade class can be used to retrieve the
* objects for all types of stereochemistry supported by the file format.
*
*
*
*
*
* @since version 2.3
*/


