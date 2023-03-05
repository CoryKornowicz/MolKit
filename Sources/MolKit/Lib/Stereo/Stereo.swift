
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
    var symmetry_classes = findSymmetry(mol)
}
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

}

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
    var symmetry_classes = findSymmetry(mol)
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
    fatalError()
}



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
    fatalError()
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
    fatalError()
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
    fatalError()
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
    fatalError()
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
func tetStereoToWedgeHash(_ mol: MKMol, _ updown: inout [MKBond: MKStereo.BondDirection], _ from: inout [MKBond: Ref]) -> Bool {
    fatalError()
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
    fatalError()
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


