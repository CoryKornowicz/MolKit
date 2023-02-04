//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 1/31/23.
//

import Foundation
import Surge

public enum MKGenericDataType: UInt {
    //! Unknown data type (default)
    case UndefinedData =      0

    //! Arbitrary key/value data, i.e., OBPairData
    case PairData      =      1

    //! Energetics data (e.g., total energy, heat of formation, etc.)
    case EnergyData    =      2

    //! Storing text comments (one per molecule, atom, bond, etc.) (for other data, e.g., author, keyword, ... use OBPairData)
    case CommentData   =      3

    //! Arbitrary information about conformers, i.e., OBConformerData
    case ConformerData =      4

    //! Bond data external to OpenBabel, i.e., OBExternalBond, OBExternalBondData
    case ExternalBondData =   5

    //! Information for generating & manipulating rotamers, i.e. OBRotamerList
    case RotamerList =        6

    //! Info. for storing bonds to atoms yet to be added, i.e. OBVirtualBond
    case VirtualBondData =    7

    //! Information on rings in a molecule, i.e., OBRingData
    case RingData =           8

    //! Information about torsion/dihedral angles, i.e., OBTorsionData and OBTorsion
    case TorsionData =        9

    //! Bond angles in a molecule, i.e., OBAngle, OBAngleData
    case AngleData =         10

    //! Residue serial numbers
    case SerialNums =        11

    //! Crystallographic unit cell data, i.e., OBUnitCell
    case UnitCell =          12

    //! Spin data, including NMR, atomic and molecular spin, etc.
    case SpinData =          13

    //! Arbitrary partial and total charges, dipole moments, etc.
    case ChargeData =        14

    //! Symmetry data -- point and space groups, transforms, etc. i.e., OBSymmetryData
    case SymmetryData =      15

    // 16 - Value unused, formerly ChiralData

    //! Atomic and molecular occupation data (e.g., for crystal structures)
    case OccupationData =    17

    //! Density (cube) data and surfaces
    case DensityData =       18

    //! Electronic levels, redox states, orbitals, etc.
    case ElectronicData =    19

    //! Vibrational modes, frequencies, etc.
    case VibrationData =     20

    //! Rotational energy information
    case RotationData =      21

    //! Nuclear transitions (e.g., decay, capture, fission, fusion)
    case NuclearData =       22

    //! Set Data (a set of OBGenericData)
    case SetData =           23

    //! Grid Data (e.g., 3D grids of data a.k.a. cubes)
    case GridData =          24

    //! Vector Data (i.e., one vector like a dipole moment)
    case VectorData =        25

    //! Matrix data (i.e., a 3x3 matrix like a rotation or quadrupole moment)
    case MatrixData =        26

    //! Stereochemistry data (see OBStereoBase)
    case StereoData =        27

    //! Density of States data (fermi energy and energy vs. density data)
    case DOSData =           28

    //! Electronic transition data (e.g., UV/Vis, excitation energies, etc.)
    case ElectronicTransitionData = 29

    // space for up to 2^14 more entries...

    //! Custom (user-defined data)
    case CustomData0 = 16384
    case CustomData1 = 16385
    case CustomData2 = 16386
    case CustomData3 = 16387
    case CustomData4 = 16388
    case CustomData5 = 16389
    case CustomData6 = 16390
    case CustomData7 = 16391
    case CustomData8 = 16392
    case CustomData9 = 16393
    case CustomData10 = 16394
    case CustomData11 = 16395
    case CustomData12 = 16396
    case CustomData13 = 16397
    case CustomData14 = 16398
    case CustomData15 = 16399
}

enum DataOrigin {
    case any                 //!< Undefined or unspecified (default)
    case fileformatInput     //!< Read from an input file
    case userInput           //!< Added by the user
    case perceived           //!< Perceived by Open Babel library methods
    case external            //!< Added by an external program
    case local                //!< Not for routine external use (e.g. in sdf or cml properties)
}


public class MKGenericData: NSObject {

    private var _attr: String
    private var _type: MKGenericDataType
    private var _source: DataOrigin

    var attr: String {
        get {
            return self._attr
        }
        set {
            self._attr = newValue
        }
    }

    var type: MKGenericDataType {
        get {
            return self._type
        }
        set {
            self._type = newValue
        }
    }
    
    var source: DataOrigin {
        get {
            return self._source
        }
        set {
            self._source = newValue
        }
    }

    init(_ attr: String?, _ type: MKGenericDataType?, _ source: DataOrigin?) {
        self._attr = attr ?? "Undefined"
        self._type = type ?? MKGenericDataType.UndefinedData
        self._source = source ?? DataOrigin.any 
    }

    static func == (lhs: MKGenericData, rhs: MKGenericData) -> Bool {
        return lhs.attr == rhs.attr && lhs.type == rhs.type && lhs.source == rhs.source
    }
    
}

public class MKCommentData: MKGenericData {
    
    var data: String = ""

    init() {
        super.init("Comment", MKGenericDataType.CommentData, DataOrigin.any)
    }

    init(_ src: MKCommentData) {
        super.init(src.attr, src.type, src.source)
        self.data = src.data
    }

    static func == (lhs: MKCommentData, rhs: MKCommentData) -> Bool {
        return lhs.data == rhs.data && lhs == rhs
    }

}

public class MKExternalBond: MKGenericData {
    
    var idx = 0
    var atom: MKAtom? = nil
    var bond: MKBond? = nil
    
    init(_ atom: MKAtom, _ bond: MKBond, _ idx: Int) {
        super.init(nil, nil, nil)
        self.atom = atom
        self.bond = bond
        self.idx = idx
    }

    static func == (lhs: MKExternalBond, rhs: MKExternalBond) -> Bool {
        return lhs.idx == rhs.idx && lhs.atom == rhs.atom && lhs.bond == rhs.bond
    }
    
}

public class MKExternalBondData: MKGenericData {
    
    var vexbonds: [MKExternalBond] = [MKExternalBond]()

    init() {
        super.init("ExternalBondData", MKGenericDataType.ExternalBondData, .perceived)
    }

    func setData(_ atom: MKAtom, _ bond: MKBond, _ idx: Int) {
        self.vexbonds.append(MKExternalBond(atom, bond, idx))
    }

}

public class MKPairData<T>: MKGenericData {

    var value: T? = nil

    init() {
        super.init("PairData", MKGenericDataType.PairData, .any)
    }
    
    init(_ t: T) {
        super.init("PairData", MKGenericDataType.PairData, .any)
        self.value = t
    }

}

public class MKSetData: MKGenericData {

    var _vdata: [MKGenericData] = [MKGenericData]()

    init() {
        super.init("SetData", MKGenericDataType.SetData, .any)
    }

    func addData(_ data: MKGenericData) {
        self._vdata.append(data)
    }

    func removeData(_ data: MKGenericData) {
        self._vdata.removeAll(where: { $0 == data })
    }

    static func == (lhs: MKSetData, rhs: MKSetData) -> Bool {
        return lhs._vdata == rhs._vdata && lhs == rhs
    }

}

public class MKVirtualBond: MKGenericData {
    
    var _begin: UInt = 0
    var _end: UInt = 0 
    var _order: UInt = 0
    var _stereo: Int = 0

    init() {
        super.init("VirtualBond", MKGenericDataType.VirtualBondData, .perceived)
    }

    init(_ begin: UInt, _ end: UInt, _ order: UInt, _ stereo: Int) {
        super.init("VirtualBond", MKGenericDataType.VirtualBondData, .perceived)
        self._begin = begin
        self._end = end
        self._order = order
        self._stereo = stereo
    }

    static func == (lhs: MKVirtualBond, rhs: MKVirtualBond) -> Bool {
        return lhs._begin == rhs._begin && lhs._end == rhs._end && lhs._order == rhs._order && lhs._stereo == rhs._stereo && lhs == rhs
    }

}

public class MKRingData: MKGenericData {

    var _vr: [MKRing] = [MKRing]()

    init() {
        super.init("RingData", MKGenericDataType.RingData, .perceived)
    }

    func addRing(_ ring: MKRing) {
        self._vr.append(ring)
    }

    func removeRing(_ ring: MKRing) {
        self._vr.removeAll(where: { $0 == ring })
    }

    func getRingIterator() -> MKIterator<MKRing> {
        return MKIterator<MKRing>(self._vr)
    }

    static func == (lhs: MKRingData, rhs: MKRingData) -> Bool {
        return lhs._vr == rhs._vr && lhs == rhs
    }

}

public enum LatticeType {
    case Undefined
    case Triclinic
    case Monoclinic
    case Orthorhombic
    case Tetragonal
    case Rhombohedral
    case Hexagonal
    case Cubic
} 

public class MKUnitCell<Scalar: ExpressibleByFloatLiteral & FloatingPoint & BinaryFloatingPoint>: MKGenericData {

    var _mOrtho = Matrix<Scalar>(rows: 3, columns: 3, repeatedValue: 0.0)
    var _mOrient = Matrix<Scalar>(rows: 3, columns: 3, repeatedValue: 0.0)
    var _offset = Vector<Scalar>.init(dimensions: 3, repeatedValue: 0.0)
    var _spaceGroupName: String = ""
    var _spaceGroup: MKSpaceGroup<Scalar>? = nil
    var _lattice: LatticeType = .Undefined
    
    public init() {
        super.init("UnitCell", MKGenericDataType.UnitCell, .fileformatInput)
    }
    
    public func setData(a: Double, b: Double, c: Double, alpha: Double, beta: Double, gamme: Double) {
        
    }
    
    public func setData(v1: Vector<Scalar>, v2: Vector<Scalar>, v3: Vector<Scalar>) {
        
    }

    public func setData(m: Matrix<Scalar>) {
        
    }
    
    public func setOffset(_ v: Vector<Scalar>) {
        
    }
    
    public func setSpaceGroup(_ spg: MKSpaceGroup<Scalar>) {
        self._spaceGroup = spg
    }
    
    func fillUnitCell(_ mol: MKMol) {
        return
    }


}
