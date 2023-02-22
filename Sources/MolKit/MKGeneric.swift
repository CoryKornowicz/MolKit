//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 1/31/23.
//

import Foundation
import OrderedCollections
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
    
    func setOrigin(_ origin: DataOrigin) {
        self.source = origin
    }
    
    func setAttribute(_ attr: String) {
        self.attr = attr
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

    func setValue(_ t: T) {
        self.value = t
    }

    func getValue() -> T? {
        return self.value
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

    func getBgn() -> UInt {
        return self._begin
    }

    func getEnd() -> UInt {
        return self._end
    }

    func getOrder() -> UInt {
        return self._order
    }

    func getStereo() -> Int {
        return self._stereo
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

    func getData() -> [MKRing] {
        return self._vr
    }

    func setData(_ data: [MKRing]) {
        self._vr = data
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

private let __spacegroups__: [String] = ["P1", "P-1", "P2", "P2(1)", "C2", "Pm", "Pc", "Cm", "Cc", "P2/m",
                                              "P2(1)/m", "C2/m", "P2/c", "P2(1)/c", "C2/c", "P222", "P222(1)",
                                              "P2(1)2(1)2", "P2(1)2(1)2(1)", "C222(1)", "C222", "F222", "I222",
                                              "I2(1)2(1)2(1)", "Pmm2", "Pmc2(1)", "Pcc2", "Pma2", "Pca2(1)", "Pnc2",
                                              "Pmn2(1)", "Pba2", "Pna2(1)", "Pnn2", "Cmm2", "Cmc2(1)", "Ccc2", "Amm2",
                                              "Abm2", "Ama2", "Aba2", "Fmm2", "Fdd2", "Imm2", "Iba2", "Ima2", "Pmmm",
                                              "Pnnn", "Pccm", "Pban", "Pmma", "Pnna", "Pmna", "Pcca", "Pbam", "Pccn",
                                              "Pbcm", "Pnnm", "Pmmn", "Pbcn", "Pbca", "Pnma", "Cmcm", "Cmca", "Cmmm",
                                              "Cccm", "Cmma", "Ccca", "Fmmm", "Fddd", "Immm", "Ibam", "Ibca", "Imma",
                                              "P4", "P4(1)", "P4(2)", "P4(3)", "I4", "I4(1)", "P-4", "I-4", "P4/m",
                                              "P4(2)/m", "P4/n", "P4(2)/n", "I4/m", "I4(1)/a", "P422", "P42(1)2",
                                              "P4(1)22", "P4(1)2(1)2", "P4(2)22", "P4(2)2(1)2", "P4(3)22", "P4(3)2(1)2",
                                              "I422", "I4(1)22", "P4mm", "P4bm", "P4(2)cm", "P4(2)nm", "P4cc", "P4nc",
                                              "P4(2)mc", "P4(2)bc", "I4mm", "I4cm", "I4(1)md", "I4(1)cd", "P-42m",
                                              "P-42c", "P-42(1)m", "P-42(1)c", "P-4m2", "P-4c2", "P-4b2", "P-4n2",
                                              "I-4m2", "I-4c2", "I-42m", "I-42d", "P4/mmm", "P4/mcc", "P4/nbm",
                                              "P4/nnc", "P4/mbm", "P4/mnc", "P4/nmm", "P4/ncc", "P4(2)/mmc",
                                              "P4(2)/mcm", "P4(2)/nbc", "P4(2)/nnm", "P4(2)/mbc", "P4(2)/mnm",
                                              "P4(2)/nmc", "P4(2)/ncm", "I4/mmm", "I4/mcm", "I4(1)/amd", "I4(1)/acd",
                                              "P3", "P3(1)", "P3(2)", "R3", "P-3", "R-3", "P312", "P321", "P3(1)12",
                                              "P3(1)21", "P3(2)12", "P3(2)21", "R32", "P3m1", "P31m", "P3c1", "P31c",
                                              "R3m", "R3c", "P-31m", "P-31c", "P-3m1", "P-3c1", "R-3m", "R-3c", "P6",
                                              "P6(1)", "P6(5)", "P6(2)", "P6(4)", "P6(3)", "P-6", "P6/m", "P6(3)/m",
                                              "P622", "P6(1)22", "P6(5)22", "P6(2)22", "P6(4)22", "P6(3)22", "P6mm",
                                              "P6cc", "P6(3)cm", "P6(3)mc", "P-6m2", "P-6c2", "P-62m", "P-62c",
                                              "P6/mmm", "P6/mcc", "P6(3)/mcm", "P6(3)/mmc", "P23", "F23", "I23",
                                              "P2(1)3", "I2(1)3", "Pm-3", "Pn-3", "Fm-3", "Fd-3", "Im-3", "Pa-3",
                                              "Ia-3", "P432", "P4(2)32", "F432", "F4(1)32", "I432", "P4(3)32",
                                              "P4(1)32", "I4(1)32", "P-43m", "F4-3m", "I-43m", "P-43n", "F-43c",
                                              "I-43d", "Pm-3m", "Pn-3n", "Pm-3n", "Pn-3m", "Fm-3m", "Fm-3c",
                                              "Fd-3m", "Fd-3c", "Im-3m", "Ia-3d"]

public class MKUnitCell: MKGenericData {

    var _mOrtho = Matrix<Double>(rows: 3, columns: 3, repeatedValue: 0.0)
    var _mOrient = Matrix<Double>(rows: 3, columns: 3, repeatedValue: 0.0)
    var _offset = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
    var _spaceGroupName: String = ""
    var _spaceGroup: MKSpaceGroup? = nil
    var _lattice: LatticeType = .Undefined
    
    public init() {
        super.init("UnitCell", MKGenericDataType.UnitCell, .fileformatInput)
    }
    
    func areDupliateAtoms(_ v1: Vector<Double>, _ v2: Vector<Double>) -> Bool {
        var dr = v2 - v1
        if dr[0] < -0.5 { dr[0] = dr[0] + 1 }
        if dr[0] > 0.5 { dr[0] = dr[0] - 1 }
        if dr[1] < -0.5 { dr[1] = dr[1] + 1 }
        if dr[1] > 0.5 { dr[1] = dr[1] - 1 }
        if dr[2] < -0.5 { dr[2] = dr[2] + 1 }
        if dr[2] > 0.5 { dr[2] = dr[2] - 1 }

        return (length(dr) < 1e-6)
    }
    
    //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#calculateOrthogonalisationMatrix">blue-obelisk:calculateOrthogonalisationMatrix</a>
    public func setData(A: Double, B: Double, C: Double, alpha: Double, beta: Double, gamma: Double) {
        self._mOrtho.fillOrtho(alpha: alpha, beta: beta, gamma: gamma, A: A, B: B, C: C)
        self._mOrient = Matrix<Double>.diagonal(rows: 3, columns: 3, repeatedValue: 1.0)
        self._spaceGroup = nil
        self._spaceGroupName = ""
        self._lattice = .Undefined
    }
    
    
    public func setData(v1: Vector<Double>, v2: Vector<Double>, v3: Vector<Double>) {
        let m: Matrix<Double> = Matrix<Double>.init([v1,v2,v3])
        self._mOrtho.fillOrtho(alpha: vector_angle(v2, v3),
                                beta: vector_angle(v1, v3), 
                                gamma: vector_angle(v1, v2), 
                                A: length(v1), B: length(v2), C: length(v3))
        self._mOrient = transpose(m) * inv(self._mOrtho)
        self._spaceGroup = nil
        self._spaceGroupName = ""
        self._lattice = .Undefined
    }

    public func setData(m: Matrix<Double>) {
        self.setData(v1: Vector(scalars: m[row: 0]), v2: Vector(scalars: m[row: 1]), v3: Vector(scalars: m[row: 2]))
    }
    
    public func setOffset(_ v: Vector<Double>) {
        self._offset = v
    }
    
    func setSpaceGroup(_ spg: MKSpaceGroup) {
        self._spaceGroup = spg
    }

    func setSpaceGroup(_ idx: Int) {
        self._spaceGroup = MolKit._SpaceGroups.getSpaceGroup(idx)
    }

    func setSpaceGroup(_ name: String) {
        self._spaceGroup = MolKit._SpaceGroups.getSpaceGroup(name)
        self._spaceGroupName = name
    }

    func setLatticeType(_ lt: LatticeType) {
        self._lattice = lt
    }
    
    func fillUnitCell(_ mol: MKMol) {

        guard let sg = self.getSpaceGroup() else {
            print("Warning SpaceGroup was not discovered in fillUnitCell")
            return
        }

        let atoms = mol.getAtomIterator()
        if mol.numAtoms() == 0 {
            print("Warning no atoms in molecule in fillUnitCell")
            return
        }

        var baseV: Vector<Double>
        var uniqueV: Vector<Double>
        var updatedCoordinate: Vector<Double>
        
        var atomsToKeep = [MKAtom]()
        var atomsToDelete = [MKAtom]()
        var coordinateSet: Set<String> = Set<String>()

        var hashString = "" 
        // For each atom, we loop through: convert the coords back to inverse space, apply the transformations and create new atoms
        for atom in atoms {
            baseV = atom.getVector()
            baseV = self.cartesianToFractional(baseV)
            baseV = self.wrapFractional(baseV)
            // "%03d,%.3f,%.3f,%.3f", atom->GetAtomicNum(), baseV.x(), baseV.y(), baseV.z());
            hashString = String(format: "%03d,%.3f,%.3f,%.3f", atom.getAtomicNum(), baseV[0], baseV[1], baseV[2])
            if coordinateSet.insert(hashString).inserted {
                atomsToKeep.append(atom)
            } else {
                atomsToDelete.append(atom)
            }
        }
        
        for at in atomsToDelete {
            mol.deleteAtom(at)
        }
        
        // Cross-check all transformations for duplicity
        for keep in atomsToKeep {
            uniqueV = keep.getVector()
            uniqueV = self.cartesianToFractional(uniqueV)
            uniqueV = self.wrapFractional(uniqueV)

            guard let transformedVectors = sg.transform(uniqueV) else { continue }
            for tv in transformedVectors {
                updatedCoordinate = self.wrapFractional(tv)
                // Check if the transformed coordinate is a duplicate of an atom
                // "%03d,%.3f,%.3f,%.3f", (*atomIter)->GetAtomicNum(), updatedCoordinate.x(),
                //  updatedCoordinate.y(), updatedCoordinate.z()
                hashString = String(format: "%03d,%.3f,%.3f,%.3f", keep.getAtomicNum(), updatedCoordinate[0], updatedCoordinate[1], updatedCoordinate[2])
                if coordinateSet.insert(hashString).inserted {
                    let newAtom = mol.newAtom()
                    newAtom.duplicate(keep)
                    newAtom.setVector(self.fractionalToCartesian(updatedCoordinate))
                }
            }
        }
        
        self.setSpaceGroup(1)
    }

    func getA() -> Double {
        return length(Vector(scalars: self._mOrtho[column: 0]))
    }

    func getB() -> Double {
        return length(Vector(scalars: self._mOrtho[column: 1]))
    }

    func getC() -> Double {
        return length(Vector(scalars: self._mOrtho[column: 2]))
    }

    func getAlpha() -> Double {
        return vector_angle(Vector(scalars: self._mOrtho[column: 1]), Vector(scalars: self._mOrtho[column: 2]))
    }

    func getBeta() -> Double {
        return vector_angle(Vector(scalars: self._mOrtho[column: 0]), Vector(scalars: self._mOrtho[column: 2]))
    }

    func getGamma() -> Double {
        return vector_angle(Vector(scalars: self._mOrtho[column: 0]), Vector(scalars: self._mOrtho[column: 1]))
    }

    func getOffset() -> Vector<Double> {
        return self._offset
    }

    //! \return the text representation of the space group for this unit cell
    func getSpaceGroup() -> MKSpaceGroup? {
        return self._spaceGroup
    }

    func getSpaceGroupName() -> String {
        return self._spaceGroupName
    }

    func getLatticeType() -> LatticeType {
        
        if self._lattice != .Undefined {
            return self._lattice
        } else if let group = self._spaceGroup {
            return self.getLatticeType(Int(group.m_id))
        }
        
        let a = self.getA()
        let b = self.getB()
        let c = self.getC()
        let alpha = self.getAlpha()
        let beta = self.getBeta()
        let gamma = self.getGamma()

        var rightAngles = 0
        if isApprox(alpha, to: 90.0, epsilon: 1.0e-3) { rightAngles += 1 }
        if isApprox(beta, to: 90.0, epsilon: 1.0e-3) { rightAngles += 1 }
        if isApprox(gamma, to: 90.0, epsilon: 1.0e-3) { rightAngles += 1 }

        switch rightAngles {
        case 3:
            if isApprox(a, to: b, epsilon: 1.0e-4) && isApprox(b, to: c, epsilon: 1.0e-4) {
                self._lattice = .Cubic
                return .Cubic
            } else if isApprox(a, to: b, epsilon: 1.0e-4) || isApprox(b, to: c, epsilon: 1.0e-4) {
                self._lattice = .Tetragonal
                return .Tetragonal
            } else {
                self._lattice = .Orthorhombic
                return .Orthorhombic
            }
        case 2:
            if isApprox(alpha, to: 120.0, epsilon: 1.0e-3) || isApprox(beta, to: 120.0, epsilon: 1.0e-3) || isApprox(gamma, to: 120.0, epsilon: 1.0e-3)
                    && isApprox(a, to: b, epsilon: 1.0e-4) || isApprox(b, to: c, epsilon: 1.0e-4) {
                self._lattice = .Hexagonal
                return .Hexagonal
            } else {
                self._lattice = .Monoclinic
                return .Monoclinic
            }
        case 1: 
            if isApprox(a, to: b, epsilon: 1.0e-4) && isApprox(b, to: c, epsilon: 1.0e-4) {
                self._lattice = .Rhombohedral
                return .Rhombohedral
            } else {
                self._lattice = .Triclinic
                return .Triclinic
            }
        default:
            break
        }

        return self._lattice
    }
    
    func getLatticeType(_ index: Int) -> LatticeType {
        
        //    1-2     Triclinic
        //    3-15    Monoclinic
        //    16-74   Orthorhombic
        //    75-142  Tetragonal
        //    143-167 Rhombohedral
        //    168-194 Hexagonal
        //    195-230 Cubic
        
        var spacegroup = index
        
        if spacegroup == 0 && self._spaceGroup != nil {
            spacegroup = Int(self._spaceGroup!.m_id)
        }
        
        if spacegroup <= 0 {
            return .Undefined
        } else if spacegroup == 1 || spacegroup == 2 {
            return .Triclinic
        } else if spacegroup >= 3 && spacegroup <= 15 {
            return .Monoclinic
        } else if spacegroup >= 16 && spacegroup <= 74 {
            return .Orthorhombic
        } else if spacegroup >= 75 && spacegroup <= 142 {
            return .Tetragonal
        } else if spacegroup >= 143 && spacegroup <= 167 {
            return .Rhombohedral
        } else if spacegroup >= 168 && spacegroup <= 194 {
            return .Hexagonal
        } else if spacegroup >= 195 && spacegroup <= 230 {
            return .Cubic
        } else {
            return .Undefined
        }
        
    }

    //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertNotionalIntoCartesianCoordinates">blue-obelisk:convertNotionalIntoCartesianCoordinates</a>
    func getCellVectors() -> [Vector<Double>] {
        
        var vectors: [Vector<Double>] = []
        
        let m = self.getCellMatrix()
        
        vectors.append(Vector.init(scalars: m[row: 0]))
        vectors.append(Vector.init(scalars: m[row: 1]))
        vectors.append(Vector.init(scalars: m[row: 2]))
        
        return vectors
    }

    func getCellMatrix() -> Matrix<Double> {
        return transpose(self._mOrient * self._mOrtho)
    }

    func getOrthoMatrix() -> Matrix<Double> {
        return self._mOrtho
    }

    func getOrientMatrix() -> Matrix<Double> {
        return self._mOrient
    }

    func getFractionalMatrix() -> Matrix<Double> {
        return inv(self._mOrtho)
    }

    func fractionalToCartesian(_ v: Vector<Double>) -> Vector<Double> {
        return self._mOrient * self._mOrient * v + self._offset
    }

    func cartesianToFractional(_ v: Vector<Double>) -> Vector<Double> {
        return inv(self._mOrtho) * inv(self._mOrient) * (v - self._offset)
    }

    func wrapCartesian(_ v: Vector<Double>) -> Vector<Double> {
        return self.fractionalToCartesian(self.wrapFractional(self.cartesianToFractional(v)))
    }

    private let LIMIT = 0.999999
    private let EPSILON = 1.0e-6
    
    func wrapFractional(_ v: Vector<Double>) -> Vector<Double> {
        var x: Double = fmod(v[0], 1.0)
        var y: Double = fmod(v[1], 1.0)
        var z: Double = fmod(v[2], 1.0)
        
        if x > LIMIT { x -= 1.0 }
        if y > LIMIT { y -= 1.0 }
        if z > LIMIT { z -= 1.0 }
        
        if x > 1 - EPSILON || x < EPSILON { x = 0.0 }
        if y > 1 - EPSILON || y < EPSILON { y = 0.0 }
        if z > 1 - EPSILON || z < EPSILON { z = 0.0 }
        
        return Vector(scalars: [x,y,z])
    }

    func unwrapCartesianNear(_ v: Vector<Double>, _ w: Vector<Double>) -> Vector<Double> {
        let bond_dir = self.minimumImageCartesian(v - w)
        return w + bond_dir
    }

    func unwrapFractionalNear(_ v: Vector<Double>, _ w: Vector<Double>) -> Vector<Double> {
        let bond_dir = self.minimumImageFractional(v - w)
        return w + bond_dir
    }

    func minimumImageCartesian(_ v: Vector<Double>) -> Vector<Double> {
        var frac = self.cartesianToFractional(v)
        frac = self.minimumImageFractional(frac)
        return self.fractionalToCartesian(frac)
    }

    func minimumImageFractional(_ v: Vector<Double>) -> Vector<Double> {
        let x = v[0] - round(v[0])
        let y = v[1] - round(v[1])
        let z = v[2] - round(v[2])
        return Vector(scalars: [x,y,z])
    }

    func getSpaceGroupNumber(_ inputName: String) -> Int {
        var name = inputName
        if name.length == 0 {
            if let group = self._spaceGroup {
                return Int(group.m_id)
            } else {
                name = self._spaceGroupName
            }
        }
        
        for (i, spacegroup) in __spacegroups__.enumerated() {
            if name == spacegroup {
                return i + 1
            }
        }
//        presumably never reached
        return 0
    }

    func getCellVolume() -> Double {
        return fabs(det(self.getCellMatrix()) ?? 0.0)
    }
    
}

//! \brief Used to hold data on conformers or geometry optimization steps
//!
//! Supplements the support for multiple coordinate sets in MKMol, e.g.,
public class MKConformerData: MKGenericData {
    var vDimension: [Int] = []
    var vEnergies: [Double] = []
    var vForces: [[Vector<Double>]] = []
    var vVelocity: [[Vector<Double>]] = []
    var vDisplace: [[Vector<Double>]] = []
    var vData: [String] = []

    public init() {
        super.init("Conformers", .ConformerData, .local)
    }
}

//! \brief Used to hold the point-group and/or space-group symmetry
//! \todo Add support for translation between symbol notations.
//!       Add symmetry perception routines.

public class MKSymmetryData: MKGenericData {

    var _spaceGroup: String = ""
    var _pointGroup: String = ""

    public init() {
        super.init("Symmetry", .SymmetryData, .local)
    }

}

//! \brief Used to hold the torsion data for a single rotatable bond
//! and all four atoms around it

public class MKTorsion {

    private var _bc: Pair<MKAtom?, MKAtom?> = (nil, nil)
    private var _ads: [Triple<MKAtom?, MKAtom?, Double?>] = []

    public init() {    }

    public init(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom, _ d: MKAtom) {
        self._ads.append(Triple(a, d, 0.0))
        self._bc = (b, c)
    }
    
    func clear() {
        self._ads.removeAll()
        self._bc = (nil, nil)
    }
    
    func empty() -> Bool {
        return self._bc.0 == nil && self._bc.1 == nil
    }
    
    /*!
     **\brief adds a new torsion to the OBTorsion object
     */
    func addTorsion(_ a: MKAtom, _ b: MKAtom, _ c: MKAtom, _ d: MKAtom) -> Bool {
        if !self.empty() && (b != self._bc.0 || c != self._bc.1) {
            return false
        }

        if self.empty() {
            self._bc.0 = b
            self._bc.1 = c
        }

        let ad = Triple(a, d, 0.0)
        self._ads.append(ad)
        
        return true
    }

    /*!
     **\brief adds a new torsion to the OBTorsion object
     */
    func addTorsion(_ atoms: Quad<MKAtom, MKAtom, MKAtom, MKAtom>) -> Bool {
        if !self.empty() && (atoms.1 != self._bc.0 || atoms.2 != self._bc.1) {
            return false
        }

        if self.empty() {
            self._bc.0 = atoms.1
            self._bc.1 = atoms.2
        }

        let ad = Triple(atoms.0, atoms.3, 0.0)
        self._ads.append(ad)

        return true
    }
    
    /*!
     **\brief Sets the angle of a torsion in OBTorsion
     **\param radians the value to assign to the torsion
     **\param index the index into the torsion of the OBTorsion
     **\return boolean success
     */
    func setAngle(_ radians: Double, _ index: Int) -> Bool {
        if index >= self._ads.count {
            return false
        }
        self._ads[index].2 = radians
        return true
    }

    /*!
     **\brief Obtains the angle of a torsion in OBTorsion
     **\param radians the value of the angle is set here
     **\param index the index into the torsion of the OBTorsion
     **\return boolean success
     */
    func getAngle(_ radians: inout Double, _ index: Int) -> Bool {
        if index >= self._ads.count {
            return false
        }
        radians = self._ads[index].2 ?? 0.0
        return true
    }

    //! Gets the bond index of the central bond
    //! \return int bond index
    func getBondIdx() -> Int {
        guard let first = self._bc.0, let second = self._bc.1 else { return 0 }
//        MARK: Potential error thrown here if bond does not exist between atom pair
        return Int(first.getBond(second)!.getIdx())
    }

    public func getSize() -> Int {
        return self._ads.count
    }

    public func getTorsions() -> [Quad<MKAtom?, MKAtom?, MKAtom?, MKAtom?>] {
        var abcd: Quad<MKAtom?, MKAtom?, MKAtom?, MKAtom?> = (nil, nil, nil, nil)

        abcd.1 = self._bc.0
        abcd.2 = self._bc.1

        var torsions: [Quad<MKAtom?, MKAtom?, MKAtom?, MKAtom?>] = []

        for ad in self._ads {
            abcd.0 = ad.0
            abcd.3 = ad.1
            torsions.append(abcd)
        }

        return torsions
    }

    //! Gets the two central atoms of ABCD torsion
    //!   \return pair<OBAtom*,OBAtom*>
    public func getBC() -> Pair<MKAtom?, MKAtom?> {
        return self._bc
    }

    //! Gets the vector of distal atoms of ABCD torsion
    //! \return vector of A,D atom pointers and a double
    public func getADs() -> [Triple<MKAtom?, MKAtom?, Double?>] {
        return self._ads
    }

    /*!
     **\brief determines if torsion has only protons on either the a or d end
     **\return boolean
     */
    public func isProtonRotor() -> Bool {
        var Aprotor: Bool = true
        var Dprotor: Bool = true

        for ad in self._ads {
            if ad.0?.getAtomicNum() != 1 {
                Aprotor = false
            }
            if ad.1?.getAtomicNum() != 1 {
                Dprotor = false
            }
        }

        return Aprotor || Dprotor
    }
}

public class MKTorsionData: MKGenericData {

    var _torsions: [MKTorsion] = []

    public init() {
        super.init("TorsionData", .TorsionData, .local)
    }

    public func clear() {
        self._torsions.removeAll()
    }

    public func getSize() -> Int {
        return self._torsions.count
    }

    public func setData(_ t: MKTorsion) {
        self._torsions.append(t)
    }
    
    /*!
    **\brief Fills a vector with the indices of the atoms in torsions (ordered abcd)
    **\param torsions reference to the vector of abcd atom sets
    **\return boolean success
    */
    func fillTorsionArray(_ tors: inout [[UInt]]) -> Bool {
        if self._torsions.isEmpty {
            return false
        }

        var tmpquads: [Quad<MKAtom?, MKAtom?, MKAtom?, MKAtom?>] = []
        var quads: [Quad<MKAtom?, MKAtom?, MKAtom?, MKAtom?>] = []

        for torsion in self._torsions {
            tmpquads = torsion.getTorsions()
            for thisQuad in tmpquads {
                quads.append(thisQuad)
            }
        }

        //fill array of torsion atoms
        
        tors.removeAll()
        tors = Array(repeating: [0, 0, 0, 0], count: quads.count)

        var ct: Int = 0

        for thisQuad in quads {
            tors[ct][0] = UInt(thisQuad.0!.getIdx() - 1)
            tors[ct][1] = UInt(thisQuad.1!.getIdx() - 1)
            tors[ct][2] = UInt(thisQuad.2!.getIdx() - 1)
            tors[ct][3] = UInt(thisQuad.3!.getIdx() - 1)
            ct += 1
        }

        return true
    }
}

public class MKAngle: Equatable {

    var _vertex: MKAtom?
    var _radians: Double = 0.0
    var _termini: Pair<MKAtom?, MKAtom?> = (nil, nil)

    public init() {
        self._vertex = nil
        self._radians = 0.0
        self._termini = (nil, nil)
    }

    public init(_ vertex: MKAtom?, _ a: MKAtom?, _ b: MKAtom?) {
        self._vertex = vertex
        self._termini.0 = a
        self._termini.1 = b
        self.sortByIndex()
    }

    public func getAtoms() -> Triple<MKAtom?, MKAtom?, MKAtom?> {
        return (self._vertex, self._termini.0, self._termini.1)
    }

    /*!
    **\brief sorts atoms in angle by order of indices
    */
    public func sortByIndex() {
        if self._termini.0!.getIdx() > self._termini.1!.getIdx() {
            let tmp = self._termini.0
            self._termini.0 = self._termini.1
            self._termini.1 = tmp
        }
    }


    public func clear() {
        self._vertex = nil
        self._radians = 0.0
        self._termini = (nil, nil)
    }

    //! Gets the OBAngle angle value
    //! \return angle in radians
    public func getAngle() -> Double {
        return self._radians
    }

    //! Sets the OBAngle to @p radians
    //! \param angle in radians
    public func setAngle(_ angle: Double) {
        self._radians = angle
    }

    public func setAtoms(_ vertex: MKAtom?, _ a: MKAtom?, _ b: MKAtom?) {
        self._vertex = vertex
        self._termini.0 = a
        self._termini.1 = b
        self.sortByIndex()
    }

    public func setAtoms(_ atoms: Triple<MKAtom?, MKAtom?, MKAtom?>) {
        self._vertex = atoms.0
        self._termini.0 = atoms.1
        self._termini.1 = atoms.2
        self.sortByIndex()
    }

    public static func == (lhs: MKAngle, rhs: MKAngle) -> Bool {
        return lhs._vertex == rhs._vertex && lhs._termini.0 == rhs._termini.0 && lhs._termini.1 == rhs._termini.1
    }

}

public class MKAngleData: MKGenericData {
    
    var _angles: [MKAngle] = []
    
    public init() {
        super.init("AngleData", .AngleData, .local)
    }
    
    //! Gets the number of angles stored
    //! \return integer count of the number of angles
    public func getSize() -> Int {
        return self._angles.count
    }
    
    public func clear() {
        self._angles.removeAll()
    }
    
    public func setData(_ angleData: MKAngle) {
        self._angles.append(angleData)
    }
    
    /*!
    **\brief Fills an array with the indices of the atoms in the angle (vertex first)
    **\param angles pointer to the pointer to an array of angles atom indices
    **\param size the current number of rows in the array
    **\return int The number of angles
    */
    public func fillAngleArray(_ angles: inout [[Int]], _ sizeInput: UInt) -> UInt {
        var size = sizeInput
        if self._angles.count > size {
            angles.removeAll()
            angles = Array(repeating: [0, 0, 0], count: self._angles.count)
            size = UInt(self._angles.count)
        }

        var angleIdx: Int = 0

        for angle in self._angles {
            angles[angleIdx][0] = angle._vertex!.getIdx()
            angles[angleIdx][1] = angle._termini.0!.getIdx()
            angles[angleIdx][2] = angle._termini.1!.getIdx()
            angleIdx += 1
        }

        return UInt(self._angles.count)
    }
    
    /*!
    **\brief Fills an array with the indices of the atoms in the angle (vertex first)
    **\param angles pointer to the pointer to an array of angles atom indices
    **\return True if successful
    */
    public func fillAngleArray(_ angles: inout [[UInt]]) -> Bool {
        if self._angles.isEmpty {
            return false
        }
        var ct: Int = 0

        angles.removeAll()
        angles = Array(repeating: [0, 0, 0], count: self._angles.count)

        for angle in self._angles {
            angles[ct][0] = UInt(angle._vertex!.getIdx() - 1)
            angles[ct][1] = UInt(angle._termini.0!.getIdx() - 1)
            angles[ct][2] = UInt(angle._termini.1!.getIdx() - 1)
            ct += 1
        }

        return true 
    }
    
}


public class MKSerialNums: MKGenericData {
    
    var _serialMap: OrderedDictionary<String, MKAtom> = [:]
    
    public init() {
        super.init("mkSerialNums", .SerialNums, .local)
    }
    
    public func setData(_ sm: OrderedDictionary<String, MKAtom>) {
        self._serialMap = sm
    }
    
}











public class MKPcharge: MKGenericData {
    
    var _partialCharge: Array<Double> = []
    
    public init() {
        super.init("MKPcharge", .ChargeData, .local)
    }
    
    func addPartialCharge(_ q: Array<Double>) {
        self._partialCharge.append(contentsOf: q)
    }
    
    func getPartialCharge() -> Array<Double> {
        return self._partialCharge
    }
    
}
