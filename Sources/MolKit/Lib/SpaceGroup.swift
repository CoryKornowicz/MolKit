

import Foundation
import OrderedCollections
import simd
import Surge

public class MKSpaceGroup<Scalar>: Hashable where Scalar: FloatingPoint, Scalar: ExpressibleByFloatLiteral {
    
    var HEXAGONAL_ORIGIN: UInt = 10
    var m_Hall: String = ""
    var m_id: UInt16 = 0
    var m_OriginAlternative: UInt16 = 0
    var m_transforms: [MKTransform3D<Scalar>] = []
    private var m_HM: String = ""
    
    init() { }

    func getHMName() -> String {
        return self.m_HM
    }
    
    func setHMName(_ name: String) {
        let index = name.firstIndex(of: ":")
        
        if let idx = index {
            let origin = name[idx...]
            if origin == "H" {
                self.m_OriginAlternative = UInt16(HEXAGONAL_ORIGIN)
            } else {
                do {
                    self.m_OriginAlternative = try UInt16(String(origin), format: .number)
                } catch {
                    self.m_OriginAlternative = 0
                }
            }
        }
        
        self.m_HM = name
    }

    func isValid() -> Bool {
        if self.m_transforms.count == 0 { return false }
        var T: OrderedDictionary<String, MKTransform3D<Scalar>> = [:]
        
        for transform in self.m_transforms {
            if T.index(forKey: transform.describeAsString()) != T.count - 1 {
                print("Duplicated transform: \(transform.describeAsString())")
                return false
            }
            T[transform.describeAsString()] = transform
        }
        
        // calculate all products and check if they are in the group
        var hasInverse: Bool = false
        for j in T {
            hasInverse = false
            for k in T {
                let s = ((j.value as! MKTransform3D<GlobalPointType>) * (k.value as! MKTransform3D<GlobalPointType>)).describeAsString()
                if (T.index(forKey: s) == nil) {
                    print("Invalid transform: \(j.key) * \(k.key) = \(s)")
                    return false
                }
                if !(hasInverse && s == "x,y,z") {
                    hasInverse = true
                }
            }
            if !hasInverse {
                print("Transform with no inverse \(j.key)")
                return false
            }
        }
        return true
    }
    
    func addTransform(_ transform: String) {
        var mat: Matrix<Scalar> = Matrix(rows: 3, columns: 3, repeatedValue: 0.0)
        var vec: Vector<Scalar> = Vector(dimensions: 3, repeatedValue: 0.0)
        
        if transform.firstIndex(of: ",") != nil {
            let s1 = transform.removeWhiteSpaceAndUnderscore()
            var j: Int
            var neg: Bool = false
            
            let _components = s1.split(separator: ",")
            for i in 0..<3 {
                neg = false
                j = 0
                let currentComponent = String(_components[i])
                while (j < currentComponent.length) {
                    switch (currentComponent[j]) {
                    case "0", ".":
                        guard let val = currentComponent.substring(toIndex: j).toScalar() as? Scalar else { break }
                        switch i {
                        case 0:
                            vec[0] = val
                        case 1:
                            vec[1] = val
                        case 2:
                            vec[2] = val
                        default:
                            break
                        }
                        j = (currentComponent.substring(toIndex: j).length + 1) - currentComponent.length - 1
                        if (neg) { vec[i] = -vec[i] }
                    case "1", "2", "3", "4", "5", "6", "7", "8", "9":
                        if (j+2 < currentComponent.length && currentComponent[j+1] == "/") {
                            guard let comp = currentComponent[j].toScalar() as? Scalar else { break }
                            guard let denom = currentComponent[j+2].toScalar() as? Scalar else { break }
                            switch (i) {
                            case 0:
                                vec[0] = (comp - 0.0) / (denom - 0.0)
                            case 1:
                                vec[1] = (comp - 0.0) / (denom - 0.0)
                            case 2:
                                vec[2] = (comp - 0.0) / (denom - 0.0)
                            default:
                                break
                            }
                            if (neg) { vec[i] = -vec[i] }
                            j+=2
                          }
                        
                    case "-":
                        neg = true
                    case "+":
                        neg = false
                    case "X", "x":
                        mat[i, 0] = (neg ? -1.0 : 1.0)
                    case "Y", "y":
                        mat[i, 1] = (neg ? -1.0 : 1.0)
                    case "Z", "z":
                        mat[i, 2] = (neg ? -1.0 : 1.0)
                    default:
                        break
                    }
                    j+=1
                }
            }
        } else if transform.firstIndex(of: " ") != nil {
            /* supposing the string is a list of at least 12 float values. If there are
                16, the last four are 0., 0., 0. and 1. and are not needed */
            mat[0,0] = transform[0].toScalar() as! Scalar
            mat[0,1] = transform[1].toScalar() as! Scalar
            mat[0,2] = transform[2].toScalar() as! Scalar
            vec[0] = transform[3].toScalar() as! Scalar
            mat[1,0] = transform[4].toScalar() as! Scalar
            mat[1,1] = transform[5].toScalar() as! Scalar
            mat[1,2] = transform[6].toScalar() as! Scalar
            vec[1] = transform[7].toScalar() as! Scalar
            mat[2,0] = transform[8].toScalar() as! Scalar
            mat[2,1] = transform[9].toScalar() as! Scalar
            mat[2,2] = transform[10].toScalar() as! Scalar
            vec[2] = transform[11].toScalar() as! Scalar
        }
        
        for i in 0..<3 {
            switch vec[i] < 0 {
            case true:
                vec[i] += 1.0
            case false:
                vec[i] -= 1.0
            }
        }
        
        let candidate: MKTransform3D<Scalar> = MKTransform3D<Scalar>(m: mat, v: vec)
        let candidateSignature = candidate.describeAsString()
        var found: Bool = false
        for transform in self.m_transforms {
            if candidateSignature == transform.describeAsString() {
                found = true
            }
        }
        if !found {
            self.m_transforms.append(candidate)
        }
    }
    
    func transform(_ v: Vector<GlobalPointType>) -> [Vector<GlobalPointType>]? {
        let prec: GlobalPointType = 2e-5
        var res: [Vector<GlobalPointType>] = []
        for m_transform in m_transforms {
            guard let transform = m_transform as? MKTransform3D<GlobalPointType> else { continue }
            var vec: Vector<GlobalPointType> = Vector<GlobalPointType>.init(dimensions: 3, repeatedValue: 0.0)
            vec = transform * v
            if vec[0] < 0 { vec[0] += 1.0 }
            if vec[0] >= 1.0 { vec[0] -= 1.0 }
            if vec[1] < 0 { vec[1] += 1.0 }
            if vec[1] >= 1.0 { vec[1] -= 1.0 }
            if vec[2] < 0 { vec[2] += 1.0 }
            if vec[2] >= 1.0 { vec[2] -= 1.0 }
            // Duplicate check
            var duplicate: Bool = false
            for re in res {
                if vec.distance(re) < prec {
                    duplicate = true
                    break
                }
            }
            if !duplicate { res.append(vec) } 
        }

        return res
    }
    
    func makeTransformIterator() -> MKIterator<MKTransform3D<Scalar>> {
        return MKIterator(self.m_transforms)
    }

    public func hash(into hasher: inout Hasher) {
        hasher.combine(self.m_id)
        hasher.combine(self.m_HM)
        hasher.combine(self.m_Hall)
        hasher.combine(self.m_OriginAlternative)
    }
    
    public static func == (lhs: MKSpaceGroup<Scalar>, rhs: MKSpaceGroup<Scalar>) -> Bool {
        
        let lhs_transforms = lhs.m_transforms
        let rhs_transforms = rhs.m_transforms
        
        if lhs_transforms.count != rhs_transforms.count { return false }
//        Collect their transformations into Set<String> and compare their sizes and then their contents
        
        let lhs_set: OrderedSet<String> = OrderedSet(lhs_transforms.map { transform in
            transform.describeAsString()
        })
        
        let rhs_set: OrderedSet<String> = OrderedSet(rhs_transforms.map { transform in
            transform.describeAsString()
        })
        
        if lhs_set.count != rhs_set.count { return false }
        
//        Iteratively compare their transforms
        for string in lhs_set {
            if let found_index = rhs_set.firstIndex(of: string) {
                if rhs_set[found_index] == rhs_set.last {
                    return false
                }
            }
        }
        return true
    }

}


extension Vector where Scalar == GlobalPointType {
    
    func distance(_ to: Vector<Scalar>) -> Scalar {
        let from_vec: SIMD3<GlobalPointType> = SIMD3<GlobalPointType>.init(x: self[0], y: self[1], z: self[2])
        let to_vec: SIMD3<GlobalPointType> = SIMD3<GlobalPointType>.init(x: to[0], y: to[1], z: to[2])
        return simd_distance(from_vec, to_vec)
    }
    
}


public enum MKSpaceGroupParseStep: Int {
    case SPACE_GROUP_ID = 0
    case SPACE_GROUP_HALL = 1
    case SPACE_GROUP_HM = 2
    case SPACE_GROUP_TRANSFORM = 3
}

public class MKSpaceGroups<Scalar>: MKGlobalDataBase where Scalar: FloatingPoint, Scalar: ExpressibleByFloatLiteral {
    
    var _sgbn: OrderedDictionary<String, MKSpaceGroup<Scalar>> = OrderedDictionary<String, MKSpaceGroup<Scalar>>()
    var _sgbi: Array<[MKSpaceGroup<Scalar>]> = Array<[MKSpaceGroup<Scalar>]>.init(repeating: [], count: 230)
    var _sgs: OrderedSet<MKSpaceGroup<Scalar>> = OrderedSet<MKSpaceGroup<Scalar>>()
    
    init() {
//             Read the file and load data
        super.init(fileName: "space-groups", subDir: "Data")
        self.readFile()
        
    }
    
    override func readFile() {
        
//        Try to load contents of file
        guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }
        
        var readState: MKSpaceGroupParseStep = .SPACE_GROUP_ID
        var spacegroup: MKSpaceGroup<Scalar>? = nil
        var HMs: String = ""
        
        filePath.foreachRow { rowContents, rowNum in
            switch readState {
            case .SPACE_GROUP_ID:
//                print("reading state : \(rowContents)")
                spacegroup = MKSpaceGroup<Scalar>()
                do {
                    spacegroup!.m_id = try UInt16(rowContents, format: .number)
                } catch {
                    print("Incorrect ID format")
                    spacegroup!.m_id = 0
                }
                readState = .SPACE_GROUP_HALL
            case .SPACE_GROUP_HALL:
//                print("reading hall : \(rowContents)")
                guard let spacegroup = spacegroup else { break }
                spacegroup.m_Hall = rowContents
                readState = .SPACE_GROUP_HM
            case .SPACE_GROUP_HM:
//                print("reading hm : \(rowContents)")
                guard let spacegroup = spacegroup else { break }
                let index = rowContents.firstIndex(of: ",")
                if let idx = index {
                    let alt = rowContents[..<idx]
                    
                    if alt.count > 0 && self._sgbn[String(alt)] == nil {
                        self._sgbn[String(alt)] = spacegroup
                    }
                    
                    let altStripped = String(alt).removeWhiteSpaceAndUnderscore()
                    
                    if (altStripped.count > 0 && self._sgbn[altStripped] == nil) {
                        self._sgbn[altStripped] = spacegroup
                    }
                    
                    spacegroup.setHMName(String(rowContents[idx...]))
                } else {
                    spacegroup.setHMName(rowContents)
                }
                readState = .SPACE_GROUP_TRANSFORM
            case .SPACE_GROUP_TRANSFORM:
//                print("reading transform : \(rowContents)")
                if rowContents.isEmpty {
                    if let spacegroup = spacegroup {
                        readState = .SPACE_GROUP_ID
                        if HMs.length > 0 {
                            self.registerSpaceGroup(spacegroup, HMs)
                        } else {
                            self.registerSpaceGroup(spacegroup)
                        }
                        HMs = ""
                    }
                    spacegroup = nil
                } else {
                    guard let spacegroup = spacegroup else { break }
                    spacegroup.addTransform(rowContents)
                }
            }
        }
        
    }
    
    override func getSize() -> Int {
        return self._sgs.count
    }
    
    
    func getSpaceGroup(_ id: Int) -> MKSpaceGroup<Scalar>? {
        return (id > 0 && id <= 230) ? self._sgbi[id - 1].first : nil
    }
    
    func getSpaceGroup(_ name: String) -> MKSpaceGroup<Scalar>? {
        
        // This needs to be more forgiving
        // First, try it without removing the white space
        var match: MKSpaceGroup<Scalar>? = self._sgbn[name] != nil ? self._sgbn[name] : nil
        if match != nil { return match }
        
        let nameTrimmed = name.removeWhiteSpaceAndUnderscore()
        match = self._sgbn[nameTrimmed] != nil ? self._sgbn[nameTrimmed] : nil
        
        if (match == nil) {
            // Try another search, e.g. Fm-3m instead of Fm3m
            var search = String(name)
            let hasMirror = name.firstIndex(of: "m") != nil || (name.firstIndex(of: "d") != nil) || (name.firstIndex(of: "n") != nil) || (name.firstIndex(of: "c") != nil)
            if (name.firstIndex(of: "4") != nil && hasMirror && name.firstIndex(of: "-") != nil) {
                search.insert("-", at: name.firstIndex(of: "4")!)
            } else if (name.firstIndex(of: "3") != nil && hasMirror && name.firstIndex(of: "-") != nil) {
                search.insert("-", at: name.firstIndex(of: "3")!)
            } else if (name.firstIndex(of: "6") != nil && hasMirror && name.firstIndex(of: "-") != nil) {
                search.insert("-", at: name.firstIndex(of: "6")!)
            }

            match = self._sgbn[search] != nil ? self._sgbn[search] : nil
        }
        
        return match
    }
    
    func findGroup(_ spg: MKSpaceGroup<Scalar>) -> MKSpaceGroup<Scalar>? {
        
        var foundGroup: MKSpaceGroup<Scalar>? = nil
        
        if spg.m_Hall.length > 0 {
            foundGroup = self._sgbn[spg.m_Hall]
            if foundGroup == nil {
                print("Unknown space group (Hall Symbol: \(spg.m_Hall)) error!")
            }
            if spg.m_transforms.count > 0 && foundGroup != spg {
                let id = spg.m_id
                if id != 3 && id != 68 {
                    print("Space group error (Hall symbol and list of transforms do not match)")
                    return foundGroup
                }
            } else {
                /* even if there is an error (this should not occur) return the found group, since
                Hall names are secure */
                return foundGroup
            }
        }
        // Identify from the HM symbol, after removing all whitespaces or underscore (which are valid separators in
        // old CIF files)
        let hmTrimmed = spg.getHMName().removeWhiteSpaceAndUnderscore()
        if let found = foundGroup {
            if hmTrimmed.length > 0 && (found == self._sgbn[found.getHMName()]) {
                if found == spg {
                    return self._sgbn[found.getHMName()]
                }
                if spg.m_transforms.count > 0 {
                    // If transforms (symmetry operations) are listed, make sure they match the tabulated ones
                    for group in self._sgbi[Int(found.m_id) - 1] {
                        if group == spg {
                            return group
                        }
                    }
                    print("Unknown space group error (H-M symbol:\(spg.getHMName())), cannot match the list of transforms")
                    return nil
                } else if spg.m_transforms.count == 0 {
                    // No transforms (symmetry operations) are listed, warn if HM symbol can match several spacegroups
                    let n: Int = self._sgbi[Int(spg.m_id) - 1].map { group in
                        if group.getHMName().removeWhiteSpaceAndUnderscore() == hmTrimmed {
                            return 1
                        } else { return 0 }
                    }.reduce(0, +)
                    
                    if n > 0 {
                        print("Ambiguous space group: HM symbol corresponds to several space groups.")
                    }
                    return found
                }
            }
            
        } else if spg.m_id > 0 && spg.m_id <= 230 {
            
            if spg.m_transforms.count > 0 {
                // If transforms (symmetry operations) are listed, make sure they match the tabulated ones
                for group in self._sgbi[Int(spg.m_id) - 1] {
                    if group == spg {
                        return group
                    }
                }
                print("Unknown space group error (H-M symbol:\(spg.getHMName())), cannot match the list of transforms")
                return nil
            } else if spg.m_transforms.count == 0 {
                if self._sgbi[Int(spg.m_id) - 1].count > 1 {
                    print("Ambiguous space group: sg number corresponds to several space groups")
                }
                return self._sgbi[Int(spg.m_id) - 1].first
            }
        }
        
        if (!spg.isValid()) {
            print("Unknown space group (HM:\(spg.getHMName()),Hall:\(spg.m_Hall)) with incomplete or wrong definition.")
            return nil
        }
        
        for group in self._sgs {
            if group == spg {
                return group
            }
        }
        print("Unknown space group error, please file a bug report")
        return nil
    }
    
    func registerSpaceGroup(_ spg: MKSpaceGroup<Scalar>, _ hmName: String) {
        self.registerSpaceGroup(spg)
        if hmName.length > 0 && self._sgbn[hmName] == nil {
            self._sgbn[hmName] = spg
        }
    }
    
    func registerSpaceGroup(_ spg: MKSpaceGroup<Scalar>) {
//        print(spg)
        self._sgs.append(spg)
        if spg.m_id > 0 && spg.m_id <= 230 {
            self._sgbi[Int(spg.m_id) - 1].append(spg)
        }
        if spg.getHMName().length > 0 {
            if spg.m_OriginAlternative != 0 {
                let name = "\(spg.getHMName()):\(spg.m_OriginAlternative)"
                if self._sgbn[name] == nil {
                    self._sgbn[name] = spg
                }
                let nameTrimmed = name.removeWhiteSpaceAndUnderscore()
                if self._sgbn[nameTrimmed] == nil && nameTrimmed.length > 0 {
                    self._sgbn[nameTrimmed] = spg
                }
            }
            if ((spg.m_OriginAlternative & 1) == 0) && self._sgbn[spg.getHMName()] == nil {
                self._sgbn[spg.getHMName()] = spg
            }
        }
        let nameTrimmed = spg.getHMName().removeWhiteSpaceAndUnderscore()
        if self._sgbn[nameTrimmed] == nil && nameTrimmed.length > 0 {
            self._sgbn[nameTrimmed] = spg
        }
        if self._sgbn[spg.m_Hall] == nil && spg.m_Hall.length > 0 {
            self._sgbn[spg.m_Hall] = spg
        }
    }
    
}
