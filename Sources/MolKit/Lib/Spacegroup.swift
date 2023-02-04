

import Foundation
import Surge



public class MKSpaceGroup<Scalar: ExpressibleByFloatLiteral & FloatingPoint & BinaryFloatingPoint>: Hashable {
    
    var HEXAGONAL_ORIGIN: Int = 10
    var m_HM: String = ""
    var m_Hall: String = ""
    var m_id: UInt16 = 0
    var m_OriginAlternative: UInt16 = 0
    var m_transforms: [MKTransform3D<Scalar>]?
    
    init() {
        
    }

    func isValid() -> Bool {
        return false
    }

    func addTransform(_ transform: String) { }
    
    func transform(_ v: Vector<Scalar>) -> [Vector<Scalar>]? {
        return nil
    }
    
    func makeTransformIterator() -> MKIterator<MKTransform3D<Scalar>>? {
        guard let transforms = self.m_transforms else { return nil }
        return MKIterator(transforms)
    }
    
    static func getSpaceGroup(_ id: Int) -> MKSpaceGroup? {
        return nil
    }
    
    static func getSpaceGroup(_ name: String) -> MKSpaceGroup? {
        return nil
    }
    
    static func findGroup(_ spg: MKSpaceGroup) -> MKSpaceGroup? {
        return nil
    }
    
    func registerSpaceGroup() { }

    public func hash(into hasher: inout Hasher) {
        
    }
    
    public static func == (lhs: MKSpaceGroup<Scalar>, rhs: MKSpaceGroup<Scalar>) -> Bool {
        
        guard let lhs_transforms = lhs.m_transforms, let rhs_transforms = rhs.m_transforms else { return false }
        
        if lhs_transforms.count != rhs_transforms.count { return false }
//        Collect their transformations into Set<String> and compare their sizes and then their contents
        
        let lhs_set: Set<String> = Set(lhs_transforms.map { transform in
            transform.describeAsString()
        })
        
        let rhs_set: Set<String> = Set(rhs_transforms.map { transform in
            transform.describeAsString()
        })
        
        if lhs_set.count != rhs_set.count { return false }
        
//        Iteratively compare their transforms
        for string in lhs_set {
            if let found_index = rhs_set.firstIndex(of: string) {
                if rhs_set[found_index] == rhs_set[rhs_set.endIndex] {
                    return false
                }
            }
        }
        return true
    }
    

}

extension MKSpaceGroup {
    
    private enum MKSpaceGroupParseStep: Int {
        case SPACE_GROUP_ID = 0
        case SPACE_GROUP_HALL = 1
        case SPACE_GROUP_HM = 2
        case SPACE_GROUP_TRANSFORM = 3
    }
    
    private class MKSpaceGroups: MKGlobalDataBase {
        
        var _filename: String = "space-groups.txt"
        var _subdir: String = "Data"
        var _sgbn: Dictionary<String, MKSpaceGroup> = Dictionary<String, MKSpaceGroup>()
        var _sgbi: Array<[MKSpaceGroups]> = Array<[MKSpaceGroups]>.init(repeating: [], count: 230)
        var _sgs: Set<MKSpaceGroup> = Set<MKSpaceGroup>()
        
        init() {
//             Read the file and load data
        }
        
        func getSize() -> Int {
            return self._sgs.count
        }
        
        func parseLine(_ line: String) {
            
        }
        
    }
}
