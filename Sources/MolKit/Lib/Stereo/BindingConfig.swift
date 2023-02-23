

import Foundation 

enum from_or_towrds: Equatable {
    case from(_ value: Ref)
    case towards(_ value: Ref)
    
    static func == (_ lhs: from_or_towrds, _ rhs: from_or_towrds) -> Bool {
        switch lhs {
        case .from(let lhsf):
            switch rhs {
            case .from(let rhsf):
                return lhsf == rhsf
            case .towards:
                return false
            }
        case .towards(let lhst):
            switch rhs {
            case .from:
                return false
            case .towards(let rhst):
                return lhst == rhst
            }
        }
    }
    
}

protocol ConfigPlanar {
    var shape: MKStereo.Shape { get }
    var refs: Refs { get }
}

protocol ConfigNonPlanar {
    var from_or_towrds: from_or_towrds { get }
    var winding: MKStereo.Winding { get }
    var view: MKStereo.View { get }
    var refs: Refs { get }
}

enum ConfigType {
    case Planar(config: ConfigPlanar)
    case NonPlanar(config: ConfigNonPlanar)
}

struct MKTetrahedralConfig {}
struct MKCisTransConfig {}
struct MKSquarePlanarConfig {}

