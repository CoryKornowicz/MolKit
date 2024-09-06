

import Foundation 

public enum from_or_towrds: Equatable {
    case from(_ value: Ref)
    case towards(_ value: Ref)
    
    var refValue: Ref {
        get {
            switch self {
            case .from(let ref):
                return ref
            case .towards(let ref):
                return ref
            }
        }
        
        set {
            switch self {
            case .from(_):
                self = .from(newValue)
            case .towards(_):
                self = .towards(newValue)
            }
        }
    }
        
    var from: from_or_towrds {
        get {
            return .from(self.refValue)
        }
        set {
            self = newValue
        }
    }
    
    var towards: from_or_towrds {
        get {
            return .towards(refValue)
        }
        set {
            self = newValue
        }
    }
    
    public static func == (_ lhs: from_or_towrds, _ rhs: from_or_towrds) -> Bool {
        switch lhs {
        case .from(let lhsRefValue):
            switch rhs {
            case .from(let rhsRefValue):
                return lhsRefValue == rhsRefValue
            case .towards(_):
                return false
            }
        case .towards(let lhsRefValue):
            switch rhs {
            case .from(_):
                return false
            case .towards(let rhsRefValue):
                return lhsRefValue == rhsRefValue
            }
        }
    }
    
    static func == (_ lhs: from_or_towrds, _ rhs: RefValue) -> Bool {
        return lhs.refValue == rhs
    }
    
    static func != (_ lhs: from_or_towrds, _ rhs: RefValue) -> Bool {
        return lhs.refValue != rhs
    }
    
}

protocol ConfigPlanar {
    init()
    var shape: MKStereo.Shape { get set }
    var refs: Refs { get set }
}

protocol ConfigNonPlanar {
    init()
    var center: Ref { get set }
    var from_or_towrds: from_or_towrds { get set }
    var winding: MKStereo.Winding { get set }
    var view: MKStereo.View { get set }
    var refs: Refs { get set }
}

enum ConfigType {
    case Planar(config: ConfigPlanar)
    case NonPlanar(config: ConfigNonPlanar)
}

struct MKTetrahedralConfig {}
struct MKCisTransConfig {}
struct MKSquarePlanarConfig {}

