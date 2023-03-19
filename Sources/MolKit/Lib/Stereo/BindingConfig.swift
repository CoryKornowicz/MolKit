

import Foundation 

enum from_or_towrds: Equatable {
    case from(_ value: Ref)
    case towards(_ value: Ref)
    
    var value: Int? {
        switch self {
        case .from(let ref):
            switch ref {
            case .NoRef:
                return nil
            case .ImplicitRef:
                return nil
            case .Ref(let intValue):
                return intValue
            }
        case .towards(let ref):
            switch ref {
            case .NoRef:
                return nil
            case .ImplicitRef:
                return nil
            case .Ref(let intValue):
                return intValue
            }
        }
    }
    
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
            case .from(var ref):
                ref = newValue
            case .towards(var ref):
                ref = newValue
            }
        }
    }
        
    var from: Ref? {
        get {
            if case .from(let value) = self {
                return value
            } else {
                return nil
            }
        }
        set {
            self = .from(newValue!)
        }
    }
    
    var towards: Ref? {
        get {
            if case .towards(let value) = self {
                return value
            } else {
                return nil
            }
        }
        set {
            self = .towards(newValue!)
        }
    }
    
    static func == (_ lhs: from_or_towrds, _ rhs: from_or_towrds) -> Bool {
        return lhs.value == rhs.value
    }
    
    static func == (_ lhs: from_or_towrds, _ rhs: RefValue) -> Bool {
        return lhs.value == rhs.intValue
    }
    
    static func != (_ lhs: from_or_towrds, _ rhs: RefValue) -> Bool {
        return lhs.value != rhs.intValue
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

