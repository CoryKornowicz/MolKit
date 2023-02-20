

import Foundation 

enum from_or_towrds {
    case from(_ value: Ref)
    case towards(_ value: Ref)
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

