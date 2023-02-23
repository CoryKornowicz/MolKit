


import Foundation 

class MKTetrahedralStereo: MKTetraNonPlanarStereo {
    
    var m_cfg: Config

    struct Config: ConfigNonPlanar, Equatable {
        
        var center: Ref
        var from_or_towrds: from_or_towrds
        var winding: MKStereo.Winding
        var view: MKStereo.View
        var specified: Bool
        var refs: Refs
        
        init() {
            center = .NoRef
            from_or_towrds = .from(.NoRef)
            winding = .Clockwise
            view = .ViewFrom
            specified = true
            refs = []
        }
        
        init(center: Ref, from_or_towrds: from_or_towrds, winding: MKStereo.Winding, view: MKStereo.View, specified: Bool, refs: Refs) {
            self.center = center
            self.from_or_towrds = from_or_towrds
            self.winding = winding
            self.view = view
            self.specified = specified
            self.refs = refs
        }
        
        static func == (lhs: MKTetrahedralStereo.Config, rhs: MKTetrahedralStereo.Config) -> Bool {
            if lhs.center == rhs.center && lhs.from_or_towrds == rhs.from_or_towrds 
            && lhs.winding == rhs.winding && lhs.view == rhs.view && lhs.specified == rhs.specified 
            && lhs.refs == rhs.refs {
                return true
            }
            return false
        }
        
    }

    override func getType() -> MKStereo.TType {
        return MKStereo.TType.Tetrahedral
    }

    func isValid() -> Bool {
        return false 
    }
    
    override init(_ mol: MKMol) {
        self.m_cfg = MKTetrahedralStereo.Config()
        super.init(mol)
    }

    func getConfig() -> MKTetrahedralStereo.Config {
        return self.m_cfg
    }

}
