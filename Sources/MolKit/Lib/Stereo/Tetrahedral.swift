


import Foundation 

class MKTetrahedralStereo: MKTetraNonPlanarStereo {
    
    private var m_cfg: Config?

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
            
            if lhs.center != rhs.center { return false }
            if lhs.refs.count != 3 || rhs.refs.count != 3 { return false }
            // return true if either is unspecified (i.e. accidental)
            if !lhs.specified || !rhs.specified { return true }
            
            // Convert both Config's refs to same from, winding and view while
            // avoiding having an ImplicitRef in the 'from' position of either
            var lhsConfig: Config 
            var rhsConfig: Config 
            
            if lhs.from_or_towrds == .from(.ImplicitRef) { // not check the towards parameter (could it ever be there?)
                lhsConfig = MKTetraNonPlanarStereo.toConfig(lhs, .from(lhs.refs[0]), lhs.winding, lhs.view)
                rhsConfig = MKTetraNonPlanarStereo.toConfig(rhs, lhsConfig.from_or_towrds, lhs.winding, lhs.view)
            } else if rhs.from_or_towrds == .from(.ImplicitRef) {
                rhsConfig = MKTetraNonPlanarStereo.toConfig(rhs, .from(rhs.refs[0]), lhs.winding, lhs.view)
                lhsConfig = MKTetraNonPlanarStereo.toConfig(lhs, rhsConfig.from_or_towrds, lhs.winding, lhs.view)
            } else {
                lhsConfig = lhs
                rhsConfig = MKTetraNonPlanarStereo.toConfig(rhs, lhsConfig.from_or_towrds, lhs.winding, lhs.view)
            }
            
            if !MKStereo.containsSameRef(lhsConfig.refs, rhsConfig.refs) {
                if MKStereo.containsRef(lhsConfig.refs, .ImplicitRef) {
                    // if both refs already contain ImplicitRef, return false
                    if MKStereo.containsRef(rhsConfig.refs, .ImplicitRef) {
                        return false
                    }
                    // example: *this       = 23H
                    //          otherConfig = 234 --> 23H
                    // for each ref in otherConfig
                    for i in 0..<rhsConfig.refs.count {
                        var found: Bool = false 
                        for ref in lhsConfig.refs {
                            if rhsConfig.refs[i] == ref {
                                found = true
                            }
                        }
                        if !found {
                            // the ref from otherConfig is not found in this config
                            rhsConfig.refs[i] = .ImplicitRef
                            break
                        }
                    }
                } else {
                    if MKStereo.containsRef(rhsConfig.refs, .ImplicitRef) {
                        // if both refs already contain ImplicitRef, return false
                        if MKStereo.containsRef(lhsConfig.refs, .ImplicitRef) {
                            return false
                        }
                        // example: *this       = 234 --> 23H
                        //          otherConfig = 23H
                        // for each ref in this config
                        for i in 0..<lhsConfig.refs.count {
                            var found: Bool = false 
                            for ref in rhsConfig.refs {
                                if lhsConfig.refs[i] == ref {
                                    found = true
                                }
                            }
                            if !found {
                                for var refIter in rhsConfig.refs {
                                    if refIter == .ImplicitRef {
                                        refIter = lhsConfig.refs[i]
                                    }
                                }
                                break
                            }
                        }
                    }
                }
            }
            
            let Ni1 = MKStereo.numInversions(lhsConfig.refs)
            let Ni2 = MKStereo.numInversions(rhsConfig.refs)
            return ((Ni1 + Ni2) % 2 == 0)
        }
    }

    override func getType() -> MKStereo.TType {
        return .Tetrahedral
    }

    override init(_ mol: MKMol) {
        super.init(mol)
        self.m_cfg = Config()
    }

    /**
     * @return True if this object is valid. This object is valid if all (center, from
     * and ref) atom ids are set.
     */
    func isValid() -> Bool {
        guard m_cfg != nil else {
            return false
        }
        if m_cfg!.center == .NoRef || m_cfg!.refs.count != 3 {
            return false
        }
        switch m_cfg!.from_or_towrds {
        case .from(let from):
            if from == .NoRef {
                return false
            }
        case .towards(let towards):
            if towards == .NoRef {
                return false
            }
        }
        return true
    }
    
    /**
     * Set the configuration using a Config struct.
     */
    func setConfig(_ cfg: Config) {
        if cfg.center == .NoRef {
            print("MKTetrahedralStereo::setConfig(): center atom id is invalid")
            m_cfg = Config()
            return
        }
        if cfg.from_or_towrds == .from(.NoRef) || cfg.from_or_towrds == .towards(.NoRef) {
            print("MKTetrahedralStereo::setConfig(): from/towards atom id is invalid")
            m_cfg = Config()
            return
        }
        if cfg.refs.count != 3 {
            print("MKTetrahedralStereo::setConfig(): number of reference atom ids is invalid, found \(cfg.refs.count), expected 3")
            m_cfg = Config()
            return
        }
    }
    
    /**
     * Get the configuration as Config struct.
     */
    func getConfig(_ winding: MKStereo.Winding = .Clockwise, _ view: MKStereo.View = .ViewFrom) -> Config {
        if !isValid() {
            return Config()
        }
        if m_cfg!.winding == .UnknownWinding {
            return MKTetraNonPlanarStereo.toConfig(m_cfg!, m_cfg!.from_or_towrds, winding, view)
        } else {
            return MKTetraNonPlanarStereo.toConfig(m_cfg!, m_cfg!.from_or_towrds, .UnknownWinding, view)
        }
    }
    
    /**
     * Get the configuration as Config struct viewing from/towards the specified id.
     */
    func getConfig(_ fromTorwards: from_or_towrds, _ winding: MKStereo.Winding = .Clockwise, _ view: MKStereo.View = .ViewFrom) -> Config {
        if !isValid() {
            return Config()
        }
        if m_cfg!.winding == .UnknownWinding {
            return MKTetraNonPlanarStereo.toConfig(m_cfg!, fromTorwards, winding, view)
        } else {
            return MKTetraNonPlanarStereo.toConfig(m_cfg!, fromTorwards, .UnknownWinding, view)
        }
    }
    
    /**
     * Compare the stereochemistry stored in the Config struct with the
     * stereochemistry specified in the Config struct from @p other.
     *
     * @copydoc Config::operator==()
     */
    static func == (lhs: MKTetrahedralStereo, rhs: MKTetrahedralStereo) -> Bool {
        if !lhs.isValid() || !rhs.isValid() {
            return false
        }
        if lhs.getConfig() == rhs.getConfig() {
            return true
        }
        return false
    }
    
    override func clone<T>(_ mol: T) -> MKGenericData? where T : MKBase {
        let data = MKTetrahedralStereo(mol as! MKMol)
        data.setConfig(m_cfg!)
        return data
    }
}

// TODO: Should extend to ExpressibleByStringLiteral
