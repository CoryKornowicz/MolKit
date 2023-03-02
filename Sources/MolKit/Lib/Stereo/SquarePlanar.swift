



import Foundation 

/**
 * @class OBSquarePlanarStereo squareplanar.h <openbabel/stereo/squareplanar.h>
 * @brief Class for handling and storing square planar stereochemistry.
 *
 * @image html squareplanar.png
 *
 * The OBSquarePlanarStereo class is used to represent square planar stereochemistry.
 * Like all OBTetraPlanarStereo subclasses, it uses the OBStereo::Shape parameters
 * to set/get the reference ids.
 *
 * This class works with reference ids only, it never uses the molecule
 * to get more information. Like all stereo classes, errors,
 * warnings or info is reported using OBMessageHandler.
 */
class MKSquarePlanarStereo: MKTetraPlanarStereo {
    
    struct Config: ConfigPlanar, Equatable {
        
        var shape: MKStereo.Shape = .ShapeU
        var center: Ref
        var specified: Bool = true
        var refs: Refs

        init() {
            center = .NoRef
            specified = true
            refs = []
        }

        init(center: Ref, shape: MKStereo.Shape = .ShapeU, refs: Refs) {
            self.center = center
            self.shape = shape
            self.specified = true
            self.refs = refs
        }

        /**
       * Equal to operator. Comparing OBSquarePlanarStereo::Config structs
       * is done using the information stored in the struct's data members
       * (i.e. center, refs and shape).
       *
       * There are a number of cases resuling in false being returned:
       * - @p center atom ids don't match
       * - One of the Refs lists does not contain 4 elements.
       * - 2 or more OBStereo::ImplicitRef values in a single Config struct
       * - (The two @p refs don't share a single common element)
       *
       * In the simplest case where both @p refs contain exactly the same elements
       * (OBStereo::ContainsSameRefs()), coould include OBStereo::ImplicitRef), both Config
       * struct are normalized to OBStereo::ShapeU starting with the same element.
       * After this normalization, there are two possible orientations to overlay the
       * shape on the double bond. From the illustration below, it can be seen only
       * @p refs[2] has to be checked in order to conclude both Config structs
       * have the same stereochemistry.
       *
         @verbatim
         1   4    1      4    1------4
          \ /     |      |           |
           C      |      |           |
          / \     |      |           |
         2   3    2------3    2------3

                  1 2 3 4     1 2 3 4
                  |   |       |   |      <- in any case, refs[0] & refs[2] remain unchanged
                  1 2 3 4     1 4 3 2
        @endverbatim
       *
       * When comparing a Config struct with explicit hydrogen(s) to one with
       * implicit hydrogen(s), both @p refs are also normalized to OBStereo::ShapeU
       * starting with the same common element. This shared element cannot be
       * OBStereo::ImplicitRef. Depending on the position of the OBStereo::ImplicitRef
       * element(s) in the @p refs, 3 cases are possible:
       *
        @verbatim

         refs[2] != OBStereo::ImplicitId:

           (analog to the case above where they contained the same elements )

           1 2 3 4
           |   |      <- refs[0] & refs[2] remain unchanged
           1 H 3 H

         else:

           1 2 3 4
           |     |    <- refs[0] & refs[3] remain unchanged
           1 H H 4

           1 2 3 4
           | |        <- refs[0] & refs[1] remain unchanged
           1 2 H H
        @endverbatim
       *
       * In each case, the orientation of the U shape is also defined since
       * there can be only one OBStereo::ImplicitRef for each side of the
       * double bond.
       *
       * @return True if both Config structs represent the stereochemistry.
       */
        static func == (lhs: MKSquarePlanarStereo.Config, rhs: MKSquarePlanarStereo.Config) -> Bool {
            if lhs.center != rhs.center { return false }
            if lhs.refs.count != 4 || rhs.refs.count != 4 { return false }
            var u1: Config?
            var u2: Config?

            if !MKStereo.containsSameRefs(lhs.refs, rhs.refs) {
                // find a ref that occurs in both
                for refIter in lhs.refs {
                    if MKStereo.containsRef(rhs.refs, refIter) {
                        u1 = MKTetraPlanarStereo.toConfig(lhs, refIter, .ShapeU) // refs[0] = u1.refs[0]
                        u2 = MKTetraPlanarStereo.toConfig(rhs, refIter, .ShapeU) // refs[0] = u2.refs[0]
                    }
                }
                // check if they actually share an id...
                if u1!.refs.isEmpty {
                    return false
                }
            } else {
                // normalize the other Config struct
                u1 = MKTetraPlanarStereo.toConfig(lhs, lhs.refs[0], .ShapeU) // refs[0] = u1.refs[0]
                u2 = MKTetraPlanarStereo.toConfig(rhs, lhs.refs[0], .ShapeU) // refs[0] = u2.refs[0]
                // both now start with the same ref
                //
                // 2 possiblilities:
                //
                //   1 2 3 4      1 2 3 4
                //   |   |        |   |      <- in any case, refs[0] & refs[2] remain unchanged
                //   1 2 3 4      1 4 3 2
                //
                return (u1!.refs[2] == u2!.refs[2])
            }

            // possibilities:
            //
            //   1 2 3 4
            //   |   |      <- refs[0] & refs[2] remain unchanged
            //   1 H 3 H
            //
            //   1 2 3 4
            //   |     |    <- refs[0] & refs[3] remain unchanged
            //   1 H H 4
            //
            //   1 2 3 4
            //   | |        <- refs[0] & refs[1] remain unchanged
            //   1 2 H H
            if u1!.refs[2] == .ImplicitRef || u2!.refs[2] == .ImplicitRef {
                if u1!.refs[3] == .ImplicitRef || u2!.refs[3] == .ImplicitRef {
                    return (u1!.refs[1] == u2!.refs[1]) // 1 2 H H
                } else {
                    return (u1!.refs[3] == u2!.refs[3]) // 1 H H 4
                }
            } else {
                return (u1!.refs[2] == u2!.refs[2]) // 1 2 3 4  &  1 H 3 4  &  1 2 3 H
            }
        }
    }
    
    private var m_cfg: Config?

    override init(_ mol: MKMol) {
        super.init(mol)
        m_cfg = Config()
    }
    

    override func getType() -> MKStereo.TType {
        return .SquarePlanar
    }

    /**
     * @return True if this object is valid. This object is valid if all (center and
     * and 4 reference) atom ids are set.
     */
    func isValid() -> Bool {
        guard m_cfg != nil else { return false }
        if m_cfg!.center == .NoRef { return false }
        if m_cfg!.refs.count != 4 { return false }
        return true
    }

    /**
     * Set the configuration using a Config struct.
     */
    func setConfig(_ cfg: Config) {
        if cfg.center == .NoRef { 
            print("MKSquarePlanarStereo::setConfig(): center id is invalid.")
            m_cfg = Config() 
            return 
        }
        if cfg.refs.count != 4 { 
            print("MKSquarePlanarStereo::setConfig(): number of reference ids is invalid found \(cfg.refs.count) expected 4.)")
            m_cfg = Config() 
            return 
        }
        // store using U Shape 
        m_cfg = MKTetraPlanarStereo.toConfig(cfg, cfg.refs[0], .ShapeU)
    }

    func getConfig(_ shape: MKStereo.Shape = .ShapeU) -> Config {
        if !isValid() { return Config() }
        return MKTetraPlanarStereo.toConfig(m_cfg!, m_cfg!.refs[0], shape)
    }

    func getConfig(_ shape: MKStereo.Shape = .ShapeU, _ start: Ref) -> Config {
        if !isValid() { return Config() }
        return MKTetraPlanarStereo.toConfig(m_cfg!, start, shape)
    }

    /**
    * Compare the stereochemistry stored in the Config struct with the
    * stereochemistry specified in the Config struct from @p other.
    *
    * @copydoc Config::operator==()
    */
    static func == (lhs: MKSquarePlanarStereo, rhs: MKSquarePlanarStereo) -> Bool {
        
        if !lhs.isValid() || !rhs.isValid() { return false }
        
        var u: Config = MKTetraPlanarStereo.toConfig(rhs.getConfig(), lhs.m_cfg!.refs[0], .ShapeU)

        var a1 = u.refs[0]
        var b1 = u.refs[2]

        if a1 == .ImplicitRef && b1 == .ImplicitRef {
            a1 = u.refs[1]
            b1 = u.refs[3]
        }

        if b1 != .ImplicitRef {
            if lhs.getTransRef(b1) == a1 {
                return true
            }
        }

        if a1 != .ImplicitRef {
            if lhs.getTransRef(a1) == b1 {
                return true
            }
        }

        return false
    }

    //! @name Query methods to compare stereochemistry.
    //@{
    /**
     * @return True if the two reference ids are placed trans configuration.
     */
    // bool IsTrans(unsigned long id1, unsigned long id2) const;
    func isTrans(_ id1: Ref, _ id2: Ref) -> Bool {
        return (getTransRef(id1) == id2)
    }
    /**
     * @return True if the two reference ids are placed in a cis configuration.
     */
    func isCis(_ id1: Ref, _ id2: Ref) -> Bool {
        guard m_cfg != nil else { return false }
        if m_cfg!.refs.count != 4 { return false }
        let cis = getCisRefs(id1)
        if cis.count != 2 { return false }
        if cis[0] == id2 || cis[1] == id2 { return true }
        return false
    }
    /**
     * @image html gettransref.png
     * Get the reference id trans from reference @p id.
     */
    func getTransRef(_ id: Ref) -> Ref {
        guard m_cfg != nil else { return .NoRef }
        if m_cfg!.refs.count != 4 { return .NoRef }

        // find id1
        for i in 0..<4 {
            if m_cfg!.refs[i] == id {
                // use it's index to compare id2 with the opposite reference id
                let j = (i > 1) ? i - 2 : i + 2
                return m_cfg!.refs[j]
            }
        }
        return .NoRef
    }
    /**
     * Get the reference id cis from reference @p id.
     */
    func getCisRefs(_ id: Ref) -> [Ref] {
        guard m_cfg != nil else { return [.NoRef] }
        if m_cfg!.refs.count != 4 { return [.NoRef] }
        var refs: [Ref] = []
        // find id
        for i in 0..<4 {
            if m_cfg!.refs[i] == id {
                // use it's index to get the left/right reference ids
                var j = (i > 0) ? i - 1 : 3
                var k = (i < 3) ? i + 1 : 0
                refs.append(m_cfg!.refs[j]) 
                refs.append(m_cfg!.refs[k])
                return refs
            }
        }
        return [.NoRef]
    }

    // The following function sits behind GetCisRef and GetTransRef
    private func getCisOrTrans(_ id: Ref, _ getcisref: Bool) -> Ref {
        if !isValid() { return .NoRef }
        if id == .ImplicitRef { return .NoRef }

        //  find id
        for i in 0..<4 {
            if m_cfg!.refs[i] == id {
                // Use its index to find the index of the cis (or trans) atom
                var j: Int
                if (getcisref) { // GetCisRef
                    j = 3 - i // Convert 0 to 3, and 3 to 0
                } else { // GetTransRef
                    j = (i > 1) ? i - 2 : i + 2
                }
                let refId = m_cfg!.refs[j]
                return refId
            }
        }
//        id is not found
        return .NoRef
    }

    override func clone<T>(_ mol: T) -> MKGenericData? where T : MKBase {
        let data = MKSquarePlanarStereo(mol as! MKMol)
        data.setConfig(m_cfg!)
        return data
    }

}




// TODO: Should extend to ExpressibleByStringLiteral
