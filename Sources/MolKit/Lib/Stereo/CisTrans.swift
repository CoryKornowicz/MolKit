

import Foundation 


/**
 * @class OBCisTransStereo cistrans.h <openbabel/stereo/cistrans.h>
 * @brief Class for handling and storing cis/trans stereochemistry.
 *
 * @image html cistrans.png
 *
 * The OBCisTransStereo class is used to represent cis/trans stereochemistry.
 * Like all OBTetraPlanarStereo subclasses, it uses the OBStereo::Shape parameters
 * to set/get the OBStereo::Ref values in the Config struct. However, since the
 * orientation of the double bond matters, all methods in the "query methods"
 * section check the bonding by actually checking the atoms in the molecule
 * provided through OBStereoBase's constructor. If OBMol::GetAtomById(id) returns
 * 0 for a single id, it will be considered a deleted or implicit hydrogen if
 * the valences confirm this.
 *
 * The use of OBStereo::Shape is illustarted in the image below.
 * @image html SPshapes.png
 *
 * @code
 * //
 *  // 3              6                    3       6
 *  //   \            /                       |       |
 *  //     0===1            =          | 0   1 |
 *  //   /           \                        |       |
 *  // 4             5                    4-------5
 * OBCisTransStereo ct;
 * ct.SetCenters(0, 1);
 * ct.SetRefs(OBStereo::MakeRefs(3, 4, 5, 6), OBStereo::ShapeU);
 *
 * ct.IsTrans(3, 5); // true
 * ct.IsTrans(3, 4); // false
 * ct.IsTrans(3, 6); // false
 *
 * ct.IsCis(3, 6); // true
 * ct.IsCis(3, 4); // false
 * ct.IsCis(3, 5); // false
 * @endcode
 *
 */
class MKCisTransStereo: MKTetraPlanarStereo {
    
    /**
     * \struct Config cistrans.h <openbabel/stereo/cistrans.h>
     * \brief Stereochemical configuration for double-bond cis/trans stereochemistry
     *
     * The config struct represents the stereochemistry in a well defined way.
     * For cis/trans stereo bonds, the following data members define the spacial
     * arrengement of the atoms.
     *
     * - OBStereo::Ref @p begin: The begin atom for the double bond.
     * - OBStereo::Ref @p end: The end atom for the double bond.
     * - OBStereo::Refs @p refs: The 4 atoms connected to the double bond.
     * - OBStereo::Shape @p shape: The shape formed by the @p refs by connecting them
     *   in the same order as they occur in @p refs.
     *
     * @image html cistrans.png
     * @image html SPshapes.png
     *
     * Only @p begin and @p end are specific for OBCisTransStereo::Config. The other
     * data members occur in all OBTetraPlanarStereo derived classes.
     */
    class Config: ConfigPlanar, Equatable {
        var begin: Ref
        var end: Ref
        var shape: MKStereo.Shape
        var refs: Refs
        var specified: Bool
        
        init(begin: Ref, end: Ref, shape: MKStereo.Shape = .ShapeU, refs: Refs) {
            self.begin = begin
            self.end = end
            self.shape = shape
            self.refs = refs
            self.specified = true
        }
        
        required init() {
            self.begin = .NoRef
            self.end = .NoRef
            self.refs = []
            self.shape = .ShapeU
            self.specified = true
        }
        
        static func == (_ lhs: MKCisTransStereo.Config, _ rhs: MKCisTransStereo.Config) -> Bool {
            if lhs.begin != rhs.begin && lhs.begin != rhs.end { return false }
            if lhs.end != rhs.begin && lhs.end != rhs.end { return false }
            if lhs.refs.count != 4 || rhs.refs.count != 4 { return false }
            
            var u1: Config?
            var u2: Config?

            if !(MKStereo.containsSameRefs(lhs.refs, rhs.refs)) {
                for i in lhs.refs {
                    if MKStereo.containsRef(rhs.refs, i) {
                        u1 = MKTetraPlanarStereo.toConfig(lhs, i, .ShapeU) // refs[0] = u1.refs[0]
                        u2 = MKTetraPlanarStereo.toConfig(rhs, i, .ShapeU) // refs[0] = u2.refs[0]
                    }
                }
                // check if they actualy share an id...
                if u1!.refs.isEmpty { return false }
            } else {
                // normalize the other Config struct
                u1 = MKTetraPlanarStereo.toConfig(lhs, lhs.refs[0], .ShapeU) // refs[0] = u1.refs[0]
                u2 = MKTetraPlanarStereo.toConfig(rhs, rhs.refs[0], .ShapeU) // refs[0] = u2.refs[0]
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
    
    private var m_cfg: Config? //!< internal configuration
    
    override func getType() -> MKStereo.TType {
        return .CisTrans
    }
    
    override init(_ mol: MKMol) {
        super.init(mol)
    }
    
    /**
     * @return True if this object is valid. This object is valid if all these
     * conditions are met:
     * - @p begin != OBStereo::NoRef
     * - @p end != OBStereo::NoRef
     * - @p refs contains 4 elements
     */
    func isValid() -> Bool {
        guard m_cfg != nil else { return false }
        if m_cfg?.begin == .NoRef || m_cfg?.end == .NoRef { return false }
        if m_cfg?.refs.count != 4 { return false }
        return true
    }
    
    /**
     * Set the configuration using a Config struct.
     */
    func setConfig(_ config: Config) {
        if config.begin == .NoRef || config.end == .NoRef || config.refs.count != 4 {
            print("MKCisTransStereo::setConfig : invalid config.")
            m_cfg = nil
            return
        }
        m_cfg = MKTetraPlanarStereo.toConfig(config, config.refs[0], .ShapeU)
    }
    /**
     * Get the configuration as Config struct.
     */
    func getConfig(_ shape: MKStereo.Shape = .ShapeU) -> Config {
        if !isValid() {
            return Config()
        }
        return MKTetraPlanarStereo.toConfig(m_cfg!, m_cfg!.refs[0], shape)
    }
    /**
     * Get the configuration as Config struct and ensure refs[0] is
     * equal to @p start.
     */
    func getConfig(_ start: Ref, _ shape: MKStereo.Shape = .ShapeU) -> Config {
        if !isValid() {
            return Config()
        }
        return MKTetraPlanarStereo.toConfig(m_cfg!, start, shape)
    }
    
    /**
     * Compare the stereochemistry stored in the Config struct with the
     * stereochemistry specified in the Config struct from @p other.
     *
     * @copydoc Config::operator==()
     */
    static func == (_ lhs: MKCisTransStereo, _ rhs: MKCisTransStereo) -> Bool {
        if !lhs.isValid() || !rhs.isValid() { return false }
        let u: Config = MKTetraPlanarStereo.toConfig(rhs.getConfig(), lhs.m_cfg!.refs[0], .ShapeU)
        var a1 = u.refs[0]
        var b1 = u.refs[2]

        if a1 == .ImplicitRef && b1 == .ImplicitRef {
            a1 = u.refs[1]
            b1 = u.refs[3]
        }

        if b1 != .ImplicitRef {
            if a1 == lhs.getTransRef(b1) {
                return true
            }
        }
        if a1 != .ImplicitRef {
            if b1 == lhs.getTransRef(a1) {
                return true
            }
        }
        return false
    }
    
    override func clone<T>(_ mol: T) -> MKGenericData? where T : MKBase {
        guard m_cfg != nil else { return nil }
        let data = MKCisTransStereo(mol as! MKMol)
        data.setConfig(m_cfg!)
        return data
    }
    
    //! @name Query methods to compare stereochemistry.
    //@{
    /**
     * Check if the two atoms for @p id1 & @p id2 are bonded to the same
     * atom. If the atoms for one of the atoms doesn't exist (anymore),
     * the valence of the begin and end atom is checked. If the exising
     * atom is bonded to the begin atom and end->GetExplicitDegree() == 2, the
     * ids are considered to be on different atoms. The reasoning behind
     * this is that hydrogens may be deleted. However, you can also use
     * OBStereo::ImplicitRef explicitly in code like:
     *
     * @code
     * //
     * //  F       F      F      F     0      3
     * //   \     /       |      |     |      |
     * //    C===C        | C  C |     | 1  2 |
     * //   /     \       |      |     |      |
     * // (H)     (H)    (H)----(H)    H------H
     * //
     * reading smiles F/C=C\F cis-difluorethene
     * OBCisTransStereo ct(mol);
     * ct.SetCenters(1, 2);
     * ct.SetRefs(OBStereo::MakeRefs(0, OBStereo::ImplicitRef, OBStereo::ImplicitRef, 3));
     * ...
     * @endcode
     *
     * @return True if @p id1 and @p id2 are bonded to the same atom
     * taking implicit hydrogens into account.
     */
    func isOnSameAtom(_ id1: Ref, _ id2: Ref) -> Bool {
        let mol = getMolecule()
        guard let begin = mol.getAtomById(m_cfg!.begin) else {
            print("MKCisTransStereo::isOnSameAtom : begin atom does not exist.")
            return false
        }
        guard let end = mol.getAtomById(m_cfg!.end) else {
            print("MKCisTransStereo::isOnSameAtom : end atom does not exist.")
            return false
        }
        let a = mol.getAtomById(id1)
        let b = mol.getAtomById(id2)
        if a != nil && b != nil {
            // both on begin atom?
            if a!.isConnected(begin) && b!.isConnected(begin) {
                return true
            }
            // both on end atom?
            if a!.isConnected(end) && b!.isConnected(end) {
                return true
            }
            return false
        } else {
            if (a != nil) {
                // b atom not found, could be a deleted hydrogen...
                if a!.isConnected(begin) {
                    // a is connected to begin. if this is the atom missing a hydrogen, return false
                    if begin.getExplicitDegree() == 2 {
                        return true
                    }
                    // check if the end atom really is missing an atom
                    if end.getExplicitDegree() != 2 {
                        print("MKCisTransStereo.isOnSameAtom : id2 is not valid and is not a missing hydrogen.")
                        return false
                    }
                    // inform user we are treating id2 as deleted hydrogen
                    print("MKCisTransStereo.isOnSameAtom: Atom with id2 doesn't exist anymore, must be a (deleted) hydrogen.")
                } else if (a!.isConnected(end)) {
                    // a is connected to end. again, if this is the atom missing a hydrogen, return false
                    if (end.getExplicitDegree() == 2) {
                        return true
                    }
                    // check if the begin atom really is missing an atom
                    if (begin.getExplicitDegree() != 2) {
                        print("MKCisTransStereo.isOnSameAtom : id2 is not valid and is not a missing hydrogen.")
                        return false
                    }
                    // inform user we are treating id2 as deleted hydrogen
                    print("MKCisTransStereo.isOnSameAtom: Atom with id2 doesn't exist anymore, must be a (deleted) hydrogen.")
                } else {
                    print("MKCisTransStereo.isOnSameAtom: Atom with id1 isn't connected to the begin or end atom.")
                    return false //error in original code** was true
                }
            } else if (b != nil) {
                // a atom not found, could be a deleted hydrogen...
                if (b!.isConnected(begin)) {
                    // b is connected to begin. if this is the atom missing a hydrogen, return false
                    if (begin.getExplicitDegree() == 2) {
                        return true
                    }
                    // check if the end atom really is missing an atom
                    if (end.getExplicitDegree() != 2) {
                        print("MKCisTransStereo.isOnSameAtom : id1 is not valid and is not a missing hydrogen.")
                        return false  //error in original code** was true
                    }
                    // inform user we are treating id1 as deleted hydrogen
                    print("MKCisTransStereo.isOnSameAtom: Atom with id1 doesn't exist anymore, must be a (deleted) hydrogen.")
                } else if (b!.isConnected(end)) {
                    // a is connected to end. again, if this is the atom missing a hydrogen, return false
                    if (end.getExplicitDegree() == 2) {
                        return true
                    }
                    // check if the begin atom really is missing an atom
                    if (begin.getExplicitDegree() != 2) {
                        print("MKCisTransStereo.isOnSameAtom : id1 is not valid and is not a missing hydrogen.")
                        return false  //error in original code** was true
                    }
                    // inform user we are treating id2 as deleted hydrogen
                    print("MKCisTransStereo.isOnSameAtom: Atom with id1 doesn't exist anymore, must be a (deleted) hydrogen.")
                }else {
                    print("MKCisTransStereo.isOnSameAtom: Atom with id1 isn't connected to the begin or end atom.")
                    return false //error in original code** was true
                }
            } else {
                    var c: MKAtom? = nil
                    var d: MKAtom? = nil
                    // no a & b, check the remaining ids which will reveal same info
                //     for (int i = 0; i < 4; ++i) {
                    for i in 0..<4 {
                        if ((m_cfg!.refs[i] == id1) || (m_cfg!.refs[i] == id2)) {
                            continue
                        }
                        if (c == nil) {
                            c = mol.getAtomById(m_cfg!.refs[i])
                        } else {
                            d = mol.getAtomById(m_cfg!.refs[i])
                        }
                    }
                    
                    if (c == nil || d == nil) {
                        print("MKCisTransStereo.isOnSameAtom: invalid stereochemistry!")
                        return false //error in original code** was true
                    }
                    if begin.getExplicitDegree() != 2 || end.getExplicitDegree() != 2 {
                        print("MKCisTransStereo.isOnSameAtom: invalid stereochemistry!")
                        return false //error in original code** was true
                    }
                    print("MKCisTransStereo.isOnSameAtom: Atoms with id1 & id2 don't exist, must be a (deleted) hydrogens.")
                    
                    return isOnSameAtom(c!.getId(), d!.getId())
            }
        }
        
        return false
    }
    /**
     * @return True if the two reference ids are placed trans configuration.
     */
    func isTrans(_ id1: Ref, _ id2: Ref) -> Bool {
        return getTransRef(id1) == id2
    }
    /**
     * @return True if the two reference ids are placed in a cis configuration.
     */
    func isCis(_ id1: Ref, _ id2: Ref) -> Bool {
        return getCisRef(id1) == id2
    }
    /**
     * @image html gettransref.png
     * Get the reference id trans from reference @p id.
     */
    func getTransRef(_ id: Ref) -> Ref {
        getCisOrTransRef(id, false)
    }
    /**
     * @image html getcisref.png
     * Get the reference id cis from reference @p id.
     */
    func getCisRef(_ id: Ref) -> Ref {
        getCisOrTransRef(id, true)
    }
    //@}
    
    // The following function sits behind GetCisRef and GetTransRef
    private func getCisOrTransRef(_ id: Ref,_ getcisref: Bool) -> Ref {
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
    
}


//  TODO: Should extend to ExpressibleByStringLiteral
