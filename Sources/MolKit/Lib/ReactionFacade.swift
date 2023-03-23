




/**
* The various roles a reaction component can have
* @sa OBReactionFacade
*/
enum MKReactionRole: UInt {
    case NO_REACTIONROLE = 0 //!< no reaction role - useful for temporarily hiding a component
    case REACTANT        = 1 //!< reactant
    case AGENT           = 2 //!< agent, a term that includes solvents and catalysts 
    case PRODUCT         = 3 //!< product
}
/**
* \class OBReactionFacade reactionfacade.h <openbabel/reactionfacade.h>
* \brief Facade to simplify manipulation of reactions stored as OBMol objects
*
* All of the information defining a reaction is stored as part of an OBMol.
* This information is:
* - a flag indicating that the OBMol represents a reaction (see OBMol::SetIsReaction(),
*   and OBMol::IsReaction())
* - #OBPairInteger data attached to each atom indicating the reaction role ("rxnrole")
* - #OBPairInteger data attached to each atom indicating the reaction component Id ("rxncomp")
*
* Everything that may need to be done with reactions can be implemented using this
* information. The purpose of this class is simply to provide convenience functions
* that cover a range of likely use cases. If your use case is not included, it should
* still be straightforward for you to manipulate the underlying information to
* achieve your goal.
*
* @since version 3.0
*/

class MKReactionFacade {

    /**
     * @name Constructor
     * @{
     */
    /**
     * The constructor requires an OBMol. Member functions retrieve or manipulate reaction
     * data stored in this molecule.
     */

    private var _d: MKReactionFacadePrivate

    init(_ mol: MKMol) {
        _d = MKReactionFacadePrivate(mol)
    }

    /**
     * @}
     * @name Low-level methods
     * These are convenience functions that set/get reaction-related data on individual atoms.
     * The actual data is stored in the underlying OBMol as OBPairInteger data with the
     * keys "rxnrole" or "rxncomp" (reaction role, reaction component).
     * The other methods of this class are implemented using these methods. Care should be
     * taken using the setter functions, as (unlike the other
     * methods) their use can result in an inconsistent state. IsValid() may be used to check
     * whether this is the case.
     * @{
     */

    /**
     * Assigns a component Id to every atom in the molecule based on connected components.
     * If @p wipe is \c false, then any atoms that already have a component Id are skipped.
     */
    func assignComponentIds(wipe: Bool = true) {
        _d.assignComponentIds(wipe)
    }
    
    /**
     * Return the reaction role of an atom.
     */
    func getRole(atom: MKAtom) -> MKReactionRole {
        return _d.getRole(atom)
    }
    
    /**
     * Return the component Id of an atom.
     */
    func getComponentId(atom: MKAtom) -> UInt {
        return _d.getComponentId(atom)
    }
    
    /**
     * Set the reaction role of an atom.
     */
    func setRole(atom: MKAtom, role: MKReactionRole) {
        _d.setRole(atom, role)
    }

    /**
     * Set the component Id of an atom.
     */
    func setComponentId(atom: MKAtom, compid: UInt) {
        _d.setComponentId(atom, compid)
    }
    
    /**
     * @}
     * @name Check validity of reaction data
     * @{
     */
    /**
     * Check whether this OBMol is in a valid state to represent a reaction.
     *
     * Note, first of all, that this does not check whether a reaction is
     * chemically valid. It simply checks whether the OBMol has the following
     * properties:
     * - (a) it is marked as a reaction
     * - (b) every atom is marked as having a reaction role
     * - (c) every atom is marked as belonging to a particular reaction component
     * - (d) every atom in a connected component is marked as belong to the same
     *       reaction component and having the same reaction role
     *
     * If (a) is not true, while it does not affect any of the methods of the
     * OBReactionFacade, the reaction may not be written appropriately by an
     * output format.
     *
     * If (b), (c) or (d) is not true, then the behavior of the toolkit is undefined
     * when handling this molecule.
     *
     * While not checked by this function, it is also a good idea to avoid reusing 
     * the same component Id for different reaction roles, e.g. having a reactant
     * with component Id 1 as well as a product with component Id 1. If you reassign
     * roles, then the duplication of component Ids may cause two components to
     * be merged into one.
     *
     */
    func isValid() -> Bool {
        return _d.isValid()
    }
    
    /**
     * @}
     * @name High-level methods
     * These are a set of methods that manipulate reactions at the component level.
     * If you use these high-level methods, and subsequently modify the underlying
     * OBMol or its reaction data, you may need to clear cached information on the
     * reaction components via ClearInternalState().
     * @{
     */
    /**
    * Add a new component to the reaction
    */
    func addComponent(mol: MKMol, role: MKReactionRole) {
        _d.addComponent(mol, role)
    }
    
    /**
    * Clear the list of found components.
    *
    * On first
    * use any of the high-level methods,
    * the number of components for each reaction role is
    * determined and each is assigned a number from 0 to the number of components.
    * This information is cached for future accesses. If you modify the underlying
    * OBMol in a way that invalidates this information (e.g. using the low-level
    * methods), you should call ClearInternalState() or create a new OBReactionFacade.
    */
    func clearInternalState() {
        _d.clearInternalState()
    }
    
    /**
    * Copy a component from the reaction into the provided OBMol.
    * Note that all
    * reaction data is copied, and so this may be useful for building up a new
    * reaction composed of just some of the components of the original reaction.
    * If copying components from multiple reactions into a single one, to avoid
    * duplicate component Ids it may be easier to copy into an empty OBMol, and
    * then use AddComponent() to add this OBMol into the new reaction.
    *
    * \return A boolean indicating whether the reaction has a component with that number.
    */
    func getComponent(mol: MKMol, role: MKReactionRole, num: UInt) -> Bool {
        return _d.getComponent(mol, role, num)
    }
    
    /**
    * Return the number of reaction components with a particular reaction role.
    */
    func numComponents(role: MKReactionRole) -> UInt {
        return UInt(_d.numComponents(role))
    }
    
    /**
    * Reassign the reaction role of a given component. If assigned to #NO_REACTIONROLE,
    * this may be used to hide the reaction component when writing to a reaction
    * output format.
    * 
    * It should be noted that reassigning a component changes the numbers assigned.
    * For example,
    * in a reaction with 2 reactants and 2 products (each numbered from 0 to 1), if you
    * reassign reactant number
    * 0 to be a product, it will now be product number 2 while former reactant number 1
    * will now be reactant number 0.
    * 
    * \return A boolean indicating whether the reaction has a component with that number.
    */
    func reassignComponent(oldrole: MKReactionRole, num: UInt, newrole: MKReactionRole) -> Bool {
        return _d.reassignComponent(oldrole, num, newrole)
    }

}


class MKReactionFacadePrivate {

    private var _mol: MKMol 
    private var _found_components: Bool = false 
    private var _unassigned_components: [UInt] = []
    private var _reactant_components: [UInt] = []
    private var _product_components: [UInt] = []
    private var _agent_components: [UInt] = []


    init(_ mol: MKMol) {
        _mol = mol
    }

    private func getId(_ idtype: String, _ atom: MKAtom) -> UInt {
        var idval: UInt = 0
        if let atomData = atom.getData(idtype) as? MKPairData<UInt> {
            idval = atomData.getValue()!
        } else {
            print("ERROR: No data for \(idtype) in atom \(atom)")
        }
        return idval
    }

    private func setId(_ idtype: String, _ atom: MKAtom, _ idval: UInt) {
        if let atomData = atom.getData(idtype) as? MKPairData<UInt> {
            atomData.setValue(idval)
        } else {
            let pi = MKPairData<UInt>()
            pi.setAttribute(idtype)
            pi.setValue(idval)
            atom.setData(pi)
        }
    }

    private func findComponents() {
        var reactant_components = Set<UInt>()
        var product_components = Set<UInt>()
        var agent_components = Set<UInt>()
        var unassigned_components = Set<UInt>()
        for atom in _mol.getAtomIterator() {
            let component = getComponentId(atom)
            switch getRole(atom) {
            case .REACTANT:
                reactant_components.insert(component)
            case .PRODUCT:
                product_components.insert(component)
            case .AGENT:
                agent_components.insert(component)
            default:
                unassigned_components.insert(component)
            }
        }
        // Convert to vector so we have random access - note: these will be in sorted order
        for component in reactant_components {
            _reactant_components.append(component)
        }
        for component in product_components {
            _product_components.append(component)
        }
        for component in agent_components {
            _agent_components.append(component)
        }
        for component in unassigned_components {
            _unassigned_components.append(component)
        }
        _found_components = true
    }

    private func getComponentIds(_ rxnrole: MKReactionRole, _ roles : (inout [UInt]) -> ()) {
        if !_found_components {
            findComponents()
        }
        switch rxnrole {
        case .NO_REACTIONROLE:
            roles(&self._unassigned_components)
        case .REACTANT:
            roles(&self._reactant_components)
        case .AGENT:
            roles(&self._agent_components)
        case .PRODUCT:
            roles(&self._product_components)
        }
    }

    private func getComponentIds(_ rxnrole: MKReactionRole) -> [UInt] {
        if !_found_components {
            findComponents()
        }
        switch rxnrole {
        case .NO_REACTIONROLE:
            return self._unassigned_components
        case .REACTANT:
            return self._reactant_components
        case .AGENT:
            return self._agent_components
        case .PRODUCT:
            return self._product_components
        }
    }

    func addComponent(_ mol: MKMol, _ rxnrole: MKReactionRole) {
        if !_found_components {
            findComponents()
        }

        var max_compid: UInt = 0
        if !_product_components.isEmpty {
            max_compid = _product_components.last!
        }
        if !_agent_components.isEmpty && _agent_components.last! > max_compid {
            max_compid = _agent_components.last!
        }
        if !_reactant_components.isEmpty && _reactant_components.last! > max_compid {
            max_compid = _reactant_components.last!
        }
        if !_unassigned_components.isEmpty && _unassigned_components.last! > max_compid {
            max_compid = _unassigned_components.last!
        }

        var new_compid = max_compid + 1
        if new_compid == 0 {
            new_compid = 1
        }

        for atom in mol.getAtomIterator() {
            setRole(atom, rxnrole)
            setComponentId(atom, new_compid)
        }
        _mol += mol

        getComponentIds(rxnrole) { $0.append(new_compid) } 
    }

    func assignComponentIds(_ wipe: Bool = true) {
        var compid: UInt = 1
        var dfs_iter = MKAtomDFSIterator(_mol)
        // This code loops each stack set, when it is empty, it briefly returns nil, 
        // then calls ++ to update the stack and ptr to the next set.
        while !dfs_iter.isEmpty() { // for each connected component
            repeat {
                guard let atom = dfs_iter.current() else { break }
                if wipe || !atom.hasData("rxncomp") {
                    setComponentId(atom, compid)
                }
            } while dfs_iter.next() != nil // next advances the stack ptr
            dfs_iter++ // updates the stack with neighbors and resets the ptr from nil
            compid += 1
        } 
    }

    func clearInternalState() {
        _found_components = false
    }

    func getComponent(_ mol: MKMol, _ rxnrole: MKReactionRole, _ num: UInt) -> Bool {
        var component_ids = getComponentIds(rxnrole)
        if num >= component_ids.count {
            return false
        }
        let componentId = component_ids[Int(num)]
        var atoms = MKBitVec()
        for atom in _mol.getAtomIterator() {
            if getRole(atom) == rxnrole && getComponentId(atom) == componentId {
                atoms.setBitOn(atom.getIdx())
            }
        }
        var atomO: [Int]? = nil
        var bondO: [Int]? = nil
        let ok = _mol.copySubstructure(mol, atoms, &atomO, &bondO)
        return ok
    }

    func getComponentId(_ atom: MKAtom) -> UInt {
        return getId("rxncomp", atom)
    }

    func getRole(_ atom: MKAtom) -> MKReactionRole {
        let rxnrole = getId("rxnrole", atom)
        switch rxnrole {
        case 0: return .NO_REACTIONROLE
        case 1: return .REACTANT
        case 2: return .AGENT
        case 3: return .PRODUCT
        default: return .NO_REACTIONROLE
        }
    }

    func isValid() -> Bool {
        if !_mol.isReaction() {
            print("The molecule is not marked as a reaction. Use setIsReaction().")
            return false
        }
        for atom in _mol.getAtomIterator() {
            if let data = atom.getData("rxncomp") {
                if let pi = data as? MKPairData<Int> {
                    if let val = pi.getValue() {
                        if val <= 0 {
                            print("Reaction component Ids should all be non-zero positive integers.")
                            return false
                        }
                    } else {
                        print("A reaction component Id has been stored using a data type that is not an MKPairData<Int>.")
                        return false
                    }
                } else {
                    print("A reaction component Id has been stored using a data type that is not an MKPairData<Int>.")
                    return false
                }
            } else {
                print("The molecule contains an atom that is missing a reaction component Id. Use setComponentId().")
                return false
            }

            if let data = atom.getData("rxnrole") {
                if let pi = data as? MKPairData<Int> {
                    if let val = pi.getValue() {
                        if val < 0 || val > 3 {
                            print("Reaction roles should be in the range 0 to 3 inclusive.")
                            return false
                        }
                    } else {
                        print("A reaction role has been stored using a data type that is not an MKPairData<Int>.")
                        return false
                    }
                } else {
                    print("Reaction role information has been stored using a data type that is not an MKPairData<Int>.")
                    return false
                }
            } else {
                print("The molecule contains an atom that is missing a reaction role information. Use setRole().")
                return false
            }
        }

        // Ensure that every atom in a particular connected component has the same component Id and rxn role
        var dfs_iter = MKAtomDFSIterator(_mol)
        while !dfs_iter.isEmpty() { // for each connected component
            guard let atom = dfs_iter.current() else { break }
            let rxncomp = getComponentId(atom)
            let rxnrole = getRole(atom)
            repeat { // for each atom in connected component
                guard let inner_atom = dfs_iter.current() else { break }
                let compid = getComponentId(inner_atom)
                if compid != rxncomp {
                    print("The molecule contains a connected component that contains atoms with different reaction component Ids. All atoms in a particular connected component should have the same value.")
                    return false
                }
                let roleid = getRole(inner_atom)
                if roleid != rxnrole {
                    print("The molecule contains a connected component that contains atoms with different reaction roles. All atoms in a particular connected component should have the same role.")
                    return false
                }
            } while dfs_iter.next() != nil // next advances the stack ptr
            dfs_iter++ // updates the stack with neighbors and resets the ptr from nil
        }
        return true
    }

    func numComponents(_ rxnrole: MKReactionRole) -> Int {
        return getComponentIds(rxnrole).count
    }

    @discardableResult
    func reassignComponent(_ oldrole: MKReactionRole, _ num: UInt, _ newrole: MKReactionRole) -> Bool {
        var componentId: UInt?
        getComponentIds(oldrole) { component_ids in
            if num >= component_ids.count {
                return
            }
            componentId = component_ids[Int(num)]
            for atom in _mol.getAtomIterator() {
                if getRole(atom) == oldrole && getComponentId(atom) == componentId {
                    setRole(atom, newrole)
                }
            }
            // remove the entry from the original component_ids and update the new one
            component_ids.remove(at: Int(num))
        }
        getComponentIds(newrole) { component_ids in
            if let cI = componentId {
                component_ids.append(cI)
            }
        }
        return true
    }

    func setComponentId(_ atom: MKAtom, _ compid: UInt) {
        setId("rxncomp", atom, compid)
    }

    func setRole(_ atom: MKAtom, _ rxnrole: MKReactionRole) {
        setId("rxnrole", atom, rxnrole.rawValue)
    }

}

