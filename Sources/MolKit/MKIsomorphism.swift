import Foundation
import Bitset


/**
* @class OBIsomorphismMapper isomorphism.h <openbabel/isomorphism.h>
* @since version 2.3
* @brief Abstract class defining interface for isomorphism (i.e. substructure) searches.
*
* The OBIsomorphism class is an abstract class which defines an interface for
* performing isomorphism (i.e. substructure) searches. It uses a OBQuery and
* tries to map this onto a queried OBMol. A single mapping is represented by
* a OBIsomorphismMapper::Mapping which is a std::map mapping query indexes to
* queried indexes. Both query and queried indexes in the map start from 0.
* Multiple mappings can be stored in a OBIsomorphismMapper::Mappings object
* which is a std::vector of OBIsomorphismMapper objects.
*
* Since this is an abstract class with pure virtual methods, this class can't
* be instantiated directly. To get a pointer to a subclass, the GetInstance()
* method can be used which also sets the query. Once an instance is obtained,
* the desired mapping function can be used to perform the mapping (i.e. MapFirst(),
* MapUnique() or MapAll()).
*
* A typical example:
* @code
* OBMol *queried;
* // ... initialize queried ...
* OBQuery *query = CompileSmilesQuery("c1ccccc1");
* OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
* OBIsomorphismMapper::Mappings maps = mapper->MapUnique(mol);
*
* std::cout << "found " << maps.size() << " unique mappings" << std::endl;
*
* delete mapper;
* delete query;
* @endcode
*
* All mapping methods take an optional mask parameter. This can be used to
* restrict the search to a part of the queried OBMol. The masked atoms in
* the OBBitVec are  indexed from 1. A special case of isomorphism search is
* an automorphism search where the query and queried molecule are the same.
* Automorphism searches can be done using the MapAll method but an additional
* FindAutomorphisms() function is provided for convenience.
*/

let red    = "\033[1;31m"
let green  = "\033[1;32m"
let yellow = "\033[1;33m"
let blue   = "\033[1;34m"
let normal = "\033[0m"

class MKIsomorphMapper {
    
    /**
     * @typedef std::vector< std::pair<unsigned int,unsigned int> > Mapping
     * Type for an individual mapping.
     */
    typealias Mapping = [Pair<UInt, UInt>]
    
    /**
     * @typedef std::vector<OBIsomorphismMapper::Mapping> Mappings
     * Type for a collection (std::vector) of Mapping objects.
     */
    typealias Mappings = [Mapping]
    
    /**
     * Constructor. OBIsomorphismMapper is an abstract class, use GetInstance() to
     * get an instance of a derived class.
     * @param query The search query.
     */
    
    var m_query: MKQuery //!< The search query.
    var m_timeout: TimeInterval = 0  //!< The timeout in seconds
    
    init(_ query: MKQuery) {
        self.m_query = query
        self.m_timeout = 60.0
    }
    
    func setTimeout(_ seconds: TimeInterval) {
        m_timeout = seconds
    }
    
    /**
     * Get a pointer to an instance of the specified @p algorithm. This pointer
     * has to be delted when the instance is no longer needed.
     * @param query The search query to be mapped.
     * @param algorithm The algorithm for the mapper.
     * @return OBIsomorphismMapper instance or 0 if there is no subclass implementing
     * the specified @p algorithm.
     */
    //   static OBIsomorphismMapper* GetInstance(OBQuery *query, const std::string &algorithm = std::string("VF2"));
    static func getInstance(_ query: MKQuery, _ algorithm: String = "VF2") -> MKIsomorphMapper? {
        // if algorithm == "VF2":
        return VF2Mapper(query)
        // return VF2 mapper as default
        // return 
    }
    
    /**
     * Find a single mapping in @p queried.
     * @param queried The molecule to search.
     * @param map Reference to the object to store the result in.
     * @param mask A mask to restrict the search to a part of the queried molecule.
     * The default empty mask will result in all atoms being considered. The mask
     * indexes start from 1 (i.e. OBAtom::GetIdx()).
     */
    open func mapFirst(_ queried: MKMol, _ map: inout Mapping, _ mask: Bitset = Bitset()) {
        fatalError("Not implemented in base class")
    }
    /**
     * Find all unique mappings in @p queried. A mapping is unique when there is no previous
     * mapping covering the same queried atoms. For two mappings, some overlap is allowed but
     * at least one atom should be different.
     * @param queried The molecule to search.
     * @param maps Reference to the object to store the results in.
     * @param mask A mask to restrict the search to a part of the queried molecule.
     * The default empty mask will result in all atoms being considered. The mask
     * indexes start from 1 (i.e. OBAtom::GetIdx()).
     */
    open func mapUnique(_ queried: MKMol, _ maps: inout Mappings, _ mask: Bitset = Bitset()) {
        fatalError("Not implemented in base class")
    }
    /**
     * Find all mappings in @p queried. This function is used by FindAutomorphisms()
     * with a query that is a copy of the queried molecule (taking the mask into
     * account).
     * @param queried The molecule to search.
     * @param maps Reference to the object to store the results in.
     * @param mask A mask to restrict the search to a part of the queried molecule.
     * @param maxMemory Memory limit for the @p maps object in bytes. Default is 300MB.
     * The default empty mask will result in all atoms being considered. The mask
     * indexes start from 1 (i.e. OBAtom::GetIdx()).
     */
    open func MapAll(_ queried: MKMol, _ maps: inout Mappings, _ mask: Bitset = Bitset(), _ maxMemory: Int = 3000000) {
        fatalError("Not implemented in base class")
    }
    
    /**
     * @class Functor isomorphism.h <openbabel/isomorphism.h>
     * @brief Functor base class to be used in combination with MapGeneric.
     * @see @ref MapGeneric
     * @since 2.3
     */
    open class Functor {
        /**
         * This function is called every time an isomorphism is discovered.
         * Returing true, will abort the mapping process. The map is passed
         * as non-const reference and may be modified (e.g. swap).
         *
         * @see @ref MapGeneric
         */
        init () {}
        
        @discardableResult
        open func operate(_ map: Mapping) -> Bool {
            fatalError("Not implemented in base class")
        }
    }
    /**
     * Find all mappings in @p queried. The functor will be called when a mapping is found.
     * @param functor The functor to handle found mappings.
     * @param queried The molecule to search.
     * @param mask A mask to restrict the search to a part of the queried molecule.
     * The default empty mask will result in all atoms being considered. The mask
     * indexes start from 1 (i.e. OBAtom::GetIdx()).
     */
    open func mapGeneric(_ functor: inout MKIsomorphMapper.Functor, _ queried: MKMol, _ mask: Bitset = Bitset()) {
        fatalError("Not implemented in base class")
    }
    
}

//@inlinable
func mapsTo(_ map: MKIsomorphMapper.Mapping, _ queryIndex: UInt, _ queriedIndex: inout UInt) -> Bool {
    for (q, qd) in map {
        if q == Int(queryIndex) {
            queriedIndex = UInt(qd)
            return true
        }
    }
    return false
}

/**
* @typedef OBIsomorphismMapper::Mapping Automorphism
* @brief A single automorphic permutation.
* @since 2.3
*/
typealias Automorphism = MKIsomorphMapper.Mapping

/**
* @typedef OBIsomorphismMapper::Mappings Automorphisms
* @brief A group of automorphic permutations.
* @since 2.3
*/
typealias Automorphisms = MKIsomorphMapper.Mappings

/**
* Find the automorphisms of a molecule by using an OBIsomorphismMapper. This
* function is a wrapper for FindAutomorphisms with a functor to store all
* automorphisms.
*
* @param mol The molecule for which to find the automorphisms.
* @param aut The result will be stored here.
* @param symmetry_classes The graph invariants to use. See OBGraphSym or use
* the FindAutomorphisms function that computes these for you below.
* @param mask A bit vector specifying which atoms to consider. An empty mask
* will consider all atoms. The bits are indexed from 1 (i.e. OBAtom::GetIdx()).
* @param maxMemory Maximum memory limit for @p aut. The number of automorphisms
* for certain graphs can be large. The default is 300MB, consider using a functor
* to process automorphisms when they are found.
*
* @since version 2.3
*/
func findAutomorphisms(_ mol: MKMol, _ maps: inout MKIsomorphMapper.Mappings, _ symmetry_classes: [UInt], _ mask: Bitset = Bitset(), _ maxMemory: Int = 3000000) -> Bool {
    maps.removeAll()
    var functor: MKIsomorphMapper.Functor = MapAllFunctor(maps, maxMemory)
    findAutomorhphisms(&functor, mol, symmetry_classes, mask)
    return !maps.isEmpty
}
/**
* Find the automorphisms of a molecule by using an OBIsomorphismMapper. This
* function will first find the graph invariants (i.e. symmetry_classes) using
* the mask. This function is a wrapper for FindAutomorphisms with a functor to
* store all automorphisms.
*
* @since version 2.3
*/
func findAutomorphisms(_ mol: MKMol, _ aut: inout MKIsomorphMapper.Mappings, _ mask: inout Bitset, _ maxMemory: Int = 3000000) -> Bool {
    // set all atoms to 1 if the mask is empty
    var queriedMask: Bitset = mask
    if queriedMask.count() == 0 {
        for i in 0..<mol.numAtoms() {
            queriedMask.add(i + 1)
        }
    }
    // get the symmetry classes
    let gs = MKGraphSym(mol, &queriedMask)
    var symClasses = [Ref]()
    gs.getSymmetry(&symClasses)
    return findAutomorphisms(mol, &aut, symClasses.map {UInt($0.intValue)}, queriedMask, maxMemory) // TODO: Check if using queriedMask here is bad. In default impl they had mask instead
}
/**
* Find the automorphisms of a molecule by using an OBIsomorphismMapper. This
* is the main implementation for finding automorphisms and uses an
* OBIsomorphismMapper::Functor to process found isomorphisms. Wrapper functions
* are provided to find and store all automorphisms but the number of automorphisms
* can be large for certain graphs making it not feasible to store all automorphisms
* in memory (RAM).
*
* @see  @ref MapGeneric
* @since version 2.3
*/
func findAutomorhphisms(_ functor: inout MKIsomorphMapper.Functor, _ mol: MKMol, _ symmetry_classes: [UInt], _ mask: Bitset = Bitset()) {
    // set all atoms to 1 if the mask is empty
    var queriedMask = mask
    if queriedMask.count() == 0 {
        for i in 0..<mol.numAtoms() {
            queriedMask.add(i + 1)
        }
    }
    // compute the connected fragments
    var visited: Bitset = Bitset()
    var fragments: [Bitset] = []

    for i in 0..<mol.numAtoms() {
        if !queriedMask.contains(i + 1) || visited.contains(i + 1) {
            continue
        }
        fragments.append(getFragment(mol.getAtom(i + 1)!, &queriedMask))
        visited |= fragments.last!
    }
    // count the symmetry classes
    var symClassCounts: [UInt] = [UInt].init(repeating: 0, count: symmetry_classes.count + 1)
    for i in 0..<symmetry_classes.count {
        if !queriedMask.contains(i + 1) {
            continue
        }
        let symClass = symmetry_classes[i]
        symClassCounts[Int(symClass)] += 1
    }
    
    for i in 0..<fragments.count {
        let query = compileAutomorphismQuery(mol, symmetry_classes, fragments[i])
        guard let mapper = MKIsomorphMapper.getInstance(query) else {
            print("ERROR: in generating MKIsomorph instance")
            fatalError()
        }
        var autFunctor: MKIsomorphMapper.Functor = AutomorphismFunctor(functor, fragments[i], UInt(mol.numAtoms()))
        mapper.mapGeneric(&autFunctor, mol, fragments[i])
    }
}

// Helper Declarations

class AutomorphismFunctor: MKIsomorphMapper.Functor {
    private var m_functor: MKIsomorphMapper.Functor
    private var m_fragment: Bitset
    private var m_indexes: [UInt] = []
    
    init(_ functor: MKIsomorphMapper.Functor, _ fragment: Bitset, _ numAtoms: UInt) {
        m_functor = functor
        m_fragment = fragment
        for j in 0..<numAtoms {
            if m_fragment.contains(Int(j) + 1) {
                m_indexes.append(j)
            }
        }
    }
    override func operate(_ map: MKIsomorphMapper.Mapping) -> Bool {
        // convert the continuous mapping map to a mapping with gaps (considering key values)
        for var iter in map {
            iter.0 = m_indexes[Int(iter.0)]
        }
        return m_functor.operate(map)
    }
}

func compileAutomorphismQuery(_ mol: MKMol, _ symmetry_classes: [UInt], _ mask: Bitset = Bitset()) -> MKQuery {
    let query: MKQuery = MKQuery()
    var offset: UInt = 0
    var indexes = [UInt]()
    for atom in mol.getAtomIterator() {
        indexes.append(UInt(atom.getIndex()) - offset)
        if !mask.contains(atom.getIndex() + 1) {
            offset += 1
            continue
        }
        query.addAtom(MKAutomorphismQueryAtom(symmetry_classes[Int(atom.getIndex())], symmetry_classes))
    }
    for bond in mol.getBondIterator() {
        if isFerroceneBond(bond) {
            continue
        }
        let beginIndex = bond.getBeginAtom().getIndex()
        let endIndex = bond.getEndAtom().getIndex()
        if !mask.contains(beginIndex + 1) || !mask.contains(endIndex + 1) {
            continue
        }
        query.addBond(MKQueryBond(query.getAtoms()[Int(indexes[beginIndex])], query.getAtoms()[Int(indexes[endIndex])], Int(bond.getBondOrder()), bond.isAromatic()))
    }

    return query;
}

/**
* @page substructure Substructure Search
* @since version 2.3
*
* Substructure searching is finding a mapping for a query to a target molecule.
* Such a mapping is also known as a graph isomorphism. A graph isomorphism maps
* the vertexes (i.e. atoms) from the query to vertexes in the molecule such that
* two vertexes adjacent in the query are also adjacent in the target molecule.
* In other words, no bonds are broken and no new bonds are formed.
*
* @section smarts SMARTS Substructure Search
* Smarts is an extension of smiles to create powerful queries. Smarts substructure
* search has been available in OpenBabel for many years. It is also used for many
* of OpenBabel's algorithms. Although smarts is only a syntax for queries, the
* implementation has it's own matching algorithm. For many purposes smarts are the
* easiest way to do substructure searches. See the OBSmartsPattern documentation
* for details on how to use smarts.
*
* @section query Queries
* Since OpenBabel version 2.3, there are some classes for representing generic
* queries. The OBQuery, OBQueryAtom and OBQueryBond class define interfaces that
* can be reimplemented to get custom behavior. The classes also contain some
* methods to access topological information which are used by the mapping
* algorithms. The default implementations allow very simple exact substructure
* searches to be performed but subclassing allows very advanced queries to be
* used (e.g. smarts).
*
* While it is possible to construct these queries manually, "compilers" are
* provided to convert a query representation to a OBQuery object. Currently,
* only two exact substructure search compilers exist. The first is
* CompileMoleculeQuery which converts an OBMol object to an OBQuery object.
* The second is CompileSmilesQuery and converts a smiles string to an OBQuery
* object.
*
* @code
* #include <openbabel/query.h>
* using namespace OpenBabel;
*
* OBMol *mol = new OBMol;
*
* // ... read molecule ...
*
* OBQuery *query;
* query = CompileMoleculeQuery(mol);
* query = CompileSmilesQuery("c1ccccc1CC(=O)O");
* @endcode
*
* @section mapping Mapping Isomorphisms
* The OBIsomorphismMapper class defined an interface for mapping queries to
* target molecules. Multiple implementations can be added but they all do the
* same. The MapFirst, MapUnique and MapAll methods are used for gettings the
* map(s).
*
* @subsection MapFirst
* This method returns the first map found. The main reason for getting only
* one map is improved performance since it is considerably faster than
* MapUnique and MapAll. However, depending on the use case a single map is
* all that is needed. For example, to check if a molecule in a database
* contains a phenyl ring, a single mapping is enough.
*
* @subsection MapUnique
* MapUnique returns all unique maps. A map is considered unique if there is
* no other map covering exactly the same atoms in the target molecule. For
* example, when a phenyl query is performed on a molecule with 2 phenyl rings,
* MapUnique will return 2 maps. These 2 maps are selected from the 24 found
* non-duplicate maps (6 atoms to start from * 2 directions (CW/CCW) * 2 rings).
*
* @subsection MapAll
* MapAll returns all non-duplicate maps. For example, when a phenyl query is
* performed on a molecule with 2 phenyl rings, MapAll will return 24 maps
* (6 atoms to start from * 2 directions (CW/CCW) * 2 rings).
*
* @subsection MapGeneric
* MapGeneric takes a functor object and calls the functor to handle found
* isomorphisms. This allows for custom mapping results to be obtained by
* filtering or avoid storing all results. To implement a custom functor,
* a simple class that inherits OBIsomorphismMapper::Functor and implements
* the required operator().
*
* @code
* #include <openbabel/isomorphism.h>
* using namespace OpenBabel;
*
* class MyCustomFunctor : public OBIsomorphismMapper::Functor
* {
*   private:
*     // store all mappings in m_data
*     std::vector<OBIsomorphismMapper::Mapping> &m_data;
*   public:
*     MyCustomFunctor(std::vector<OBIsomorphismMapper::Mapping> &data) : m_data(data) {}
*     bool operator()(OBIsomorphismMapper::Mapping &map)
*     {
*       // code to handle found isomorphism...
*       // examples: - store the mapping
*       //           - filter mappings
*       //           - use the found map in some way
*
*       m_data.push_back(map);
*
*       // continue mapping
*       return false;
*     }
* }
* @endcode
*
* @section automorphisms Automorphisms
* The automorphisms of a graph or molecule are a group of isomorphism mappings
* of the molecule onto itself (i.e. the query and target are the same). The
* automorphisms make it easy to take symmetry into account. See FindAutomorphisms
* for detials.
*
*
*/

class MapAllFunctor: MKIsomorphMapper.Functor {

    var m_maps: MKIsomorphMapper.Mappings = []
    var m_memory: Int = 0
    var m_maxMemory: Int = 0 

    init(_ maps: MKIsomorphMapper.Mappings, _ maxMemory: Int) {
        m_maps = maps
        m_memory = 0
        m_maxMemory = maxMemory
    }

    override func operate(_ map: MKIsomorphMapper.Mapping) -> Bool {
        m_maps.append(map)
        m_memory += 2 * sizeof(map)
        if m_memory > m_maxMemory {
            print("ERROR: memory limit exceeded...")
            return true
        }
        // continue mapping
        return false
    }

}

class VF2Mapper: MKIsomorphMapper {

    struct Candidate: Equatable {
        var queryAtom: MKQueryAtom?
        var queriedAtom: MKAtom?

        init () {}

        init(_ queryAtom: MKQueryAtom, _ queriedAtom: MKAtom) {
            self.queryAtom = queryAtom
            self.queriedAtom = queriedAtom
        }

        static func == (lhs: VF2Mapper.Candidate, rhs: VF2Mapper.Candidate) -> Bool {
            return lhs.queryAtom == rhs.queryAtom && lhs.queriedAtom == rhs.queriedAtom
        }

    }

    struct State {

        var abort: Bool 
        var functor: MKIsomorphMapper.Functor
        var query: MKQuery                       // the query
        var queried: MKMol                       // the queried molecule
        var queriedMask: Bitset                // the queriedMask
        var queryPath: [UInt]                    // the path in the query
        var queriedPath: [UInt]                  // the path in the queried molecule
        var mapping: [MKAtom?]
        // the terminal sets
        var queryPathBits: Bitset
        var queriedPathBits: Bitset
        var queryDepths: [UInt]
        var queriedDepths: [UInt]

        init(_ functor: inout MKIsomorphMapper.Functor, _ query: MKQuery, _ queried: MKMol, _ mask: Bitset) {
            self.abort = false
            self.functor = functor
            self.query = query
            self.queried = queried
            self.queriedMask = mask
            self.queryPath = []
            self.queriedPath = []
            self.mapping = [MKAtom?](repeating: nil, count: query.numAtoms())
            self.queryPathBits = Bitset()
            self.queriedPathBits = Bitset()
            self.queryDepths = [UInt](repeating: 0, count: query.numAtoms())
            self.queriedDepths = [UInt](repeating: 0, count: queried.numAtoms())
        }

    }

    var m_startTime: Date

    override init(_ query: MKQuery) {
        m_startTime = Date()
        super.init(query)
    }

    func isInTerminalSet(_ depths: [UInt], _ path: Bitset, _ i: UInt) -> Bool {
        if depths[Int(i)] == 0 {
            return false
        }
        if path.contains(Int(i)) {
            return false
        }
        return true
    }

    /**
    * Check bonds around newly mapped atom.
    */
    func checkBonds(_ state: inout State, _ queryAtom: MKQueryAtom) -> Bool {
        let qbonds = queryAtom.getBonds()
        for qbond in qbonds {
            let beginIndex = qbond.getBegin().getIndex()
            let endIndex = qbond.getEnd().getIndex()
            let begin = state.mapping[Int(beginIndex)]
            let end = state.mapping[Int(endIndex)]
            if begin == nil || end == nil {
                continue
            }
            let bond = state.queried.getBond(begin!, end!)
            if bond == nil {
                return false
            }
            if !qbond.matches(bond!) {
                return false
            }
        }
        return true
    }
    
    /**
    * Check if the current state is a full mapping of the query.
    */
    func checkForMap(_ state: inout State) -> Bool {
        if state.queryPath.count != state.query.numAtoms() {
            return false
        }
        var map = MKIsomorphMapper.Mapping()
        map.reserveCapacity(state.queryPath.count)
        for k in 0..<state.queryPath.count {
            map.append((state.queryPath[k], state.queriedPath[k]))
        }
        return state.functor.operate(map)
    }

    /**
    * Match the candidate atoms and bonds.
    */
    func matchCandidate(_ state: inout State, _ queryAtom: MKQueryAtom, _ queriedAtom: MKAtom) -> Bool {
        if !queryAtom.matches(queriedAtom) {
            return false
        }
        // add the neighbors to the paths
        state.queryPath.append(UInt(queryAtom.getIndex()))
        state.queriedPath.append(UInt(queriedAtom.getIndex()))
        // update the terminal sets
        state.queryPathBits.add(queryAtom.getIndex())
        state.queriedPathBits.add(queriedAtom.getIndex())
        // update mapping
        state.mapping[Int(queryAtom.getIndex())] = queriedAtom
        //
        // update queryDepths
        //
        if (state.queryDepths[Int(queryAtom.getIndex())] == 0) {
            state.queryDepths[Int(queryAtom.getIndex())] = UInt(state.queryPath.count)
        }
        let queryNbrs = queryAtom.getNbrs()
        for i in 0..<queryNbrs.count {
            let index = queryNbrs[i].getIndex()
            if (state.queryDepths[Int(index)] == 0) {
                state.queryDepths[Int(index)] = UInt(state.queryPath.count)
            }
        }
        //
        // update queriedDepths
        //
        if (state.queriedDepths[Int(queriedAtom.getIndex())] == 0) {
            state.queriedDepths[Int(queriedAtom.getIndex())] = UInt(state.queriedPath.count)
        }
        guard let nbors = queriedAtom.getNbrAtomIterator() else { return false }
//        TODO: Throw error here
        for nbor in nbors {
            let index = nbor.getIndex()
            if !state.queriedMask.contains(Int(index + 1)) {
                continue
            }
            if state.queriedDepths[Int(index)] == 0 {
                state.queriedDepths[Int(index)] = UInt(state.queriedPath.count)
            }
        }
        // check if the bonds match
        if (!checkBonds(&state, queryAtom)) {
            backtrack(state: &state)
            return false
        }
        //
        // Feasibility rules for the VF2 algorithm:
        //
        //  size( T1(s) ) <= size( T2(s) )
        //
        //  size( N1 - M1(s) - T1(s) ) <= size( N2 - M2(s) - T2(s) )
        //
        // compute T1(s) size
        var numT1 = 0
        for i in 0..<state.query.numAtoms() {
            if isInTerminalSet(state.queryDepths, state.queryPathBits, UInt(i)) {
                numT1 += 1
            }
        }
        // compute T2(s) size
        var numT2 = 0
        for i in 0..<state.queried.numAtoms() {
            if isInTerminalSet(state.queriedDepths, state.queriedPathBits, UInt(i)) {
                numT2 += 1
            }
        }
        // T1(s) > T()
        if numT1 > numT2 {
            backtrack(state: &state)
            return false
        }
        //  N1 - M1(s) - T1(s) > N2 - M2(s) - T2(s)
        if (state.query.numAtoms() - state.queryPath.count - numT1) > (state.queried.numAtoms() - state.queriedPath.count - numT2) {
            backtrack(state: &state)
            return false
        }
        // Check if there is a mapping found
        state.abort = checkForMap(&state)
        return true
    }

    func nextCandidate(_ state: inout State, _ lastCandidate: inout Candidate) -> Candidate {
        var lastQueryAtom = lastCandidate.queryAtom?.getIndex() ?? 0
        var lastQueriedAtom = (lastCandidate.queriedAtom != nil) ? lastCandidate.queriedAtom!.getIndex() + 1 : 0

        let querySize = state.query.numAtoms()
        let queriedSize = state.queried.numAtoms()

        let queryTerminalSize = state.queryDepths.count - state.queryDepths.filter({ $0 == 0 }).count
        let queriedTerminalSize = state.queriedDepths.count - state.queriedDepths.filter({ $0 == 0 }).count

        let mappingSize = state.queryPath.count

        if (queryTerminalSize > mappingSize && queriedTerminalSize > mappingSize) {
            while (lastQueryAtom < querySize && (state.queryPathBits.contains(lastQueryAtom) || (state.queryDepths[lastQueryAtom] == 0))) {
                lastQueryAtom += 1
                lastQueriedAtom = 0 
            }
        } else {
            while(lastQueryAtom < querySize && state.queryPathBits.contains(lastQueryAtom)) {
                lastQueryAtom += 1
                lastQueriedAtom = 0
            }
        }

        if (queryTerminalSize > mappingSize && queriedTerminalSize > mappingSize) {
            while (lastQueriedAtom < queriedSize && (state.queriedPathBits.contains(lastQueriedAtom) || (state.queriedDepths[lastQueriedAtom] == 0))) {
                lastQueriedAtom += 1
            }
        } else {
            while(lastQueriedAtom < queriedSize && state.queriedPathBits[lastQueriedAtom]) {
                lastQueriedAtom += 1
            }
        }

        if (lastQueryAtom < querySize && lastQueriedAtom < queriedSize) {
            return Candidate(state.query.getAtoms()[lastQueryAtom], state.queried.getAtom(lastQueriedAtom + 1)!)
        }
        
        return Candidate()
    }

    /**
    * The depth-first isomorphism algorithm.
    */
    func mapNext(_ state: inout State, _ queryAtom: MKQueryAtom, _ queriedAtom: MKAtom) {
        if (Date().timeIntervalSince(m_startTime) > m_timeout) {
            return
        }
        if (state.abort) {
            return
        }
        var candidate = Candidate() 
        while (!state.abort) {
            candidate = nextCandidate(&state, &candidate)
            if (candidate.queryAtom == nil) {
                return
            }
            if (matchCandidate(&state, candidate.queryAtom!, candidate.queriedAtom!)) {
                mapNext(&state, candidate.queryAtom!, candidate.queriedAtom!)
                backtrack(state: &state)
            }
        }
    }
   
    func backtrack(state: inout State) {
        // remove last atoms from the mapping
        if (state.queryPath.count > 0) {
            state.mapping[Int(state.queryPath.last!)] = nil
            state.queryPathBits.remove(Int(state.queryPath.last!))
            state.queryPath.removeLast()
        }
        if state.queriedPath.count > 0 {
            state.queriedPathBits.remove(Int(state.queriedPath.last!))
            state.queriedPath.removeLast()
        }
        // restore queryDepths and queriedDepths
        let depth: UInt = UInt(state.queryPath.count + 1)
        state.queryDepths.replace(depth, with: 0)
        state.queriedDepths.replace(depth, with: 0) // O(n)  n = # vertices in the queried
    }

    class MapFirstFunctor : MKIsomorphMapper.Functor {
        private var m_map: MKIsomorphMapper.Mapping
    
        init(_ map: MKIsomorphMapper.Mapping) {
            self.m_map = map
        }
        
        override func operate(_ map: MKIsomorphMapper.Mapping) -> Bool {
            self.m_map = map
            // stop mapping
            return true
        }
    }

    /**
    * Get the first mappings of the query on the queried molecule.
    * @param queried The queried molecule.
    * @return The mapping.
    */
    override func mapFirst(_ queried: MKMol, _ map: inout MKIsomorphMapper.Mapping, _ mask: Bitset = Bitset()) {
        var functor: MKIsomorphMapper.Functor = MapFirstFunctor(map)
        mapGeneric(&functor, queried, mask)
    }

    class MapUniqueFunctor: MKIsomorphMapper.Functor  {
          
            private var m_maps: MKIsomorphMapper.Mappings
            
            init(_ maps: MKIsomorphMapper.Mappings) {
                self.m_maps = maps
            }
        
            override func operate(_ map: MKIsomorphMapper.Mapping) -> Bool {
                // get the values from the map
                var values: [UInt] = []
                for (_, value) in map {
                    values.append(value)
                }
                values.sort()
                // print_vector("values ", values);

                var isUnique = true
                for k in 0..<m_maps.count {
                    var kValues: [UInt] = []
                    for (_, value) in m_maps[k] {
                        kValues.append(value)
                    }
                    kValues.sort()

                    //  print_vector("kValues", kValues);
                    if (values == kValues) {
                        isUnique = false
                    }
                }

                if (isUnique) {
                    m_maps.append(map)
                }

                // continue mapping
                return false
            }
        }

    /**
    * Get all unique mappings of the query on the queried molecule.
    * @param queried The queried molecule.
    * @return The unique mappings
    */
    override func mapUnique(_ queried: MKMol, _ maps: inout MKIsomorphMapper.Mappings, _ mask: Bitset = Bitset()) {
        maps.removeAll()
        var functor: MKIsomorphMapper.Functor = MapUniqueFunctor(maps)
        mapGeneric(&functor, queried, mask)
    }

    /**
    * Get all mappings of the query on the queried molecule. Duplicates are
    * ignored but unlinke MapUnique, multiple mappings of the query on the same
    * part of the queried structure are allowed. This makes it possible to use
    * MapAll for finding the automorphism group.
    * @param queried The queried molecule.
    * @return The mappings.
    */
    override func MapAll(_ queried: MKMol, _ maps: inout MKIsomorphMapper.Mappings, _ mask: Bitset = Bitset(), _ maxMemory: Int = 3000000) {
        maps.removeAll()
        var functor: MKIsomorphMapper.Functor = MapAllFunctor(maps, maxMemory)
        mapGeneric(&functor, queried, mask)
        // TODO: Maybe add debugging here if it is needed
        //     if (DEBUG)
        //       for (unsigned int i =0; i < maps.size(); ++i) {
        //         cout << "mapping:" << endl;
        //         for (Mapping::iterator it = maps[i].begin(); it != maps[i].end(); ++it)
        //           cout << "    " << it->first << " -> " << it->second << endl;
        //       }
    }

    override func mapGeneric(_ functor: inout MKIsomorphMapper.Functor, _ queried: MKMol, _ mask: Bitset = Bitset()) {
        m_startTime = Date()
        if (m_query.numAtoms() == 0) {
            return
        }
        // set all atoms to 1 if the mask is empty
        let queriedMask = mask
        if (queriedMask.count() == 0) {
            for i in 0..<queried.numAtoms() {
                queriedMask.add(i + 1)
            }
        }
        let queryAtom = m_query.getAtoms()[0]
        for i in 0..<queried.numAtoms() {
            if (!queriedMask.contains(Int(i + 1))) {
                continue
            }
            var state = State(&functor, m_query, queried, queriedMask)
            guard let queriedAtom = queried.getAtom(i+1) else { continue }
            if (!queryAtom.matches(queriedAtom)) {
                continue
            }
            if (m_query.numAtoms() > 1) {
                if (matchCandidate(&state, queryAtom, queriedAtom)) {
                    mapNext(&state, queryAtom, queriedAtom)
                }
            } else {
                var map: MKIsomorphMapper.Mapping = []
                map.append((UInt(queryAtom.getIndex()), UInt(queriedAtom.getIndex())))
                functor.operate(map)
            }
        }

        if (Date().timeIntervalSince(m_startTime) > m_timeout) {
            print("ERROR: Time limi exceeded...")
//            TODO: Maybe throw catchable error here to gracely exit from
            return
        }
    }
}

class MKAutomorphismQueryAtom: MKQueryAtom {

    var symClass: UInt = 0
    var symClasses: [UInt] = []

    init(_ symClass: UInt, _ symClasses: [UInt]) {
        super.init()
        self.symClass = symClass
        self.symClasses = symClasses
    }

    override func matches(_ atom: MKAtom) -> Bool {
        return symClasses[atom.getIndex()] == symClass
    }
    
}
