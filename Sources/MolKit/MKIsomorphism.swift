import Foundation


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
    typealias Mapping = [(UInt, UInt)]
    
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
    
    private var m_query: MKQuery //!< The search query.
    private var m_timeout: UInt = 0  //!< The timeout in seconds
    
    init(_ query: MKQuery) {
        self.m_query = query
        self.m_timeout = 60
    }
    
    func setTimeout(_ seconds: UInt) {
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
    open func mapFirst(_ queried: MKMol, _ map: inout Mapping, _ mask: MKBitVec = MKBitVec()) {
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
    open func mapUnique(_ queried: MKMol, _ maps: inout Mappings, _ mask: MKBitVec = MKBitVec()) {
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
    open func MapAll(_ queried: MKMol, _ maps: inout Mappings, _ mask: MKBitVec = MKBitVec(), _ maxMemory: Int = 3000000) {
        fatalError("Not implemented in base class")
    }
    
    /**
     * @class Functor isomorphism.h <openbabel/isomorphism.h>
     * @brief Functor base class to be used in combination with MapGeneric.
     * @see @ref MapGeneric
     * @since 2.3
     */
    class Functor {
        /**
         * This function is called every time an isomorphism is discovered.
         * Returing true, will abort the mapping process. The map is passed
         * as non-const reference and may be modified (e.g. swap).
         *
         * @see @ref MapGeneric
         */
        init () {}
        
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
    open func mapGeneric(_ functor: MKIsomorphMapper.Functor, _ queried: MKMol, _ mask: MKBitVec = MKBitVec()) {
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
func findAutomorphisms(_ mol: MKMol, _ aut: inout MKIsomorphMapper.Mappings, _ symmetry_classes: [UInt], _ mask: MKBitVec = MKBitVec(), _ maxMemory: Int = 3000000) -> Bool {
    fatalError()
}
/**
* Find the automorphisms of a molecule by using an OBIsomorphismMapper. This
* function will first find the graph invariants (i.e. symmetry_classes) using
* the mask. This function is a wrapper for FindAutomorphisms with a functor to
* store all automorphisms.
*
* @since version 2.3
*/
func findAutomorphisms(_ mol: MKMol, _ aut: inout MKIsomorphMapper.Mappings, _ mask: MKBitVec = MKBitVec(), _ maxMemory: Int = 3000000) -> Bool {
    fatalError()

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
func findAutomorhphisms(_ functor: inout MKIsomorphMapper.Functor, _ mol: MKMol, _ symmetry_classes: [UInt], _ mask: MKBitVec = MKBitVec()) {
    fatalError()

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
        var queriedMask: MKBitVec                // the queriedMask
        var queryPath: [UInt]                    // the path in the query
        var queriedPath: [UInt]                  // the path in the queried molecule
        var mapping: [MKAtom?]
        // the terminal sets
        var queryPathBits: MKBitVec
        var queriedPathBits: MKBitVec
        var queryDepths: [UInt]
        var queriedDepths: [UInt]

        init(_ functor: inout MKIsomorphMapper.Functor, _ query: MKQuery, _ queried: MKMol, _ mask: MKBitVec) {
            self.abort = false
            self.functor = functor
            self.query = query
            self.queried = queried
            self.queriedMask = mask
            self.queryPath = []
            self.queriedPath = []
            self.mapping = [MKAtom?](repeating: nil, count: query.numAtoms())
            self.queryPathBits = MKBitVec()
            self.queriedPathBits = MKBitVec()
            self.queryDepths = [UInt](repeating: 0, count: query.numAtoms())
            self.queriedDepths = [UInt](repeating: 0, count: queried.numAtoms())
        }

    }

    var m_startTime: Darwin.time_t?

    override init(_ query: MKQuery) {
        super.init(query)
        m_startTime = Darwin.time(nil)
    }

    func isInTerminalSet(_ depths: [UInt], _ path: MKBitVec, _ i: UInt) -> Bool {
        if depths[Int(i)] == 0 {
            return false
        }
        if path.bitIsSet(Int(i)) {
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
    func matchCandidate(_ state: inout State, _ queryAtom: MKQueryAtom, _ queryBond: MKQueryBond) -> Bool {
    //     if (!queryAtom->Matches(queriedAtom))
    //       return false;

    //     // add the neighbors to the paths
    //     state.queryPath.push_back(queryAtom->GetIndex());
    //     state.queriedPath.push_back(queriedAtom->GetIndex());
    //     // update the terminal sets
    //     state.queryPathBits.SetBitOn(queryAtom->GetIndex());
    //     state.queriedPathBits.SetBitOn(queriedAtom->GetIndex());
    //     // update mapping
    //     state.mapping[queryAtom->GetIndex()] = queriedAtom;

    //     //
    //     // update queryDepths
    //     //
    //     if (!state.queryDepths[queryAtom->GetIndex()])
    //       state.queryDepths[queryAtom->GetIndex()] = state.queryPath.size();

    //     std::vector<OBQueryAtom*> queryNbrs = queryAtom->GetNbrs();
    //     for (unsigned int i = 0; i < queryNbrs.size(); ++i) {
    //       unsigned int index = queryNbrs[i]->GetIndex();
    //       if (!state.queryDepths[index])
    //         state.queryDepths[index] = state.queryPath.size();
    //     }

    //     //
    //     // update queriedDepths
    //     //
    //     if (!state.queriedDepths[queriedAtom->GetIndex()])
    //       state.queriedDepths[queriedAtom->GetIndex()] = state.queriedPath.size();

    //     FOR_NBORS_OF_ATOM (nbr, queriedAtom) {
    //       unsigned int index = nbr->GetIndex();
    //       // skip atoms not in the mask
    //       if (!state.queriedMask.BitIsSet(index + 1))
    //         continue;
    //       if (!state.queriedDepths[index])
    //         state.queriedDepths[index] = state.queriedPath.size();
    //     }

    //     // check if the bonds match
    //     if (!checkBonds(state, queryAtom)) {
    //       Backtrack(state);
    //       return false;
    //     }

    //     //
    //     // Feasibility rules for the VF2 algorithm:
    //     //
    //     //  size( T1(s) ) <= size( T2(s) )
    //     //
    //     //  size( N1 - M1(s) - T1(s) ) <= size( N2 - M2(s) - T2(s) )
    //     //

    //     // compute T1(s) size
    //     unsigned int numT1 = 0;
    //     for (unsigned int i = 0; i < state.query->NumAtoms(); ++i)
    //       if (isInTerminalSet(state.queryDepths, state.queryPathBits, i))
    //           numT1++;
    //     // compute T2(s) size
    //     unsigned int numT2 = 0;
    //     for (unsigned int i = 0; i < state.queried->NumAtoms(); ++i)
    //       if (isInTerminalSet(state.queriedDepths, state.queriedPathBits, i))
    //           numT2++;

    //     // T1(s) > T()
    //     if (numT1 > numT2) {
    //       Backtrack(state);
    //       return false;
    //     }
    //     //  N1 - M1(s) - T1(s) > N2 - M2(s) - T2(s)
    //     if ((state.query->NumAtoms() - state.queryPath.size() - numT1) > (state.queried->NumAtoms() - state.queriedPath.size() - numT2)) {
    //       Backtrack(state);
    //       return false;
    //     }

    //     // Check if there is a mapping found
    //     state.abort = checkForMap(state);

    //     return true
    //   }
    }

    func nextCandidate(_ state: inout State, _ lastCandidate: inout Candidate) -> Candidate {
        // std::size_t lastQueryAtom = lastCandidate.queryAtom ? lastCandidate.queryAtom->GetIndex() : 0;
        // std::size_t lastQueriedAtom = lastCandidate.queriedAtom ? lastCandidate.queriedAtom->GetIndex() + 1 : 0;

        // std::size_t querySize = state.query->NumAtoms();
        // std::size_t queriedSize = state.queried->NumAtoms();

        // std::size_t queryTerminalSize = state.queryDepths.size() - std::count(state.queryDepths.begin(), state.queryDepths.end(), 0);
        // std::size_t queriedTerminalSize = state.queriedDepths.size() - std::count(state.queriedDepths.begin(), state.queriedDepths.end(), 0);

        // std::size_t mappingSize = state.queryPath.size();

        // if (queryTerminalSize > mappingSize && queriedTerminalSize > mappingSize) {
        //     while (lastQueryAtom < querySize && (state.queryPathBits.BitIsSet(lastQueryAtom) || !state.queryDepths[lastQueryAtom])) {
        //     lastQueryAtom++;
        //     lastQueriedAtom = 0;
        //     }
        // } else {
        //     while(lastQueryAtom < querySize && state.queryPathBits.BitIsSet(lastQueryAtom)) {
        //     lastQueryAtom++;
        //     lastQueriedAtom = 0;
        //     }
        // }

        // if (queryTerminalSize > mappingSize && queriedTerminalSize > mappingSize) {
        //     while (lastQueriedAtom < queriedSize && (state.queriedPathBits.BitIsSet(lastQueriedAtom) || !state.queriedDepths[lastQueriedAtom]))
        //     lastQueriedAtom++;
        // } else {
        //     while(lastQueriedAtom < queriedSize && state.queriedPathBits[lastQueriedAtom])
        //     lastQueriedAtom++;
        // }

        // if (lastQueryAtom < querySize && lastQueriedAtom < queriedSize)
        //     return Candidate(state.query->GetAtoms()[lastQueryAtom], state.queried->GetAtom(lastQueriedAtom + 1));

        // return Candidate();
        
    }

    /**
       * The depth-first isomorphism algorithm.
       */
    //   void MapNext(State &state, OBQueryAtom *queryAtom, OBAtom *queriedAtom)
    //   {
    //     if (time(nullptr) - m_startTime > m_timeout)
    //       return;
    //     if (state.abort)
    //       return;

    //     Candidate candidate;
    //     while (!state.abort) {
    //       candidate = NextCandidate(state, candidate);

    //       if (!candidate.queryAtom)
    //         return;

    //       if (DEBUG)
    //         cout << yellow << "candidate: " << candidate.queryAtom->GetIndex() << " -> " << candidate.queriedAtom->GetIndex() << normal << endl;


    //       if (matchCandidate(state, candidate.queryAtom, candidate.queriedAtom)) {
    //         MapNext(state, candidate.queryAtom, candidate.queriedAtom);
    //         Backtrack(state);
    //       }
    //     }

    //   }

    //   void Backtrack(State &state)
    //   {
    //     if (DEBUG)
    //       cout << red << "backtrack... " << normal << state.queryPath.size()-1 << endl;
    //     // remove last atoms from the mapping
    //     if (state.queryPath.size()) {
    //       state.mapping[state.queryPath.back()] = nullptr;
    //       state.queryPathBits.SetBitOff(state.queryPath.back());
    //       state.queryPath.pop_back();
    //     }
    //     if (state.queriedPath.size()) {
    //       state.queriedPathBits.SetBitOff(state.queriedPath.back());
    //       state.queriedPath.pop_back();
    //     }
    //     // restore queryDepths and queriedDepths
    //     unsigned int depth = state.queryPath.size() + 1;
    //     std::replace(state.queryDepths.begin(), state.queryDepths.end(), depth, static_cast<unsigned int>(0));
    //     std::replace(state.queriedDepths.begin(), state.queriedDepths.end(), depth, static_cast<unsigned int>(0)); // O(n)  n = # vertices in the queried
    //   }

    //   /**
    //    * Get the first mappings of the query on the queried molecule.
    //    * @param queried The queried molecule.
    //    * @return The mapping.
    //    */
    //   void MapFirst(const OBMol *queried, Mapping &map, const OBBitVec &mask) override
    //   {
    //     class MapFirstFunctor : public Functor
    //     {
    //       private:
    //         Mapping &m_map;
    //       public:
    //         MapFirstFunctor(Mapping &map) : m_map(map)
    //         {
    //         }
    //         bool operator()(Mapping &map) override
    //         {
    //           m_map = map;
    //           // stop mapping
    //           return true;
    //         }
    //     };

    //     MapFirstFunctor functor(map);
    //     MapGeneric(functor, queried, mask);
    //   }

    //   /**
    //    * Get all unique mappings of the query on the queried molecule.
    //    * @param queried The queried molecule.
    //    * @return The unique mappings
    //    */
    //   void MapUnique(const OBMol *queried, Mappings &maps, const OBBitVec &mask) override
    //   {
    //     class MapUniqueFunctor : public OBIsomorphismMapper::Functor
    //     {
    //       private:
    //         OBIsomorphismMapper::Mappings &m_maps;
    //       public:
    //         MapUniqueFunctor(OBIsomorphismMapper::Mappings &maps) : m_maps(maps)
    //         {
    //         }
    //         bool operator()(OBIsomorphismMapper::Mapping &map) override
    //         {
    //           // get the values from the map
    //           std::vector<unsigned int> values;
    //           for (OBIsomorphismMapper::Mapping::const_iterator it = map.begin(); it != map.end(); ++it)
    //             values.push_back(it->second);
    //           std::sort(values.begin(), values.end());
    //            // print_vector("values ", values);

    //           bool isUnique = true;
    //           for (unsigned int k = 0; k < m_maps.size(); ++k) {
    //             std::vector<unsigned int> kValues;
    //             for (OBIsomorphismMapper::Mapping::iterator it = m_maps[k].begin(); it != m_maps[k].end(); ++it)
    //               kValues.push_back(it->second);
    //             std::sort(kValues.begin(), kValues.end());

    //           //  print_vector("kValues", kValues);
    //             if (values == kValues)
    //               isUnique = false;
    //           }

    //           if (isUnique)
    //             m_maps.push_back(map);

    //           // continue mapping
    //           return false;
    //         }
    //     };


    //     maps.clear();
    //     MapUniqueFunctor functor(maps);
    //     MapGeneric(functor, queried, mask);

    //     if (DEBUG)
    //       for (unsigned int i =0; i < maps.size(); ++i) {
    //         cout << "mapping:" << endl;
    //         for (Mapping::iterator it = maps[i].begin(); it != maps[i].end(); ++it)
    //           cout << "    " << it->first << " -> " << it->second << endl;
    //       }
    //   }

    //   /**
    //    * Get all mappings of the query on the queried molecule. Duplicates are
    //    * ignored but unlinke MapUnique, multiple mappings of the query on the same
    //    * part of the queried structure are allowed. This makes it possible to use
    //    * MapAll for finding the automorphism group.
    //    * @param queried The queried molecule.
    //    * @return The mappings.
    //    */
    //   void MapAll(const OBMol *queried, Mappings &maps, const OBBitVec &mask, std::size_t maxMemory) override
    //   {
    //     maps.clear();
    //     MapAllFunctor functor(maps, maxMemory);
    //     MapGeneric(functor, queried, mask);

    //     if (DEBUG)
    //       for (unsigned int i =0; i < maps.size(); ++i) {
    //         cout << "mapping:" << endl;
    //         for (Mapping::iterator it = maps[i].begin(); it != maps[i].end(); ++it)
    //           cout << "    " << it->first << " -> " << it->second << endl;
    //       }

    //   }

    //   void MapGeneric(Functor &functor, const OBMol *queried, const OBBitVec &mask) override
    //   {
    //     m_startTime = time(nullptr);
    //     if(m_query->NumAtoms() == 0) return;
    //     // set all atoms to 1 if the mask is empty
    //     OBBitVec queriedMask = mask;
    //     if (!queriedMask.CountBits())
    //       for (unsigned int i = 0; i < queried->NumAtoms(); ++i)
    //         queriedMask.SetBitOn(i + 1);

    //     OBQueryAtom *queryAtom = m_query->GetAtoms()[0];
    //     for (unsigned int i = 0; i < queried->NumAtoms(); ++i) {
    //       if (!queriedMask.BitIsSet(i + 1)) {
    //         continue;
    //       }
    //       State state(functor, m_query, queried, queriedMask);
    //       OBAtom *queriedAtom = queried->GetAtom(i+1);
    //       if (!queryAtom->Matches(queriedAtom)) {
    //         continue;
    //       }
    //       if (DEBUG)
    //         std::cout << blue << "START: 0 -> " << queriedAtom->GetIndex() << normal << std::endl;

    //       if (m_query->NumAtoms() > 1) {
    //         if (matchCandidate(state, queryAtom, queriedAtom))
    //           MapNext(state, queryAtom, queriedAtom);
    //       } else {
    //         Mapping map;
    //         map.push_back(std::make_pair(queryAtom->GetIndex(), queriedAtom->GetIndex()));
    //         functor(map);
    //       }
    //     }

    //     if (time(nullptr) - m_startTime > m_timeout)
    //       obErrorLog.ThrowError(__FUNCTION__, "time limit exceeded...", obError);

    //   }
}

class MKAutomorphismQueryAtom: MKQueryAtom {

    var symClass: UInt 
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
