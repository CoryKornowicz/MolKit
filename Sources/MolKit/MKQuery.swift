import Foundation
import Bitset
/**
* @class OBQueryAtom query.h <openbabel/query.h>
* @brief Atom in an OBQuery
*
* The OBQueryAtom class defines an interface for query atoms. The class provides
* some general methods and properties to access the topology information. The Matches
* method can be reimplemented in subclasses to get custom matching behavior.
*
* The default Matches implementation only checks the atomic number.
*
* See @ref substructure for more information.
*
* @sa OBQuery OBQueryBond OBIsomorphismMapper
* @since version 2.3
*/
class MKQueryAtom: Equatable {
    private var m_atomicNum: Int
    private var m_isInRing: Bool = false
    private var m_isAromatic: Bool = false
    var m_index: Int = 0
    var m_bonds: [MKQueryBond] = []
    var m_nbrs: [MKQueryAtom] = []

    /**
    * Constructor.
    * @param atomicNum The atomic number for this query atom.
    * @param isInRing Specify whether the query atom is in a ring. Default is false.
    * @param isAromatic Specify whether the query atom is aromatic. Default is false.
    */
    init(_ atomicNum: Int = 6, _ isInRing: Bool = false, _ isAromatic: Bool = false) {
        self.m_atomicNum = atomicNum
        self.m_isInRing = isInRing
        self.m_isAromatic = isAromatic
    }

    /**
    * Get the index for this query atom. Atoms are indexed starting from 0.
    * This method is used by OBIsomorphismMapper implementations.
    */
    func getIndex() -> Int {
        return m_index
    }

    /**
    * Get the query bonds for this atom.
    * This method is used by OBIsomorphismMapper implementations.
    */
    func getBonds() -> [MKQueryBond] {
        return m_bonds
    }

    /**
    * Get the neighbor query atoms.
    * This method is used by OBIsomorphismMapper implementations.
    */
    func getNbrs() -> [MKQueryAtom] {
        return m_nbrs
    }

    /**
    * This is the match method to verify if an OBQueryAtom and OBAtom class match.
    * The default implementation only checks if the atomic numbers match. Reimplement
    * this method in a subclass for more advances matching.
    * This method is used by OBIsomorphismMapper implementations.
    * @param atom The OBAtom object to compare this OBQueryAtom with.
    */
    open func matches(_ atom: MKAtom) -> Bool {
        if (atom.getAtomicNum() != m_atomicNum) { return false }
        if (atom.isAromatic() != m_isAromatic) { return false }
        if m_isInRing { 
          if !atom.isInRing() {
            return false
          }
        }
        return true
    }
    
    static func == (_ lhs: MKQueryAtom, _ rhs: MKQueryAtom) -> Bool {
        return lhs.m_index == rhs.m_index && lhs.m_atomicNum == rhs.m_atomicNum
    }
}

/**
* @class OBQueryBond query.h <openbabel/query.h>
* @brief Bond in an OBQuery
*
* The OBQueryBond class defines an interface for query bonds. The class provides
* some general methods and properties to access the topology information. The Matches
* method can be reimplemented in subclasses to get custom matching behavior.
*
* The default Matches implementation only checks if the bonds are both aromatic,
* otherwise the bond orders are compared.
*
* See @ref substructure for more information.
*
* @sa OBQuery OBQueryAtom OBIsomorphismMapper
* @since version 2.3
*/
class MKQueryBond {

    var m_index: Int = 0
    private var m_begin: MKQueryAtom
    private var m_end: MKQueryAtom
    private var m_order: Int = 1
    private var m_aromatic: Bool = false

    init(_ begin: MKQueryAtom, _ end: MKQueryAtom, _ order: Int = 1, _ aromatic: Bool = false) {
        self.m_begin = begin
        self.m_end = end
        self.m_order = order
        self.m_aromatic = aromatic

        m_begin.m_bonds.append(self)
        m_end.m_bonds.append(self)
        m_begin.m_nbrs.append(m_end)
        m_end.m_nbrs.append(m_begin)
    }


    /**
    * Get the index for this query bonds. Query bonds are indexed starting from 0.
    */
    func getIndex() -> Int {
        return m_index
    }

    /**
    * Get the begin atom.
    */
    func getBegin() -> MKQueryAtom {
        return m_begin
    }

    /**
    * Get the end atom.
    */
    func getEnd() -> MKQueryAtom {
        return m_end
    }

    /**
    * This is the match method to verify if an OBQueryBond and OBBond class match.
    * The default implementation checks if both bonds are aromatic and compares the
    * bond orders otherwise. Reimplement this method in a subclass for more
    * advances matching.
    * This method is used by OBIsomorphismMapper implementations.
    * @param bond The OBBond object to compare this OBQueryBond with.
    */
    open func matches(_ bond: MKBond) -> Bool {
        if m_aromatic {
            return bond.isAromatic()
        }
        return bond.getBondOrder() == m_order
    }

}


/**
* @class OBQuery query.h <openbabel/query.h>
* @brief A substructure query
*
* See @ref substructure for more information.
* @since version 2.3
*/
class MKQuery {
    private var m_atoms: [MKQueryAtom] = []
    private var m_bonds: [MKQueryBond] = []

    /**
    * Get the number of query atoms.
    */
    func numAtoms() -> Int {
        return m_atoms.count
    }

    /**
    * Get the number of query bonds.
    */
    func numBonds() -> Int {
        return m_bonds.count
    }

    /**
    * Get the query atoms.
    */
    func getAtoms() -> [MKQueryAtom] {
        return m_atoms
    }

    /**
    * Get the query bonds.
    */
    func getBonds() -> [MKQueryBond] {
        return m_bonds
    }

    /**
    * @return The query bond between @p begin and @p end. If there is no
    * bond between @p begin and @p end, this function returns 0.
    */
    func getBond(_ begin: MKQueryAtom, _ end: MKQueryAtom) -> MKQueryBond? {
        for i in 0..<begin.m_bonds.count {
            if begin.getNbrs()[i] == end {
                return begin.getBonds()[i]
            }
        }
        return nil
    }

    /**
     * Add a query atom to the query. This function steals the pointer.
     */
    func addAtom(_ atom: MKQueryAtom) {
        atom.m_index = m_atoms.count
        m_atoms.append(atom)
    }
    /**
     * Add a query atom to the query. This function steals the pointer.
     */
    func addBond(_ bond: MKQueryBond) {
        bond.m_index = m_bonds.count
        m_bonds.append(bond)
    }
    
}

/**
* Create an OBQuery object from an OBMol object. 
* @param mol The query molecule.
* @param mask The mask specifying the atoms to use. Indexed from 1 (i.e. OBAtom::GetIdx()).
* @return A pointer to an OBQuery object for the smiles string. This pointer should be deleted.
* @since version 2.3
*/
func compileMoleculeQuery(_ mol: MKMol, _ mask: Bitset = Bitset()) -> MKQuery {
    // set all atoms to 1 if the mask is empty
    var mask2 = mask
    if mask2.count() == 0 {
        for i in 0..<mol.numAtoms() {
            mask2.add(i + 1)
        }
    }
    var query = MKQuery()
    var offset = 0
    var indexes: [Int] = []
    for atom in mol.getAtomIterator() {
        indexes.append(atom.getIndex() - offset)
        if !mask2.contains(atom.getIndex() + 1) {
            offset += 1
            continue
        }
        query.addAtom(MKQueryAtom(atom.getAtomicNum(), atom.isInRing(), atom.isAromatic()))
    }
    for bond in mol.getBondIterator() {
        let beginIndex = bond.getBeginAtom().getIndex()
        let endIndex = bond.getEndAtom().getIndex()
        if !mask2.contains(beginIndex + 1) || !mask2.contains(endIndex + 1) {
            continue
        }
        query.addBond(MKQueryBond(query.getAtoms()[indexes[beginIndex]], query.getAtoms()[indexes[endIndex]], Int(bond.getBondOrder()), bond.isAromatic()))
    }
    return query
}

/**
* Create an OBQuery object from a smiles string. 
* @param smiles The query smiles string.
* @param mask The mask specifying the atoms to use. Indexed from 1 (i.e. OBAtom::GetIdx()).
* @return A pointer to an OBQuery object for the smiles string. This pointer should be deleted.
* @since version 2.3
*/
func compileSmilesQuery(_ smiles: String, _ mask: Bitset = Bitset()) -> MKQuery {
    let conv = MKConversion()
    conv.setInFormat("smi")
    var mol: MKMol = MKMol()
    conv.readString(&mol, smiles)
    return compileMoleculeQuery(mol, mask)
}
