//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 1/31/23.
//

import Foundation
import Bitset

//! Residue types.
public enum MKResidueProperty: Int {
    case AMINO = 0
    case AMINO_NUCLEO = 1
    case COENZYME = 2
    case ION = 3
    case NUCLEO = 4
    case PROTEIN = 5
    case PURINE = 6
    case PYRIMIDINE = 7
    case SOLVENT = 8
    case WATER = 9
}

//! Residue names (index into Residue[] array)
public enum MKResidueIndex: Int, Comparable, Equatable {
    case ALA   =  0
    case GLY   =  1
    case LEU   =  2
    case SER   =  3
    case VAL   =  4
    case THR   =  5
    case LYS   =  6
    case ASP   =  7
    case ILE   =  8
    case ASN   =  9
    case GLU   = 10
    case PRO   = 11
    case ARG   = 12
    case PHE   = 13
    case GLN   = 14
    case TYR   = 15
    case HIS   = 16
    case CYS   = 17
    case MET   = 18
    case TRP   = 19
    case ASX   = 20
    case GLX   = 21
    case PCA   = 22
    case HYP   = 23
    case A     = 24
    case C     = 25
    case G     = 26
    case T     = 27
    case U     = 28
    case UPLUS = 29
    case I     = 30
    case _1MA  = 32
    case _5MC  = 54 // MARK: This will be an issue when converting from raw value to enum case <-> vice versa | 32
    case OMC   = 33
    case _1MG  = 34
    case _2MG  = 35
    case M2G   = 36
    case _7MG  = 37
    case OMG   = 38
    case YG    = 39
    case H2U   = 40
    case _5MU  = 41
    case PSU   = 42
    case UNK   = 43
    case ACE   = 44
    case FOR   = 45
    case HOH   = 46
    case DOD   = 47
    case SO4   = 48
    case PO4   = 49
    case NAD   = 50
    case COA   = 51
    case NAP   = 52
    case NDP   = 53
    
    public static func < (lhs: MKResidueIndex, rhs: MKResidueIndex) -> Bool {
        return lhs.rawValue < rhs.rawValue
    }
    
    public static func == (lhs: MKResidueIndex, rhs: MKResidueIndex) -> Bool {
        return lhs.rawValue == rhs.rawValue
    }
}

//! Residue atom properties
public enum MKResidueAtomProperty: Int  {
    case ALPHA_CARBON     = 0
    case AMINO_BACKBONE   = 1
    case BACKBONE         = 2
    case CYSTEINE_SULPHUR = 3
    case LIGAND           = 4
    case NUCLEIC_BACKBONE = 5
    case SHAPELY_BACKBONE = 6
    case SHAPELY_SPECIAL  = 7
    case SIDECHAIN        = 8
    case SUGAR_PHOSPHATE  = 9
}

//! Residue property definitions
public enum MKAminoAcidProperty: Int {

    case ACIDIC      =  0
    case ACYCLIC     =  1
    case ALIPHATIC   =  2
    case AROMATIC    =  3
    case BASIC       =  4
    case BURIED      =  5
    case CHARGED     =  6
    case CYCLIC      =  7
    case HYDROPHOBIC =  8
    case LARGE       =  9
    case MEDIUM      = 10
    case NEGATIVE    = 11
    case NEUTRAL     = 12
    case POLAR       = 13
    case POSITIVE    = 14
    case SMALL       = 15
    case SURFACE     = 16
}

public let MAXSETNO = 40
public let MAXELEM = 29
public let MAXRES = 55 // MARK: Potential Error since MKResiudeIndex is +1 more since the values all need to be unique

public let AA_ALA = (1<<1)
public let AA_GLY = (1<<2)
public let AA_LEU = (1<<3)
public let AA_SER = (1<<4)
public let AA_VAL = (1<<5)
public let AA_THR = (1<<6)
public let AA_LYS = (1<<7)
public let AA_ASP = (1<<8)
public let AA_ILE = (1<<9)
public let AA_ASN = (1<<10)
public let AA_GLU = (1<<11)
public let AA_PRO = (1<<12)
public let AA_ARG = (1<<13)
public let AA_PHE = (1<<14)
public let AA_GLN = (1<<15)
public let AA_TYR = (1<<16)
public let AA_HIS = (1<<17)
public let AA_CYS = (1<<18)
public let AA_MET = (1<<19)
public let AA_TRP = (1<<20)

/////////////////////////////////////////////////////////////////////////////
// Amino Acid Property Definitions
/////////////////////////////////////////////////////////////////////////////

public func is_acidic(_ x: Int) -> Bool {
    return (x & (AA_ASP | AA_GLU)) != 0
}

public func is_acyclic(_ x: Int) -> Bool {
    return (x & (AA_ALA | AA_GLY | AA_LEU | 
                 AA_SER | AA_VAL | AA_THR | 
                 AA_LYS | AA_ASP | AA_ILE | 
                 AA_ASN | AA_GLU | AA_GLN | 
                 AA_CYS | AA_MET)) != 0
}

public func is_aliphatic(_ x: Int) -> Bool {
    return (x & (AA_ALA | AA_GLY | AA_ILE | 
                 AA_LEU | AA_VAL)) != 0
}

public func is_aromatic(_ x: Int) -> Bool {
    return (x & (AA_HIS | AA_PHE | AA_TRP | AA_TYR)) != 0
}

public func is_basic(_ x: Int) -> Bool {
    return (x & (AA_ARG | AA_HIS | AA_LYS)) != 0
}

public func is_buried(_ x: Int) -> Bool {
    return (x & (AA_ALA | AA_CYS | AA_ILE | 
                 AA_LEU | AA_MET | AA_PHE | 
                 AA_TRP | AA_VAL)) != 0
}

public func is_charged(_ x: Int) -> Bool {
    return (x & (AA_ASP | AA_GLU | AA_ARG | 
                 AA_HIS | AA_LYS)) != 0
}

public func is_cyclic(_ x: Int) -> Bool {
    return (x & (AA_HIS | AA_PHE | AA_PRO | 
                 AA_TRP | AA_TYR)) != 0
}

public func is_hydrophobic(_ x: Int) -> Bool {
    return (x & (AA_ALA | AA_LEU | AA_VAL | 
                 AA_ILE | AA_PRO | AA_PHE | 
                 AA_MET | AA_TRP)) != 0
}

public func is_large(_ x: Int) -> Bool {
    return (x & (AA_ARG | AA_PHE | AA_GLN | 
                 AA_TYR | AA_HIS | AA_LEU | 
                 AA_LYS | AA_ILE | AA_GLU | 
                 AA_MET | AA_TRP)) != 0
}

public func is_medium(_ x: Int) -> Bool {
    return (x & (AA_VAL | AA_THR | AA_ASP | 
                 AA_ASN | AA_PRO | AA_CYS)) != 0
}

public func is_negative(_ x: Int) -> Bool {
    return (x & (AA_ASP | AA_GLU)) != 0
}

public func is_neutral(_ x: Int) -> Bool {
    return (x & (AA_ALA | AA_GLY | AA_LEU | 
                 AA_SER | AA_VAL | AA_THR | 
                 AA_PHE | AA_GLN | AA_TYR | 
                 AA_HIS | AA_CYS | AA_MET | 
                 AA_TRP | AA_ILE | AA_ASN | 
                 AA_PRO)) != 0
}

public func is_polar(_ x: Int) -> Bool {
    return (x & (AA_ASP | AA_ILE | AA_ASN | 
                 AA_GLU | AA_SER | AA_THR | 
                 AA_ARG | AA_GLN | AA_CYS | 
                 AA_HIS)) != 0
}

public func is_positive(_ x: Int) -> Bool {
    return (x & (AA_ARG | AA_HIS | AA_LYS)) != 0
}

public func is_small(_ x: Int) -> Bool {
    return (x & (AA_ALA | AA_GLY | AA_SER)) != 0
}

public func is_surface(_ x: Int) -> Bool {
    return (x & (AA_THR | AA_LYS | AA_ASP | 
                 AA_ILE | AA_ASN | AA_GLU | 
                 AA_PRO | AA_ARG | AA_GLY | 
                 AA_SER | AA_GLN | AA_TYR | 
                 AA_HIS)) != 0
}

public class MKResidue: MKBase {

    private let Residue =  [
    /*===============*/
    /*  Amino Acids  */
    /*===============*/

    /* Ordered by Cumulative Frequency in Brookhaven *
     * Protein Databank, December 1991               */

    "ALA", /* 8.4% */     "GLY", /* 8.3% */
    "LEU", /* 8.0% */     "SER", /* 7.5% */
    "VAL", /* 7.1% */     "THR", /* 6.4% */
    "LYS", /* 5.8% */     "ASP", /* 5.5% */
    "ILE", /* 5.2% */     "ASN", /* 4.9% */
    "GLU", /* 4.9% */     "PRO", /* 4.4% */
    "ARG", /* 3.8% */     "PHE", /* 3.7% */
    "GLN", /* 3.5% */     "TYR", /* 3.5% */
    "HIS", /* 2.3% */     "CYS", /* 2.0% */
    "MET", /* 1.8% */     "TRP", /* 1.4% */

    "ASX", "GLX", "PCA", "HYP",

    /*===================*/
    /*  DNA Nucleotides  */
    /*===================*/
    "  A", "  C", "  G", "  T",

    /*===================*/
    /*  RNA Nucleotides  */
    /*===================*/
    "  U", " +U", "  I", "1MA",
    "5MC", "OMC", "1MG", "2MG",
    "M2G", "7MG", "OMG", " YG",
    "H2U", "5MU", "PSU",

    /*=================*/
    /*  Miscellaneous  */
    /*=================*/
    "UNK", "ACE", "FOR", "HOH",
    "DOD", "SO4", "PO4", "NAD",
    "COA", "NAP", "NDP" 
    ]

    // MARK: Need to Type this as OrderedDictionary
    private let ElemDesc = [
    [ " ", "N", " ", " " ],  /* 0*/
    [ " ", "C", "A", " " ],  /* 1*/
    [ " ", "C", " ", " " ],  /* 2*/
    [ " ", "O", " ", " " ],  /* 3*/   /* 0-3   Amino Acid Backbone    */
    [ " ", "C", "\\", " " ], /* 4*/
    [ " ", "O", "T", " " ],  /* 5*/
    [ " ", "S", " ", " " ],  /* 6*/
    [ " ", "P", " ", " " ],  /* 7*/   /* 4-7   Shapely Amino Backbone */
    [ " ", "O", "1", "P" ],  /* 8*/
    [ " ", "O", "2", "P" ],  /* 9*/
    [ " ", "O", "5", "*" ],  /*10*/
    [ " ", "C", "5", "*" ],  /*11*/
    [ " ", "C", "4", "*" ],  /*12*/
    [ " ", "O", "4", "*" ],  /*13*/
    [ " ", "C", "3", "*" ],  /*14*/
    [ " ", "O", "3", "*" ],  /*15*/
    [ " ", "C", "2", "*" ],  /*16*/
    [ " ", "O", "2", "*" ],  /*17*/
    [ " ", "C", "1", "*" ],  /*18*/   /* 7-18  Nucleic Acid Backbone  */
    [ " ", "C", "A", "2" ],  /*19*/   /* 19    Shapely Special        */
    [ " ", "S", "G", " " ],  /*20*/   /* 20    Cysteine Sulphur       */
    [ " ", "N", "1", " " ],  /*21*/
    [ " ", "N", "2", " " ],  /*22*/
    [ " ", "N", "3", " " ],  /*23*/
    [ " ", "N", "4", " " ],  /*24*/
    [ " ", "N", "6", " " ],  /*25*/
    [ " ", "O", "2", " " ],  /*26*/
    [ " ", "O", "4", " " ],  /*27*/
    [ " ", "O", "6", " " ]   /*28*/   /* 21-28 Nucleic Acid H-Bonding */
    ]
    
    var _idx: Int = 0
    var _aakey: Int = 0
    var _chain: String = ""
    var _reskey: MKResidueIndex = MKResidueIndex.UNK
    var _resnum: String = ""
    var _resname: String = ""
    var _insertioncode: String = ""

    var _hetatm: [Bool] = []
    var _atomid: [String] = []
    var _atoms: [MKAtom] = []
    var _sernum: [Int] = []

    public override init() {
        super.init()
    }

    func addAtom(_ atom: MKAtom) { 
        atom.setResidue(self)
        self._atoms.append(atom)
        self._atomid.append("")
        self._hetatm.append(false)
        self._sernum.append(0)
    }

    func insertAtom(_ atom: MKAtom) { 
        self.addAtom(atom)
    }

    func removeAtom(_ atom: MKAtom) {
        
        if let idx = self._atoms.firstIndex(of: atom) {
            self._atoms[idx].setResidue(nil)
            self._atoms.remove(at: idx)
            self._atomid.remove(at: idx)
            self._hetatm.remove(at: idx)
            self._sernum.remove(at: idx)
        }
    }

    deinit {
        for atom: MKAtom in self._atoms {
            atom.setResidue(nil)
        }
        self._atoms = []
    }

    func copyData(_ src: MKResidue) {
        if src != self {
            _chain = src._chain
            _aakey = src._aakey
            _reskey = src._reskey
            _resnum = src._resnum
            _resname = src._resname
            _atomid = src._atomid
            _hetatm = src._hetatm
            _sernum = src._sernum
            _insertioncode = src._insertioncode
        }
    }

    override func clear() {
        for atom in self._atoms {
            atom.setResidue(nil)
        }
        self._atoms = []
        self._atomid = []
        self._hetatm = []
        self._sernum = []

        self._chain = ""
        self._reskey = MKResidueIndex.UNK
        self._resnum = ""
        self._resname = ""
        self._idx = 0
        self._aakey = 0
        self._insertioncode = ""
        super.clear()
    }

    //! \brief Set the name of this residue (e.g., "ALA"). Use 3-char PDB standard names.
    //! http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/part_79.html
    //! MODRES records for modified residues:
    //! http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/part_36.html
    func setName(_ resName: String) {
        self._resname = resName 
        MKResidue.setResidueKeys(self._resname, &self._reskey, &self._aakey)
     }

    //! Set the residue number (in the sequence)
    func setNum(_ resNum: Int) { 
        self._resnum = String(resNum)
    }
    func setNum(_ resNum: String) {
        self._resnum = resNum
    }

    //! Set the chain ID for this residue
    func setChain(_ chain: String) {
        self._chain = chain
    }

    //! Set the chain number for this residue
    func setChainNum(_ chainNum: Int) {
        self._chain = String(chainNum)
    }

    //! Set the internal index of this residue in the parent OBMol.
    //! Intended mostly for internal use
    func setIdx(_ idx: Int) { }

    //! Set  PDB insertion code information for this residue. This allows
    //! consecutive residues to have the same number. Some communities
    //! that work in a well-conserved structural world use this, e.g.
    //! for immunoglobulins.
    func setInsertionCode(_ insCode: String) { 
        self._insertioncode = insCode
    }

    //! Set the character code ID for an ATOM record for the supplied atom
    //! This does nothing if the supplied atom is not found in the residue
    func setAtomID(_ atom: MKAtom, _ id: String) {
        if let idx = self._atoms.firstIndex(of: atom) {
            self._atomid[idx] = id
        }
    }

    func setHetAtom(_ atom: MKAtom, _ hetatm: Bool) {
        if let idx = self._atoms.firstIndex(of: atom) {
            self._hetatm[idx] = hetatm
        }
    }

    //! Set the atomic serial number for a given atom (see OBSerialNums)
    func setSerialNum(_ atom: MKAtom, _ sernum: Int) {
        if let idx = self._atoms.firstIndex(of: atom) {
            self._sernum[idx] = sernum
        }
    }

    //! \return The residue name
    func getName() -> String {
        return self._resname
    }

    //! \return The residue number (in the sequence)
    func getNum() -> Int {
        do {
            let num = try Int(self._resnum, format: .number)
            return num
        } catch {
            return -1
        }
    }

    func getNumString() -> String {
        return self._resnum
    }

    //! \return The number of atoms in this residue
    func getNumAtoms() -> Int {
        return self._atoms.count
    }

    //! \return The number of heavy atoms in this residue
    func getNumHvyAtoms() -> Int {
        return self._atoms.map({ $0.getAtomicNum() != 1 ? 1 : 0 }).reduce(0, +)
    }

    //! \return The ID of the chain which includes this residue
    func getChain() -> String {
        return self._chain
    }

    //! \return The number of the chain which includes this residue
    func getChainNum() -> Int {
        if self._chain.isNumber {
            return Int(self._chain)!
        } else {
            do {
                let num = try Int(self._chain, format: .number) + 1
                return num
            } catch {
                return -1
            }
        }
    }

    //! \return The internal index of this residue in the parent OBMol
    func getIdx() -> Int {
        return self._idx
    }

    //! \return The residue key (i.e., an entry in the MKResidueIndex namespace)
    func getResKey() -> MKResidueIndex {
        return self._reskey
    }

    //! \return a vector of all atoms in this residue
    func getAtoms() -> [MKAtom] {
        return self._atoms
    }

    //! \return all bonds in this residue. @p exterior includes bonds to atoms
    //!  outside this residue (default is true)
    func getBonds(_ exterior: Bool = true) -> [MKBond] {
        let idxs: Bitset = Bitset()
        let bonds: [MKBond] = [MKBond]()
        for atom in self._atoms {
            guard let bonds = atom.getBondIterator() else { continue }
            for bond in bonds {
                if !idxs.contains(Int(bond.getIdx())) {
                    if !exterior {
                        if bond.getNbrAtom(atom).getResidue() == self {
                            bonds.append(bond)
                        }
                    } else {
                        bonds.append(bond)
                    }
                    idxs.add(Int(bond.getIdx()))
                }
            }
        }
        return bonds
    }

    //! \return the atom ID (character code) for the supplied atom or ""
    //!  if the atom is not found in this residue
    func getAtomID(_ atom: MKAtom) -> String {
        for i in 0..<self._atoms.count {
            if self._atoms[i] === atom {
                return self._atomid[i]
            }
        }
        return ""
    }

    //! \return the serial number of the supplied atom (uses OBSerialNums)
    func getSerialNum(_ atom: MKAtom) -> Int {
        for i in 0..<self._atoms.count {
            if self._atoms[i] === atom {
                return self._sernum[i]
            }
        }
        return 0
    }

    //! \return The Insertion Code (i.e., an extra position motivated by a
    //! multiple sequence alignment against a template with defined numbers)
    func getInsertionCode() -> String {
        return self._insertioncode
    }

    func getAtomIterator() -> MKIterator<MKAtom> {
        return MKIterator<MKAtom>(_atoms)
    }
    
///////////////////////////////////////////////////////////////////////////////
// MKResidue: Information Functions
///////////////////////////////////////////////////////////////////////////////

    //! \return Whether this residue has the supplied amino acid property
    //!  defined from the OBAminoAcidProperty namespace
    func getAminoAcidProperty(_ prop: Int) -> Bool {
        switch MKAminoAcidProperty(rawValue: prop) {
            case .ACIDIC:
                return is_acidic(_aakey)
            case .ACYCLIC:
                return is_acyclic(_aakey)
            case .ALIPHATIC:
                return is_aliphatic(_aakey)
            case .AROMATIC:
                return is_aromatic(_aakey)
            case .BASIC:
                return is_basic(_aakey)
            case .BURIED:
                return is_buried(_aakey)
            case .CHARGED:
                return is_charged(_aakey)
            case .CYCLIC:
                return is_cyclic(_aakey)
            case .HYDROPHOBIC:
                return is_hydrophobic(_aakey)
            case .LARGE:
                return is_large(_aakey)
            case .MEDIUM:
                return is_medium(_aakey)
            case .NEGATIVE:
                return is_negative(_aakey)
            case .NEUTRAL:
                return is_neutral(_aakey)
            case .POLAR:
                return is_polar(_aakey)
            case .POSITIVE:
                return is_positive(_aakey)
            case .SMALL:
                return is_small(_aakey)
            case .SURFACE:
                return is_surface(_aakey)
            default:
                return false
        }
    }

    //! \return Whether atom @p a has the supplied residue atom property
    //!  defined from the MKResidueAtomProperty namespace
    func getAtomProperty(_ a: MKAtom, _ prop: MKResidueAtomProperty) -> Bool {
        let atomID = MKResidue.getAtomIDNumber(String(self.getAtomID(a)))
        switch prop {
            case .ALPHA_CARBON: 
                return atomID == 1
            case .AMINO_BACKBONE: 
                return atomID <= 3
            case .BACKBONE:
                return atomID <= 18
            case .CYSTEINE_SULPHUR:
                return atomID == 20
            case .LIGAND:
                return a.isHetAtom() && !self.getResidueProperty(MKResidueProperty.SOLVENT)
            case .NUCLEIC_BACKBONE:
                return (atomID >= 7) && (atomID <= 18)
            case .SHAPELY_BACKBONE:
                return atomID <= 7
            case .SHAPELY_SPECIAL:
                return atomID == 19
            case .SIDECHAIN:
                return self.getResidueProperty(MKResidueProperty.AMINO_NUCLEO) && (atomID > 18)
            case .SUGAR_PHOSPHATE:
                return atomID == 7
        }
    }

    //! \return Whether this residue has the supplied property
    //!  defined from the OBResidueProperty namespace
    func getResidueProperty(_ prop: MKResidueProperty) -> Bool {

        switch prop {
            case .AMINO:
            return self._reskey <= MKResidueIndex.HYP
            case .AMINO_NUCLEO:
                return self._reskey <= MKResidueIndex.PSU
            case .COENZYME:
                return self._reskey >= MKResidueIndex.NAD && self._reskey <= MKResidueIndex.NDP
            case .ION:
                return self._reskey == MKResidueIndex.SO4 || self._reskey == MKResidueIndex.PO4
            case .NUCLEO:
                return self._reskey >= MKResidueIndex.A && self._reskey <= MKResidueIndex.PSU
            case .PROTEIN:
                return self._reskey <= MKResidueIndex.HYP || (self._reskey >= MKResidueIndex.UNK && self._reskey <= MKResidueIndex.FOR)
            case .PURINE:
                return self._reskey == MKResidueIndex.A || self._reskey == MKResidueIndex.G
            case .PYRIMIDINE:
                return self._reskey == MKResidueIndex.C || self._reskey == MKResidueIndex.T
            case .SOLVENT:
                return self._reskey >= MKResidueIndex.HOH && self._reskey <= MKResidueIndex.PO4
            case .WATER:
                return self._reskey == MKResidueIndex.HOH || self._reskey == MKResidueIndex.DOD
        }
    }

    public func isHetAtom(_ atom: MKAtom) -> Bool {
        for i in 0..<self._atoms.count {
            if self._atoms[i] == atom {
                return self._hetatm[i]
            }
        }
        return false
    }

    public func isResidueType(_ type: MKResidueIndex) -> Bool { 
        return self._reskey == type
    }


    private static func getAtomIDNumber(_ atomID: String) -> Int {
        // TODO: Replace with real warning, if any 
        // Also, this needs to be shrunk from all of these absolute comparisons
        if atomID.length < 4 {
            print("atomID string should contain at least 4 characters")
            return MAXELEM
        }

        let ch1 = atomID[0].uppercased()
        let ch2 = atomID[1].uppercased()
        let ch3 = atomID[2].uppercased()
        let ch4 = atomID[3].uppercased()

        if ch1 == " " {
            switch ch2 {
                case "C":
                    switch ch3 {
                        case "A":
                            if ch4 == " " {
                                return 1
                            } else if ch4 == "2" {
                                return 19
                            }
                            break
                        case " ":
                            if ch4 == " " {
                                return 2
                            }
                            break
                        case "'":
                            if ch4 == " " {
                                return 4
                            }
                            break
                        case "1":
                            if ch4 == "*" {
                                return 18
                            }
                            break
                        case "2":
                            if ch4 == "*" {
                                return 16
                            }
                            break
                        case "3":
                            if ch4 == "*" {
                                return 14
                            }
                            break
                        case "4":
                            if ch4 == "*" {
                                return 12
                            }
                            break
                        case "5":
                            if ch4 == "*" {
                                return 11
                            }
                            break
                        default:
                            break
                    }
                    break
                case "N":
                    if ch4 == " " {
                        switch ch3 {
                            case " ":
                                return 0
                            case "1":
                                return 21
                            case "2":
                                return 22
                            case "3":
                                return 23
                            case "4":
                                return 24
                            case "6":
                                return 25
                            default:
                                break
                        }
                    }
                    break
                case "O":
                    switch ch3 {
                        case " ":
                            if ch4 == " " {
                                return 3
                            }
                            break
                        case "T":
                            if ch4 == " " {
                                return 5
                            }
                            break
                        case "1":
                            if ch4 == "P" {
                                return 8
                            }
                            break
                        case "2":
                            if ch4 == "P" {
                                return 9
                            } else if ch4 == "*" {
                                return 17
                            } else if ch4 == " " {
                                return 26
                            }
                            break
                        case "3":
                            if ch4 == "*" {
                                return 15
                            }
                            break
                        case "4":
                            if ch4 == "*" {
                                return 13
                            } else if ch4 == " " {
                                return 27
                            }
                            break
                        case "5":
                            if ch4 == "*" {
                                return 10
                            }
                            break
                        case "6":
                            if ch4 == " " {
                                return 28
                            }
                            break
                        default:
                            break
                    }
                case "P":
                    if ch3 == " " && ch4 == " " {
                        return 7
                    }
                    break
                case "S":
                    if ch4 == " " {
                        if ch3 == " " {
                            return 6
                        } else if ch3 == "G" {
                            return 20
                        }
                    }
                    break
                default:
                    break
            }
        }
        return MAXELEM
    }

    private static func getResidueNumber(_ resname: String) -> Int {
        if resname.length < 3 {
            print("resname string should contain at least 3 characters")
            return MAXRES
        }

        let ch1 = resname[0].uppercased()
        let ch2 = resname[1].uppercased()
        let ch3 = resname[2].uppercased()

        switch (ch1) {
            case " ":
                if ch2 == " " {
                    switch ch3 {
                        case "A": return 24
                        case "C": return 25
                        case "G": return 26
                        case "T": return 27
                        case "U": return 28
                        case "I": return 30
                        default: break 
                    }
                } else if ch2 == "+" {
                    if ch3 == "U" { return 29 }
                } else if ch2 == "Y" {
                    if ch3 == "G" { return 39 }
                }
            case "0": 
                if ch2 == "M" {
                    if ch3 == "C" { return 33 }
                    else if ch3 == "G" { return 38 }
                }
            case "1":
                if ch2 == "M" {
                    if ch3 == "A" { return 31 }
                    else if ch3 == "G" { return 34 }
                }
            case "2":
                if ch2 == "M" {
                    if ch3 == "G" { return 35 }
                }
            case "5":
                if ch2 == "M" {
                    if ch3 == "C" { return 32 }
                    else if ch3 == "U" { return 41 }
                }
            case "7":
                if ch2 == "M" {
                    if ch3 == "G" { return 37 }
                }
            case "A":
                if ch2 == "L" {
                    if ch3 == "A" { return 0 }
                } else if ch2 == "S" {
                    if ch3 == "P" { return 7 }
                    else if ch3 == "N" { return 9 }
                    else if ch3 == "X" { return 20 }
                } else if ch2 == "R" {
                    if ch3 == "G" { return 12 }
                } else if ch2 == "C" {
                    if ch3 == "E" { return 44 }
                } else if ch2 == "D" {
                    if ch3 == "E" { return 24 } /* "ADE" -> "  A" */
                }
            case "C":
                if ch2 == "Y" {
                    if ch3 == "S" { return 17 }
                    else if ch3 == "H" { return 17 } /* "CYH" -> "CYS" */
                    else if ch3 == "T" { return 25 } /* "CYT" -> "  C" */
                } else if ch2 == "O" {
                    if ch3 == "A" { return 51 }
                } else if ch2 == "P" {
                    if ch3 == "R" { return 11 } /* "CPR" -> "PRO" */
                } else if ch2 == "S" {
                    if ch3 == "H" { return 17 } /* "CSH" -> "CYS" */
                    else if ch3 == "M" { return 17 } /* "CSM" -> "CYS" */
                }
            case "D":
                if ch2 == "O" {
                    if ch3 == "D" { return 47 }
                } else if ch2 == "2" {
                    if ch3 == "O" { return 47 } /* "D2O" -> "DOD" */
                }
            case "F":
                if ch2 == "O" {
                    if ch3 == "R" { return 45 }
                }
            case "G":
                if ch2 == "L" {
                    if ch3 == "Y" { return 1 }
                    else if ch3 == "U" { return 10 }
                    else if ch3 == "N" { return 14 }
                    else if ch3 == "X" { return 21 }
                } else if ch2 == "U" {
                    if ch3 == "A" { return 26 } /* "GUA" -> "  G" */
                }
            case "H":
                if ch2 == "I" {
                    if ch3 == "S" { return 16 }
                } else if ch2 == "O" {
                    if ch3 == "H" { return 46 }
                } else if ch2 == "Y" {
                    if ch3 == "P" { return 23 }
                } else if ch2 == "2" {
                    if ch3 == "O" { return 46 } /* "H20" -> "HOH" */
                    else if ch3 == "U" { return 40 }
                }
            case "I":
                if ch2 == "L" {
                    if ch3 == "E" { return 8 }
                }
            case "L":
                if ch2 == "E" {
                    if ch3 == "U" { return 2 }
                } else if ch2 == "Y" {
                    if ch3 == "S" { return 6 }
                }
            case "M":
                if ch2 == "E" {
                    if ch3 == "T" { return 18 }
                } else if ch2 == "2" {
                    if ch3 == "G" { return 36 }
                }
            case "N":
                if ch2 == "A" {
                    if ch3 == "D" { return 50 }
                    else if ch3 == "P" { return 52 }
                } else if ch2 == "D" {
                    if ch3 == "P" { return 53 }
                }
            case "P":
                if ch2 == "R" {
                    if ch3 == "O" { return 11 }
                } else if ch2 == "H" {
                    if ch3 == "E" { return 13 }
                } else if ch2 == "C" {
                    if ch3 == "A" { return 22 }
                } else if ch2 == "O" {
                    if ch3 == "4" { return 49 }
                } else if ch2 == "S" {
                    if ch3 == "U" { return 42 }
                }
            case "S":
                if ch2 == "E" {
                    if ch3 == "R" { return 3 }
                } else if ch2 == "O" {
                    if ch3 == "4" { return 48 }
                    else if ch3 == "L" { return 46 } /* "SOL" -> "HOH" */
                } else if ch2 == "U" {
                    if ch3 == "L" { return 48 }      /* "SUL" -> "SO4" */
                }
            case "T":
                if ch2 == "H" {
                    if ch3 == "R" { return 5 }
                    else if ch3 == "Y" { return 27 } /* "THY" -> "  T" */
                } else if ch2 == "Y" {
                    if ch3 == "R" { return 15 }
                } else if ch2 == "R" {
                    if ch3 == "P" { return 19 }
                    else if ch3 == "Y" { return 19 } /* "TRY" -> "TRP" */
                } else if ch2 == "I" {
                    if ch3 == "P" { return 46 }      /* "TIP" -> "HOH" */
                }
            case "U":
                if ch2 == "N" {
                    if ch3 == "K" { return 43 }
                } else if ch2 == "R" {
                    if ch3 == "A" { return 28 } /* "URA" -> "  U" */
                    else if ch3 == "I" { return 28 } /* "URI" -> "  U" */
                }
            case "V":
                if ch2 == "A" {
                    if ch3 == "L" { return 4 }
                }
            case "W":
                if ch2 == "A" {
                    if ch3 == "T" { return 46 } /* "WAT" -> "HOH" */
                }
            default:
                break
        }
        return MKResidueIndex.UNK.rawValue
    }

    private static func setResidueKeys(_ residue: String, _ reskey: inout MKResidueIndex, _ aakey: inout Int) {
        reskey = MKResidueIndex(rawValue: MKResidue.getResidueNumber(residue))!
        switch reskey {
            case MKResidueIndex.ALA: aakey = AA_ALA
            case MKResidueIndex.ARG: aakey = AA_ARG
            case MKResidueIndex.ASN: aakey = AA_ASN
            case MKResidueIndex.ASP: aakey = AA_ASP
            case MKResidueIndex.CYS: aakey = AA_CYS
            case MKResidueIndex.GLN: aakey = AA_GLN
            case MKResidueIndex.GLU: aakey = AA_GLU
            case MKResidueIndex.GLY: aakey = AA_GLY
            case MKResidueIndex.HIS: aakey = AA_HIS
            case MKResidueIndex.ILE: aakey = AA_ILE
            case MKResidueIndex.LEU: aakey = AA_LEU
            case MKResidueIndex.LYS: aakey = AA_LYS
            case MKResidueIndex.MET: aakey = AA_MET
            case MKResidueIndex.PHE: aakey = AA_PHE
            case MKResidueIndex.PRO: aakey = AA_PRO
            case MKResidueIndex.SER: aakey = AA_SER
            case MKResidueIndex.THR: aakey = AA_THR
            case MKResidueIndex.TRP: aakey = AA_TRP
            case MKResidueIndex.TYR: aakey = AA_TYR
            case MKResidueIndex.VAL: aakey = AA_VAL
            default: aakey = 0
        }
    }
}
