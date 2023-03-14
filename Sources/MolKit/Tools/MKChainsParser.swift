
import Foundation

protocol _ByteCodeLabel {
    
}

private let MaxMonoAtom = 20
private let MaxMonoBond = 20

//! Chemical graph matching virtual machine
struct ByteCode {
    var type: Int
    init(type: Int) {
        self.type = type
        self._wrappedValue = nil
    }
    private var _wrappedValue: _ByteCodeLabel?
    //!< Eval - push current neighbors onto the stack
    var eval: MonOpStruct? {
        get {
            return _wrappedValue as? MonOpStruct
        }
        set {
            if newValue != nil {
                _wrappedValue = newValue
            }
        }
    }
    //!< Count - test the number of eval bonds
    var count: BinOpStruct? {
        get {
            return _wrappedValue as? BinOpStruct
        }
        set {
            if newValue != nil {
                _wrappedValue = newValue
            }
        }
    }
    //!< Element - test the element of current atom
    var elem: BinOpStruct? {
        get {
            return _wrappedValue as? BinOpStruct
        }
        set {
            if newValue != nil {
                _wrappedValue = newValue
            }
        }
    }
    //!< Ident - test the atom for backbone identity
    var ident: BinOpStruct? {
        get {
            return _wrappedValue as? BinOpStruct
        }
        set {
            if newValue != nil {
                _wrappedValue = newValue
            }
        }
    }
    //!< Local - test whether the atom has been visited
    var local: BinOpStruct? {
        get {
            return _wrappedValue as? BinOpStruct
        }
        set {
            if newValue != nil {
                _wrappedValue = newValue
            }
        }
    }
    //!< Assign - assign residue name, atom name and bond type to output
    var assign: AssignStruct? {
    get {
        return _wrappedValue as? AssignStruct
    }
    set {
        if newValue != nil {
            _wrappedValue = newValue
        }
    }
}
}

struct MonoAtomType: _ByteCodeLabel {
    var atomid: Int
    var elem: Int
    var bcount: Int
    var index: Int
}

struct MonoBondType: _ByteCodeLabel {
    var src: Int
    var dst: Int
    var index: Int
    var flag: Int
}

struct MonOpStruct: _ByteCodeLabel {
    var type: Int
    var next: _ByteCodeLabel?
}

struct BinOpStruct: _ByteCodeLabel {
    var type: Int
    var value: Int
    var tcond: _ByteCodeLabel?
    var fcond: _ByteCodeLabel?
}

//! Output array -- residue id, atom id, bond flags, etc.
struct AssignStruct: _ByteCodeLabel {
    var type: Int
    var resid: Int
    var atomid: [Int]?
    var bflags: [Int]?
}

struct StackType {
    var atom: Int
    var bond: Int
    var prev: Int
}

// static MonoAtomType MonoAtom[MaxMonoAtom];
private var MonoAtom = [MonoAtomType]()
// static MonoBondType MonoBond[MaxMonoBond];
private var MonoBond = [MonoBondType]()
// static int MonoAtomCount;
private var MonoAtomCount = 0
// static int MonoBondCount;
private var MonoBondCount = 0
// static StackType Stack[STACKSIZE];
private var Stackk = [StackType]()
// static int StackPtr;
private var StackPtr = 0
// static int  AtomIndex;
private var AtomIndex = 0
// static int  BondIndex;
private var BondIndex = 0
// static bool StrictFlag = false;
private var StrictFlag = false

//! The first available index for actual residues
//! 0, 1, 2 reserved for UNK, HOH, UNL
private let RESIDMIN =        4
//! The maximum number of residue IDs for this code
private let RESIDMAX =        64

//! An index of the residue names perceived during a run
//! 0, 1, and 2 reserved for UNK, HOH, LIG
private let ChainsResName = [
  /*0*/ "UNK",
  /*1*/ "HOH",
  /*2*/ "UNL",
  /*3*/ "ACE"]

private let ATOMMINAMINO =    4
private let ATOMMINNUCLEIC =  50
private let MAXPEPTIDE =      11
private let MAXNUCLEIC =      15
//! The number of amino acids recognized by this code
//! Currently: ILE, VAL, ALA, ASN, ASP, ARG, CYS, GLN, GLU
//!  GLY, HIS, HYP, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR
private let AMINOMAX =        21
//! The number of nucleic acids recognized by this code
//! Currently A, C, T, G, U, I
private let NUCLEOMAX =       6
private let STACKSIZE =       20

private let AI_N =            0
private let AI_CA =           1
private let AI_C =            2
private let AI_O =            3
private let AI_OXT =          37

private let AI_P =            38
private let AI_O1P =          39
private let AI_O2P =          40
private let AI_O5 =           41
private let AI_C5 =           42
private let AI_C4 =           43
private let AI_O4 =           44
private let AI_C3 =           45
private let AI_O3 =           46
private let AI_C2 =           47
private let AI_O2 =           48
private let AI_C1 =           49

private let BitN =            0x0001
private let BitNTer =         0x0002
private let BitNPro =         0x0004
private let BitNPT =          0x0008
private let BitCA =           0x0010
private let BitCAGly =        0x0020
private let BitC =            0x0100
private let BitCTer =         0x0200
private let BitCOXT =         0x0400
private let BitO =            0x1000
private let BitOXT =          0x2000

private let BitNAll =         0x000F
private let BitCAAll =        0x0030
private let BitCAll =         0x0700
private let BitOAll =         0x3000

private let BitP =            0x0001
private let BitPTer =         0x0002
private let BitOP =           0x0004
private let BitO5 =           0x0008
private let BitO5Ter =        0x0010
private let BitC5 =           0x0020
private let BitC4 =           0x0040
private let BitO4 =           0x0080
private let BitC3 =           0x0100
private let BitO3 =           0x0200
private let BitO3Ter =        0x0400
private let BitC2RNA =        0x0800
private let BitC2DNA =        0x1000
private let BitO2 =           0x2000
private let BitC1 =           0x4000

private let BitPAll =         0x0003
private let Bit05All =        0x0018
private let BitO3All =        0x0600
private let BitC2All =        0x1800

private let BC_ASSIGN =       0x01
private let BC_COUNT =        0x02
private let BC_ELEM =         0x03
private let BC_EVAL =         0x04
private let BC_IDENT =        0x05
private let BC_LOCAL =        0x06

private let BF_SINGLE =       0x01
private let BF_DOUBLE =       0x02
private let BF_TRIPLE =       0x04
private let BF_AROMATIC =     0x08


//////////////////////////////////////////////////////////////////////////////
// Global Variables / Tables
//////////////////////////////////////////////////////////////////////////////

//! The number of PDB atom type names recognized by this code
private let ATOMMAX =       68

//! PDB atom types (i.e., columns 13-16 of a PDB file)
//!  index numbers from this array are used in the pseudo-SMILES format
//!  for side-chains in the AminoAcids[] & Nucleotides[] global arrays below
let ChainsAtomName: [[String]] = [
    /*  0 */  [ " ", "N", " ", " " ],
    /*  1 */  [ " ", "C", "A", " " ],
    /*  2 */  [ " ", "C", " ", " " ],
    /*  3 */  [ " ", "O", " ", " " ],
    /*  4 */  [ " ", "C", "B", " " ],
    /*  5 */  [ " ", "S", "G", " " ],
    /*  6 */  [ " ", "O", "G", " " ],
    /*  7 */  [ " ", "C", "G", " " ],
    /*  8 */  [ " ", "O", "G", "1" ],
    /*  9 */  [ " ", "C", "G", "1" ],
    /* 10 */  [ " ", "C", "G", "2" ],
    /* 11 */  [ " ", "C", "D", " " ],
    /* 12 */  [ " ", "O", "D", " " ],
    /* 13 */  [ " ", "S", "D", " " ],
    /* 14 */  [ " ", "C", "D", "1" ],
    /* 15 */  [ " ", "O", "D", "1" ],
    /* 16 */  [ " ", "N", "D", "1" ],
    /* 17 */  [ " ", "C", "D", "2" ],
    /* 18 */  [ " ", "O", "D", "2" ],
    /* 19 */  [ " ", "N", "D", "2" ],
    /* 20 */  [ " ", "C", "E", " " ],
    /* 21 */  [ " ", "N", "E", " " ],
    /* 22 */  [ " ", "C", "E", "1" ],
    /* 23 */  [ " ", "O", "E", "1" ],
    /* 24 */  [ " ", "N", "E", "1" ],
    /* 25 */  [ " ", "C", "E", "2" ],
    /* 26 */  [ " ", "O", "E", "2" ],
    /* 27 */  [ " ", "N", "E", "2" ],
    /* 28 */  [ " ", "C", "E", "3" ],
    /* 29 */  [ " ", "C", "Z", " " ],
    /* 30 */  [ " ", "N", "Z", " " ],
    /* 31 */  [ " ", "C", "Z", "2" ],
    /* 32 */  [ " ", "C", "Z", "3" ],
    /* 33 */  [ " ", "O", "H", " " ],
    /* 34 */  [ " ", "N", "H", "1" ],
    /* 35 */  [ " ", "N", "H", "2" ],
    /* 36 */  [ " ", "C", "H", "2" ],
    /* 37 */  [ " ", "O", "X", "T" ],
    /* 38 */  [ " ", "P", " ", " " ],
    /* 39 */  [ " ", "O", "1", "P" ],
    /* 40 */  [ " ", "O", "2", "P" ],
    /* 41 */  [ " ", "O", "5", "*" ],
    /* 42 */  [ " ", "C", "5", "*" ],
    /* 43 */  [ " ", "C", "4", "*" ],
    /* 44 */  [ " ", "O", "4", "*" ],
    /* 45 */  [ " ", "C", "3", "*" ],
    /* 46 */  [ " ", "O", "3", "*" ],
    /* 47 */  [ " ", "C", "2", "*" ],
    /* 48 */  [ " ", "O", "2", "*" ],
    /* 49 */  [ " ", "C", "1", "*" ],
    /* 50 */  [ " ", "N", "9", " " ],
    /* 51 */  [ " ", "C", "8", " " ],
    /* 52 */  [ " ", "N", "7", " " ],
    /* 53 */  [ " ", "C", "5", " " ],
    /* 54 */  [ " ", "C", "6", " " ],
    /* 55 */  [ " ", "O", "6", " " ],
    /* 56 */  [ " ", "N", "6", " " ],
    /* 57 */  [ " ", "N", "1", " " ],
    /* 58 */  [ " ", "C", "2", " " ],
    /* 59 */  [ " ", "O", "2", " " ],
    /* 60 */  [ " ", "N", "2", " " ],
    /* 61 */  [ " ", "N", "3", " " ],
    /* 62 */  [ " ", "C", "4", " " ],
    /* 63 */  [ " ", "O", "4", " " ],
    /* 64 */  [ " ", "N", "4", " " ],
    /* 65 */  [ " ", "C", "5", " " ],
    /* 66 */  [ " ", "C", "5", "M" ],
    /* 67 */  [ " ", "C", "6", " " ]
]

//! Definition of side chains, associating overall residue name with
//!  the pseudo-SMILES pattern
struct ResidType {
    var name: String //!< Residue name, standardized by PDB
    var data: String  //!< pseudo-SMILES definition of side-chain
    init(_ name: String, _ data: String) {
        self.name = name
        self.data = data
    }
}

/**
   * Side chains for recognized amino acids using a pseudo-SMARTS syntax
   * for branching and bonds. Numbers indicate atom types defined by
   * OpenBabel::ChainsAtomName global array above.
   */
let AminoAcids: [ResidType] = [
    ResidType( "ILE", "1-4(-9-14)-10"                        ),
    ResidType( "VAL", "1-4(-9)-10"                           ),
    ResidType( "ALA", "1-4"                                  ),
    ResidType( "ASN", "1-4-7(=15)-19"                        ),
    ResidType( "ASP", "1-4-7(=15)-18"                        ),
    ResidType( "ARG", "1-4-7-11-21-29(=34)-35"               ),
    ResidType( "CYS", "1-4-5"                                ),
    ResidType( "GLN", "1-4-7-11(=23)-27"                     ),
    ResidType( "GLU", "1-4-7-11(=23)-26"                     ),
    ResidType( "GLY", "1"                                    ),
    ResidType( "HIS", "1-4-7^16~22^27^17~7"                  ),
    ResidType( "HYP", "1-4-7(-12)-11-0"                      ),
    ResidType( "LEU", "1-4-7(-14)-17"                        ),
    ResidType( "LYS", "1-4-7-11-20-30"                       ),
    ResidType( "MET", "1-4-7-13-20"                          ),
    ResidType( "PHE", "1-4-7~14^22~29^25~17^7"               ),
    ResidType( "PRO", "1-4-7-11-0"                           ),
    ResidType( "SER", "1-4-6"                                ),
    ResidType( "THR", "1-4(-8)-10"                           ),
    ResidType( "TRP", "1-4-7~14^24^25~17(^7)^28~32^36~31^25" ),
    ResidType( "TYR", "1-4-7~14^22~29(-33)^25~17^7"          )
]

// Other possible amino acid templates (less common)
/* Pyroglutamate (PCA):        1-4-7-11(=" OE ")-0  PDB Example: 1CEL */
/* Amino-N-Butyric Acid (ABA): 1-4-7                PDB Example: 1BBO */
/* Selenic Acid (SEC):         1-4-"SEG "(-15)-18   PDB Example: 1GP1 */

/**
* Side chains for recognized nucleotides using a pseudo-SMARTS syntax
* for branching and bonds. Numbers indicate atom types defined by
* OpenBabel::ChainsAtomName global array above.
*/
let Nucleotides: [ResidType] = [
    ResidType( "  A", "49-50-51-52-53-54(-56)-57-58-61-62(-53)-50"      ),
    ResidType( "  C", "49-57-58(-59)-61-62(-64)-65-67-57"               ),
    ResidType( "  G", "49-50-51-52-53-54(-55)-57-58(-60)-61-62(-53)-50" ),
    ResidType( "  T", "49-57-58(-59)-61-62(-63)-65(-66)-67-57"          ),
    ResidType( "  U", "49-57-58(-59)-61-62(-63)-65-67-57"               ),
    ResidType( "  I", "49-50-51-52-53-54(-55)-57-58-61-62(-53)-50"      )
]
//////////////////////////////////////////////////////////////////////////////
// Structure / Type Definitions
//////////////////////////////////////////////////////////////////////////////

//! Structure template for atomic patterns in residues for OBChainsParser
struct StructureTemplate {
    var flag: Int        //!< binary flag representing this atom type
    var elem: Int        //!< atomic number of this element
    var count: Int       //!< expected valence for this atom type
    var n1: Int          //!< mask 1 used by ConstrainBackbone() and MatchConstraint()
    var n2: Int          //!< mask 2 used by ConstrainBackbone() and MatchConstraint()
    var n3: Int          //!< mask 3 used by ConstrainBackbone() and MatchConstraint()
    var n4: Int          //!< mask 4 used by ConstrainBackbone() and MatchConstraint()
    init(_ flag: Int, _ elem: Int, _ count: Int, _ n1: Int, _ n2: Int, _ n3: Int, _ n4: Int) {
        self.flag = flag
        self.elem = elem
        self.count = count
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.n4 = n4
    }
}

/**
 * Generic template for peptide residue backbone. \n
 * col 1: bitmask \n
 * col 2: element number \n
 * col 3: neighbour count \n
 * col 4-7: 1-4 bitmasks for neighbour atoms (-6 means carbon)
 */
let Peptide: [StructureTemplate] = [
    /* N     */    StructureTemplate(  0x0001, 7, 2, 0x0030, 0x0100,      0, 0 ), //!< N
    /* NTer  */    StructureTemplate(  0x0002, 7, 1, 0x0030,      0,      0, 0 ), //!< NTre
    /* NPro  */    StructureTemplate(  0x0004, 7, 3, 0x0030, 0x0100,     -6, 0 ), //!< NPro
    /* NPT   */    StructureTemplate(  0x0008, 7, 2, 0x0030,     -6,      0, 0 ), //!< NPT
    /* CA    */    StructureTemplate(  0x0010, 6, 3, 0x000F, 0x0700,     -6, 0 ), //!< CA
    /* CAGly */    StructureTemplate(  0x0020, 6, 2, 0x0003, 0x0700,      0, 0 ), //!< CAGly
    /* C     */    StructureTemplate(  0x0100, 6, 3, 0x0030, 0x1000, 0x0005, 0 ), //!< C
    /* CTer  */    StructureTemplate(  0x0200, 6, 2, 0x0030, 0x1000,      0, 0 ), //!< CTer
    /* COXT  */    StructureTemplate(  0x0400, 6, 3, 0x0030, 0x1000, 0x2000, 0 ), //!< COXT
    /* O     */    StructureTemplate(  0x1000, 8, 1, 0x0700,      0,      0, 0 ), //!< O
    /* OXT   */    StructureTemplate(  0x2000, 8, 1, 0x0400,      0,      0, 0 )  //!< OXT
]

//! Generic template for peptide nucleotide backbone
let Nucleotide: [StructureTemplate] = [
    /* P     */    StructureTemplate(  0x0001, 15, 4, 0x0004, 0x0004, 0x0008, 0x0200 ),
    /* PTer  */    StructureTemplate(  0x0002, 15, 3, 0x0004, 0x0004, 0x0008,      0 ),
    /* OP    */    StructureTemplate(  0x0004,  8, 1, 0x0003,      0,      0,      0 ),
    /* O5    */    StructureTemplate(  0x0008,  8, 2, 0x0020, 0x0003,      0,      0 ),
    /* O5Ter */    StructureTemplate(  0x0010,  8, 1, 0x0020,      0,      0,      0 ),
    /* C5    */    StructureTemplate(  0x0020,  6, 2, 0x0018, 0x0040,      0,      0 ),
    /* C4    */    StructureTemplate(  0x0040,  6, 3, 0x0020, 0x0080, 0x0100,      0 ),
    /* O4    */    StructureTemplate(  0x0080,  8, 2, 0x0040, 0x4000,      0,      0 ),
    /* C3    */    StructureTemplate(  0x0100,  6, 3, 0x0040, 0x0600, 0x1800,      0 ),
    /* O3    */    StructureTemplate(  0x0200,  8, 2, 0x0100, 0x0001,      0,      0 ),
    /* O3Ter */    StructureTemplate(  0x0400,  8, 1, 0x0100,      0,      0,      0 ),
    /* C2RNA */    StructureTemplate(  0x0800,  6, 3, 0x0100, 0x4000, 0x2000,      0 ),
    /* C2DNA */    StructureTemplate(  0x1000,  6, 2, 0x0100, 0x4000,      0,      0 ),
    /* O2    */    StructureTemplate(  0x2000,  8, 1, 0x0800,      0,      0,      0 ),
    /* C1    */    StructureTemplate(  0x4000,  6, 3, 0x0080, 0x1800,     -7,      0 )
]

typealias Template = StructureTemplate

private func allocateByteCode(_ type: Int) -> ByteCode {
    var bytecode = ByteCode(type: type)
    return bytecode
}





/** @class OBChainsParser chains.h <openbabel/chains.h>
    @brief Perceives peptide or nucleotide chains and residues in an OBMol

    Perceive peptide or nucleotide chains and residues from atom connectivity.
    Based on original RasMol code by Roger Sayle and modified by Joe Corkery.
    For more on Roger's original talk, see:
    http://www.daylight.com/meetings/mug96/sayle/sayle.html
 */
class MKChainsParser {

    var PDecisionTree: [Template] = []
    var NDecisionTree: [Template] = []

    var ResMonoAtom: [Int] = []
    var ResMonoBond: [Int] = []

    var bitmasks: [UInt16] = []
    var visits: [Bool] = []  //!< mark visits to prevent looping
    var resids: [UInt8] = []
    var flags: [UInt8] = []
    var hetflags: [Bool] = []
    var atomids: [Int] = []
    var resnos: [Int] = []
    var sernos: [Int] = []   //!< array of residue serial numbers
    var hcounts: [Int] = []
    var chains: [Int] = []

    init() {}

    /**
       * Perceive macromolecular (peptide and nucleotide) residues and chains
       * @param mol The molecule to parse and update
       * @param nukeSingleResidue If only one residue is found, clear information
       * default = false  -- single residue files should still be recognized.
       */
    func perceiveChains(_ mol: MKMol) {
        
    }
    
    //! @name Step 1: Determine hetero atoms
      //@{
      /**
       * Determine HETATOM records for all atoms with a heavy valance of 0.
       * This includes HOH, Cl, Fe, ...
       *
       * Sets resids[i] & hetflags[i] for these atoms.
       * @todo add ions (Cl, Fe, ...)
       */
    private func determineAtoms(_ mol: MKMol) -> Bool  {
        fatalError()
    }
    
    //! @name Step 2: Determine connected chains
      //@{
      /**
       * Determine connected chains (e.g., subunits). Chains will be labeled A, B, C, ...
       * Ligands also get assigned a chain label. The chain for ligands will later be
       * replaced by ' '. The residue numbers will also be updated in this process to
       * make sure all ligands, HOH, ions, etc. have a unique residue number in the ' '
       * chain.
       *
       * Sets chains[i] for all atoms. (through RecurseChain())
       */
    private func determineConnectedChains(_ mol: MKMol) -> Bool {
        fatalError()
    }
    /**
       * Perform the actual work for DetermineConnectedChains(). Set chains[i]
       * to @p c for all atoms of the recursed chain.
       * @param mol The molecule.
       * @param i Index for the current atom. (RecurseChain() will be called for all neighbours)
       * @param c The chain which we are recusring. ('A' + count)
       * @return The number of heavy atoms in the recursed chain.
       */
    private func recurseChain(_ mol: MKMol, _ i: Int, _ c: Int) -> Int {
        fatalError()
    }

    //! @name Step 3: Determine peptide backbone
      //@{
      /**
       * Walk a peptide "backbone" atom sequence, from one residue to the next. This
       * function will look for N-CA-C-O sequences and mark these atoms.
       *
       * Sets bitmaks[i] for these atoms. (through ConstrainBackbone())
       * Sets resnos[i] for these atoms. (through TracePeptideChain())
       */
    private func determinePeptideBackbone(_ mol: MKMol) -> Bool {
        fatalError()
    }

    /**
       * First the bitmasks[i] will be OR-ed with Template::flag for all atoms based on
       * on Template::element and Template::count.
       *
       * Next, the bitmasks[i] are iteratively resolved by matching the
       * constraints in OpenBabel::Peptide or OpenBabel::Nucleotide.
       * @param mol The molecule.
       * @param templ OpenBabel::Peptide or OpenBabel::Nucleotide
       * @param tmax Number of entries in @p templ
       */
    //   void ConstrainBackbone(OBMol &mol, Template *templ, int tmax);
    
    //   /**
    //    * @return True if the bitmasks[i] for @p atom matches @p mask.
    //    */
    //   bool MatchConstraint(OBAtom *atom, int mask);

    //   /**
    //    * @return True if atom @p na and @p nb match the Template::n1 and
    //    * Template::n2.
    //    */
    //   bool Match2Constraints(Template *templ, OBAtom *na, OBAtom *nb);

    //   /**
    //    * @return True if atom @p na, @p nb and @p nc match the Template::n1,
    //    * Template::n2 and Template::n3.
    //    */
    //   bool Match3Constraints(Template *templ, OBAtom *na, OBAtom *nb, OBAtom *nc);

    //   /**
    //    * @return True if atom @p na, @p nb, @p nc and @p nd match the Template::n1,
    //    * Template::n2, Template::n3 and Template::n4.
    //    */
    //   bool Match4Constraints(Template *templ, OBAtom *na, OBAtom *nb, OBAtom *nc, OBAtom *nd);

    //   /**
    //    * Now we have the constrained bitmaks[i], trace N-CA-C-O-... and set
    //    * resnos[i] and atomids[i] for each N-CA-C-O sequence.
    //    *
    //    * Also adds BF_DOUBLE to flags[b] for< each carbonyl bond in N-CA-C=O.
    //    * @param mol The molecule.
    //    * @param i Index for the current atom. (TracePeptideChain() will be called for all neighbours)
    //    * @param r The residue number which we are tracing.
    //    */
    //   void TracePeptideChain(OBMol &mol, unsigned int i, int r);
    //   //@}

    //   //! @name Step 4: Determine peptide side chains
    //   //@{
    //   /**
    //    * Look for atoms with atomids[i] CA and identify their side chain.
    //    *
    //    * Sets resnos[i] and resids[i] for all identified residues (including the N-CA-C-O).
    //    * (through IdentifyResidue() and AssignResidue())
    //    */
    //   bool  DeterminePeptideSidechains(OBMol &);

    //   /**
    //    * Identify a residue based on the @p tree ByteCode.
    //    *
    //    * Sets resnos[i] for all sidechain atoms to the residue number of
    //    * the seed CA/C1 atom.
    //    * @param tree Bytecode for the residues. (OBChainsParser::PDecisionTree or OBChainsParser::NDecisionTree)
    //    * @param mol The molecule.
    //    * @param seed Atom index for the CA (peptides) or C1 (nucleotides) atom.
    //    * @param resno The residue number for this residue.
    //    * @return The resids[i] for the identified residue.
    //    */
    //   int IdentifyResidue(void *tree, OBMol &mol, unsigned int seed, int resno); // ByteCode *

    //   /**
    //    * Set resids[i] for all atoms where resids[i] = @p r and chains[i] = @p c.
    //    * @param mol The molecule.
    //    * @param r The residue number.
    //    * @param c The chain number.
    //    * @param i The residue id (resids[i] returned by IdentifyResidue())
    //    */
    //   void  AssignResidue(OBMol &mol, int r, int c, int i);
    //   //@}

    //   //! @name Step 5: Assign hydrogens
    //   //@{
    //   /**
    //    * Assign the resids[i], resnos[i], ... for all hydrogens based on the
    //    * atom they are bound to.
    //    */
    //   bool  DetermineHydrogens(OBMol &);
    //   //@}

    //   //! @name Step 6: Set the residue information
    //   //@{
    //   /**
    //    * Convert the private data vectors to OBResidue objects and add them to @p mol.
    //    * @param mol The molecule to parse and update
    //    * @param nukeSingleResidue If only one residue is found, clear information
    //    * default = false  -- single residue files should still be recognized.
    //    */
    //   void  SetResidueInformation(OBMol &, bool nukeSingleResidue);
    //   //@}

    //   //! @name Nucleic acids (analog to peptides)
    //   //@{
    //   /**
    //    * Walk a nucleic "backbone" atom sequence, from one residue to the next. This
    //    * function will look for ribose-5-P sequences and mark these atoms.
    //    *
    //    * Sets bitmaks[i] for these atoms. (through ConstrainBackbone())
    //    * Sets resnos[i] for these atoms. (through TraceNucleicChain())
    //    */
    //   bool  DetermineNucleicBackbone(OBMol &);

    //   /**
    //    * Now we have the constrained bitmaks[i], trace nucleic backbone and set
    //    * resnos[i] and atomids[i] for each ribose-5-P sequence.
    //    * @param mol The molecule.
    //    * @param i Index for the current atom. (TraceNucleicChain() will be called for all neighbours)
    //    * @param r The residue number which we are tracing.
    //    */
    //   void  TraceNucleicChain(OBMol &, unsigned int i, int r);

    //   /**
    //    * Look for atoms with atomids[i] C1 and identify their side chain.
    //    *
    //    * Sets resnos[i] and resids[i] for all identified residues.
    //    * (through IdentifyResidue() and AssignResidue())
    //    */
    //   bool  DetermineNucleicSidechains(OBMol &);
    //   //@}

    //   /**
    //    * Set up the chain perception to operate on the supplied molecule
    //    * by resizing and initializing the private data vectors.
    //    */
    //   void  SetupMol(OBMol &);
    
    //   /**
    //    * Delete all residues in @p mol
    //    */
    //   void  ClearResidueInformation(OBMol &mol);
    
    //   /**
    //    * Clear all private data vectors
    //    */
    //   void CleanupMol();
    
    //   /**
    //    * Construct and add ByteCode to the @p tree for a single residue.
    //    * @param tree Bytecode for the residues. (OBChainsParser::PDecisionTree or OBChainsParser::NDecisionTree)
    //    * @param resid The residue id.
    //    * @param smiles The pseudo-smiles string (from OpenBabel::AminoAcids or OpenBabel::Nucleotides)
    //    */
    //   void  DefineMonomer(void **tree, int resid, const char *smiles); // ByteCode **
    
    //   /**
    //    * @param ptr Element id (from OpenBabel::ChainsAtomName)
    //    * @return The element number.
    //    */
    //   int   IdentifyElement(char *ptr);
    
       /**
        * Parse a pseudo smiles from OpenBabel::AminoAcids or OpenBabel::Nucleotides.
        * @param smiles The pseudo-smiles string.
        * @param prev The previous position (used for recursing, use -1 to start).
        */
    private func parseSmiles(_ smiles: String, _ prev: Int) -> String {
        fatalError()
    }
    
    
    
    //   /**
    //    * Debugging function.
    //    */
    //   void DumpState();


    
    

}
