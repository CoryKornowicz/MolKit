
import Foundation

protocol _ByteCodeLabel {
    
}

private let MaxMonoAtom = 20
private let MaxMonoBond = 20

//! Chemical graph matching virtual machine
struct ByteCode: _ByteCodeLabel {
    
    var type: Int
    
    init(type: Int) {
        self.type = type
    }
    
    private var _wrappedEval: MonOpStruct?
    //!< Eval - push current neighbors onto the stack
    var eval: MonOpStruct {
        mutating get {
            if _wrappedEval == nil {
                _wrappedEval = MonOpStruct(type: type)
            }
            return _wrappedEval!
        }
        set {
            _wrappedEval = newValue
        }
    }

    private var _wrappedCount: BinOpStruct?
    //!< Count - test the number of eval bonds
    var count: BinOpStruct {
        mutating get {
            if _wrappedCount == nil {
                _wrappedCount = BinOpStruct(type: type)
            }
            return _wrappedCount!
        }
        set {
            _wrappedCount = newValue
        }
    }


    //!< Element - test the element of current atom
    private var _wrappedElement: BinOpStruct?
    var elem: BinOpStruct {
        mutating get {
            if _wrappedElement == nil {
                _wrappedElement = BinOpStruct(type: type)
            }
            return _wrappedElement!
        }
        set {
            _wrappedElement = newValue
        }
    }

    //!< Ident - test the atom for backbone identity
    private var _wrappedIdent: BinOpStruct?
    var ident: BinOpStruct {
        mutating get {
            if _wrappedIdent == nil {
                _wrappedIdent = BinOpStruct(type: type)
            }
            return _wrappedIdent!
        }
        set {
            _wrappedIdent = newValue
        }
    }
    
    //!< Local - test whether the atom has been visited
    private var _wrappedLocal: BinOpStruct?
    var local: BinOpStruct {
        mutating get {
            if _wrappedLocal == nil {
                _wrappedLocal = BinOpStruct(type: type)
            }
            return _wrappedLocal!
        }
        set {
            _wrappedLocal = newValue
        }
    }

    //!< Assign - assign residue name, atom name and bond type to output
    private var _wrappedAssign: AssignStruct?
    var assign: AssignStruct {
        mutating get {
            if _wrappedAssign == nil {
                _wrappedAssign = AssignStruct(type: type)
            }
            return _wrappedAssign!
        }
        set {
            _wrappedAssign = newValue
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
    var value: Int?
    var tcond: _ByteCodeLabel?
    var fcond: _ByteCodeLabel?
}

//! Output array -- residue id, atom id, bond flags, etc.
struct AssignStruct: _ByteCodeLabel {
    var type: Int
    var resid: Int?
    var atomid: [Int]?
    var bflags: [Int]?
}

struct StackType {
    var atom: Int?
    var bond: Int?
    var prev: Int?
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
private var ChainsResName: [String] = [String].init(repeating: "", count: RESIDMAX)

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
    let bytecode = ByteCode(type: type)
    return bytecode
}

fileprivate var ChainStack: [StackType] = [StackType].init(repeating: StackType(), count: STACKSIZE)

private func generateByteCodes(_ node: inout ByteCode?, _ resid: Int, _ curr: Int, _ prev: Int, _ bond: Int) {

    var done: Bool = false
    var found: Bool = false
    var ptr: ByteCode
    var neighbour: [StackType] = [StackType].init(repeating: StackType(), count: 4)

    if curr != prev {
        if MonoAtom[curr].atomid < ATOMMINAMINO {
            found = false 
            while node != nil && (node?.type == BC_IDENT) {
                if node?.ident.value == MonoAtom[curr].atomid {
                    if let tcond = node?.ident.tcond {
                        node = tcond as? ByteCode
                        found = true
                        break 
                    }
                } else {
                    if let fcond = node?.ident.fcond {
                        node = fcond as? ByteCode
                    }
                }
            }
            if !found {
                // TODO: this looks sloppy, can it be better?
                ptr = allocateByteCode(BC_IDENT)
                ptr.ident.tcond = nil
                ptr.ident.fcond = node
                node = ptr
                node = ptr.ident.tcond as? ByteCode
                ptr.ident.value = MonoAtom[curr].atomid
            }
            MonoBond[bond].index = BondIndex
            BondIndex += 1
            done = true
        } else if MonoAtom[curr].index != -1 {
            while node != nil && (node?.type == BC_IDENT) {
                node = node?.ident.fcond as? ByteCode
            }
            found = false
            while node != nil && node?.type == BC_LOCAL {
                if node?.local.value == MonoAtom[curr].index {
                    if let tcond = node?.local.tcond {
                        node = tcond as? ByteCode
                        found = true
                        break
                    }
                } else {
                    if let fcond = node?.local.fcond {
                        node = fcond as? ByteCode
                    }
                }
            }
            if !found {
                // TODO: this looks sloppy, can it be better?
                ptr = allocateByteCode(BC_LOCAL)
                ptr.local.tcond = nil
                ptr.local.fcond = node
                node = ptr
                node = ptr.local.tcond as? ByteCode
                ptr.local.value = MonoAtom[curr].atomid
            }
            MonoBond[bond].index = BondIndex
            BondIndex += 1
            done = true            
        } else {
            while node != nil && (node?.type == BC_IDENT) {
                node = node?.ident.fcond as? ByteCode
            }
            while node != nil && (node?.type == BC_LOCAL) {
                node = node?.local.fcond as? ByteCode
            }
            found = false
            while node != nil && (node?.type == BC_ELEM) {
                if node?.elem.value == MonoAtom[curr].elem {
                    if let tcond = node?.elem.tcond {
                        node = tcond as? ByteCode
                        found = true
                        break
                    }
                } else {
                    if let fcond = node?.elem.fcond {
                        node = fcond as? ByteCode
                    }
                }
            }
            if !found {
                ptr = allocateByteCode(BC_ELEM)
                ptr.elem.tcond = nil
                ptr.elem.fcond = node
                node = ptr
                node = ptr.elem.tcond as? ByteCode
                ptr.elem.value = MonoAtom[curr].elem
            }
            MonoAtom[curr].index = AtomIndex
            AtomIndex += 1
            MonoBond[bond].index = BondIndex
            BondIndex += 1
            done = false
        }
    } else {
        MonoAtom[curr].index = AtomIndex
        AtomIndex += 1
        done = false
    }   
    
    var count = 0 
    if !done {
        for i in 0..<MonoBondCount {
            if MonoBond[i].src == curr {
                if MonoBond[i].dst != prev {
                    neighbour[count].atom = MonoBond[i].dst
                    neighbour[count].bond = i
                    count += 1
                }
            } else if MonoBond[i].dst == curr {
                if MonoBond[i].src != prev {
                    neighbour[count].atom = MonoBond[i].src
                    neighbour[count].bond = i
                    count += 1
                }
            }
        }
        if node != nil && node?.type == BC_EVAL {
            found = false 
            node = node?.eval.next as? ByteCode
            while node != nil && node?.type == BC_COUNT {
                if node?.count.value == count {
                    if let tcond = node?.count.tcond {
                        node = tcond as? ByteCode
                        found = true
                        break
                    }
                } else {
                    if let fcond = node?.count.fcond {
                        node = fcond as? ByteCode
                    }
                }
            }
            if !found {
                ptr = allocateByteCode(BC_COUNT)
                ptr.count.tcond = nil
                ptr.count.fcond = node
                node = ptr
                node = ptr.count.tcond as? ByteCode
                ptr.count.value = count
            }
        } else if ((count != 0) || StrictFlag || (StackPtr != 0)) {
            ptr = allocateByteCode(BC_EVAL)
            ptr.eval.next = node
            node = ptr
            node = ptr.eval.next as? ByteCode

            ptr = allocateByteCode(BC_COUNT)
            ptr.count.tcond = nil
            ptr.count.fcond = node
            node = ptr
            node = ptr.count.tcond as? ByteCode
            ptr.count.value = count
        }
    }
    
    if count == 1 {
        generateByteCodes(&node, resid, neighbour[0].atom!, curr, neighbour[0].bond!)
    } else if count == 2 {
        let original = ChainStack[StackPtr]
        StackPtr += 1
        ChainStack[StackPtr-1] = neighbour[0]
        ChainStack[StackPtr-1].prev = curr
        generateByteCodes(&node, resid, neighbour[1].atom!, curr, neighbour[1].bond!)
        ChainStack[StackPtr-1] = neighbour[1]
        ChainStack[StackPtr-1].prev = curr
        generateByteCodes(&node, resid, neighbour[0].atom!, curr, neighbour[0].bond!)
        StackPtr -= 1
        ChainStack[StackPtr] = original
    } else if count != 0 {
        print("Maximum Monomer Fanout Exceeded!")
        print("Residue \(ChainsResName[resid]) atom \(curr)")
        print("Previous = \(prev) Fanout = \(count)")
        fatalError()
    } else if (StackPtr != 0) {
        StackPtr -= 1
        generateByteCodes(&node, resid, ChainStack[StackPtr].atom!, ChainStack[StackPtr].prev!, ChainStack[StackPtr].bond!)
        StackPtr += 1
    } else if node == nil {
        ptr = allocateByteCode(BC_ASSIGN)
        ptr.assign.resid = resid
        ptr.assign.atomid = [Int](repeating: 0, count: AtomIndex)
        for i in 0..<MonoAtomCount {
            if MonoAtom[i].index != -1 {
                ptr.assign.atomid![MonoAtom[i].index] = MonoAtom[i].atomid
            }
        }
        if BondIndex != 0 {
            ptr.assign.bflags = [Int](repeating: 0, count: BondIndex)
            for i in 0..<MonoBondCount {
                if MonoBond[i].index != -1 {
                    ptr.assign.bflags![MonoBond[i].index] = MonoBond[i].flag
                }
            }
        }
        node = ptr
    } else if node?.type == BC_ASSIGN {
        if node?.assign.resid != resid {
            print("Duplicated Monomer Specification!")
            print("Residue \(ChainsResName[resid]) matches residue \(ChainsResName[node!.assign.resid!])")
        }
    }
    if curr != prev {
        if !done {
            MonoAtom[curr].index = -1
            AtomIndex -= 1
        }
        MonoBond[bond].index = -1
        BondIndex -= 1
    }
}



/** @class OBChainsParser chains.h <openbabel/chains.h>
    @brief Perceives peptide or nucleotide chains and residues in an OBMol

    Perceive peptide or nucleotide chains and residues from atom connectivity.
    Based on original RasMol code by Roger Sayle and modified by Joe Corkery.
    For more on Roger's original talk, see:
    http://www.daylight.com/meetings/mug96/sayle/sayle.html
 */
class MKChainsParser {

    var PDecisionTree: ByteCode? = nil 
    var NDecisionTree: ByteCode? = nil

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
    var chains: [Character] = []

    init() {
        // moved here since they are not allowed at the top level
        ChainsResName[0] = "UNK"
        ChainsResName[1] = "HOH"
        ChainsResName[2] = "UNL"
        ChainsResName[3] = "ACE"
        
        
        var res = RESIDMIN
        PDecisionTree = nil
        for i in 0..<AMINOMAX {
            ChainsResName[res] = AminoAcids[i].name
            defineMonomer(&PDecisionTree, res, AminoAcids[i].data)
            res += 1
        }
        NDecisionTree = nil
        for i in 0..<NUCLEOMAX {
            ChainsResName[res] = Nucleotides[i].name
            defineMonomer(&NDecisionTree, res, Nucleotides[i].data)
            res += 1
        }
    }

    /**
       * Perceive macromolecular (peptide and nucleotide) residues and chains
       * @param mol The molecule to parse and update
       * @param nukeSingleResidue If only one residue is found, clear information
       * default = false  -- single residue files should still be recognized.
       */
    @discardableResult
    func perceiveChains(_ mol: MKMol, _ nukeSingleResidue: Bool = false) -> Bool {

        var result: Bool = true
        var idx: Int = 0
        
        setupMol(mol)
        clearResidueInformation(mol)
        result = determineHetAtoms(mol)          && result
        result = determineConnectedChains(mol)   && result
        result = determinePeptideBackbone(mol)   && result
        result = determinePeptideSidechains(mol) && result
        result = determineNucleicBackbone(mol)   && result
        result = determineNucleicSidechains(mol) && result
        result = determineHydrogens(mol)         && result

        // Partially identified residues
        // example: CSD in 1LWF (CYS with two Os on the S)
        var changed: Bool = false
        var invalidResidues: [(Character, Int)] = []
        repeat {
            changed = false
            for atom in mol.getAtomIterator() {
                idx = atom.getIdx() - 1
                if resids[idx] == 0 { // UNK
                    for nbr in atom.getNbrAtomIterator()! {
                        let idx2 = nbr.getIdx() - 1
                        if resids[idx2] != 0 { // !UNK
                            if atomids[idx2] == AI_N || atomids[idx2] == AI_C {
                                // bound to backbone-N/C
                                hetflags[idx] = true
                                resids[idx] = 3 // ACE
                                atomids[idx] = -1
                            } else {
                                resnos[idx] = resnos[idx2]
                                resids[idx] = resids[idx2]
                                changed = true
                                var addResidue = true
                                for i in 0..<invalidResidues.count {
                                    if invalidResidues[i].0 == chains[idx2] &&
                                        invalidResidues[i].1 == resnos[idx2] {
                                        addResidue = false
                                    }
                                }
                                if addResidue {
                                    invalidResidues.append((chains[idx2], resnos[idx2]))
                                }
                            }
                        }
                    }
                }
            }
        } while changed
        
        for i in 0..<invalidResidues.count {
            for atom in mol.getAtomIterator() {
                idx = atom.getIdx() - 1
                if ( (invalidResidues[i].0 == chains[idx]) &&
                    (invalidResidues[i].1 == resnos[idx]) ) {
                    hetflags[idx] = true
                    resids[idx] = 0 // UNK
                    atomids[idx] = -1
                }
            }
        }
        invalidResidues.removeAll()

        // number the element in the ' ' chain (water, ions, ligands, ...)
        var resno: Int = 1
        for atom in mol.getAtomIterator() {
            idx = atom.getIdx() - 1
            if atom.getAtomicNum() == MKElements.Hydrogen.atomicNum {
                chains[idx] = " "
                resnos[idx] = resno
                resno += 1
            } else {
                if resids[idx] != 0 { // UNK
                    continue
                }
                if hetflags[idx] {
                    continue
                }
                let chain = chains[idx]
                for b in mol.getAtomIterator() {
                    let idx2 = b.getIdx() - 1
                    if chains[idx2] == chain && !hetflags[idx2] {
                        hetflags[idx2] = true
                        chains[idx2] = " "
                        resnos[idx2] = resno
                        resids[idx2] = 2 // unknown ligand
                    }
                }
                resno += 1
            }
        }
        setResidueInformation(mol, nukeSingleResidue)
        cleanupMol()
        mol.setChainsPerceived()
        return result
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
    @discardableResult
    private func determineHetAtoms(_ mol: MKMol) -> Bool  {
        // find un-connected atoms (e.g., HOH oxygen atoms)
        for atom in mol.getAtomIterator() where atom.getAtomicNum() == MKElements.Hydrogen.atomicNum || atom.getHeavyDegree() != 0 {
             let idx = atom.getIdx() - 1 
             if atom.getAtomicNum() == MKElements.Oxygen.atomicNum {
                resids[idx] = 1 
                hetflags[idx] = true
             }
        }
        return true
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
    @discardableResult
    private func determineConnectedChains(_ mol: MKMol) -> Bool {
        var count = 0
        var resno = 1
        let numAtoms = mol.numAtoms()

        for atom in mol.getAtomIterator() {
            let idx = atom.getIdx() - 1
            if !hetflags[idx] && chains[idx] == " " && atom.getAtomicNum() != MKElements.Hydrogen.atomicNum {
                // recurse chain
                var char: Character = "A"
                let size = recurseChain(mol, idx, char.offset(by: count))
                // size = number of heavy atoms in residue chain
                if size < 4 { // small ligand, probably
                    if size == 1 && atom.getAtomicNum() == MKElements.Oxygen.atomicNum {
                        resids[idx] = 1 // HOH
                    } else {
                        resids[idx] = 2 // Unknown ligand
                    }
                    for i in 0..<numAtoms where chains[i] == char.offset(by: count) {
                        hetflags[i] = true
                        resnos[i] = resno
                        chains[i] = " "
                    }
                    resno += 1
                } else {
                    count += 1 // number of connected chains
                    if count > 26 { // out of chain IDs
                        break
                    }
                }
            }
        }
        return true
    }

    /**
       * Perform the actual work for DetermineConnectedChains(). Set chains[i]
       * to @p c for all atoms of the recursed chain.
       * @param mol The molecule.
       * @param i Index for the current atom. (RecurseChain() will be called for all neighbours)
       * @param c The chain which we are recusring. ('A' + count)
       * @return The number of heavy atoms in the recursed chain.
       */
    private func recurseChain(_ mol: MKMol, _ i: Int, _ c: Character) -> Int {
        guard let atom = mol.getAtom(i + 1) else { return 0}
        // ignore hydrogens 
        if atom.getAtomicNum() == MKElements.Hydrogen.atomicNum {
            return 0
        }
        var result = 1 
        chains[i] = c
        // recurse till we have all atoms for this chain 
        for nbr in atom.getNbrAtomIterator()! {
            let idx = nbr.getIdx() - 1
            if chains[idx] == " " {
                result += recurseChain(mol, idx, c)
            }
        }
        // and return how many we found 
        return result
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
        constrainBackbond(mol, Peptide, MAXPEPTIDE)
        
        let numAtoms = mol.numAtoms()
        
        // Cyclic peptides have no NTer (1SKI)
        
        var foundNTer = false
        for i in 0..<numAtoms {
            if (Int(bitmasks[i]) & BitNTer) != 0 {
                foundNTer = true
                break
            }
        }
        if !foundNTer {
            for i in 0..<numAtoms {
                if (Int(bitmasks[i]) & BitNAll) != 0 {
                    bitmasks[i] |= UInt16(BitNTer)
                }
            }
        }

        /* Order Peptide Backbone */
        for i in 0..<numAtoms where atomids[i] == -1 {
            if (Int(bitmasks[i]) & BitNTer) != 0 {
                atomids[i] = AI_N
                tracePeptideChain(mol, i, 1)
            } else if (((Int(bitmasks[i]) & BitNPT) != 0) && ((Int(bitmasks[i]) & BitN) != 0)) {
                atomids[i] = AI_N
                tracePeptideChain(mol, i, 1)
            } 
        }
        
        /* Carbonyl Double Bond */
        for bond in mol.getBondIterator() {
            if (atomids[bond.getBeginAtomIdx() - 1] == 2 && atomids[bond.getEndAtomIdx() - 1] == 3) ||
                (atomids[bond.getBeginAtomIdx() - 1] == 3 && atomids[bond.getEndAtomIdx() - 1] == 2) {
                flags[Int(bond.getIdx())] |= UInt8(BF_DOUBLE)
            }
        }
        
        return true
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
    func constrainBackbond(_ mol: MKMol, _ templ: [Template], _ tmax: Int) {

        var neighhbour: [MKAtom?] = [nil, nil, nil, nil, nil, nil]
        var change: Bool = false
        var result: Bool = false 
        var idx: Int 
        var count: Int 

        var na: MKAtom?
        var nb: MKAtom?
        var nc: MKAtom?
        var nd: MKAtom?


        // first pass 

        for atom in mol.getAtomIterator() {
            idx = atom.getIdx() - 1
            bitmasks[idx] = 0
            for i in 0..<tmax {
                if templ[i].elem == atom.getAtomicNum() || templ[i].count == atom.getHeavyDegree() {
                    bitmasks[idx] |= UInt16(templ[i].flag)
                }
            }
        }
        
        // second pass 

        repeat {
            change = false
            for atom in mol.getAtomIterator() {
                idx = atom.getIdx() - 1
                if bitmasks[idx] != 0 { // determine neighbors 
                    count = 0 
                    for nbr in atom.getNbrAtomIterator()! where nbr.getAtomicNum() != MKElements.Hydrogen.atomicNum {
                        neighhbour[count] = nbr
                        count += 1
                    }
                    if count >= 1 {
                        na = neighhbour[0]
                    }
                    if count >= 2 {
                        nb = neighhbour[1]
                    }
                    if count >= 3 {
                        nc = neighhbour[2]
                    }
                    if count >= 4 {
                        nd = neighhbour[3]
                    }
                    for i in 0..<tmax where templ[i].flag & Int(bitmasks[idx]) != 0 {
                        let pep = templ[i]
                        result = true
                        if count == 4 {
                            if na != nil, nb != nil, nc != nil, nd != nil {
                                result = match4Constraints(pep, na!, nb!, nc!, nd!)
                            }
                        } else if count == 3 {
                            if na != nil, nb != nil, nc != nil {
                                result = match3Constraints(pep, na!, nb!, nc!)
                            }
                        } else if count == 2 {
                            if na != nil, nb != nil {
                                result = match2Constraints(pep, na!, nb!)
                            }
                        } else if count == 1 {
                            if na != nil {
                                result = match1Constraint(na!, pep.n1)
                            }
                        }
                        if result == false {
                            bitmasks[idx] &= ~UInt16(pep.flag)
                            change = true
                        }
                    }
                }
            }
        } while change
    }
    
    /**
    * @return True if the bitmasks[i] for @p atom matches @p mask.
    */
    func match1Constraint(_ atom: MKAtom, _ mask: Int) -> Bool {
        if mask < 0 {
            return atom.getAtomicNum() == -mask
        } else {
            return Int(bitmasks[atom.getIdx() - 1]) & mask != 0
        }
    }
    
    /**
    * @return True if atom @p na and @p nb match the Template::n1 and
    * Template::n2.
    */
    func match2Constraints(_ templ: Template, _ na: MKAtom, _ nb: MKAtom) -> Bool {
        if match1Constraint(na, templ.n2) {
            if match1Constraint(nb, templ.n1) { return true }
        } 
        if match1Constraint(nb, templ.n2) {
            if match1Constraint(na, templ.n1) { return true }
        }
        return false
    }

    /**
    * @return True if atom @p na, @p nb and @p nc match the Template::n1,
    * Template::n2 and Template::n3.
    */
    func match3Constraints(_ templ: Template, _ na: MKAtom, _ nb: MKAtom, _ nc: MKAtom) -> Bool {
        if match1Constraint(na, templ.n3) {
            if match2Constraints(templ, nb, nc) { return true }
        }
        if match1Constraint(nb, templ.n3) {
            if match2Constraints(templ, na, nc) { return true }
        }
        if match1Constraint(nc, templ.n3) {
            if match2Constraints(templ, na, nb) { return true }
        }
        return false
    }

    /**
    * @return True if atom @p na, @p nb, @p nc and @p nd match the Template::n1,
    * Template::n2, Template::n3 and Template::n4.
    */
    func match4Constraints(_ templ: Template, _ na: MKAtom, _ nb: MKAtom, _ nc: MKAtom, _ nd: MKAtom) -> Bool {
        if match1Constraint(na, templ.n4) {
            if match3Constraints(templ, nb, nc, nd) { return true }
        }
        if match1Constraint(nb, templ.n4) {
            if match3Constraints(templ, na, nc, nd) { return true }
        }
        if match1Constraint(nc, templ.n4) {
            if match3Constraints(templ, na, nb, nd) { return true }
        }
        if match1Constraint(nd, templ.n4) {
            if match3Constraints(templ, na, nb, nc) { return true }
        }
        return false
    }

      /**
       * Now we have the constrained bitmaks[i], trace N-CA-C-O-... and set
       * resnos[i] and atomids[i] for each N-CA-C-O sequence.
       *
       * Also adds BF_DOUBLE to flags[b] for< each carbonyl bond in N-CA-C=O.
       * @param mol The molecule.
       * @param i Index for the current atom. (TracePeptideChain() will be called for all neighbours)
       * @param r The residue number which we are tracing.
       */
      //@}
    func tracePeptideChain(_ mol: MKMol, _ i: Int, _ r: Int ) {
        var neighbor: [Int] = [0,0,0,0]
        var na: Int = 0
        var nb: Int = 0
        var nc: Int = 0
        var count: Int = 0 

        //determine neighbors 
        guard let atom = mol.getAtom(i+1) else { return }
        let idx = atom.getIdx() - 1 
        if visits[i] {
            return 
        }
        visits[i] = true

        for nbr in atom.getNbrAtomIterator()! where nbr.getAtomicNum() != MKElements.Hydrogen.atomicNum {
            neighbor[count] = nbr.getIdx() - 1
            count += 1
        }

        resnos[idx] = r 
        if count >= 1 { 
            na = neighbor[0]
        }
        if count >= 2 { 
            nb = neighbor[1]
        }
        if count >= 3 { 
            nc = neighbor[2]
        }

        switch atomids[i] {
            case AI_N: 
                for j in 0..<count {
                    if (Int(bitmasks[neighbor[j]]) & BitCAll) != 0 {
                        atomids[neighbor[j]] = AI_CA
                        if !visits[neighbor[j]] {
                            tracePeptideChain(mol, neighbor[j], r)
                        }
                    } 
                }
            case AI_CA: 
            if count == 3 {
                if bitmasks[na] & BitNAll != 0 {
                    na = nc 
                } else if bitmasks[nb] & BitNAll != 0 {
                    nb = nc 
                } 
                var j: Int = 0
                var k: Int = 0 
                if bitmasks[na] & BitC != 0 {
                    j = na
                    k = nb
                } else if bitmasks[nb] & BitC != 0 {
                    j = nb
                    k = na
                } else if bitmasks[na] & BitCAll != 0 {
                    j = na
                    k = nb
                } else if bitmasks[nb] & BitCAll != 0 {
                    j = nb
                    k = na
                } 

                atomids[j] = AI_C
                atomids[k] = 0 

                if !visits[j] {
                    tracePeptideChain(mol, j, r)
                }
            } else if count == 2 {
                if bitmasks[na] & BitCAll != 0 {
                    atomids[na] = AI_C
                    if !visits[na] {
                        tracePeptideChain(mol, na, r)
                    }
                } else if bitmasks[nb] & BitCAll != 0 {
                    atomids[nb] = AI_C
                    if !visits[nb] {
                        tracePeptideChain(mol, nb, r)
                    }
                }
            }
            case AI_C: 
                var k = AI_O 
                for j in 0..<count {
                    if bitmasks[neighbor[j]] & BitNAll != 0 {
                        atomids[neighbor[j]] = AI_N
                        if !visits[neighbor[j]] {
                            tracePeptideChain(mol, neighbor[j], r + 1)
                        }
                    } else if bitmasks[neighbor[j]] & BitOAll != 0 {
                        atomids[neighbor[j]] = k
                        resnos[neighbor[j]] = r
                        k = AI_OXT // OXT
                    }
                }
            default: break 
        }
    }

    //! @name Step 4: Determine peptide side chains
    //@{
    /**
    * Look for atoms with atomids[i] CA and identify their side chain.
    *
    * Sets resnos[i] and resids[i] for all identified residues (including the N-CA-C-O).
    * (through IdentifyResidue() and AssignResidue())
    */
    @discardableResult
    func determinePeptideSidechains(_ mol: MKMol) -> Bool {
        var resid = 0
        let max = mol.numAtoms()
        for i in 0..<max where atomids[i] == AI_CA {
            resid = identifyResidue(&PDecisionTree, mol, i, resnos[i])
            assignResidue(mol, resnos[i], chains[i], resid)
        }
        return true
    }

      /**
       * Identify a residue based on the @p tree ByteCode.
       *
       * Sets resnos[i] for all sidechain atoms to the residue number of
       * the seed CA/C1 atom.
       * @param tree Bytecode for the residues. (OBChainsParser::PDecisionTree or OBChainsParser::NDecisionTree)
       * @param mol The molecule.
       * @param seed Atom index for the CA (peptides) or C1 (nucleotides) atom.
       * @param resno The residue number for this residue.
       * @return The resids[i] for the identified residue.
       */
    func identifyResidue(_ tree: inout ByteCode?, _ mol: MKMol, _ seed: Int, _ resno: Int) -> Int {
        var ptr: ByteCode? = tree 
        var AtomCount = 0
        var BondCount = 0
        var curr = 0
        var prev = 0
        var bond = 0
        var bcount = 0

        Stackk[0].atom = seed
        Stackk[0].prev = seed
        StackPtr = 0

        ResMonoAtom[0] = seed
        AtomCount = 1
        BondCount = 0
        
        while ptr != nil {
            switch ptr?.type {
            case BC_IDENT:
                curr = Stackk[StackPtr - 1].atom!
                if atomids[curr] == ptr?.ident.value {
                    bond = Stackk[StackPtr - 1].bond!
                    ResMonoBond[BondCount] = bond
                    BondCount += 1
                    ptr = ptr?.ident.tcond as? ByteCode
                    StackPtr -= 1
                } else {
                    ptr = ptr?.ident.fcond as? ByteCode
                }
            case BC_LOCAL:
                curr = Stackk[StackPtr - 1].atom!
                if curr == ResMonoAtom[ptr!.local.value!] {
                    bond = Stackk[StackPtr - 1].bond!
                    ResMonoBond[BondCount] = bond
                    BondCount += 1
                    ptr = ptr?.local.tcond as? ByteCode
                    StackPtr -= 1
                } else {
                    ptr = ptr?.local.fcond as? ByteCode
                }
            case BC_ELEM:
                curr = Stackk[StackPtr - 1].atom! 
                if mol.getAtom(curr + 1)?.getAtomicNum() == ptr?.elem.value {
                    bond = Stackk[StackPtr - 1].bond!
                    ResMonoAtom[AtomCount] = curr
                    AtomCount += 1
                    ResMonoBond[BondCount] = bond
                    resnos[curr] = resno
                    BondCount += 1
                    ptr = ptr?.elem.tcond as? ByteCode
                    StackPtr -= 1
                } else {
                    ptr = ptr?.elem.fcond as? ByteCode
                }
            case BC_EVAL:
                curr = Stackk[StackPtr].atom!
                prev = Stackk[StackPtr].prev!

                if let atom = mol.getAtom(curr + 1) {
                    for nbr in atom.getNbrAtomIterator()! where nbr.getAtomicNum() != MKElements.Hydrogen.atomicNum {  
                        let j = nbr.getIdx() - 1

                        if !((curr == prev) && bitmasks[j] != 0) && j != prev {
                            Stackk[StackPtr].prev = curr
                            Stackk[StackPtr].atom = j
                            if let b = mol.getBond(atom, nbr) {
                                Stackk[StackPtr].bond = Int(b.getIdx())
                            } else {
                                fatalError("Bond not found when it should be accessible")
                            }
                            StackPtr += 1
                            bcount += 1
                        }
                    }
                    ptr = ptr?.eval.next as? ByteCode
                } else {
                    fatalError("ATOM accessed out of bounds") 
                } // TODO: fix potential index issue here

            case BC_COUNT:
                if bcount == ptr?.count.value {
                    ptr = ptr?.count.tcond as? ByteCode
                } else {
                    ptr = ptr?.count.fcond as? ByteCode
                }
            case BC_ASSIGN:
                for i in 0..<AtomCount {
                    if (bitmasks[ResMonoAtom[i]] == 0) {
                        let j = ptr?.assign.atomid![i]
                        atomids[ResMonoAtom[i]] = j!
                    }
                }
                for i in 0..<BondCount {
                    let j = ptr?.assign.bflags![i]
                    flags[ResMonoBond[i]] = UInt8(j!)
                }
                return (ptr?.assign.resid)!
            default: // illegal instruction
                return 0
            } // switch
        } // while loop through atoms
        return 0
    }

    /**
    * Set resids[i] for all atoms where resids[i] = @p r and chains[i] = @p c.
    * @param mol The molecule.
    * @param r The residue number.
    * @param c The chain number.
    * @param i The residue id (resids[i] returned by IdentifyResidue())
    */
    //@}
    func assignResidue(_ mol: MKMol, _ r: Int, _ c: Character, _ i: Int) {
        let max = mol.numAtoms()
        for j in 0..<max where resnos[j] == r && chains[j] == c && !hetflags[j]{
            resids[j] = UInt8(i)
        }
    }

    //! @name Step 5: Assign hydrogens
    //@{
    /**
    * Assign the resids[i], resnos[i], ... for all hydrogens based on the
    * atom they are bound to.
    */
    //@}
    func determineHydrogens(_ mol: MKMol) -> Bool {
        let max = mol.numAtoms()
        for i in 0..<max {
            hcounts[i] = 0
        }
        /* First Pass */
        for atom in mol.getAtomIterator() {
            if atom.getAtomicNum() == MKElements.Hydrogen.atomicNum {
                if let nbr = atom.getNbrAtomIterator()!.next() {
                    let idx = atom.getIdx() - 1
                    let sidx = nbr.getIdx() - 1

                    hcounts[idx] = hcounts[sidx] + 1
                    hetflags[idx] = hetflags[sidx]
                    atomids[idx] = atomids[sidx]
                    resids[idx] = resids[sidx]
                    resnos[idx] = resnos[sidx]
                    chains[idx] = chains[sidx]
                }
            }
        }
        /* Second Pass */
        for atom in mol.getAtomIterator() {
            if atom.getAtomicNum() == MKElements.Hydrogen.atomicNum {
                if let nbr = atom.getNbrAtomIterator()!.next() {
                    if hcounts[nbr.getIdx() - 1] == 1 {
                        hcounts[atom.getIdx() - 1] = 0
                    }
                }
            }
        }
        return true
    }

    //! @name Step 6: Set the residue information
    //@{
    /**
    * Convert the private data vectors to OBResidue objects and add them to @p mol.
    * @param mol The molecule to parse and update
    * @param nukeSingleResidue If only one residue is found, clear information
    * default = false  -- single residue files should still be recognized.
    */
    //@}
    func setResidueInformation(_ mol: MKMol, _ nukeSingleResidue: Bool = false) {

        var buffer = ""
        var atomid: String
        var name: String
        var symbol = ""

        let numAtoms = mol.numAtoms()
        var resmap = [Character: [Int: MKResidue]]()
        //DumpState();

        // correct serine OG
        for i in 0..<numAtoms where resids[i] == RESIDMIN + 17 { // serine
            if atomids[i] == -1 {
                guard let atom = mol.getAtom(i + 1) else {
                    fatalError("Atom not found when it should be accessible")
                }
                for nbr in atom.getNbrAtomIterator()! {
                    if atomids[nbr.getIdx() - 1] == 4 { // CB
                        atomids[i] = 6 // OG
                    }
                }
            }
        }
        for i in 0..<numAtoms {
            guard let atom = mol.getAtom(i + 1) else { // TODO: fix potential atom index issue
                fatalError("Atom not found when it should be accessible")
            }

            if atomids[i] == -1 {
                symbol = MKElements.getSymbol(atom.getAtomicNum())
                // TODO: Maybe enforce symbol sizes here
                if symbol.count > 1 {
                    buffer = String(symbol.prefix(1)) + String(symbol.suffix(1).uppercased())
                } else {
                    buffer = " " + symbol
                }
                buffer += "  "
                buffer += "  "
                assert(buffer.length == 4, "Incorrect amount of characters in buffer")
                 
            } else if atom.getAtomicNum() == MKElements.Hydrogen.atomicNum {
                // TODO: Fix this because it is not moving the characters around currently 
            //     if (hcounts[i]) {
                //     snprintf(buffer, BUFF_SIZE, "H%.2s%c", ChainsAtomName[atomids[i]]+2, hcounts[i]+'0');
                //     if (buffer[1] == ' ') {
                //         buffer[1] = buffer[3];
                //         buffer[2] = '\0';
                //     }
                //     else if (buffer[2] == ' ') {
                //         buffer[2] = buffer[3];
                //         buffer[3] = '\0';
                //     }
            //     } else {
            //     snprintf(buffer, BUFF_SIZE, "H%.2s", ChainsAtomName[atomids[i]]+2);
            //     }
                if hcounts[i] != 0 {
                    buffer = "H" + ChainsAtomName[atomids[i]].suffix(2).reduce("", +) + String(hcounts[i])
                    if buffer.length == 3 {
                        buffer += " "
                    }
                } else {
                    buffer = "H" + ChainsAtomName[atomids[i]].suffix(2).reduce("", +)
                    buffer += "  "
                    buffer += "  "
                }
                assert(buffer.length == 4, "Incorrect amount of characters in buffer")
            } else {
                buffer = ChainsAtomName[atomids[i]].reduce("", +)
                assert(buffer.length == 4, "Incorrect amount of characters in buffer")
            }

            if buffer[3] == " " {
                buffer = String(buffer.prefix(3))
            }
            
            atomid = buffer[0] == " " ? String(buffer.suffix(3)) : buffer

            if let residue = resmap[chains[i]]?[resnos[i]] {
                residue.addAtom(atom)
                residue.setAtomID(atom, atomid)
                residue.setHetAtom(atom, hetflags[i])
                residue.setSerialNum(atom, sernos[i])
            } else {
                name = ChainsResName[Int(resids[i])]
                let residue = MKResidue()
                residue.setName(name)
                residue.setNum(resnos[i])
                residue.setChain(String(chains[i]))
                residue.addAtom(atom)
                residue.setAtomID(atom, atomid)
                residue.setHetAtom(atom, hetflags[i])
                residue.setSerialNum(atom, sernos[i])
                if resmap[chains[i]] == nil {
                    resmap[chains[i]] = [Int: MKResidue]()
                }
                resmap[chains[i]]![resnos[i]] = residue
            }

        }

        if mol.numResidues() == 1 && nukeSingleResidue {
            mol.deleteResidue(mol.getResidue(0)!)
        } else if mol.numResidues() == 1 && mol.getResidue(0)!.getName() == "UNK" {
            mol.deleteResidue(mol.getResidue(0)!)
        }
    }

    //! @name Nucleic acids (analog to peptides)
    //@{
    /**
    * Walk a nucleic "backbone" atom sequence, from one residue to the next. This
    * function will look for ribose-5-P sequences and mark these atoms.
    *
    * Sets bitmaks[i] for these atoms. (through ConstrainBackbone())
    * Sets resnos[i] for these atoms. (through TraceNucleicChain())
    */
    func determineNucleicBackbone(_ mol: MKMol) -> Bool {
        constrainBackbond(mol, Nucleotide, MAXNUCLEIC)
        let max = mol.numAtoms()
        // Order Nucleic backbone 
        for i in 0..<max where atomids[i] == -1 {
            if( bitmasks[i] & BitPTer ) != 0 {
              atomids[i] = AI_P
              traceNucleicChain(mol, i, 1)
            } else if( bitmasks[i] & BitO5Ter ) != 0 {
              atomids[i] = AI_O5
              traceNucleicChain(mol, i, 1)
            }
        }
        return true
    }

    /**
    * Now we have the constrained bitmaks[i], trace nucleic backbone and set
    * resnos[i] and atomids[i] for each ribose-5-P sequence.
    * @param mol The molecule.
    * @param i Index for the current atom. (TraceNucleicChain() will be called for all neighbours)
    * @param r The residue number which we are tracing.
    */
    func traceNucleicChain(_ mol: MKMol, _ i: Int, _ r: Int) {
        var neighbor: [Int] = [0, 0, 0, 0]
        if visits[i] {
            return
        }
        visits[i] = true
        var count = 0
        guard let atom = mol.getAtom(i + 1) else {
            fatalError("No atom at \(i + 1)")
        }

        for nbr in atom.getNbrAtomIterator()! {
            if nbr.getAtomicNum() != MKElements.Hydrogen.atomicNum {
                neighbor[count] = nbr.getIdx() - 1
                count += 1
            }
        }
        resnos[i] = r 
        switch atomids[i] {
        case AI_P: 
            var k = AI_O1P /* O1P */
            for j in 0..<count {
                if (bitmasks[neighbor[j]] & BitO5) != 0 {
                    atomids[neighbor[j]] = AI_O5
                    if !visits[neighbor[j]] {
                        traceNucleicChain(mol, neighbor[j], r)
                    }
                } else if (bitmasks[neighbor[j]] & BitOP) != 0 {
                    atomids[neighbor[j]] = k
                    resnos[neighbor[j]] = r
                    k = AI_O2P  /* O2P */
                }
            }
        case AI_O5:
            for j in 0..<count where (bitmasks[neighbor[j]] & BitC5) != 0 {
                atomids[neighbor[j]] = AI_C5
                if !visits[neighbor[j]] {
                    traceNucleicChain(mol, neighbor[j], r)
                }
            }
        case AI_C5:
            for j in 0..<count where (bitmasks[neighbor[j]] & BitC4) != 0 {
                atomids[neighbor[j]] = AI_C4
                if !visits[neighbor[j]] {
                    traceNucleicChain(mol, neighbor[j], r)
                }
            }
        case AI_C4:
            for j in 0..<count {
                if (bitmasks[neighbor[j]] & BitC3) != 0 {
                    atomids[neighbor[j]] = AI_C3
                    if !visits[neighbor[j]] {
                        traceNucleicChain(mol, neighbor[j], r)
                    }
                } else if (bitmasks[neighbor[j]] & BitO4) != 0 {
                    atomids[neighbor[j]] = AI_O4
                    resnos[neighbor[j]] = r
                }
            }
        case AI_C3:
            for j in 0..<count {
                if (bitmasks[neighbor[j]] & BitO3All) != 0 {
                    atomids[neighbor[j]] = AI_O3
                    if !visits[neighbor[j]] {
                        traceNucleicChain(mol, neighbor[j], r)
                    }
                } else if (bitmasks[neighbor[j]] & BitC2All) != 0 {
                    atomids[neighbor[j]] = AI_C2
                    if !visits[neighbor[j]] {
                        traceNucleicChain(mol, neighbor[j], r)
                    }
                }
            }
        case AI_O3:
            for j in 0..<count where (bitmasks[neighbor[j]] & BitP) != 0 {
                atomids[neighbor[j]] = AI_P
                if !visits[neighbor[j]] {
                    traceNucleicChain(mol, neighbor[j], r + 1)
                }
            }
        case AI_C2:
            for j in 0..<count {
                if (bitmasks[neighbor[j]] & BitC1) != 0 {
                    atomids[neighbor[j]] = AI_C1
                    resnos[neighbor[j]] = r
                } else if (bitmasks[neighbor[j]] & BitO2) != 0 {
                    atomids[neighbor[j]] = AI_O2
                    resnos[neighbor[j]] = r
                }
            }
        default:
            fatalError("Unknown atomid \(atomids[i])")
        }
    }

    /**
    * Look for atoms with atomids[i] C1 and identify their side chain.
    *
    * Sets resnos[i] and resids[i] for all identified residues.
    * (through IdentifyResidue() and AssignResidue())
    */
    //@}
    func determineNucleicSidechains(_ mol: MKMol) -> Bool {
        for i in 0..<mol.numAtoms() {
            if atomids[i] == 49 {
                let resid = identifyResidue(&NDecisionTree, mol, i, resnos[i])
                assignResidue(mol, resnos[i], chains[i], resid)
            }
        }
        return true
    }

    /**
    * Set up the chain perception to operate on the supplied molecule
    * by resizing and initializing the private data vectors.
    */
    func setupMol(_ mol: MKMol) {
        cleanupMol()
        
        let asize = mol.numAtoms()
        let bsize = mol.numBonds()
        bitmasks = [UInt16](repeating: 0, count: asize)
        visits = [Bool](repeating: false, count: asize)
        resids = [UInt8](repeating: 0, count: asize)
        flags = [UInt8](repeating: 0, count: bsize)
        hetflags = [Bool](repeating: false, count: asize)
        atomids = [Int](repeating: 0, count: asize)
        resnos = [Int](repeating: 0, count: asize)
        sernos = [Int](repeating: 0, count: asize)
        hcounts = [Int](repeating: 0, count: asize)
        chains = [Character](repeating: " ", count: asize)

        for i in 0..<asize {
            atomids[i] = -1
        }
    }
    
    /**
    * Delete all residues in @p mol
    */
    func clearResidueInformation(_ mol: MKMol) {
        if mol.numResidues() == 0 { return }
        // TODO: Make sure this works
        for residue in mol.getResidueIterator() {
            mol.deleteResidue(residue)
        }
    }
    /**
    * Clear all private data vectors
    */
    func cleanupMol() {
        bitmasks.removeAll()
        visits.removeAll()
        atomids.removeAll()
        resnos.removeAll()
        resids.removeAll()
        chains.removeAll()
        hcounts.removeAll()
        sernos.removeAll()
        hetflags.removeAll()
        flags.removeAll()
    }
    
    /**
    * Construct and add ByteCode to the @p tree for a single residue.
    * @param tree Bytecode for the residues. (OBChainsParser::PDecisionTree or OBChainsParser::NDecisionTree)
    * @param resid The residue id.
    * @param smiles The pseudo-smiles string (from OpenBabel::AminoAcids or OpenBabel::Nucleotides)
    */
    func defineMonomer(_ tree: inout ByteCode?, _ resid: Int, _ smiles: String) {
        MonoAtomCount = 0
        MonoBondCount = 0
        var prev = -1
        parseSmiles(smiles, &prev)
        for i in 0..<MonoBondCount {
            MonoBond[i].index = -1
        }
        for i in 0..<MonoAtomCount {
            MonoAtom[i].index = -1
        }
        AtomIndex = 0
        BondIndex = 0
        StackPtr = 0
        generateByteCodes(&tree, resid, 0, 0, 0)
    }

    /**
    * @param ptr Element id (from OpenBabel::ChainsAtomName)
    * @return The element number.
     */
    func identifyElement(_ ptr: [Character]) -> Int {
        
        let ch = ptr[1].uppercased()
        switch ptr[0].uppercased() {
        case " ":
            switch ch {
            case "B": return 5
            case "C": return 6
            case "D": return 1
            case "F": return 9
            case "H": return 1
            case "I": return 53
            case "K": return 19
            case "L": return 1
            case "N": return 7
            case "O": return 8
            case "P": return 15
            case "S": return 16
            case "U": return 92
            case "V": return 23
            case "W": return 74
            case "Y": return 39
            default: return 0
            }
        case "A":
            switch ch {
            case "C": return 89
            case "G": return 47
            case "L": return 13
            case "M": return 95
            case "R": return 18
            case "S": return 33
            case "T": return 85
            case "U": return 79
            default: return 0
            }
        case "B":
            switch ch {
            case "A": return 56
            case "E": return 4
            case "I": return 83
            case "K": return 97
            case "R": return 35
            case " ": return 5
            default: return 0
            }
        case "C":
            switch ch {
            case "A": return 20
            case "D": return 48
            case "E": return 58
            case "F": return 98
            case "L": return 17
            case "M": return 96
            case "O": return 27
            case "R": return 24
            case "S": return 55
            case "U": return 29
            case " ": return 6
            default: return 0
            }
            
        case "D":
            if ch == "Y" {
                return 66
            } else if ch == " " {
                return 1
            }
        case "E":
            if ch == "R" {
                return 68
            } else if ch == "S" {
                return 99
            } else if ch == "U" {
                return 63
            }
        case "F":
            if ch == "E" {
                return 26
            } else if ch == "M" {
                return 100
            } else if ch == "R" {
                return 87
            } else if ch == "F" {
                return 9
            }
        case "G":
            if ch == "A" {
                return 31
            } else if ch == "D" {
                return 64
            } else if ch == "E" {
                return 32
            }
        case "H":
            if ch == "E" {
                return 2
            } else if ch == "F" {
                return 72
            } else if ch == "G" {
                return 80
            } else if ch == "O" {
                return 67
            } else if ch == " " {
                return 1
            }
        case "I":
            if ch == "N" {
                return 49
            } else if ch == "R" {
                return 77
            } else if ch == " " {
                return 53
            }
        case "K":
            if ch == "R" {
                return 36
            } else if ch == " " {
                return 19
            }
        case "L":
            if ch == "A" {
                return 57
            } else if ch == "I" {
                return 3
            } else if ch == "R" || ch == "W" {
                return 103
            } else if ch == "U" {
                return 71
            } else if ch == " " {
                return 1
            }
        case "M":
            if ch == "D" {
                return 101
            } else if ch == "G" {
                return 12
            } else if ch == "N" {
                return 25
            } else if ch == "O" {
                return 42
            }
        case "N":
            if ch == "A" {
                return 11
            } else if ch == "B" {
                return 41
            } else if ch == "D" {
                return 60
            } else if ch == "E" {
                return 10
            } else if ch == "I" {
                return 28
            } else if ch == "O" {
                return 102
            } else if ch == "P" {
                return 93
            } else if ch == " " {
                return 7
            }
        case "O":
            if ch == "S" {
                return 76
            } else if ch == " " {
                return 8
            }
        case "P":
            if ch == "A" {
                return 91
            } else if ch == "B" {
                return 82
            } else if ch == "D" {
                return 46
            } else if ch == "M" {
                return 61
            } else if ch == "O" {
                return 84
            } else if ch == "R" {
                return 59
            } else if ch == "T" {
                return 78
            } else if ch == "U" {
                return 94
            } else if ch == " " {
                return 15
            }
        case "R":
            if ch == "A" {
                return 88
            } else if ch == "B" {
                return 37
            } else if ch == "E" {
                return 75
            } else if ch == "H" {
                return 45
            } else if ch == "N" {
                return 86
            } else if ch == "U" {
                return 44
            }
        case "S":
            if ch == "B" {
                return 51
            } else if ch == "C" {
                return 21
            } else if ch == "E" {
                return 34
            } else if ch == "I" {
                return 14
            } else if ch == "M" {
                return 62
            } else if ch == "N" {
                return 50
            } else if ch == "R" {
                return 38
            } else if ch == " " {
                return 16
            }
        case "T":
            if ch == "A" {
                return 73
            } else if ch == "B" {
                return 65
            } else if ch == "C" {
                return 43
            } else if ch == "E" {
                return 52
            } else if ch == "H" {
                return 90
            } else if ch == "I" {
                return 22
            } else if ch == "L" {
                return 81
            } else if ch == "M" {
                return 69
            }
        case "U":
            if ch == " " {
                return 92
            }
        case "V":
            if ch == " " {
                return 23
            }
        case "W":
            if ch == " " {
                return 74
            }
        case "X":
            if ch == "E" {
                return 54
            }
        case "Y":
            if ch == "B" {
                return 70
            } else if ch == " " {
                return 39
            }
        case "Z":
            if ch == "N" {
                return 30
            } else if ch == "R" {
                return 40
            }
        default:
            break
        }
        
        if ptr[0].isNumber {
            if ptr[0] == "H" || ptr[0] == "D" {
                return 1
            }
        } 
        return 0 
    }
    
    /**
     * Parse a pseudo smiles from OpenBabel::AminoAcids or OpenBabel::Nucleotides.
     * @param smiles The pseudo-smiles string.
     * @param prev The previous position (used for recursing, use -1 to start).
     */
    @discardableResult
    private func parseSmiles(_ smiles: String, _ prev: inout Int) -> String {
        
        var ptr = Iterator<Character>([Character](smiles))
        
        var type: Int = 0
        var ch: Character?
        ch = ptr.next()
        while ch != nil {
            switch ch {
            case "-": 
                type = BF_SINGLE
            case "=": 
                type = BF_DOUBLE
            case "#": 
                type = BF_TRIPLE
            case "^": 
                type = BF_SINGLE|BF_AROMATIC
            case "~": 
                type = BF_DOUBLE|BF_AROMATIC

            case ")": 
                return String(ptr)
            case ".": 
                prev = -1
            case "(": 
                ptr = Iterator<Character>([Character](parseSmiles(String(ptr.constructToEnd()), &prev)))
            default: 
                guard var atomid = ch!.wholeNumberValue else {
                    fatalError("Invalid character in smiles")
                }
                while ptr[0]!.isNumber {
                    atomid = atomid * 10 + (ptr[0]?.wholeNumberValue!)!
                    ptr.ignore()
                }
                var next = 0
                while next < MonoAtomCount {
                    if MonoAtom[next].atomid == atomid {
                        break
                    }
                    next += 1
                }
                if next == MonoAtomCount {
                    let name = ChainsAtomName[atomid]
                    MonoAtom[next].elem = identifyElement(name.map(Character.init))
                    MonoAtom[next].atomid = atomid
                    MonoAtom[next].bcount = 0
                    MonoAtomCount += 1
                }
                if prev != -1 {
                    MonoBond[MonoBondCount].flag = type
                    MonoBond[MonoBondCount].src = prev
                    MonoBond[MonoBondCount].dst = next
                    MonoBondCount += 1
                    MonoAtom[prev].bcount += 1
                    MonoAtom[next].bcount += 1
                }
                prev = next
            }
            ch = ptr.next()
        }
        ptr.unget()
        return String(ptr.constructToEnd())
    }
    
    //   /**
    //    * Debugging function.
    //    */
    //   void DumpState();
    
}
