//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/21/23.
//

import Foundation

//Daylight SMARTS parser


private let ATOMPOOL =          1
private let BONDPOOL =          1

private let AE_ANDHI =          1
private let AE_ANDLO =          2
private let AE_OR =             3
private let AE_RECUR =          4
private let AE_NOT =            5
private let AE_TRUE =           6
private let AE_FALSE =          7
private let AE_AROMATIC =       8
private let AE_ALIPHATIC =      9
private let AE_CYCLIC =         10
private let AE_ACYCLIC =        11
private let AE_MASS =           12
private let AE_ELEM =           13
private let AE_AROMELEM =       14
private let AE_ALIPHELEM =      15
private let AE_HCOUNT =         16
private let AE_CHARGE =         17
private let AE_CONNECT =        18
private let AE_DEGREE =         19
private let AE_IMPLICIT =       20
private let AE_RINGS =          21
private let AE_SIZE =           22
private let AE_VALENCE =        23
private let AE_CHIRAL =         24
private let AE_HYB =            25
private let AE_RINGCONNECT =    26

private let AL_CLOCKWISE =      1
private let AL_ANTICLOCKWISE =  2
private let AL_UNSPECIFIED =    0
/* Each BondExpr node is identified by an integer type field
 selected from the list of BE_ codes below.  BE_ANDHI,
 BE_ANDLO and BE_OR are binary nodes of type BondExpr.bin,
 BE_NOT is unary node of type BondExpr.mon, and the remaining
 code are all leaf nodes.  */
private let BE_ANDHI =          1
private let BE_ANDLO =          2
private let BE_OR =             3
private let BE_NOT =            4
private let BE_ANY =            5
private let BE_DEFAULT =        6
private let BE_SINGLE =         7
private let BE_DOUBLE =         8
private let BE_TRIPLE =         9
private let BE_QUAD =           10
private let BE_AROM =           11
private let BE_RING =           12
private let BE_UP =             13
private let BE_DOWN =           14
private let BE_UPUNSPEC =       15
private let BE_DOWNUNSPEC =     16

private let SmartsImplicitRef = -9999 // Used as a placeholder when recording atom nbrs for chiral atoms

//! \union _AtomExpr parsmart.h <openbabel/parsmart.h>
//! \brief An internal (SMARTS parser) atomic expression

protocol _AtomExprProtocol {
    var type: Int { get set }
}

struct _AtomExprLeaf: _AtomExprProtocol {
    var type: Int
    var value: Int
}

struct _AtomExprRecur: _AtomExprProtocol {
    var type: Int
    var recur: Pattern?
}

struct _AtomExprMon: _AtomExprProtocol {
    var type: Int
    var arg: _AtomExpr
}

struct _AtomExprBin: _AtomExprProtocol {
    var type: Int
    var lft: _AtomExpr
    var rgt: _AtomExpr
}

indirect enum _AtomExpr {
    case leaf(_AtomExprLeaf)
    case recur(_AtomExprRecur)
    case mon(_AtomExprMon)
    case bin(_AtomExprBin)
}

//! \union _BondExpr parsmart.h <openbabel/parsmart.h>
//! \brief An internal (SMARTS parser) bond expression

protocol _BondExprProtocol {
    var type: Int { get set }
}

struct _BondExprLeaf: _BondExprProtocol {
    var type: Int
}

struct _BondExprMon: _BondExprProtocol {
    var type: Int
    var arg: _BondExpr
}

struct _BondExprBin: _BondExprProtocol {
    var type: Int
    var lft: _BondExpr
    var rgt: _BondExpr
}

indirect enum _BondExpr {
    case leaf(_BondExprLeaf)
    case mon(_BondExprMon)
    case bin(_BondExprBin)
}

//! \struct BondSpec parsmart.h <openbabel/parsmart.h>
//! \brief An internal (SMARTS parser) bond specification
struct BondSpec {
    var expr: _BondExpr
    var src, dst: Int
    var visit: Int?
    var grow: Bool?
}

//! \struct AtomSpec parsmart.h <openbabel/parsmart.h>
//! \brief An internal (SMARTS parser) atom specification
struct AtomSpec {
    var expr: _AtomExpr
    var visit: Int?
    var part: Int
    var chiral_flag: Int?
    var vb: Int
    var nbrs: [Int]?
}

//! \struct Pattern parsmart.h <openbabel/parsmart.h>
//! \brief A SMARTS parser internal pattern
struct Pattern {
    var aalloc, acount: Int
    var balloc, bcount: Int
    var isChiral: Bool
    var atom: [AtomSpec]
    var bond: [BondSpec]
    var parts: Int
    var hasExplicitH: Bool
}

//! \struct ParseState parsmart.h <openbabel/parsmart.h>
//! \brief A SMARTS parser internal state
struct ParseState {
    var closord: [_BondExpr]
    var closure: [Int]
    var closindex: Int
}


/*=============================*/
/*  Standard Utility Routines  */
/*=============================*/


//! \brief SMARTS (SMiles ARbitrary Target Specification) substructure searching
class MKSmartsPattern {
    
    enum MatchType {
        case All
        case Single
        case AllUnique
    }
    
    private var _mlist: [[Int]] = []
    private var _pat: Pattern?
    private var _str: String = ""
    // OBSmartsPrivate                *_d;        //!< Internal data storage for future expansion
    // OB_DEPRECATED std::vector<bool> _growbond; //!< \deprecated (Not used)

    private var _buffer: String = "" 
    private var LexPtr: Int = 0
    private var MainPtr: Int = 0
    
    
    //! \name Initialization Methods

    //**********************************
    //********Pattern Matching**********
    //**********************************

    init() {
        _pat = nil
    }
    
    func initialize(_ pattern: String) -> Bool {
        _pat = self.parseSMARTSRecord(pattern)
        _str = pattern
        return _pat != nil 
    }
    
    //! \name Pattern Properties
    func getSMARTS() -> String {
        return self._str
    }
    
    //! \return If the SMARTS pattern is an empty expression (e.g., invalid)
    func isEmpty() -> Bool {
        return _pat == nil
    }
    
    //! \return If the SMARTS pattern is a valid expression
    func isValid() -> Bool {
        return _pat != nil
    }
    
    func numAtoms() -> Int {
        return _pat?.acount ?? 0
    }
    
    //! \return the number of bonds in the SMARTS pattern
    func numBonds() -> Int {
        return self._pat?.bcount ?? 0
    }
    
    //! Access the bond @p idx in the internal pattern
    //! \param src The index of the beginning atom
    //! \param dst The index of the end atom
    //! \param ord The bond order of this bond
    //! \param idx The index of the bond in the SMARTS pattern
    func getBond(_ src: Int, _ dst: Int, _ ord: Int, _ idx: Int) {
        
    }
    
    func getAtomicNum(_ idx: Int) -> Int {
        
    }
    
    func getCharge(_ idx: Int) -> Int {
        
    }
    
    //! \name Matching methods (SMARTS on a specific OBMol)
    //@{
    //! Perform SMARTS matching for the pattern specified using Init().
    //! \param mol The molecule to use for matching
    //! \param single Whether only a single match is required (faster). Default is false.
    //! \return Whether matches occurred
    func match(_ mol: MKMol, _ single: Bool = false) -> Bool {
        
    }
    
    //! \name Matching methods (SMARTS on a specific OBMol)
    //@{
    //! Perform SMARTS matching for the pattern specified using Init().
    //! This version is (more) thread safe.
    //! \param mol The molecule to use for matching
    //! \param mlist The resulting match list
    //! \param mtype The match type to use. Default is All.
    //! \return Whether matches occurred
    func match(_ mol: MKMol, _ mlist: [[Int]], _ mtype: MatchType = .All) -> Bool {
        
    }
    
    //! \name Matching methods (SMARTS on a specific OBMol)
    //@{
    //! Thread safe check for any SMARTS match
    //! \param mol The molecule to use for matching
    //! \return Whether there exists any match
    func hasMatch(_ mol: MKMol) -> Bool {
        
    }
    
    func restrictedMatch(_ mol: MKMol, _ pairs: [Pair<Int, Int>], _ single: Bool = false) -> Bool {}
    func restrictedMatch(_ mol: MKMol, _ bv: MKBitVec, _ single: Bool = false) -> Bool {}
    
    //! \return the number of non-unique SMARTS matches
    //! To get the number of unique SMARTS matches, query GetUMapList()->size()
    func numMatches() -> Int {
        return self._mlist.count
    }
    
    //! \return the entire list of non-unique matches for this pattern
    //! \see GetUMapList()
    func getMapList() -> [[Int]] { return self._mlist }
    
    //! \return the entire list of unique matches for this pattern
        /**
           A unique match is defined as one which does not cover the
           identical atoms that a previous match has covered.

           For instance, the pattern [OD1]~C~[OD1] describes a
           carboxylate group. This pattern will match both atom number
           permutations of the carboxylate, and if GetMapList() is called, both
           matches will be returned. If GetUMapList() is called only unique
           matches of the pattern will be returned.
        **/
    func getUMapList() -> [[Int]] { }
    
    /*=========================*/
    /*  SMARTS Syntax Parsing  */
    /*=========================*/

    private func parseSMARTSPattern() -> Pattern {}
    private func parseSMARTSPart(_ pat: Pattern, _ idx: Int) -> Pattern {}
    // private func SMARTSError(_ pat: Pattern) -> Pattern {}
    private func parseSMARTSError(_ pat: Pattern, _ expr: _BondExpr) -> Pattern {}
    
    private func parseSimpleAtomPrimitive() -> _AtomExpr {}
    private func parseComplexAtomPrimitive() -> _AtomExpr {}
    private func parseAtomExpr(_ level: Int) -> _AtomExpr {}
    
    private func parseBondPrimitive() -> _BondExpr {}
    private func parseBondExpr(_ level: Int) -> _BondExpr {}
    
    private func parseSMARTSString(_ ptr: String) -> Pattern {
        
    }

    private func parseSMARTSRecord(_ ptr: String) -> Pattern {
        // Trim whitespace
        return self.parseSMARTSString(ptr.trimmingCharacters(in: .whitespacesAndNewlines))
    }
    
    //! \return the vector binding of the atom @p idx in the internal pattern
    func getVectorBinding(_ idx: Int) -> Int {
        return _pat?.atom[idx].vb ?? 0
    }
    
    private func SMARTSParser(_ pat: Pattern, _ stat: ParseState, _ prev: Int, _ part: Int) -> Pattern {}


}


/*================================*/
/*  Atom Expression Manipulation  */
/*================================*/

func buildAtomPred(_ type: Int) -> _AtomExpr {
    let res = _AtomExprLeaf(type: type, value: 0)
    return .leaf(res)
}

func buildAtomLeaf(_ type: Int, _ value: Int) -> _AtomExpr {
    let res = _AtomExprLeaf(type: type, value: value)
    return .leaf(res)
}

func buildAtomNot(_ expr: _AtomExpr) -> _AtomExpr {
    let res = _AtomExprMon(type: AE_NOT, arg: expr)
    return .mon(res)
}

func buildAtomBinary(_ type: Int, _ lhs: _AtomExpr, _ rhs: _AtomExpr) -> _AtomExpr {
    let res = _AtomExprBin(type: type, lft: lhs, rgt: rhs)
    return .bin(res)
}

func buildAtomRecurs(_ pat: Pattern) -> _AtomExpr {
    let res = _AtomExprRecur(type: AE_RECUR, recur: pat)
    return .recur(res)
}

func generateElement(_ elem: Int) -> _AtomExpr {
    return buildAtomLeaf(AE_ELEM, elem)
}

func generateAromElem(_ elem: Int, _ flag: Bool) -> _AtomExpr {
    // MARK: sneakily add aliphatic elements here as well (for now)
    return flag ? buildAtomLeaf(AE_AROMELEM, elem) : buildAtomLeaf(AE_ALIPHELEM, elem)
}

/*================================*/
/*  Bond Expression Manipulation  */
/*================================*/

/**
   * Check if two BondExpr objects are the same. This is used for ring closures
   * to identify invalid SMARTS like:
   *
   *   C-1CCCCC#1
   *   C=1CCCCC:1
   *
   * However, the SMARTS below are valid and the bond expression next to the the
   * second closure digit is used.
   *
   *   C1CCCCC#1
   *   C1CCCCC=1
   */
func equivalentBonExpr(_ lhs: _BondExpr, _ rhs: _BondExpr) -> Bool {
    if (lhs as! _BondExprProtocol ).type != (rhs as! _BondExprProtocol ).type {
        return false
    }
    
    switch lhs {
    case .leaf(let _expr):
        break
    case .mon(let expr1):
        switch rhs {
        case .leaf(let _expr):
            break
        case .mon(let expr2):
            switch expr1.type {
            case BE_NOT:
                return equivalentBonExpr(expr1.arg, expr2.arg)
            default: break
            }
        case .bin(_):
            break
        }
    case .bin(let expr1):
        switch rhs {
        case .leaf(let _expr):
            break
        case .mon(_):
            break
        case .bin(let expr2):
            switch expr1.type {
            case BE_ANDHI, BE_ANDLO, BE_OR:
                return equivalentBonExpr(expr1.lft, expr2.lft) &&
                       equivalentBonExpr(expr1.rgt, expr2.rgt)
            default: break
            }
        }
    }
    return true
}


func buildBondLeaf(_ type: Int) -> _BondExpr {
    let res = _BondExprLeaf(type: type)
    return .leaf(res)
}

func buildBondNot(_ expr: _BondExpr) -> _BondExpr {
    let res = _BondExprMon(type: BE_NOT, arg: expr)
    return .mon(res)
}

func buildBondBin(_ type: Int, _ lhs: _BondExpr, _ rhs: _BondExpr) -> _BondExpr {
    let res = _BondExprBin(type: type, lft: lhs, rgt: rhs)
    return .bin(res)
}

func generateDefaultBond() -> _BondExpr {
    return buildBondLeaf(BE_DEFAULT)
}

  /*===============================*/
  /*  SMARTS Pattern Manipulation  */
  /*===============================*/

func createAtom(_ pat: inout Pattern, _ expr: _AtomExpr, _ part: Int, _ vb: Int) -> Int {
    let tmp = AtomSpec(expr: expr, part: part, vb: vb)
    pat.atom.append(tmp)
    pat.acount += 1
    return pat.acount
}

func createBond(_ pat: inout Pattern, _ expr: _BondExpr, _ src: Int, _ dst: Int) -> Int {
    let tmp = BondSpec(expr: expr, src: src, dst: dst)
    pat.bond.append(tmp)
    pat.bcount += 1
    return pat.bcount
}



class MKSmartsMatcher {
	
    //recursive smarts cache
    var rscache: [Pair<Pattern, [Bool]>] = []
    // list of fragment patterns (e.g., (*).(*)
    var fragments: [Pattern] = []

    init() {

    }

    private func evalAtomExpr(_ expr: _AtomExpr, _ atom: MKAtom) -> Bool {
        return false 
    }

    private func evalBondExpr(_ expr: _BondExpr, _ bond: MKBond) -> Bool {
        return false 
    }

    private func setupAtomMatchTable(_ ttab: [[Bool]], _ pat: Pattern, _ mol: MKMol) {
        return  
    }

    private func fastSingleMatch(_ mol: MKMol, _ pat: Pattern, _ mlist: [[Int]]) {

    }

    public func match(_ mol: MKMol, _ pattern: MKSmartsPattern, _ mlist: [[Int]], _ single: Bool = false) -> Bool {
        return false
    }

}


  //! \brief Internal class: performs fast, exhaustive matching used to find
  //! just a single match in match() using recursion and explicit stack handling.
class MKSSMatch {

    private var _uatoms: Bool
    private var _mol: MKMol 
    private var _pat: Pattern
    private var _map: [Int]

    init(mol: MKMol, pat: Pattern) {
       
    }

    func match(_ v: [[Int]], _ bidx: Int = -1) {

    }
}

public func smartsLexReplace(_ s: String, _ p: [Pair<String, String>]) {

}
