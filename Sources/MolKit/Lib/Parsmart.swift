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
    var aalloc: Int = 0
    var acount: Int = 0
    var balloc: Int = 0
    var bcount: Int = 0
    var isChiral: Bool = false
    var atom: [AtomSpec] = []
    var bond: [BondSpec] = []
    var parts: Int = 1
    var hasExplicitH: Bool = false
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
    private var LexPtr: LexicalParser = LexicalParser()
    
    
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
        return 0
    }
    
    func getCharge(_ idx: Int) -> Int {
        return 0
    }
    
    //! \name Matching methods (SMARTS on a specific OBMol)
    //@{
    //! Perform SMARTS matching for the pattern specified using Init().
    //! \param mol The molecule to use for matching
    //! \param single Whether only a single match is required (faster). Default is false.
    //! \return Whether matches occurred
    func match(_ mol: MKMol, _ single: Bool = false) -> Bool {
        return false
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
        return false
    }
    
    //! \name Matching methods (SMARTS on a specific OBMol)
    //@{
    //! Thread safe check for any SMARTS match
    //! \param mol The molecule to use for matching
    //! \return Whether there exists any match
    func hasMatch(_ mol: MKMol) -> Bool {
        return false
    }
    
    func restrictedMatch(_ mol: MKMol, _ pairs: [Pair<Int, Int>], _ single: Bool = false) -> Bool {
        return false
    }
    
    func restrictedMatch(_ mol: MKMol, _ bv: MKBitVec, _ single: Bool = false) -> Bool {
        return false
    }
    
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
    func getUMapList() -> [[Int]] {
        if self._mlist.isEmpty || self._mlist.count == 1 {
            return self._mlist
        }
        
        var ok: Bool
        var bv: MKBitVec = MKBitVec()
        var vbv: [MKBitVec] = []
        var mlist: [[Int]] = []
        
        for i in self._mlist {
            ok = true
            bv.clear()
            bv.fromVecInt(i)
            for j in 0..<vbv.count {
                if vbv[j] == bv && ok {
                    ok = false
                }
            }
            if ok {
                mlist.append(i)
                vbv.append(bv)
            }
        }
        
        self._mlist = mlist
        return self._mlist
    }
    
    /*=========================*/
    /*  SMARTS Syntax Parsing  */
    /*=========================*/

    private func parseSMARTSPattern() -> Pattern? {
        
        var result: Pattern = Pattern()
        
        while self.LexPtr == "(" {
            
            self.LexPtr.inc()
            
//            TODO: reconcile passing references and values with return types
            guard var result = self.parseSMARTSPart(&result, result.parts) else { return nil }
            result.parts += 1
            
            if self.LexPtr != ")" {
                print("SMARTS Error")
//                MARK: install throw here
            }
            self.LexPtr.inc()
            
            if self.LexPtr.empty() || self.LexPtr == ")" {
                return result
            }
            
            if self.LexPtr != "." {
                print("SMARTS Error")
//                MARK: install throw here
            }
            // Here's where we'd handle fragments
            //        cerr << " conjunction " << LexPtr[0] << endl;
            self.LexPtr.inc()
        }
        
        return self.parseSMARTSPart(&result, 0)
    }
    
    private func parseSMARTSPart(_ pat: inout Pattern, _ idx: Int) -> Pattern? {
        return nil
    }
    
    // private func SMARTSError(_ pat: Pattern) -> Pattern {}
    // private func parseSMARTSError(_ pat: Pattern, _ expr: _BondExpr) -> Pattern {}
    
    private func parseSimpleAtomPrimitive() -> _AtomExpr? {
        switch self.LexPtr.next() {
        case "*":
            return buildAtomPred(AE_TRUE)
        case "A":
            return buildAtomPred(AE_ALIPHATIC)
        case "B":
            if self.LexPtr == "r" {
                self.LexPtr.inc()
                return generateElement(35)
            }
            return generateElement(5)
        case "C":
            if self.LexPtr == "l" {
                self.LexPtr.inc()
                return generateElement(17)
            }
            return generateAromElem(6, false)
        case "F":
            return generateElement(9)
        case "I":
            return generateElement(53)
        case "N":
            return generateAromElem(7, false)
        case "O":
            return generateAromElem(8, false)
        case "P":
            return generateAromElem(15, false)
        case "S":
            return generateAromElem(16, false)
        case "a":
            if self.LexPtr == "s" {
                self.LexPtr.inc()
                return generateAromElem(33, true)
            }
            return buildAtomPred(AE_AROMATIC)
        case "c":
            return generateAromElem(6, true)
        case "n":
            return generateAromElem(7, true)
        case "o":
            return generateAromElem(8, true)
        case "p":
            return generateAromElem(15, true)
        case "s":
            if self.LexPtr == "e" {
                self.LexPtr.inc()
                return generateAromElem(34, true)
            }
            return generateAromElem(16, true)
        default: break
        }

        self.LexPtr.dec()
        return nil
    }
    
    private func parseComplexAtomPrimitive() -> _AtomExpr? {
        var pat: Pattern? 
        var index: Int

        switch self.LexPtr.next() {
        case "#":
            if !self.LexPtr.cur().isNumber { return nil }
            index = 0
            while self.LexPtr.cur().isNumber {
                index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
            }
            if index > 255 {
                self.LexPtr.dec()
                return nil
            }
            return generateElement(index)
        case "$":
            if self.LexPtr != "(" { return nil }
            self.LexPtr.inc()
            pat = parseSMARTSPattern()
            if self.LexPtr != ")" { return nil }
            self.LexPtr.inc()
            if pat != nil {
                return buildAtomRecurs(pat!)
            } else { return nil }
        case "*":
            return buildAtomPred(AE_TRUE)
        case "+":
            if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
            } else {
                index = 1
                while self.LexPtr.cur() == "+" {
                    self.LexPtr.inc()
                    index += 1
                }
            }
            return buildAtomLeaf(AE_CHARGE, index)
        case "-":
            if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
            } else {
                index = 1
                while self.LexPtr.cur() == "-" {
                    self.LexPtr.inc()
                    index += 1
                }
            }
            return buildAtomLeaf(AE_CHARGE, -index)
        case "@":
            if self.LexPtr.cur() == "?" {
                self.LexPtr.inc()
                return buildAtomLeaf(AE_CHIRAL,AL_UNSPECIFIED)
            } else if self.LexPtr.cur() != "@" {
                return buildAtomLeaf(AE_CHIRAL,AL_ANTICLOCKWISE)
            } else {
                self.LexPtr.inc()
                return buildAtomLeaf(AE_CHIRAL,AL_CLOCKWISE)
            }
        case "^":
            if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
                return buildAtomLeaf(AE_HYB, index)
            } else {
                return buildAtomLeaf(AE_HYB, 1)
            }
        case "0", "1", "2", "3", "4", "5", "6", "7", "8", "9":
            index = self.LexPtr.prev().wholeNumberValue!
            while self.LexPtr.cur().isNumber {
                index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
            }
            return buildAtomLeaf(AE_MASS, index)
        case "A":
            switch self.LexPtr.next() {
            case "c":  return  generateElement(89)
            case "g":  return  generateElement(47)
            case "l":  return  generateElement(13)
            case "m":  return  generateElement(95)
            case "r":  return  generateElement(18)
            case "s":  return  generateElement(33)
            case "t":  return  generateElement(85)
            case "u":  return  generateElement(79)
            default: break
            }
            self.LexPtr.dec()
            return buildAtomPred(AE_ALIPHATIC)
        case "B":
            switch self.LexPtr.next() {
            case "a":  return  generateElement(56)
            case "e":  return  generateElement(4)
            case "i":  return  generateElement(83)
            case "k":  return  generateElement(97)
            case "r":  return  generateElement(35)
            default: break
            }
            self.LexPtr.dec()
            return generateElement(5)
        case "C":
            switch self.LexPtr.next() {
            case "a":  return  generateElement(20)
            case "d":  return  generateElement(48)
            case "e":  return  generateElement(58)
            case "f":  return  generateElement(98)
            case "l":  return  generateElement(17)
            case "m":  return  generateElement(96)
            case "o":  return  generateElement(27)
            case "r":  return  generateElement(24)
            case "s":  return  generateElement(55)
            case "u":  return  generateElement(29)
            default: break
            }
            self.LexPtr.dec()
            return generateAromElem(6, false)
        case "D":
            if self.LexPtr.cur() == "y" {
                self.LexPtr.inc()
                return generateElement(66)
            } else if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
                return buildAtomLeaf(AE_DEGREE, index)
            }
            return buildAtomLeaf(AE_DEGREE, 1)
        case "E":
            if self.LexPtr.cur() == "r" {
                self.LexPtr.inc()
                return generateElement(68)
            } else if self.LexPtr.cur() == "s" {
                self.LexPtr.inc()
                return generateElement(99)
            } else if self.LexPtr.cur() == "u" {
                self.LexPtr.inc()
                return generateElement(63)
            }
            break
        case "F":
            if self.LexPtr.cur() == "e" {
                self.LexPtr.inc()
                return generateElement(26)
            } else if self.LexPtr.cur() == "m" {
                self.LexPtr.inc()
                return generateElement(100)
            } else if self.LexPtr.cur() == "r" {
                self.LexPtr.inc()
                return generateElement(87)
            }
            return generateElement(9)
        case "G":
            if self.LexPtr.cur() == "a" {
                self.LexPtr.inc()
                return generateElement(31)
            } else if self.LexPtr.cur() == "d" {
                self.LexPtr.inc()
                return generateElement(64)
            } else if self.LexPtr.cur() == "e" {
                self.LexPtr.inc()
                return generateElement(32)
            }
            break
        case "H":
            if self.LexPtr.cur() == "e" {
                self.LexPtr.inc()
                return generateElement(2)
            } else if self.LexPtr.cur() == "f" {
                self.LexPtr.inc()
                return generateElement(72)
            } else if self.LexPtr.cur() == "g" {
                self.LexPtr.inc()
                return generateElement(80)
            } else if self.LexPtr.cur() == "o" {
                self.LexPtr.inc()
                return generateElement(67)
            } else if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
                return buildAtomLeaf(AE_HCOUNT, index)
            }
            return buildAtomLeaf(AE_HCOUNT, 1)
        case "I":
            if self.LexPtr.cur() == "n" {
                self.LexPtr.inc()
                return generateElement(49)
            } else if self.LexPtr.cur() == "r" {
                self.LexPtr.inc()
                return generateElement(77)
            }
            return generateElement(53)
        case "K":
            if self.LexPtr.cur() == "r" {
                self.LexPtr.inc()
                return generateElement(36)
            }
            return generateElement(19)
        case "L":
            if self.LexPtr.cur() == "a" {
                self.LexPtr.inc()
                return generateElement(57)
            } else if self.LexPtr.cur() == "i" {
                self.LexPtr.inc()
                return generateElement(3)
            } else if self.LexPtr.cur() == "r" {
                self.LexPtr.inc()
                return generateElement(103)
            } else if self.LexPtr.cur() == "u" {
                self.LexPtr.inc()
                return generateElement(71)
            }
            break
        case "M":
            if self.LexPtr.cur() == "d" {
                self.LexPtr.inc()
                return generateElement(101)
            } else if self.LexPtr.cur() == "g" {
                self.LexPtr.inc()
                return generateElement(12)
            } else if self.LexPtr.cur() == "n" {
                self.LexPtr.inc()
                return generateElement(25)
            } else if self.LexPtr.cur() == "o" {
                self.LexPtr.inc()
                return generateElement(42)
            }
            break
        case "N":
            switch self.LexPtr.next() {
            case "a":  return  generateElement(11)
            case "b":  return  generateElement(41)
            case "d":  return  generateElement(60)
            case "e":  return  generateElement(10)
            case "i":  return  generateElement(28)
            case "o":  return  generateElement(102)
            case "p":  return  generateElement(93)
            default: break
            }
            self.LexPtr.dec()
            return generateAromElem(7, false)
        case "O":
            if self.LexPtr.cur() == "s" {
                self.LexPtr.inc()
                return generateElement(76)
            }
            return generateAromElem(8, false)
        case "P":
            switch self.LexPtr.next() {
            case "a":  return  generateElement(91)
            case "b":  return  generateElement(82)
            case "d":  return  generateElement(46)
            case "m":  return  generateElement(61)
            case "o":  return  generateElement(84)
            case "r":  return  generateElement(59)
            case "t":  return  generateElement(78)
            case "u":  return  generateElement(94)
            default: break
            }
            self.LexPtr.dec()
            return generateElement(15)
        case "R":
            switch self.LexPtr.next() {
            case "a":  return  generateElement(88)
            case "b":  return  generateElement(37)
            case "e":  return  generateElement(75)
            case "h":  return  generateElement(45)
            case "n":  return  generateElement(86)
            case "u":  return  generateElement(44)
            default: break
            }
            self.LexPtr.dec()
            if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
                if index == 0 {
                    return buildAtomPred(AE_CYCLIC)
                }
                return buildAtomLeaf(AE_RINGS, index)
            }
            return buildAtomPred(AE_CYCLIC)
        case "S":
            switch self.LexPtr.next() {
            case "b":  return  generateElement(51)
            case "c":  return  generateElement(21)
            case "e":  return  generateElement(34)
            case "i":  return  generateElement(14)
            case "m":  return  generateElement(62)
            case "n":  return  generateElement(50)
            case "r":  return  generateElement(38)
            default: break
            }
            self.LexPtr.dec()
            return generateAromElem(16, false)
        case "T":
            switch self.LexPtr.next() {
            case "a":  return  generateElement(73)
            case "b":  return  generateElement(65)
            case "c":  return  generateElement(43)
            case "e":  return  generateElement(52)
            case "h":  return  generateElement(90)
            case "i":  return  generateElement(22)
            case "l":  return  generateElement(81)
            case "m":  return  generateElement(69)
            default: break
            }
            self.LexPtr.dec()
            break
        case "U": return  generateElement(92)
        case "V": return  generateElement(23)
        case "W": return  generateElement(74)
        case "X": 
            if self.LexPtr.cur() == "e" {
                self.LexPtr.inc()
                return generateElement(54)
            } else if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
                if index == 0 { // default to 1 (if no number present)
                    index = 1 
                }
                return buildAtomLeaf(AE_CONNECT, index)
            }
            return buildAtomLeaf(AE_CONNECT, 1)
        case "Y": 
            if self.LexPtr.cur() == "b" {
                self.LexPtr.inc()
                return generateElement(70)
            }
            return  generateElement(39)
        case "Z":
            if self.LexPtr.cur() == "n" {
                self.LexPtr.inc()
                return generateElement(30)
            } else if self.LexPtr.cur() == "r" {
                self.LexPtr.inc()
                return generateElement(40)
            }
            break
        case "a":
            if self.LexPtr.cur() == "s" {
                self.LexPtr.inc()
                return generateAromElem(33, true)
            }
            return buildAtomPred(AE_AROMATIC)
        case "c": return generateAromElem(6, true)
        case "h":
            if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
            } else {
                index = 1 
            }
            return buildAtomLeaf(AE_IMPLICIT, index)
        case "n": return generateAromElem(7, true)
        case "o": return generateAromElem(8, true)
        case "p": return generateAromElem(15, true)
        case "r":
            if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
                if index == 0 {
                    return buildAtomPred(AE_CYCLIC)
                }
                return buildAtomLeaf(AE_SIZE, index)
            }
            return buildAtomPred(AE_CYCLIC)
        case "s":
            if self.LexPtr.cur() == "e" {
                self.LexPtr.inc()
                return generateAromElem(34, true)
            }
            return generateAromElem(16, true)
        case "v":
            if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
                return buildAtomLeaf(AE_VALENCE, index)
            }
            return buildAtomLeaf(AE_VALENCE, 1)
        case "x":
            if self.LexPtr.cur().isNumber {
                index = 0
                while self.LexPtr.cur().isNumber {
                    index = index * 10 + (self.LexPtr.next()?.wholeNumberValue)!
                }
                return buildAtomLeaf(AE_RINGCONNECT,index)
            }
            return buildAtomPred(AE_CYCLIC)
        default: break
        }
        self.LexPtr.dec()
        return nil
    }
    
    private func parseAtomExpr(_ level: Int) -> _AtomExpr? {
        return nil
    }
    
    private func parseBondPrimitive() -> _BondExpr? {
        switch self.LexPtr.next() {
        case "-": return buildBondLeaf(BE_SINGLE)
        case "=": return buildBondLeaf(BE_DOUBLE)
        case "#": return buildBondLeaf(BE_TRIPLE)
        case "$": return buildBondLeaf(BE_QUAD)
        case ":": return buildBondLeaf(BE_AROM)
        case "@": return buildBondLeaf(BE_RING)
        case "~": return buildBondLeaf(BE_ANY)
        // return BuildBondLeaf(*LexPtr == '?' ? BE_UPUNSPEC : BE_UP);
        case "/": return buildBondLeaf(BE_SINGLE)
        // return BuildBondLeaf(*LexPtr == '?' ? BE_DOWNUNSPEC : BE_DOWN);
        case "\\": return buildBondLeaf(BE_SINGLE)
        default: break
        }
        self.LexPtr.dec()
        return nil
    }
    
    private func parseBondExpr(_ level: Int) -> _BondExpr? {
        return nil
    }
    
    private func parseSMARTSString(_ ptr: String) -> Pattern? {
        
        self.LexPtr.setLex(ptr)
        
        return self.parseSMARTSPattern()
    }

    private func parseSMARTSRecord(_ ptr: String) -> Pattern? {
        // Trim whitespace
        return self.parseSMARTSString(ptr.trimmingCharacters(in: .whitespacesAndNewlines))
    }
    
    //! \return the vector binding of the atom @p idx in the internal pattern
    func getVectorBinding(_ idx: Int) -> Int {
        return _pat?.atom[idx].vb ?? 0
    }
    
    private func SMARTSParser(_ pat: Pattern, _ stat: ParseState, _ prev: Int, _ part: Int) -> Pattern? {
        return nil 
    }


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

func markGrowBonds(_ pat: inout Pattern) {
    let bv = MKBitVec()
    
    for i in 0..<pat.bcount {
        pat.bond[i].grow = (bv[pat.bond[i].src] && bv[pat.bond[i].dst]) ? false : true
        
        bv.setBitOn(UInt32(pat.bond[i].src))
        bv.setBitOn(UInt32(pat.bond[i].dst))
    }
}

func getChiralFlag(_ expr: _AtomExpr) -> Int {

    switch expr {
    case .leaf(let _expr1):
        if (expr as! _AtomExprProtocol).type == AE_CHIRAL {
            return _expr1.value
        }
    case .mon(let expr1):
        // Treat [!@] as [@@], and [!@@] as [@]
        if (expr as! _AtomExprProtocol).type == AE_NOT {
            var tmp1 = getChiralFlag(expr1.arg)
            if tmp1 == AL_ANTICLOCKWISE { return AL_CLOCKWISE }
            if tmp1 == AL_CLOCKWISE { return AL_ANTICLOCKWISE }
        }
    case .bin(let expr1):
        if (expr as! _AtomExprProtocol).type == AE_OR {
            var tmp1 = getChiralFlag(expr1.lft)
            var tmp2 = getChiralFlag(expr1.rgt)
            if tmp1 == 0 || tmp2 == 0 { return 0 }
            if tmp1 == tmp2 { return tmp1 }
        }
        if (expr as! _AtomExprProtocol).type == AE_ANDHI ||
           (expr as! _AtomExprProtocol).type == AE_ANDLO {
            var tmp1 = getChiralFlag(expr1.lft)
            var tmp2 = getChiralFlag(expr1.rgt)
            if tmp1 == 0 { return 0 }
            if tmp2 == 0 { return 0 }
            if tmp1 == tmp2 { return tmp1 }
        }
    case .recur(_):
        break
    }
    return 0

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

    private var _uatoms: Bool = false
    private var _mol: MKMol 
    private var _pat: Pattern
    private var _map: [Int] = []

    init(mol: MKMol, pat: Pattern) {
        self._mol = mol
        self._pat = pat
    }

    func match(_ v: [[Int]], _ bidx: Int = -1) {

    }
}

public func smartsLexReplace(_ s: String, _ p: [Pair<String, String>]) {

}


// Helper class to act more like a C++ char buffer with pointer indexing

class LexicalParser: IteratorProtocol, Equatable {
   
    private var _lexCharacters: [Character] = []
    private var _index: Int = 0
    
    typealias Element = Character
    
    init() { }
    
    func setLex(_ string: String) {
        self._lexCharacters = Array<Character>(string)
    }
    
    public func next() -> Character? {
        defer {
            self._index += 1
        }
        return self._index >= _lexCharacters.count ? "\0" : _lexCharacters[self._index]
    }
    
    public subscript (_ idx: Int) -> Character? {
        if idx < self._lexCharacters.count {
            return _lexCharacters[idx]
        }
        return nil
    }
    
    public func inc() {
        self._index += 1
    }
    
    public func dec() {
        self._index -= 1
    }
    
    public func cur() -> Character {
        self._lexCharacters[self._index]
    }

    public func prev() -> Character {
        if self._index > 0 {
            return self._lexCharacters[self._index - 1]
        } else {
            return self._lexCharacters.first!
        }   
    }
    
    static public func += (_ lhs: LexicalParser, _ rhs: Int) {
        lhs._index += rhs
        if lhs._index < 0 { lhs._index = 0 }
    }
    
    static public func -= (_ lhs: LexicalParser, _ rhs: Int) {
        lhs._index -= rhs
        if lhs._index < 0 { lhs._index = 0 }
    }
    
    func empty() -> Bool {
        return _lexCharacters.count > 0 || self._index >= _lexCharacters.count || self._index == -1
    }
    
    static func == (lhs: LexicalParser, rhs: LexicalParser) -> Bool {
        return lhs._lexCharacters == rhs._lexCharacters
    }
    
    static func == (lhs: LexicalParser, rhs: Character) -> Bool {
        return lhs._lexCharacters[lhs._index] == rhs
    }
    
    static func != (lhs: LexicalParser, rhs: Character) -> Bool {
        return lhs._lexCharacters[lhs._index] != rhs
    }
    
    
}

