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

protocol _MKAtomExpr {}

protocol _AtomExprProtocol {
    var type: Int { get set }
}

struct _AtomExprLeaf: _AtomExprProtocol, _MKAtomExpr {
    var type: Int
    var value: Int
}

struct _AtomExprRecur: _AtomExprProtocol, _MKAtomExpr {
    var type: Int
    var recur: Pattern?
}

struct _AtomExprMon: _AtomExprProtocol, _MKAtomExpr {
    var type: Int
    var arg: _MKAtomExpr
}

struct _AtomExprBin: _AtomExprProtocol, _MKAtomExpr {
    var type: Int
    var lft: _MKAtomExpr
    var rgt: _MKAtomExpr
}

indirect enum _AtomExpr: _MKAtomExpr {
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
struct BondSpec: Equatable {
    
    var expr: _BondExpr
    var src, dst: Int
    var visit: Int?
    var grow: Bool?
    
    static func == (lhs: BondSpec, rhs: BondSpec) -> Bool {
        lhs.src == rhs.src && lhs.dst == rhs.dst && lhs.grow == rhs.grow && lhs.visit == rhs.visit
    }
    
}

//! \struct AtomSpec parsmart.h <openbabel/parsmart.h>
//! \brief An internal (SMARTS parser) atom specification
struct AtomSpec: Equatable {

    var expr: _MKAtomExpr
    var visit: Int?
    var part: Int
    var chiral_flag: Int?
    var vb: Int
    var nbrs: [Int]?
    
    static func == (lhs: AtomSpec, rhs: AtomSpec) -> Bool {
        lhs.part == rhs.part && lhs.visit == rhs.visit && lhs.chiral_flag == rhs.chiral_flag && lhs.vb == rhs.vb && lhs.nbrs == rhs.nbrs
        // add _expr comparators if truly needed
    }
}

//! \struct Pattern parsmart.h <openbabel/parsmart.h>
//! \brief A SMARTS parser internal pattern
struct Pattern: Equatable {
    
    
    var aalloc: Int = 0
    var acount: Int  {
        get {
            return atom.count
        }
    }
    var balloc: Int = 0
    var bcount: Int {
        get {
            return bond.count
        }
    }
    var isChiral: Bool = false
    var atom: [AtomSpec] = []
    var bond: [BondSpec] = []
    var parts: Int = 1
    var hasExplicitH: Bool = false
    
    static func == (lhs: Pattern, rhs: Pattern) -> Bool {
        if lhs.acount != rhs.acount || lhs.bcount != rhs.bcount {
            return false
        }
        if lhs.isChiral != rhs.isChiral {
            return false
        }
        if lhs.parts != rhs.parts || lhs.hasExplicitH != rhs.hasExplicitH {
            return false
        }
        if lhs.atom != rhs.atom {
            return false
        }
        if lhs.bond != rhs.bond {
            return false
        }
        return true
    }
    
}

//! \struct ParseState parsmart.h <openbabel/parsmart.h>
//! \brief A SMARTS parser internal state
struct ParseState {
    var closord: [_BondExpr?] = []
    var closure: [Int] = []
    var closindex: Int = 0
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
    func getBond(_ src: inout Int, _ dst: inout Int, _ ord: inout Int, _ idx: Int) {
        guard let pat = _pat else { return }
        src = pat.bond[idx].src
        dst = pat.bond[idx].dst
        ord = getExprOrder(pat.bond[idx].expr)
    }
    
    func getAtomicNum(_ idx: Int) -> Int? {
        //        TODO: throw error here on nil
        return getExprAtomicNum((_pat?.atom[idx].expr)!)
    }
    
    func getCharge(_ idx: Int) -> Int? {
        //        TODO: throw error here on nil
        return getExprCharge((_pat?.atom[idx].expr)!)
    }
    
    //! \name Matching methods (SMARTS on a specific OBMol)
    //@{
    //! Perform SMARTS matching for the pattern specified using Init().
    //! \param mol The molecule to use for matching
    //! \param single Whether only a single match is required (faster). Default is false.
    //! \return Whether matches occurred
    func match(_ mol: MKMol, _ single: Bool = false) -> Bool {
        let matcher = MKSmartsMatcher()
        guard var _pat = _pat else { return false }
        if _pat.hasExplicitH { //The SMARTS pattern contains [H]
            //Do matching on a copy of mol with explicit hydrogens
            var tmol = mol
            tmol.addHydrogens(false, false)
            return matcher.match(tmol, &_pat, &_mlist, single)
        }
        return matcher.match(mol, &_pat, &_mlist, single)
    }
    
    //! \name Matching methods (SMARTS on a specific OBMol)
    //@{
    //! Perform SMARTS matching for the pattern specified using Init().
    //! This version is (more) thread safe.
    //! \param mol The molecule to use for matching
    //! \param mlist The resulting match list
    //! \param mtype The match type to use. Default is All.
    //! \return Whether matches occurred
    func match(_ mol: MKMol, _ mlist: inout [[Int]], _ mtype: MatchType = .All) -> Bool {
        let matcher = MKSmartsMatcher()
        mlist.removeAll()
        guard var _pat = _pat else { return false }
        
        if _pat.hasExplicitH { //The SMARTS pattern contains [H]
            //Do matching on a copy of mol with explicit hydrogens
            let tmol = mol
            tmol.addHydrogens(false, false)
            if !matcher.match(tmol, &_pat, &mlist, mtype == .Single) { return false }
        } else if !matcher.match(mol, &_pat, &mlist, mtype == .Single) { return false }
        
        if mtype == .AllUnique && mlist.count > 1 {
            //uniquify
            var ok = true
            var bv = MKBitVec()
            var vbv: [MKBitVec] = []
            var ulist: [[Int]] = []
            
            for i in 0..<mlist.count {
                ok = true
                bv.clear()
                bv.fromVecInt(mlist[i])
                for j in 0..<vbv.count {
                    if ok && bv == vbv[j] {
                        ok = false
                    }
                }
                if ok {
                    ulist.append(mlist[i])
                    vbv.append(bv)
                }
            }
            mlist = ulist
        }
        return true
    }
    
    //! \name Matching methods (SMARTS on a specific OBMol)
    //@{
    //! Thread safe check for any SMARTS match
    //! \param mol The molecule to use for matching
    //! \return Whether there exists any match
    func hasMatch(_ mol: MKMol) -> Bool {
        var dummy: [[Int]] = []
        return match(mol, &dummy, .Single)
    }
    
    func restrictedMatch(_ mol: MKMol, _ pr: [Pair<Int, Int>], _ single: Bool = false) -> Bool {
        var ok = false
        var mlist: [[Int]] = []
        guard var _pat = _pat else { return false }
        let matcher = MKSmartsMatcher()
        matcher.match(mol, &_pat, &mlist)
        
        _mlist.removeAll()
        if mlist.isEmpty { return false }
        
        let mlistIterator = mlist.makeIterator()
        let prIterator = pr.makeIterator()
        
        for i in mlistIterator {
            ok = true
            for j in prIterator {
                if ok && (i[j.0] != j.1) {
                    ok = false
                }
            }
            if ok {
                _mlist.append(i)
            }
            if single && !_mlist.isEmpty {
                return true
            }
        }
        return _mlist.isEmpty ? false : true
    }
    
    func restrictedMatch(_ mol: MKMol, _ vres: MKBitVec, _ single: Bool = false) -> Bool {
        var ok = false
        var mlist: [[Int]] = []
        guard var _pat = _pat else { return false }
        let matcher = MKSmartsMatcher()
        matcher.match(mol, &_pat, &mlist)
        
        _mlist.removeAll()
        if mlist.isEmpty { return false }
        
        let mlistIterator = mlist.makeIterator()
        
        for i in mlistIterator {
            ok = true
            for j in i.makeIterator() {
                if !vres[j] {
                    ok = false
                    break
                }
            }
            if !ok {
                continue
            }
            _mlist.append(i)
            if single && !_mlist.isEmpty {
                return true
            }
        }
        return _mlist.isEmpty ? false : true
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
        let bv: MKBitVec = MKBitVec()
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
        
        var result: Pattern? = Pattern()
        
        while self.LexPtr == "(" {
            
            self.LexPtr.inc()
            
            //            TODO: reconcile passing references and values with return types
            if result == nil { return nil }
            guard var result = self.parseSMARTSPart(&result, result!.parts) else { return nil }
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
    
    private func parseSMARTSPart(_ result: inout Pattern?, _ part: Int) -> Pattern? {
        
        var stat: ParseState = ParseState()
        stat.closure.reserveCapacity(100)
        var flag: Bool
        
        for i in 0..<100 {
            stat.closure[i] =  -1
        }
        var prev = -1
        var result: Pattern? = SMARTSParser(&result, &stat, &prev, part)
        
        flag = false
        for i in 0..<100 {
            if stat.closure[i] != -1 {
                flag = true
            }
        }
        
        if var result = result {
            if flag {
                // TODO: throw error here
                print("SMARTS Error")
                return result
            } else {
                markGrowBonds(&result)
                result.isChiral = false
                for i in 0..<result.acount {
                    result.atom[i].chiral_flag = getChiralFlag(result.atom[i].expr)
                    if (result.atom[i].chiral_flag != nil) {
                        result.isChiral = true
                    }
                }
                return result
            }
        } else {
            return nil
        }
    }
    
    // private func SMARTSError(_ pat: Pattern) -> Pattern {}
    // private func parseSMARTSError(_ pat: Pattern, _ expr: _BondExpr) -> Pattern {}
    
    private func parseSimpleAtomPrimitive() -> _MKAtomExpr? {
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
    
    private func parseComplexAtomPrimitive() -> _MKAtomExpr? {
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
    
    private func parseAtomExpr(_ level: Int) -> _MKAtomExpr? {
        var expr1: _MKAtomExpr?
        var expr2: _MKAtomExpr?
        var prev: Character
        
        switch level {
        case 0: /* Low Precedence Conjunction */
            expr1 = parseAtomExpr(1) ?? nil
            if expr1 != nil  {
                repeat {
                    self.LexPtr.inc()
                    expr2 = parseAtomExpr(1) ?? nil
                    if expr2 != nil {
                        //                        MARK: probably need better error handling here but they should not be nil
                        expr1 = buildAtomBinary(AE_ANDLO, expr1!, expr2!)
                    } else { return nil}
                } while self.LexPtr.cur() == ";"
            } else { return nil }
            return expr1
        case 1: /* Disjunction */
            expr1 = parseAtomExpr(2) ?? nil
            if expr1 != nil  {
                repeat {
                    self.LexPtr.inc()
                    expr2 = parseAtomExpr(2) ?? nil
                    if expr2 != nil {
                        //                        MARK: probably need better error handling here but they should not be nil
                        expr1 = buildAtomBinary(AE_OR, expr1!, expr2!)
                    } else { return nil}
                } while self.LexPtr.cur() == ","
            } else { return nil }
            return expr1
        case 2: /* High Precedence Conjunction */
            expr1 = parseAtomExpr(3) ?? nil
            if expr1 != nil  {
                repeat {
                    if self.LexPtr.cur() == "&" {
                        self.LexPtr.inc()
                    }
                    prev = self.LexPtr.cur()
                    expr2 = parseAtomExpr(3) ?? nil
                    if expr2 == nil {
                        if prev != self.LexPtr.cur() {
                            return nil
                        } else { return expr1 }
                    } else {
                        expr1 = buildAtomBinary(AE_ANDHI, expr1!, expr2!)
                    }
                } while ((self.LexPtr.cur() != "]") &&
                         (self.LexPtr.cur() != ";") &&
                         (self.LexPtr.cur() != ",") && !self.LexPtr.empty())
            } else { return nil }
            return expr1
        case 3: /* Negation or Primitive */
            if self.LexPtr.cur() == "!" {
                self.LexPtr.inc()
                expr1 = parseAtomExpr(3) ?? nil
                if expr1 != nil {
                    buildAtomNot(expr1!)
                } else { return nil }
            }
            return parseComplexAtomPrimitive()
        default: break
        }
        
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
        var expr1: _BondExpr?
        var expr2: _BondExpr?
        var prev: Character
        
        switch level {
        case 0: /* Low Precedence Conjunction */
            expr1 = parseBondExpr(1) ?? nil
            if expr1 != nil  {
                repeat {
                    self.LexPtr.inc()
                    expr2 = parseBondExpr(1) ?? nil
                    if expr2 != nil {
                        //                        MARK: probably need better error handling here but they should not be nil
                        expr1 = buildBondBin(BE_ANDLO, expr1!, expr2!)
                    } else { return nil}
                } while self.LexPtr.cur() == ";"
            } else { return nil }
            return expr1
        case 1: /* Disjunction */
            expr1 = parseBondExpr(2) ?? nil
            if expr1 != nil  {
                repeat {
                    self.LexPtr.inc()
                    expr2 = parseBondExpr(2) ?? nil
                    if expr2 != nil {
                        //                        MARK: probably need better error handling here but they should not be nil
                        expr1 = buildBondBin(BE_OR, expr1!, expr2!)
                    } else { return nil}
                } while self.LexPtr.cur() == ","
            } else { return nil }
            return expr1
        case 2: /* High Precedence Conjunction */
            expr1 = parseBondExpr(3) ?? nil
            if expr1 != nil  {
                repeat {
                    if self.LexPtr.cur() == "&" {
                        self.LexPtr.inc()
                    }
                    prev = self.LexPtr.cur()
                    expr2 = parseBondExpr(3) ?? nil
                    if expr2 == nil {
                        if prev != self.LexPtr.cur() {
                            return nil
                        } else { return expr1 }
                    } else {
                        expr1 = buildBondBin(BE_ANDHI, expr1!, expr2!)
                    }
                } while ((self.LexPtr.cur() != "]") &&
                         (self.LexPtr.cur() != ";") &&
                         (self.LexPtr.cur() != ",") && !self.LexPtr.empty())
            } else { return nil }
            return expr1
        case 3: /* Negation or Primitive */
            if self.LexPtr.cur() == "!" {
                self.LexPtr.inc()
                expr1 = parseBondExpr(3) ?? nil
                if expr1 != nil {
                    buildBondNot(expr1!)
                } else { return nil }
            }
            return parseBondPrimitive()
        default: break
        }
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
    
    func getVectorBinding() -> Int {
        var vb = 0
        self.LexPtr.inc()
        if self.LexPtr.cur().isNumber {
            vb = 0
            while self.LexPtr.cur().isNumber {
                vb = vb * 10 + (self.LexPtr.next()?.wholeNumberValue)!
            }
        }
        return vb
    }
    
    private func SMARTSParser(_ pat: inout Pattern?, _ stat: inout ParseState, _ prev: inout Int, _ part: Int) -> Pattern? {
        var vb = 0
        var aexpr: _MKAtomExpr? = nil
        var bexpr: _BondExpr? = nil
        var index: Int = 0
        while !self.LexPtr.empty() {
            switch self.LexPtr.next() {
            case ".":
                print("ParseSMARTSError")
                return nil
            case "-", "=", "#", "$", ":", "~", "@", "/", "\\", "!":
                self.LexPtr.dec()
                if prev == -1 || (bexpr != nil) {
                    print("ParseSMARTSError")
                    return nil
                }
                bexpr = parseBondExpr(0) ?? nil
                if bexpr == nil {
                    print("ParseSMARTSError")
                    return nil
                }
            case "(":
                if bexpr != nil {
                    self.LexPtr.dec()
                    print("ParseSMARTSError")
                    return nil
                }
                if prev == -1 {
                    index = pat?.acount ?? -1
                    prev = -1
                    pat = SMARTSParser(&pat, &stat, &prev, part) ?? nil
                    if pat == nil { return nil }
                    if index == pat!.acount {
                        print("ParseSMARTSError")
                        return nil
                    }
                    prev = index
                } else {
                    pat = SMARTSParser(&pat, &stat, &prev, part)
                    if pat == nil { return nil }
                }
                
                if self.LexPtr.cur() != ")" {
                    print("ParseSMARTSError")
                    return nil
                }
                self.LexPtr.inc()
            case ")":
                self.LexPtr.dec()
                if prev == -1 || (bexpr != nil) {
                    print("ParseSMARTSError")
                    return nil
                }
                return pat
            case "%":
                if prev == -1 {
                    print("ParseSMARTSError")
                    return nil
                }
                if self.LexPtr[0].isNumber && self.LexPtr[1].isNumber {
                    index = 10 * (self.LexPtr.next()?.wholeNumberValue)! + (self.LexPtr.next()?.wholeNumberValue)!
                    self.LexPtr += 2
                } else {
                    print("ParseSMARTSError")
                    return nil
                }
                if stat.closure[index] == -1 {
                    stat.closord[index] = bexpr
                    stat.closure[index] = prev
                } else if stat.closure[index] != prev {
                    if bexpr == nil {
                        if stat.closord[index] == nil {
                            bexpr = generateDefaultBond()
                        } else {
                            bexpr = stat.closord[index]
                        }
                    } else if (stat.closord[index] != nil) && !equivalentBonExpr(bexpr!, stat.closord[index]!) {
                        print("ParseSMARTSError")
                        return nil
                    }
                    createBond(&pat!, bexpr!, prev, stat.closure[index])
                    stat.closure[index] = -1
                    bexpr = nil
                } else {
                    print("ParseSMARTSError")
                    return nil
                }
                // case 0 - 9
            case "0", "1", "2", "3", "4", "5", "6", "7", "8", "9":
                self.LexPtr.dec()
                if prev == -1 {
                    print("ParseSMARTSError")
                    return nil
                }
                index = self.LexPtr.next()?.wholeNumberValue ?? -1
                if stat.closure[index] == -1 {
                    // Ring Opening
                    stat.closord[index] = bexpr
                    stat.closure[index] = prev
                    pat!.atom[prev].nbrs?.append(-1) // Store the BC idx as a -ve
                    bexpr = nil
                } else if stat.closure[index] != prev {
                    if bexpr == nil {
                        if stat.closord[index] == nil {
                            bexpr = generateDefaultBond()
                        } else {
                            bexpr = stat.closord[index]
                        }
                    } else if (stat.closord[index] != nil) && !equivalentBonExpr(bexpr!, stat.closord[index]!) {
                        print("ParseSMARTSError")
                        return nil
                    }
                    createBond(&pat!, bexpr!, prev, stat.closure[index])
                    
                    pat!.atom[prev].nbrs?.append(stat.closure[index])
                    for nbr_idx in pat!.atom[stat.closure[index]].nbrs!.makeIterator() {
                        if pat!.atom[stat.closure[index]].nbrs?[nbr_idx] == -index {
                            pat!.atom[stat.closure[index]].nbrs?[nbr_idx] = prev
                        }
                    }
                    stat.closure[index] = -1
                    bexpr = nil
                } else {
                    print("ParseSMARTSError")
                    return nil
                }
            case "[":
                // shortcut for '[H]' primitive (PR#1463791)
                if self.LexPtr.cur() == "H" && self.LexPtr[1] == "]" {
                    aexpr = generateElement(1)
                    self.LexPtr.inc()
                    pat!.hasExplicitH = true
                } else {
                    aexpr = parseAtomExpr(0)
                }
                vb = self.LexPtr.cur() == ":" ? getVectorBinding() : 0
                if aexpr == nil || self.LexPtr.cur() == "]" {
                    print("ParseSMARTSError")
                    return nil
                }
                index = createAtom(&pat!, aexpr!, part, vb)
                if prev != -1 {
                    if bexpr == nil {
                        bexpr = generateDefaultBond()
                    }
                    createBond(&pat!, bexpr!, prev, index)
                    pat?.atom[index].nbrs?.append(prev)
                    pat?.atom[prev].nbrs?.append(index)
                    bexpr = nil
                }
                if self.LexPtr[-1] == "H" && self.LexPtr[-2] == "@" { // i.e. [C@H] or [C@@H]
                    pat?.atom[index].nbrs?.append(SmartsImplicitRef)
                }
                prev = index
                self.LexPtr.inc()
            default:
                self.LexPtr.dec()
                aexpr = parseSimpleAtomPrimitive()
                if aexpr == nil {
                    print("ParseSMARTSError")
                    return nil
                }
                index = createAtom(&pat!, aexpr!, part)
                if prev != -1 {
                    if bexpr == nil {
                        bexpr = generateDefaultBond()
                    }
                    createBond(&pat!, bexpr!, prev, index)
                    pat?.atom[index].nbrs?.append(prev)
                    pat?.atom[prev].nbrs?.append(index)
                    bexpr = nil
                }
                prev = index
            }
        }
        
        if (prev == -1) || (bexpr != nil) {
            print("ParseSMARTSError")
            return nil
        }
        
        return pat
    }
    
}


/*================================*/
/*  Atom Expression Manipulation  */
/*================================*/

func buildAtomPred(_ type: Int) -> _MKAtomExpr {
    let res = _AtomExprLeaf(type: type, value: 0)
    return _AtomExpr.leaf(res)
}

func buildAtomLeaf(_ type: Int, _ value: Int) -> _MKAtomExpr {
    let res = _AtomExprLeaf(type: type, value: value)
    return _AtomExpr.leaf(res)
}

func buildAtomNot(_ expr: _MKAtomExpr) -> _MKAtomExpr {
    let res = _AtomExprMon(type: AE_NOT, arg: expr)
    return _AtomExpr.mon(res)
}

func buildAtomBinary(_ type: Int, _ lhs: _MKAtomExpr, _ rhs: _MKAtomExpr) -> _MKAtomExpr {
    let res = _AtomExprBin(type: type, lft: lhs, rgt: rhs)
    return _AtomExpr.bin(res)
}

func buildAtomRecurs(_ pat: Pattern) -> _MKAtomExpr {
    let res = _AtomExprRecur(type: AE_RECUR, recur: pat)
    return _AtomExpr.recur(res)
}

func generateElement(_ elem: Int) -> _MKAtomExpr {
    return buildAtomLeaf(AE_ELEM, elem)
}

func generateAromElem(_ elem: Int, _ flag: Bool) -> _MKAtomExpr {
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

func createAtom(_ pat: inout Pattern, _ expr: _MKAtomExpr, _ part: Int, _ vb: Int = 0) -> Int {
    let tmp = AtomSpec(expr: expr, part: part, vb: vb)
    pat.atom.append(tmp)
    return pat.acount
}

func createBond(_ pat: inout Pattern, _ expr: _BondExpr, _ src: Int, _ dst: Int) -> Int {
    let tmp = BondSpec(expr: expr, src: src, dst: dst)
    pat.bond.append(tmp)
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

func getChiralFlag(_ expr: _MKAtomExpr) -> Int {

    switch (expr as! _AtomExpr) {
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

func getExprCharge(_ expr: _MKAtomExpr) -> Int {
    switch expr as! _AtomExpr {
    case .leaf(let _expr1):
        if (expr as! _AtomExprProtocol).type == AE_CHARGE {
            return _expr1.value
        }
    case .mon(let expr1):
        break
    case .bin(let expr1):
        if (expr as! _AtomExprProtocol).type == AE_OR {
            var tmp1 = getExprCharge(expr1.lft)
            if tmp1 == 0  { return 0 }
            var tmp2 = getExprCharge(expr1.rgt)
            if tmp2 == 0 { return 0 }
            if tmp1 == tmp2 { return tmp1 }
        }
        if (expr as! _AtomExprProtocol).type == AE_ANDHI ||
           (expr as! _AtomExprProtocol).type == AE_ANDLO {
            var tmp1 = getExprCharge(expr1.lft)
            var tmp2 = getExprCharge(expr1.rgt)
            if tmp1 == 0 { return tmp2 }
            if tmp2 == 0 { return tmp1 }
            if tmp1 == tmp2 { return tmp1 }
        }
    case .recur(_):
        break
    }
    return 0
}

func getExprAtomicNum(_ expr: _MKAtomExpr) -> Int {
    switch expr as! _AtomExpr {
    case .leaf(let _expr1):
        if (expr as! _AtomExprProtocol).type == AE_ELEM ||
           (expr as! _AtomExprProtocol).type == AE_AROMELEM || 
           (expr as! _AtomExprProtocol).type == AE_ALIPHELEM {
            return _expr1.value
        }
    case .mon(let expr1):
        break
    case .bin(let expr1):
        if (expr as! _AtomExprProtocol).type == AE_OR {
            var tmp1 = getExprAtomicNum(expr1.lft)
            if tmp1 == 0  { return 0 }
            var tmp2 = getExprAtomicNum(expr1.rgt)
            if tmp2 == 0 { return 0 }
            if tmp1 == tmp2 { return tmp1 }
        }
        if (expr as! _AtomExprProtocol).type == AE_ANDHI ||
           (expr as! _AtomExprProtocol).type == AE_ANDLO {
            var tmp1 = getExprAtomicNum(expr1.lft)
            var tmp2 = getExprAtomicNum(expr1.rgt)
            if tmp1 == 0 { return tmp2 }
            if tmp2 == 0 { return tmp1 }
            if tmp1 == tmp2 { return tmp1 }
        }
    case .recur(_):
        break
    }
    return 0
}

func getExprOrder(_ expr: _BondExpr) -> Int { 
    switch (expr as! _BondExprProtocol).type {
    case BE_SINGLE: return 1
    case BE_DOUBLE: return 2
    case BE_TRIPLE: return 3
    case BE_QUAD: return 4
    case BE_AROM: return 5

    case BE_UP, BE_DOWN, BE_UPUNSPEC, BE_DOWNUNSPEC: return 1
    case BE_ANDHI, BE_ANDLO:
    switch expr{
        case .bin(let expr1):
            var tmp1 = getExprOrder(expr1.lft)
            var tmp2 = getExprOrder(expr1.rgt)
            if tmp1 == 0 { return tmp2 }
            if tmp2 == 0 { return tmp1 }
            if tmp1 == tmp2 { return tmp1 }
        default: break
    }
    case BE_OR:
    switch expr{
        case .bin(let expr1):
            var tmp1 = getExprOrder(expr1.lft)
            if tmp1 == 0  { return 0 }
            var tmp2 = getExprOrder(expr1.rgt)
            if tmp2 == 0 { return 0 }
            if tmp1 == tmp2 { return tmp1 }
        default: break
    }
    default: break
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

    func evalAtomExpr(_ expr: _MKAtomExpr, _ atom: MKAtom) -> Bool {
        switch (expr as! _AtomExprProtocol).type {
        case AE_TRUE: return true
        case AE_FALSE: return false
        case AE_AROMATIC: return atom.isAromatic()
        case AE_ALIPHATIC: return !atom.isAromatic()
        case AE_CYCLIC: return atom.isInRing()
        case AE_ACYCLIC: return !atom.isInRing()
        default: break
        }

        switch (expr as! _AtomExpr) {
        case .leaf(var expr):
            switch (expr as! _AtomExprProtocol).type {
            case AE_MASS: 
                return expr.value == atom.getIsotope()
            case AE_ELEM:
                return expr.value == atom.getAtomicNum()
            case AE_AROMELEM:
                return atom.isAromatic() && expr.value == atom.getAtomicNum()
            case AE_ALIPHELEM:
                return !atom.isAromatic() && expr.value == atom.getAtomicNum()
            case AE_HCOUNT:
                return expr.value == (atom.getImplicitHCount() + atom.explicitHydrogenCount())
            case AE_CHARGE:
                return expr.value == atom.getFormalCharge()
            case AE_CONNECT:
                return expr.value == atom.getTotalDegree()
            case AE_DEGREE:
                return expr.value == atom.getExplicitDegree()
            case AE_IMPLICIT:
                return expr.value == atom.getImplicitHCount()
            case AE_RINGS:
                return expr.value == atom.memberOfRingCount()
            case AE_SIZE:
                return atom.isInRingSize(expr.value)
            case AE_VALENCE:
                return expr.value == atom.getTotalValence()
            case AE_CHIRAL:
            // always return true (i.e. accept the match) and check later
                return true
            case AE_HYB:
                return expr.value == atom.getHyb()
            case AE_RINGCONNECT:
                return expr.value == atom.countRingBonds()
            default: break
            }
        case .mon(var expr):
            if (expr as _AtomExprProtocol).type == AE_NOT {
                return !evalAtomExpr(expr.arg, atom)
            }
        case .bin(var expr):
            if (expr as _AtomExprProtocol).type == AE_ANDHI ||
               (expr as _AtomExprProtocol).type == AE_ANDLO {
                if !evalAtomExpr(expr.lft, atom) { return false }
                expr = expr.rgt as! _AtomExprBin
            } else if (expr as _AtomExprProtocol).type == AE_OR {
                if evalAtomExpr(expr.lft, atom) { return true }
                expr = expr.rgt as! _AtomExprBin
            }
        case .recur(var expr):
            //see if pattern has been matched
            for i in 0..<rscache.count {
                if rscache[i].0 == expr.recur {
                    return rscache[i].1[atom.getIdx()]
                }
            }
            //perceive and match pattern
            var vb: [Bool] = []
            var mlist : [[Int]] = []
            guard let par = atom.getParent() else { return false }
            if match(par, &expr.recur!, &mlist) {
                for j in mlist {
                    vb[j[0]] = true
                }
            }
            //cache result
            rscache.append(Pair(expr.recur!, vb))
            return vb[atom.getIdx()]
        }
        return false
    }

    func evalBondExpr(_ expr: inout _BondExpr, _ bond: MKBond) -> Bool {
        repeat {
            switch (expr as! _BondExprProtocol).type {
            case BE_ANDHI, BE_ANDLO:
                // get expr to _BondExprBin type
                switch expr {
                case .bin(var expr1):
                    if !evalBondExpr(&expr1.lft, bond) {
                        return false
                    }
                    expr = expr1.rgt
                default: break
                }
            case BE_OR:
                // get expr to _BondExprBin type 
                switch expr {
                case .bin(var expr1):
                    if evalBondExpr(&expr1.lft, bond) {
                        return true
                    }
                    expr = expr1.rgt
                default: break
                }
            case BE_NOT:
                // get expr to _BondExprMon type 
                switch expr {
                case .mon(var expr1):
                    return !evalBondExpr(&expr1.arg, bond)
                default: break
                }

            case BE_ANY: return true
            case BE_DEFAULT:
                return bond.getBondOrder() == 1 || bond.isAromatic()
            case BE_SINGLE: return bond.getBondOrder() == 1 && !bond.isAromatic()
            case BE_DOUBLE: return bond.getBondOrder() == 2 && !bond.isAromatic()
            case BE_TRIPLE: return bond.getBondOrder() == 3
            case BE_QUAD: return bond.getBondOrder() == 4
            case BE_AROM: return bond.isAromatic()
            case BE_RING: return bond.isInRing()
            //case BE_UP:
            //  return bond->IsUp();
            //case BE_DOWN:
            //  return bond->IsDown();
            //case BE_UPUNSPEC: // up or unspecified (i.e., not down)
            //  return !bond->IsDown();
            //case BE_DOWNUNSPEC: // down or unspecified (i.e., not up)
            //  return !bond->IsUp();
            default: return false 
            }
        } while true 
    }
    
    private func setupAtomMatchTable(_ ttab: inout [[Bool]], _ pat: Pattern, _ mol: MKMol) {
        
        ttab.reserveCapacity(pat.acount)
        for i in 0..<pat.acount {
            ttab[i].reserveCapacity(mol.numAtoms()+1)
        }
        
        for i in 0..<pat.acount {
            for atom in mol.getAtomIterator() {
                if evalAtomExpr(pat.atom[0].expr, atom) {
                    ttab[i][atom.getIdx()] = true
                }
            }
        }
    }

    private func fastSingleMatch(_ mol: MKMol, _ pat: inout Pattern, _ mlist: inout [[Int]]) {
        var bv = MKBitVec(UInt32(mol.numAtoms() + 1))
        var map: [Int] = []
        map.reserveCapacity(pat.acount)
        var vif: [Bool] = []
        var atom, a1, nbr: MKAtom
        var a1Iter: [MKIterator<MKAtom>?] = []

        if (pat.bcount != 0) {
            vif.reserveCapacity(pat.bcount)
        }

        var bcount: Int = 0 
        for atom in mol.getAtomIterator() {
            if evalAtomExpr(pat.atom[0].expr, atom) {
                
                map[0] = atom.getIdx()
                if (pat.bcount != 0) {
                    vif[0] = false
                }
                bv.clear() 
                bv.setBitOn(UInt32(atom.getIdx()))

                repeat {
                    //***entire pattern matched***
                    if bcount == pat.bcount { //save full match here
                        mlist.append(map)
                        bcount -= 1
                        return //found a single match
                    }
                    
                    //***match the next bond***
                    if !(pat.bond[bcount].grow ?? false) {
                        if !vif[bcount] {
                            if let bond = mol.getBond(map[pat.bond[bcount].src], map[pat.bond[bcount].dst])  {
                                if evalBondExpr(&pat.bond[bcount].expr, bond) {
                                    vif[bcount] = true
                                    bcount += 1
                                    if bcount < pat.bcount {
                                        vif[bcount] = false
                                    }
                                }
                            }
                            else { 
                                bcount -= 1
                            }
                        } else { //bond must have already been visited - backtrack
                            bcount -= 1
                        }
                    } else { //need to map atom and check bond
                        a1 = mol.getAtom(map[pat.bond[bcount].src])!
                        if !vif[bcount] { //figure out which nbr atom we are mapping
                            if a1Iter[bcount] == nil {
                                a1Iter[bcount] = a1.getNbrAtomIterator() ?? nil
                            }
                            nbr = a1Iter[bcount]!.next()!
                        } else {
                            bv.setBitOff(UInt32(map[pat.bond[bcount].dst]))
                            nbr = a1Iter[bcount]!.next()!
                        }

                        for nbr in a1Iter[bcount]! {
                            if !bv[nbr.getIdx()] {
                                if evalAtomExpr(pat.atom[pat.bond[bcount].dst].expr, nbr) {
                                    if let bond = mol.getBond(a1.getIdx(), nbr.getIdx()) {
                                        if evalBondExpr(&pat.bond[bcount].expr, bond) {
                                            bv.setBitOn(UInt32(nbr.getIdx()))
                                            map[pat.bond[bcount].dst] = nbr.getIdx()
                                            vif[bcount] = true
                                            bcount += 1
                                            if bcount < pat.bcount {
                                                vif[bcount] = false
                                            }
                                            break
                                        }
                                    }
                                }
                            }
                        }
                        if a1Iter[bcount]!.next() == nil {
                            bcount -= 1
                        }
                    }
                } while bcount >= 0 
            }
        }

    }

    public func match(_ mol: MKMol, _ pat: inout Pattern, _ mlist: inout [[Int]], _ single: Bool = false) -> Bool {
        
        mlist.removeAll()
        if single && !pat.isChiral {
            // perform a fast single match (only works for non-chiral SMARTS)
            fastSingleMatch(mol, &pat, &mlist)
        } else {
            // perform normal match (chirality ignored and checked below)
            var ssm: MKSSMatch = MKSSMatch(mol, pat)
            ssm.match(&mlist)
        }

        if pat.isChiral {
            var tmp: [[Int]] = []
            // iterate over the atom mappings
            for m in mlist {
                
                var allStereoCentersMatch = true 

                // for each pattern atom
                for j in 0..<pat.acount {
                    // skip non-chiral pattern atoms
                    if (pat.atom[j].chiral_flag == nil) { continue }
                    // ignore @? in smarts, parse like any other smarts
                    if pat.atom[j].chiral_flag! == AL_UNSPECIFIED { continue }
                    
                    // use the mapping the get the chiral atom in the molecule being queried
                    guard let center = mol.getAtom(m[j]) else { continue }

                    // get the OBTetrahedralStereo::Config from the molecule

                    var stereo: MKStereoFacade = MKStereoFacade(mol)
                    var ts: MKTetrahedralStereo? = stereo.getTetrahedralStereo(center.getId().rawValue) 
                    if ts == nil {
                        allStereoCentersMatch = false
                        break
                    } else if !ts!.getConfig().specified {
                        // no stereochemistry specified in molecule for the atom
                        // corresponding to the chiral pattern atom using the current
                        // mapping --> no match
                        allStereoCentersMatch = false
                        break
                    }

                    guard let nbrs: [Int] = pat.atom[j].nbrs else { continue }

                    if nbrs.count != 4 { // 3 nbrs currently not supported. Other values are errors.
                        print("ERROR: SMARTS chiral atom has \(nbrs.count) neighbors, only works with 4")
                        continue
                    }

                    // construct a OBTetrahedralStereo::Config using the smarts pattern
                    var smartsConfig = MKTetrahedralStereo.Config()
                    smartsConfig.center = center.getId().ref
                    if nbrs[0] == SmartsImplicitRef {
                        smartsConfig.from_or_towrds = .from(.ImplicitRef)
                    } else {
                        guard let ma = mol.getAtom(m[nbrs[0]])?.getId() else { continue }
                        smartsConfig.from_or_towrds = .from(ma.ref)
                    }
                    
                    var firstref: Ref
                    if nbrs[1] == SmartsImplicitRef {
                        firstref = .ImplicitRef
                    } else {
                        guard let ma = mol.getAtom(m[nbrs[1]])?.getId() else { continue }
                        firstref = ma.ref
                    }

                    guard var ra2 = mol.getAtom(m[nbrs[2]]) else { continue }
                    guard var ra3 = mol.getAtom(m[nbrs[3]]) else { continue }

                    smartsConfig.refs = MKStereo.makeRefs(firstref, ra2.getId().ref, ra3.getId().ref)
                    
                    smartsConfig.view = MKStereo.View.ViewFrom
                    
                    switch pat.atom[j].chiral_flag! {
                    case AL_CLOCKWISE:
                        smartsConfig.winding = .Clockwise
                    case AL_ANTICLOCKWISE:
                        smartsConfig.winding = .AntiClockwise
                    default:
                        smartsConfig.specified = false
                    }
                    
                    // and save the match if the two configurations are the same
                    if ts?.getConfig() != smartsConfig {
                        allStereoCentersMatch = false
                    }

                    // don't waste time checking more stereocenters using this mapping if one didn't match
                    if !allStereoCentersMatch { break }
                }
                // if all the atoms in the molecule match the stereochemistry specified
                // in the smarts pattern, save this mapping as a match
                if allStereoCentersMatch {
                    tmp.append(m)
                }
            }
            mlist = tmp
        }
        return !mlist.isEmpty
    }

}

  //*******************************************************************
  //  The OBSSMatch class performs exhaustive matching using recursion
  //  Explicit stack handling is used to find just a single match in
  //  match()
  //*******************************************************************

  //! \brief Internal class: performs fast, exhaustive matching used to find
  //! just a single match in match() using recursion and explicit stack handling.
class MKSSMatch {

    private var _uatoms: [Bool] = []
    private var _mol: MKMol 
    private var _pat: Pattern
    private var _map: [Int] = []

    init(_ mol: MKMol, _ pat: Pattern) {
        self._mol = mol
        self._pat = pat

        if !mol.isEmpty() {
            _uatoms.reserveCapacity(mol.numAtoms()+1)
        }
    }

    func match(_ mlist: inout [[Int]], _ bidx: Int = -1) {
        let matcher = MKSmartsMatcher()
        if bidx == -1 {
            for atom in _mol.getAtomIterator() {
                if matcher.evalAtomExpr(_pat.atom[0].expr, atom) {
                    _map[0] = atom.getIdx()
                    _uatoms[atom.getIdx()] = true
                    match(&mlist, 0)
                    _map[0] = 0
                    _uatoms[atom.getIdx()] = false
                }
            }
            return
        }
        if bidx == _pat.bcount { //save full match here
            mlist.append(_map)
            return
        }
        let grow = _pat.bond[bidx].grow != nil ? _pat.bond[bidx].grow! : false
        if grow { //match the next bond
            let src = _pat.bond[bidx].src
            let dst = _pat.bond[bidx].dst
            
            if _map[src] <= 0 || _map[src] > _mol.numAtoms() {
                return
            }
            
            var aexpr = _pat.atom[dst].expr
            var bexpr = _pat.bond[bidx].expr
            
            guard let atom = _mol.getAtom(_map[src]) else { return }
            if let nbratoms = atom.getNbrAtomIterator() {
                for nbr in nbratoms {
//                    MARK: could bond be nil?
                    if !_uatoms[nbr.getIdx()] && matcher.evalAtomExpr(aexpr, nbr) &&
                        matcher.evalBondExpr(&bexpr, atom.getBond(nbr)!) {
                        _map[dst] = nbr.getIdx()
                        _uatoms[nbr.getIdx()] = true
                        match(&mlist, bidx+1)
                        _map[dst] = 0
                        _uatoms[nbr.getIdx()] = false
                    }
                }
            }
            
        } else { //just check bond here
            guard let bond = _mol.getBond(_map[_pat.bond[bidx].src],
                                          _map[_pat.bond[bidx].src]) else { return }
            if matcher.evalBondExpr(&_pat.bond[bidx].expr, bond) {
                match(&mlist, bidx+1)
            }
        }
    }
}

public func smartsLexReplace(_ s: String, _ p: [Pair<String, String>]) {
// ???
// size_t j,pos;
//     std::string token,repstr;
//     std::vector<std::pair<std::string,std::string> >::iterator i;

//     for (pos = 0,pos = s.find("$",pos);pos < s.size();pos = s.find("$",pos))
//       //for (pos = 0,pos = s.find("$",pos);pos != std::string::npos;pos = s.find("$",pos))
//       {
//         pos++;
//         for (j = pos;j < s.size();++j)
//           if (!isalpha(s[j]) && !isdigit(s[j]) && s[j] != '_')
//             break;
//         if (pos == j)
//           continue;

//         token = s.substr(pos,j-pos);
//         for (i = vlex.begin();i != vlex.end();++i)
//           if (token == i->first)
//             {
//               repstr = "(" + i->second + ")";
//               s.replace(pos,j-pos,repstr);
//               j = 0;
//             }
//         pos = j;
//       }
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
    
    public subscript (_ idx: Int) -> Character {
        if _index + idx < self._lexCharacters.count {
            return _lexCharacters[_index + idx]
        }
        return _lexCharacters.last != nil ? _lexCharacters.last! : "\0"
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
        if lhs._index > lhs._lexCharacters.count { lhs._index = lhs._lexCharacters.count }
    }
    
    static public func -= (_ lhs: LexicalParser, _ rhs: Int) {
        lhs._index -= rhs
        if lhs._index < 0 { lhs._index = 0 }
        if lhs._index > lhs._lexCharacters.count { lhs._index = lhs._lexCharacters.count }
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

