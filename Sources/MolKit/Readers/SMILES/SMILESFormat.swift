//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/24/23.
//

import Foundation
import Collections
import Bitset


// This code uses the old OpenEye SMILES parser
// but replaces the SMILES export with Craig James canonical smiles
// (For regular SMILES, the canonical order is not computed and ignored)

let IMPLICIT_CIS_RING_SIZE = 8

// some constant variables
let BondUpChar: Character = "\\"
let BondDownChar: Character = "/"

// This function return true for sulfur and nitrogen
// (I'm not sure that is the right approach, longterm)
// TODO: Oxygen should be able to have lone pair as well, no?
func canHaveLonePair(_ elem: Int) -> Bool {
    switch elem {
    case MKElements.Nitrogen.atomicNum, MKElements.Sulfur.atomicNum:
        return true
    default: return false
    }
}

//Base class for SMIFormat and CANSIFormat with most of the functionality
public class SMIBaseFormat: MKMoleculeFormat {
    
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    ///
    /////////////////////////////////////////////////////////////////
    /* Lines starting with # are ignored. Whitespace at the start (including
     blank lines) terminate the input unless -e option is used.
     Valid SMILES reactions such as [C]=O.O>[Fe]>O=C=O.[H][H] with non-null
     reactant and product are accepted and the reactant, product and
     possibly the agent molecules are output when using the Convert interface
     (babel commandline). With the OBConversion functions Read, ReadString
     and ReadFile all SMILES reactions give an error when read with this format.
     */
    override public func readMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        let pmol = pOb.castAndClear(true) as! MKMol
        let ifs = pConv.getInStream()
        var ln: String = ""
        var title: String
        var smiles: String = ""
        var pos: String.Index
        
        //Ignore lines that start with #
        while ifs != nil && ifs!.peek() == "#" {
            if !ifs!.getLine(&ln) { return false }
        }
        //Get title
        if ifs!.getLine(&ln) {
            pos = ln.firstIndex(of: " ") ?? ln.firstIndex(of: "\t") ?? ln.endIndex
            if pos != ln.endIndex {
                smiles = String(ln[ln.startIndex..<pos])
                title = String(ln[ln.index(after: pos)..<ln.endIndex])
                title = title.trimmingCharacters(in: .whitespacesAndNewlines)
                pmol.setTitle(title)
            } else {
                let temp_smiles = ln.trimmingCharacters(in: .whitespacesAndNewlines)
                guard !temp_smiles.isEmpty else { return false }
                smiles = temp_smiles
            }
        }
        
        pmol.setDimension(0)
        let sp = MKSmilesParser(preserve_aromaticity: pConv.isOption("a", .INOPTIONS))
        if !pConv.isOption("S", .INOPTIONS) {
            pmol.setChiralityPerceived()
        }
        
        return sp.smiToMol(pmol, smiles) //normal return
    }
    
    override public func writeMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        var pmol: MKMol = pOb as! MKMol
        // Define some references so we can use the old parameter names
        let ofs = pConv.getOutStream()

        // Inchified SMILES? If so, then replace mol with the new 'normalised' one
        if pConv.isOption("I") {
            let success = getInchifiedSMILESMolecules(&pmol, false)
            if !success {
                do {
                    try ofs?.writeNewLine()
                } catch {
                    print("ERROR: unable to write to buffer")
                }
                print("Cannot generate Universal NSMILES for this molecule")
                return false
            }
        }

        // Title only option?
        if pConv.isOption("t") {
            do {
                try ofs?.write(data: pmol.getTitle().data(using: .utf8)!)
                try ofs?.writeNewLine()
            } catch {
                print("ERROR: unable to write to buffer")
            }
            return true
        }

        // Option 'x' needs "SMILES Atom Order" to be set
        // FIXME: When we support features of CXN extended SMILES
        //        we can rewrite this
        if pConv.isOption("x") {
            pConv.addOption("O")
        }

        var buffer = ""

        // If there is data attached called "SMILES_Fragment", then it's
        // an ascii OBBitVec, representing the atoms of a fragment.  The
        // SMILES generated will only include these fragment atoms.

        var fragatoms = Bitset()

        if let dp = pmol.getData("SMILES_Fragment") as? MKPairData<String> {
            fragatoms.initialize(fromString: dp.getValue()!)
        } else if let ppF = pConv.isOption("F") {
            fragatoms.initialize(fromString: ppF)
        } else {
        // If no "SMILES_Fragment" data, fill the entire OBBitVec
        // with 1's so that the SMILES will be for the whole molecule.
            for a in pmol.getAtomIterator() {
                fragatoms.add(a.getIdx())
            }
        }

        if pmol.numAtoms() > 0 || pmol.isReaction() {
            createCansmiString(&pmol, &buffer, &fragatoms, pConv)
        }

        var writenewline = false
    
        if !pConv.isOption("smilesonly") {
            if !pConv.isOption("n") {
                buffer += "\t"
                buffer += pmol.getTitle()
            }
            if pConv.isOption("x") && pmol.hasData("SMILES Atom Order") {
                guard let canorder = (pmol.getData("SMILES Atom Order") as? MKPairData<String>)?.getValue() else {
                    print("ERROR: Could not find atom order in mol data")
                    return false
                }
                let vs = canorder.components(separatedBy: " ")
                buffer += "\t"
                var tmp = ""
                for i in 0..<vs.count {
                    let idx = Int(vs[i])!
                    guard let atom = pmol.getAtom(idx) else {
                        print("ERROR: Could not find atom in mol at \(idx)")
                        return false
                    }
                    if i > 0 {
                        buffer += ","
                    }
                    tmp = String(format: "%.4f", atom.getX())
                    buffer += tmp
                    buffer += ","
                    tmp = String(format: "%.4f", atom.getY())
                    buffer += tmp
                }
            }
        }

        if !pConv.isOption("nonewline") {
            writenewline = true
        }

        do {
            try ofs?.write(data: buffer.data(using: .utf8)!)
            if writenewline {
                try ofs?.writeNewLine()
            }
        } catch {
            print("ERROR: unable to write to buffer")
        }
        return true 

    }
    
    override func getMIMEType() -> String {
        return "chemical/x-daylight-smiles"
    }
    
    override func targetClassDescription() -> String {
        return MKMol.classDescription()
    }
    
    override func specificationURL() -> String {
        return "http://www.daylight.com/smiles/"
    }
    
    override func skipObjects(_ n: Int, _ pConv: MKConversion) -> Int {
        if n == 0 { return 1 } // already points after the current line
        let ifs = pConv.getInStream()
        if ifs!.isEOF { return -1 }
        var i = 0
        var tmp: String = ""
        while i < n {
            if ifs!.isEOF { return -1 }
            if ifs!.getLine(&tmp) {
                i += 1
            } else {
                break
            }
        }
        return ifs!.isEOF ? -1 : 1
    }
    
    private func getInchifiedSMILESMolecules(_ mol: inout MKMol, _ useFixedHRecMet: Bool) -> Bool {
        let MolConv = MKConversion()
        guard let pInChIFormat = MKConversion.findFormat("InChI") else {
            print("InChI format not available")
            return false
        }
        let newstream = OutputStringStream()
        MolConv.setOutStream(newstream)
        if useFixedHRecMet {
            MolConv.addOption("w", .OUTOPTIONS)
            MolConv.addOption("X", .OUTOPTIONS, "RecMet FixedH")
        } else {
            MolConv.addOption("w", .OUTOPTIONS)
        }
        var success = pInChIFormat.writeMolecule(mol, MolConv)
        if !success { return false }
        let inchi = newstream.string
        if inchi.count == 0 { return false }
        let vs = inchi.components(separatedBy: " ")
        MolConv.setInFormat(pInChIFormat)
        success = MolConv.readString(&mol, vs[0])
        mol.deleteData("inchi") // Tidy up this side-effect
        return success
    }
    
}

public class SMIFormat: SMIBaseFormat {
    
    override public init() {
        super.init()
        MKConversion.registerFormat("smi", self, "chemical/x-daylight-smiles")
        MKConversion.registerFormat("smiles", self, "chemical/x-daylight-smiles")
        MKConversion.registerOptionParam("n", self)
        MKConversion.registerOptionParam("t", self)
        MKConversion.registerOptionParam("r", self)
        MKConversion.registerOptionParam("a", self)
        MKConversion.registerOptionParam("h", self)
        MKConversion.registerOptionParam("x", self)
        MKConversion.registerOptionParam("C", self)    // "anti-canonical" form (random order)
    }
    
    required init(_ id: String, _ isDefault: Bool) {
        fatalError("init(_:_:) has not been implemented")
    }
    
    override func description() -> String? {
        return """
SMILES format
A linear text format which can describe the connectivity and chirality of a molecule
Open Babel implements the `OpenSMILES specification <http://opensmiles.org>`_.

It also implements an extension to this specification for radicals.

Note that the ``l <atomno>`` option, used to specify a \"last\" atom, is
intended for the generation of SMILES strings to which additional atoms
will be concatenated. If the atom specified has an explicit H within a bracket
(e.g. ``[nH]`` or ``[C@@H]``) the output will have the H removed along with any
associated stereo symbols.

  The :ref:`Canonical_SMILES_format` produces a canonical representation
  of the molecule in SMILES format. This is the same as the ``c`` option
  below but may be more convenient to use.

Write Options e.g. -xt
  a  Output atomclass like [C:2], if available
  c  Output in canonical form
  U  Universal SMILES
  I  Inchified SMILES
  h  Output explicit hydrogens as such
  i  Do not include isotopic or chiral markings
  k  Create Kekule SMILES instead of aromatic
  n  No molecule name
  r  Radicals lower case eg ethyl is Cc
  t  Molecule name only
  x  append X/Y coordinates in canonical-SMILES order
  C  'anti-canonical' random order (mostly for testing)
  o  <ordering> Output in user-specified order
     Ordering should be specified like 4-2-1-3 for a 4-atom molecule.
     This gives canonical labels 1,2,3,4 to atoms 4,2,1,3 respectively,
     so that atom 4 will be visited first and the remaining atoms
     visited in a depth-first manner following the lowest canonical labels.
  O  Store the SMILES atom order as a space-separated string
     The string is stored as an OBPairData wth the name
     'SMILES Atom Order'.
  F  <atom numbers> Generate SMILES for a fragment
     The atom numbers should be specified like \"1 2 4 7\".
  R  Do not reuse bond closure symbols
  f  <atomno> Specify the first atom
     This atom will be used to begin the SMILES string.
  l  <atomno> Specify the last atom
     The output will be rearranged so that any additional
     SMILES added to the end will be attached to this atom.
  T  <max seconds> Specify the canonicalization timeout
     Canonicalization can take a while for symmetric molecules and a
     timeout is used. The default is 5 seconds.

Read Options e.g. -aa
  a  Preserve aromaticity present in the SMILES
     This option should only be used if reading aromatic SMILES
     generated by the same version of Open Babel. Any other
     use will lead to undefined behavior. The advantage of this
     option is that it avoids aromaticity perception, thus speeding
     up reading SMILES.
  S  Clean stereochemistry
     By default, stereochemistry is accepted as given. If you wish
     to clean up stereochemistry (e.g. by removing tetrahedral
     stereochemistry where two of the substituents are identical)
     then specifying this option will reperceive stereocenters.
"""
    }
}

class CANSMIFormat: SMIBaseFormat {
    
    override init() {
        super.init()
        MKConversion.registerFormat("can", self, "chemical/x-daylight-cansmiles")
    }
    
    required init(_ id: String, _ isDefault: Bool) {
        fatalError("init(_:_:) has not been implemented")
    }
    
    override func writeMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        //The "c" option sets us to use canonical ordering
        pConv.addOption("c", .OUTOPTIONS)
        return super.writeMolecule(pOb, pConv)
    }
    
    override func description() -> String? {
        return """
        Canonical SMILES format
        A canonical form of the SMILES linear text format
        The SMILES format is a linear text format which can describe the
           "connectivity "
        and chirality of a molecule. Canonical SMILES gives a single
           "'canonical' form for any particular molecule.
        
        .. seealso::
        
          The \"regular\" :ref:`SMILES_format` gives faster
          output, since no canonical numbering is performed.
        
        Write Options e.g. -xt
          a  Output atomclass like [C:2], if available
          h  Output explicit hydrogens as such
          i  Do not include isotopic or chiral markings
          n  No molecule name
          r  Radicals lower case eg ethyl is Cc
          t  Molecule name only
          F  <atom numbers> Generate Canonical SMILES for a fragment
             The atom numbers should be specified like \"1 2 4 7\".
          f  <atomno> Specify the first atom
             This atom will be used to begin the SMILES string.
          l  <atomno> Specify the last atom
             The output will be rearranged so that any additional
             SMILES added to the end will be attached to this atom.
             See the :ref:`SMILES_format` for more information.
        """
    }
    
}

class MKSmilesParser {
    // simple structs to make code more readable
    
    // see _extbond
    struct ExternalBond {
        var digit: Int
        var prev: Int
        var order: Int
        var updown: Character
    }
    
    // see _rclose
    struct RingClosureBond {
        var digit: Int
        var prev: Int
        var order: Int
        var updown: Character
        var numConnections: Int
    }
    
    struct StereoRingBond {
        var atoms: Array<MKAtom> = []
        var updown: Array<Character> = []
    }
    
    var _updown: Character
    var _order: Int
    var _prev: Int
    var _rxnrole: Int
    // const char *_ptr; ??
    var _ptr: LexicalParser = LexicalParser()
    var _preserve_aromaticity: Bool
    var _vprev: Array<Int>
    var _rclose: Array<RingClosureBond>
    var _extbond: Array<ExternalBond>
    var _path: Array<Int>
    var _avisit: Array<Bool>
    var _bvisit: Array<Bool>
    var _hcount: Array<Int>
    var posDouble: Array<Int>  //for extension: lc atoms as conjugated double bonds
    
    var _stereorbond: OrderedDictionary<MKBond, StereoRingBond> // Remember info on the stereo ring closure bonds
    
    // stereochimistry
    var chiralWatch: Bool = false // set when a tetrahedral atom is read
    var _tetrahedralMap: OrderedDictionary<MKAtom, MKTetrahedralStereo.Config>  // map of tetrahedral atoms and their data
    var _upDownMap: OrderedDictionary<MKBond, Character>  // store the '/' & '\' as they occurred in smiles
    var _chiralLonePair: OrderedDictionary<UInt, Character> // for atoms with potential chiral lone pairs, remember when the l.p. was encountered
    var squarePlanarWatch: Bool = false  // set when a square planar atom is read
    var _squarePlanarMap: OrderedDictionary<MKAtom, MKSquarePlanarStereo.Config>
    
    init(preserve_aromaticity: Bool = false) {
        self._preserve_aromaticity = preserve_aromaticity
        self._rxnrole = 1
        self._updown = " "
        self._order = 0
        self._prev = 0
        self._vprev = Array<Int>()
        self._rclose = Array<RingClosureBond>()
        self._extbond = Array<ExternalBond>()
        self._path = Array<Int>()
        self._avisit = Array<Bool>()
        self._bvisit = Array<Bool>()
        self._hcount = Array<Int>()
        self.posDouble = Array<Int>()
        self._stereorbond = OrderedDictionary<MKBond, StereoRingBond>()
        self._tetrahedralMap = OrderedDictionary<MKAtom, MKTetrahedralStereo.Config>()
        self._upDownMap = OrderedDictionary<MKBond, Character>()
        self._chiralLonePair = OrderedDictionary<UInt, Character>()
        self._squarePlanarMap = OrderedDictionary<MKAtom, MKSquarePlanarStereo.Config>()
    }
    
    
    func smiToMol(_ mol: MKMol, _ s: String) -> Bool {
        _vprev.removeAll()
        _rclose.removeAll()
        _prev = 0
        chiralWatch = false
        squarePlanarWatch = false
        // We allow the empty reaction (">>") but not the empty molecule ("")
        if !parseSmiles(mol, s) || (mol.isReaction() && mol.numAtoms() == 0) {
            mol.clear()
            return false
        }
        _tetrahedralMap.removeAll()        
        _squarePlanarMap.removeAll()
        mol.setAutomaticFormalCharge(false)
        return true
    }
    
    func parseSmiles(_ mol: MKMol, _ smiles: String) -> Bool {
        mol.setAromaticPerceived() // Turn off perception until the end of this function
        mol.beginModify()
        
        _ptr = LexicalParser()
        _ptr.setLex(smiles)
        
        while !_ptr.empty() {
            switch _ptr.cur() {
            case "\r":
                return false
            case "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "%":
                if _prev == 0 {
                    return false
                }
                if !parseRingBond(mol) {
                    return false
                }
            case "&":
                if _prev == 0 {
                    return false
                }
                if !parseExternalBond(mol) {
                    return false
                }
            case ".":
                _prev = 0
            case ">":
                _prev = 0
                _rxnrole += 1
                if _rxnrole == 2 {
                    mol.setIsReaction()
                    for atom in mol.getAtomIterator() {
                        let pi = MKPairData<Int>()
                        pi.setAttribute("rxnrole")
                        pi.setValue(1)
                        atom.setData(pi)
                    }
                } else if _rxnrole == 4 {
                    print("Too many greater-than signs in SMILES string")
                    print("Title is \(mol.getTitle())")
                    return false
                }
            case "(":
                _vprev.append(_prev)
            case ")":
                if _vprev.isEmpty { //CM
                    return false
                }
                _prev = _vprev.last!
                _vprev.removeLast()
            case "[":
                if !parseComplex(mol) {
                    mol.endModify()
                    mol.clear()
                    return false
                }
            case "-":
                if _prev == 0 {
                    return false
                }
                _order = 1

            case "=":
                if _prev == 0 {
                    return false
                }
                _order = 2
            case "#":
                if _prev == 0 {
                    return false
                }
                _order = 3
            case "$":
                if _prev == 0 {
                    return false
                }
                _order = 4
            case ":":
                if _prev == 0 {
                    return false
                }
                _order = 0 // no-op
            case "/":
                if _prev == 0 {
                    return false
                }
                _order = 1
                _updown = BondDownChar
            case "\\":
                if _prev == 0 {
                    return false
                }
                _order = 1
                _updown = BondUpChar
            default:
                if !parseSimple(mol) {
                    mol.endModify()
                    mol.clear()
                    return false
                }
            } // end switch
            _ptr.inc() // advance to next character
        } // end ofr _ptr
        // place dummy atoms for each unfilled external bond
        if !_extbond.isEmpty {
            capExternalBonds(mol)
        }
           // Check to see if we've balanced out all ring closures
           // They are removed from _rclose when matched
        if !_rclose.isEmpty {
            mol.endModify()
            mol.clear()
            print("Invalid SMILES string: \(_rclose.count) unmatched ring bonds.")
            return false // invalid SMILES since rings aren't properly closed
        }
           // Check to see if we've the right number of '>' for reactions
        if _rxnrole > 1 && _rxnrole != 3 {
            mol.endModify()
            print("Invalid reaction SMILES string: only a single '>' sign found (two required to be valid).")
            return false // invalid SMILES since rings aren't properly closed
        }
        if mol.isReaction() {
            let facade = MKReactionFacade(mol)
            facade.assignComponentIds()
        }
        // Apply the SMILES valence model
        for atom in mol.getAtomIterator() {
            let idx = atom.getIdx()
            let hcount = _hcount[idx - 1]
            if hcount == -1 { // Apply SMILES implicit valence model
                var bosum = 0
                for bond in atom.getBondIterator()! {
                    bosum += Int(bond.getBondOrder())
                }
                let impval = smilesValence(atom.getAtomicNum(), bosum)
                var imph = impval - bosum
                if imph > 0 && atom.isAromatic() {
                    imph -= 1
                }
                atom.setImplicitHCount(UInt(imph))
            } else { // valence is explicit e.g. [CH3]
                atom.setImplicitHCount(UInt(hcount))
            }
        }
        mol.endModify(false)
        // Unset any aromatic bonds that *are not* in rings where the two aromatic atoms *are* in a ring
        // This is rather subtle, but it's correct and reduces the burden of kekulization
        
        for bond in mol.getBondIterator() {
            if bond.isAromatic() && !bond.isInRing() {
                if bond.getBeginAtom().isInRing() && bond.getEndAtom().isInRing() {
                    bond.setAromatic(false)
                }
            }
        }
        // TODO: Only Kekulize if the molecule has a lower case atom
        let ok = MKKekulize(mol)
        
        if !ok {
            print("Failed to kekulize aromatic SMILES")
            print("Title is \(mol.getTitle())")
            return false // Should we return false for a kekulization failure?
        }
        // Add the data stored inside the _tetrahedralMap to the atoms now after
        // endModify so they don't get lost.
        if !_tetrahedralMap.isEmpty {
            for ChiralSearch in _tetrahedralMap {
                let atom = ChiralSearch.key
                var ts = ChiralSearch.value
                if ts.refs.count != 3 { continue }
                if ts.refs[2] == .NoRef {
                    // This happens where there is chiral lone pair or where there simply aren't enough connections
                    // around a chiral atom. We handle the case where there is a S with a chiral lone pair.
                    // All other cases are ignored, and raise a warning. (Note that S can be chiral even without
                    // a lone pair, think of C[S@](=X)(=Y)Cl.
                    
                    // We have remembered where to insert the lone pair in the _chiralLonePair map
                    if let chiralLonePair = _chiralLonePair.first(where: { $0.key == atom.getIdx() }) {
                        if canHaveLonePair(atom.getAtomicNum()) {
                            ts.refs[2] = ts.refs[1]
                            ts.refs[1] = ts.refs[0]
                            if chiralLonePair.value.wholeNumberValue! == 0 { // Insert in the 'from' position
                                // TODO: probably need to guard this statement
                                ts.refs[0] = ts.from_or_towrds.refValue
                                ts.from_or_towrds = .from(.ImplicitRef)
                            } else { // Insert in the refs[0] position
                                ts.refs[0] = .ImplicitRef
                            }
                        }
                    } else {
                        print("Ignoring stereochemistry. Not enough connections to this atom. \(mol.getTitle())")
                        continue
                    }
                }
                let obts = MKTetrahedralStereo(mol)
                obts.setConfig(ts)
                mol.setData(obts)
            }
        }
        // Add the data stored inside the _squarePlanarMap to the atoms now after end
        // modify so they don't get lost.
        
        if !_squarePlanarMap.isEmpty {
            for ChiralSearch in _squarePlanarMap {
                let sp = ChiralSearch.value
                if sp.refs.count != 4 { continue }
                let obsp = MKSquarePlanarStereo(mol)
                obsp.setConfig(sp)
                mol.setData(obsp)
            }
        }
        
        if !_preserve_aromaticity {
            mol.setAromaticPerceived(false)
        }
        
        createCisTrans(mol)
        return true
    }
    
    func isUp(_ bond: MKBond) -> Bool {
        if _upDownMap[bond] == BondUpChar {
            return true
        }
        return false
    }
    
    func isDown(_ bond: MKBond) -> Bool {
        if _upDownMap[bond] == BondDownChar {
            return true
        }
        return false
    }
    
    
    // Ring Closure bonds appear twice (at opening and closure).
    // If involved in cis/trans stereo, then the stereo may be
    // specified at either end or indeed both. Although Open Babel
    // will only write out SMILES with the stereo at one end (the end
    // on the double bond), it must handle all cases when reading.
    
    // For example:
    //
    //         C
    //        /|
    //   C = C |
    //  /     \|
    // C       N
    //
    // Can be written as:
    // (a) C/C=C/1\NC1 -- preferred
    // (b) C/C=C1\NC\1
    // (c) C/C=C/1\NC\1
    //  or indeed by replacing the "\N" with "N".
    
    // If the stereo chemistry for a ring closure is inconsistently specified,
    // it is ignored. In that case, if a stereo symbol does not exist for its
    // partner bond on the double bond (e.g. (b) below), then the stereo is unspecified.
    
    // (a) C/C=C/1NC\1 -- specified stereo
    // (b) C/C=C/1NC/1  -- ignore ring closure stereo => treated as C/C=C1NC1  => CC=C1NC1
    // (c) C/C=C/1\NC/1 -- ignore ring closure stereo => treated as C/C=C1\NC1 => C/C=C/1\NC1
    
    // The ring closure bond is either up or down with respect
    // to the double bond. Our task here is to figure out which it is,
    // based on the contents of _stereorbond.
    func setRingClosureStereo(_ rcstereo: StereoRingBond, _ dbl_bond: MKBond) -> Character? {
        var found = false // We have found the answer
        var updown = true // The answer
        if rcstereo.updown[0] == BondUpChar || rcstereo.updown[0] == BondDownChar { // Is there a stereo symbol at the opening?
            let on_dbl_bond = rcstereo.atoms[0] == dbl_bond.getBeginAtom() || rcstereo.atoms[0] == dbl_bond.getEndAtom()
            updown = (rcstereo.updown[0]==BondUpChar) != on_dbl_bond
            found = true
        }
        if rcstereo.updown[1] == BondUpChar || rcstereo.updown[1] == BondDownChar { // Is there a stereo symbol at the closing?
            let on_dbl_bond = rcstereo.atoms[1] == dbl_bond.getBeginAtom() || rcstereo.atoms[1] == dbl_bond.getEndAtom()
            let new_updown = (rcstereo.updown[1]==BondUpChar) != on_dbl_bond
            if !found {
                updown = new_updown
                found = true
            } else if new_updown != updown {
                print("Ignoring the cis/trans stereochemistry specified for the ring closure as it is inconsistent.")
                found = false
            }
        }
        if !found {
            return nil
        } else {
            return updown ? "1" : "2"
        }
    }
    
    func createCisTrans(_ mol: MKMol) {
        // Create a vector of CisTransStereo objects for the molecule
        for dbi in mol.getBondIterator() {
            // Not a double bond?
            if dbi.getBondOrder() != 2 || dbi.isAromatic() { // Maybe we shouldn't check for aromaticity here? Could aromatic bonds be cis/trans?
                continue
            }
            
            // Find the single bonds around the atoms connected by the double bond.
            
            let a1 = dbi.getBeginAtom()
            let a2 = dbi.getEndAtom()
            // Check that both atoms on the double bond have at least one
            // other neighbor, but not more than two other neighbors;
            // Note: In theory, we could relax the second requirement but we would
            //       need to change the data structure we use to store cis/trans
            //       stereo to only store 2 refs instead of 4
            
            let v1 = a1.getExplicitValence()
            let v2 = a2.getExplicitValence()
            if v1 < 2 || v1 > 3 || v2 < 2 || v2 > 3 {
                continue
            }
            
            var dbl_bond_atoms = [MKAtom]()
            dbl_bond_atoms.append(a1)
            dbl_bond_atoms.append(a2)
            
            var bond_stereo = [Bool](repeating: true, count: 2)      // Store the stereo of the chosen bonds at each end of the dbl bond
            var stereo_bond = [MKBond?](repeating: nil, count: 2) // These are the chosen stereo bonds
            var other_bond = [MKBond?](repeating: nil, count: 2)  // These are the 'other' bonds at each end
            
            for i in 0..<2 { // Loop over each end of the double bond in turn
                for b in dbl_bond_atoms[i].getBondIterator()! {
                    if b == dbi {
                        continue
                    }
                    if !(isUp(b) || isDown(b)) {
                        other_bond[i] = b  // Use this for the 'other' bond
                        continue
                    }
                    
                    var found: Bool = false
                    var stereo: Bool  = false
                    if let sb_it = _stereorbond[b] {
                        let bc_result = setRingClosureStereo(sb_it, dbi)
                        if bc_result != nil {
                            stereo = bc_result == "1" ? true : false
                        } else {
                            found = false
                        }
                    } else { // Not a ring closure
                        // True/False for "up/down if moved to before the double bond C"
                        // stereo = !(IsUp(b) ^ (b->GetNbrAtomIdx(dbl_bond_atoms[i]) < dbl_bond_atoms[i]->GetIdx()))
                        stereo = !boolean_xor(isUp(b), (b.getNbrAtomIdx(dbl_bond_atoms[i]) < dbl_bond_atoms[i].getIdx()))
                    }
                    
                    if (!found) { // This cannot be used as the stereo bond
                        other_bond[i] = b // Use this for the 'other' bond
                        continue
                    }
                    
                    if (stereo_bond[i] == nil) { // This is a first stereo bond
                        stereo_bond[i] = b // Use this for the 'stereo' bond
                        bond_stereo[i] = stereo
                    } else { // This is a second stereo bond
                        if (stereo != bond_stereo[i]) { // Verify that the other stereo bond (on the same atom) has opposite stereo
                            other_bond[i] = b // Use this for the 'other' bond
                        }
                        else  {
                            print("Error in cis/trans stereochemistry specified for the double bond")
                            stereo_bond[i] = nil
                        }
                    }
                }
            }
            
            if stereo_bond[0] == nil || stereo_bond[1] == nil {
                continue
            }
            // other_bond will contain NULLs if there are bonds to implicit hydrogens
            let second: Ref = (other_bond[0] == nil) ? .ImplicitRef : other_bond[0]!.getNbrAtom(a1).getId()
            let fourth: Ref = (other_bond[1] == nil) ? .ImplicitRef : other_bond[1]!.getNbrAtom(a2).getId()
            let ct = MKCisTransStereo(mol)
            let cfg = MKCisTransStereo.Config()
            cfg.begin = a1.getId()
            cfg.end = a2.getId()
            // If bond_stereo[0]==bond_stereo[1], this means cis for stereo_bond[0] and stereo_bond[1].
            if bond_stereo[0] == bond_stereo[1] {
                cfg.refs = MKStereo.makeRefs(stereo_bond[0]!.getNbrAtom(a1).getId(), second, fourth, stereo_bond[1]!.getNbrAtom(a2).getId())
            } else {
                cfg.refs = MKStereo.makeRefs(stereo_bond[0]!.getNbrAtom(a1).getId(), second, stereo_bond[1]!.getNbrAtom(a2).getId(), fourth)
            }
            ct.setConfig(cfg)
            // add the data to the atom
            mol.setData(ct)
        }
    }
    
    
    func insertTetrahedralRef(_ mol: MKMol, id: Ref) {
        
//        print(id)
        
        guard let ChiralSearchIndex = _tetrahedralMap.keys.firstIndex(where: { $0 == mol.getAtom(_prev) }) else {
            return
        }
        
        let ChiralSearchKey = _tetrahedralMap.keys[ChiralSearchIndex]
        
        if _tetrahedralMap[ChiralSearchKey] != nil {
            let insertpos = numConnections(ChiralSearchKey, id == .ImplicitRef) - 2 // -1 indicates "from"
            if insertpos > 2 {
                return
            }
            if insertpos < 0 {
                if _tetrahedralMap[ChiralSearchKey]!.from_or_towrds.from != .NoRef {
                    MKLogger.throwError(errorMsg:
                                            "Warning: Overwriting previous from reference id \(_tetrahedralMap[ChiralSearchKey]!.from_or_towrds.from) with \(id)"
                    )
                }
                _tetrahedralMap[ChiralSearchKey]!.from_or_towrds = .from(id)
            } else {
                if _tetrahedralMap[ChiralSearchKey]!.refs[insertpos] != .NoRef {
                    MKLogger.throwError(errorMsg:
                                            "Warning: Overwriting previously set reference id \(_tetrahedralMap[ChiralSearchKey]?.refs[insertpos]) with \(id)."
                    )
                }
                _tetrahedralMap[ChiralSearchKey]!.refs[insertpos] = id
            }
        }
    }
    
    func insertSquarePlanarRef(_ mol: MKMol, id: Ref) {
        
        guard let ChiralSearchIndex = _squarePlanarMap.keys.firstIndex(where: { $0 == mol.getAtom(_prev)}) else {
            return
        }
        
        let ChiralSearchKey = _squarePlanarMap.keys[ChiralSearchIndex]
        
        if _squarePlanarMap[ChiralSearchKey] != nil {
            
            let insertpos = numConnections(ChiralSearchKey, id == .ImplicitRef) - 1
            switch insertpos {
            case -1:
                if _squarePlanarMap[ChiralSearchKey]!.refs[0] != .NoRef {
                    MKLogger.throwError(errorMsg: "Warning: Overwriting previous from reference id \(String(describing: _squarePlanarMap[ChiralSearchKey]?.refs[0])) with \(id)")
                }
                _squarePlanarMap[ChiralSearchKey]!.refs[0] = id
            case 0, 1, 2, 3:
                if _squarePlanarMap[ChiralSearchKey]!.refs[insertpos] != .NoRef {
                    MKLogger.throwError(errorMsg: "Warning: Overwriting previously set reference id \(String(describing: _squarePlanarMap[ChiralSearchKey]?.refs[insertpos])) with \(id)")
                }
                _squarePlanarMap[ChiralSearchKey]!.refs[insertpos] = id
            default:
                MKLogger.throwError(errorMsg: "Warning: Square planar stereo specified for atom with more than 4 connections.")
            }
        }
    }
    
    // NumConnections finds the number of connections already made to
    // a particular atom. This is used to figure out the correct position
    // to insert an atom ID into atom4refs
    func numConnections(_ atom: MKAtom, _ isImplicitRef: Bool = false) -> Int {
        var val = atom.getExplicitDegree()
        // The implicit H is not included in "val" so we need to adjust by 1
        if isImplicitRef {
            return val + 1
        }
        
        let idx = atom.getIdx()
        
        if idx - 1 < _hcount.count && _hcount[idx - 1] > 0 {
            val += _hcount[idx - 1]
        }
        
        for bond in _rclose {
            if bond.prev == idx {
                val += 1
            }
        }
        return val
    }
    
    @discardableResult
    func parseSimple(_ mol: MKMol) -> Bool {
        
        var element: Int = 0
        var arom: Bool = false
        
        switch _ptr.cur() {
        case "*":
            element = 0
            arom = false
        case "C":
            _ptr.inc()
            if _ptr.cur() == "l" {
                element = 17
            } else {
                element = 6
                _ptr.dec()
            }
        case "N":
            element = 7
        case "O":
            element = 8
        case "S":
            element = 16
        case "P":
            element = 15
        case "F":
            element = 9
        case "I":
            element = 53
        case "B":
            _ptr.inc()
            if _ptr.cur() == "r" {
                element = 35
            } else {
                element = 5
                _ptr.dec()
            }
            // aromatics
        case "b":
            arom = true
            element = 5
        case "c":
            arom = true
            element = 6
        case "n":
            arom = true
            element = 7
        case "o":
            arom = true
            element = 8
        case "p":
            arom = true
            element = 15
        case "s":
            arom = true
            element = 16
        default:
            print("SMILES string contains a character \(_ptr.cur()) which is invalid")
            return false
        }
        
        let atom = mol.newAtom()
        atom.setAtomicNum(element)
        if _rxnrole > 1 { // Quick test for reaction
            // Set reaction role
            let pi = MKPairData<Int>()
            pi.setAttribute("rxnrole")
            pi.setValue(_rxnrole)
            atom.setData(pi)
        }
        
        if arom {
            atom.setAromatic()
        }
        
        if _prev != 0, let prevatom = mol.getAtom(_prev) { // need to add bond
            if arom && prevatom.isAromatic() && _order == 0 {
                mol.addBond(_prev, mol.numAtoms(), 1, Int(OB_AROMATIC_BOND)) // this will be kekulized later
            } else {
                mol.addBond(_prev, mol.numAtoms(), _order == 0 ? 1 : _order)
            }
            // store up/down
            if _updown == BondUpChar || _updown == BondDownChar {
                _upDownMap[mol.getBond(_prev, mol.numAtoms())!] = _updown
            }
            insertTetrahedralRef(mol, id: mol.getLastAtom()!.getId()) // .Ref(mol.numAtoms() - 1)
            insertSquarePlanarRef(mol, id: mol.getLastAtom()!.getId()) // .Ref(mol.numAtoms() - 1)
        }
        
        // set values
        _prev = mol.numAtoms()
        _order = 0 // the default is that no bond symbol has been seen
        _updown = " "
        
        _hcount.append(-1) // implicit hydrogen count
        
        return true
    }
    
    @discardableResult
    func parseComplex(_ mol: MKMol) -> Bool {
        var element = 0
        var arom = false
        var isotope = 0
        var size = 0
        // Parse isotope information
        // - we parse anything with 1 to 4 letters
        // - any bigger and we risk overflowing the short int used to
        //   store the element information (max 65536)
        _ptr.inc()
        while _ptr.cur().isNumber && size < 5 {
            isotope *= 10
            isotope += _ptr.cur().wholeNumberValue!
            _ptr.inc()
            size += 1
        }
        
        if size == 5 {
            return false
        }
        
        // Parse element information
        switch _ptr.cur() {
        case "*":
            element = 0
            break
        case "C":
            _ptr.inc()
            switch _ptr.cur() {
            case "a":
                element = 20
            case "d":
                element = 48
            case "e":
                element = 58
            case "f":
                element = 98
            case "l":
                element = 17
            case "m":
                element = 96
            case "n":
                element = 112
            case "o":
                element = 27
            case "r":
                element = 24
            case "s":
                element = 55
            case "u":
                element = 29
            default:
                element = 6
                _ptr.dec()
            }
        case "N":
            _ptr.inc()
            switch _ptr.cur() {
            case "a":
                element = 11
            case "b":
                element = 41
            case "d":
                element = 60
            case "e":
                element = 10
            case "h":
                element = 113
            case "i":
                element = 28
            case "o":
                element = 102
            case "p":
                element = 93
            default:
                element = 7
                _ptr.dec()
            }
        case "O":
            _ptr.inc()
            switch _ptr.cur() {
            case "g":
                element = 118
            case "s":
                element = 76
            default:
                element = 8
                _ptr.dec()
            }
        case "P":
            _ptr.inc()
            switch _ptr.cur() {
            case "a":
                element = 91
            case "b":
                element = 82
            case "d":
                element = 46
            case "m":
                element = 61
            case "o":
                element = 84
            case "r":
                element = 59
            case "t":
                element = 78
            case "u":
                element = 94
            default:
                element = 15
                _ptr.dec()
            }
        case "S":
            _ptr.inc()
            switch _ptr.cur() {
            case "b":
                element = 51
            case "c":
                element = 21
            case "e":
                element = 34
            case "g":
                element = 106
            case "i":
                element = 14
            case "m":
                element = 62
            case "n":
                element = 50
            case "r":
                element = 38
            default:
                element = 16
                _ptr.dec()
            }
        case "B":
            _ptr.inc()
            switch _ptr.cur() {
            case "a":
                element = 56
            case "e":
                element = 4
            case "h":
                element = 107
            case "i":
                element = 83
            case "k":
                element = 97
            case "r":
                element = 35
            default:
                element = 5
                _ptr.dec()
            }
        case "F":
            _ptr.inc()
            switch _ptr.cur() {
            case "e":
                element = 26
            case "l":
                element = 114
            case "m":
                element = 100
            case "r":
                element = 87
            default:
                element = 9
                _ptr.dec()
            }
        case "I":
            _ptr.inc()
            switch _ptr.cur() {
            case "n":
                element = 49
            case "r":
                element = 77
            default:
                element = 53
                _ptr.dec()
            }
        case "A":
            _ptr.inc()
            switch _ptr.cur() {
            case "c":
                element = 89
            case "g":
                element = 47
            case "l":
                element = 13
            case "m":
                element = 95
            case "r":
                element = 18
            case "s":
                element = 33
            case "t":
                element = 85
            case "u":
                element = 79
            default:
                return false
            }
        case "D":
            _ptr.inc()
            switch _ptr.cur() {
            case "b":
                element = 105
            case "s":
                element = 110
            case "y":
                element = 66
            default:
                return false
            }
        case "E":
            _ptr.inc()
            switch _ptr.cur() {
            case "r":
                element = 68
            case "s":
                element = 99
            case "u":
                element = 63
            default:
                return false
            }
        case "G":
            _ptr.inc()
            switch _ptr.cur() {
            case "a":
                element = 31
            case "d":
                element = 64
            case "e":
                element = 32
            default:
                return false
            }
        case "H":
            _ptr.inc()
            switch _ptr.cur() {
            case "e":
                element =  2
            case "f":
                element = 72
            case "g":
                element = 80
            case "o":
                element = 67
            case "s":
                element = 108
            default:
                element = 1
                _ptr.dec()
            }
        case "K":
            _ptr.inc()
            if _ptr.cur() == "r" {
                element = 36
            } else {
                element = 19
                _ptr.dec()
            }
        case "L":
            _ptr.inc()
            switch _ptr.cur() {
            case "a":
                element =  57
            case "i":
                element =   3
            case "r":
                element = 103
            case "u":
                element =  71
            case "v":
                element = 116
            default:
                return false
            }
        case "M":
            _ptr.inc()
            switch _ptr.cur() {
            case "c":
                element = 115
            case "d":
                element = 101
            case "g":
                element =  12
            case "n":
                element =  25
            case "o":
                element =  42
            case "t":
                element =  109
            default:
                return false
            }
        case "R":
            _ptr.inc()
            switch _ptr.cur() {
            case "a":
                element = 88
            case "b":
                element = 37
            case "e":
                element = 75
            case "f":
                element = 104
            case "g":
                element = 111
            case "h":
                element = 45
            case "n":
                element = 86
            case "u":
                element = 44
            default:
                return false
            }
        case "T":
            _ptr.inc()
            switch _ptr.cur() {
            case "a":
                element = 73
            case "b":
                element = 65
            case "c":
                element = 43
            case "e":
                element = 52
            case "h":
                element = 90
            case "i":
                element = 22
            case "l":
                element = 81
            case "m":
                element = 69
            case "s":
                element = 117
            default:
                return false
            }
        case "U":  element = 92
        case "V":  element = 23
        case "W":  element = 74
        case "X":
            _ptr.inc()
            if _ptr.cur() == "e" {
                element = 54
            } else {
                return false
            }
        case "Y":
            _ptr.inc()
            if _ptr.cur() == "b" {
                element = 70
            } else {
                element = 39
                _ptr.dec()
            }
        case "Z":
            _ptr.inc()
            switch _ptr.cur() {
            case "n":
                element = 30
            case "r":
                element = 40
            default:
                return false
            }
        case "a":
            _ptr.inc()
            if _ptr.cur() == "s" {
                arom = true
                element = 33
            } else {
                return false
            }
        case "b":
            _ptr.inc()
            if _ptr.cur() == "i" {
                arom = true
                element = 83
            } else {
                arom = true
                element = 5
                _ptr.dec()
            }
        case "c":
            arom = true
            element = 6
        case "g":
            _ptr.inc()
            if _ptr.cur() == "e" {
                arom = true
                element = 32
            } else {
                return false
            }
        case "n":
            arom = true
            element = 7
        case "o":
            arom = true
            element = 8
        case "p":
            arom = true
            element = 15
        case "s":
            arom = true
            _ptr.inc()
            switch _ptr.cur() {
            case "e":
                element = 34
            case "i":
                element = 14
            case "n":
                element = 50
            case "b":
                element = 51
            default:
                element = 16
                _ptr.dec()
            }
        case "t":
            _ptr.inc()
            if _ptr.cur() == "e" {
                arom = true
                element = 52
            } else {
                return false
            }
        case "#":
            // Only support three digits for this extension
            if let firstNum = _ptr[1], firstNum == "1" || firstNum == "2",
               let secondNum = _ptr[2], secondNum  >= "0" && secondNum <= "9",
               let thirdNum = _ptr[3], thirdNum >= "0" && thirdNum <= "9" {
                
                element = (firstNum.wholeNumberValue!) * 100 + (secondNum.wholeNumberValue!) * 10 + (thirdNum.wholeNumberValue!)
                if element > 255 {
                    let err = "Element number must be <= 255)"
                    print(err)
                    return false
                }
                _ptr += 3
                break
            }
            fallthrough
        default:
            let err = "SMILES string contains a character \(_ptr.cur()) which is invalid"
            print(err)
            return false
        }
        
        
        //handle hydrogen count, stereochemistry, and charge
        
        let atom: MKAtom = mol.newAtom()
        var hcount = 0
        var charge = 0
        var rad: Int = 0
        var clval: Int = 0
        _ptr.inc()
        lex_repeat: while !_ptr.empty() && _ptr.cur() != "]" {
            switch _ptr.cur() {
            case "@":
                _ptr.inc()
                if _ptr.cur() == "S" && _ptr[1] == "P" { // @SP1/2/3
                    // square planar atom found
                    squarePlanarWatch = true
                    if _squarePlanarMap[atom] == nil {
                        _squarePlanarMap[atom] = MKSquarePlanarStereo.Config()
                    }
                    _squarePlanarMap[atom]!.refs = [.NoRef, .NoRef, .NoRef, .NoRef]
                    _squarePlanarMap[atom]!.center = atom.getId()
                    _ptr += 2
                    switch _ptr.cur() {
                    case "1":
                        _squarePlanarMap[atom]!.shape = .ShapeU
                    case "2":
                        _squarePlanarMap[atom]!.shape = .Shape4
                    case "3":
                        _squarePlanarMap[atom]!.shape = .ShapeZ
                    default:
                        let err = "Square planar stereochemistry must be one of SP1, SP2 or SP3"
                        print(err)
                        return false
                    }
                } else {
                    // tetrahedral atom found
                    chiralWatch = true
                    if _tetrahedralMap[atom] == nil {
                        _tetrahedralMap[atom] = MKTetrahedralStereo.Config()
                    }
                    _tetrahedralMap[atom]!.refs = [.NoRef, .NoRef, .NoRef]
                    _tetrahedralMap[atom]!.center = atom.getId()
                    if _ptr.cur() == "@" {
                        _tetrahedralMap[atom]!.winding = .Clockwise
                    } else if _ptr.cur() == "?" {
                        _tetrahedralMap[atom]!.specified = false
                    } else {
                        _tetrahedralMap[atom]!.winding = .AntiClockwise
                        _ptr.dec()
                    }
                }
                break
            case "-":
                if charge != 0 {
                    let err = "Charge can only be specified once"
                    print(err)
                    return false
                }
                while _ptr.inc().cur() == "-" {
                    charge -= 1 // handle [O--]
                }
                if charge == 0 {
                    while _ptr.cur().isNumber {
                        charge = charge * 10 - (_ptr.cur().wholeNumberValue!)
                        _ptr.inc()
                    }
                    if charge == 0 { // handle [Cl-]
                        charge = -1
                    }
                } else {
                    charge -= 1 // finish handling [Ca++]
                }
                _ptr.dec()
                break
            case "+":
                if charge != 0 {
                    let err = "Charge can only be specified once"
                    print(err)
                    return false
                }
                while _ptr.inc().cur() == "+" {
                    charge += 1 // handle [Ca++]
                }
                if charge == 0 {
                    while _ptr.cur().isNumber {
                        charge = charge * 10 + (_ptr.cur().wholeNumberValue!)
                        _ptr.inc()
                    }
                    if charge == 0 { // handle [Na+]
                        charge = 1
                    }
                } else {
                    charge += 1 // finish handling [Ca++]
                }
                _ptr.dec()
                break
            case "H":
                _ptr.inc()
                if _ptr.cur().isNumber {
                    hcount = _ptr.cur().wholeNumberValue!
                } else {
                    hcount = 1
                    _ptr.dec()
                }
                break
            case ".":
                rad = 2
                if _ptr.inc().cur() == "." {
                    rad = 3
                } else {
                    _ptr.dec()
                }
                break
            case ":":
                if !_ptr.inc().cur().isNumber {
                    let err = "The atom class following : must be a number"
                    print(err)
                    return false
                }
                while _ptr.cur().isNumber && clval < 100000000 {
                    clval = clval * 10 + (_ptr.cur().wholeNumberValue!)
                    _ptr.inc()
                }
                _ptr.dec()
                let atomclass = MKPairData<Int>()
                atomclass.setAttribute("Atom Class")
                atomclass.setValue(clval)
                atomclass.setOrigin(.fileformatInput)
                atom.setData(atomclass)
                break
            default:
                return false
            }
            _ptr.inc()
        }
        
        if _ptr.empty() || _ptr.cur() != "]" {
            let err = "Expected a closing ]"
            print(err)
            return false
        }
        
        if charge != 0 {
            atom.setFormalCharge(charge)
            if abs(charge) > 10 || (element != 0 && charge > element) { // if the charge is +/- 10 or more than the number of electrons
                let err = "Atom \(atom.getId().intValue) had an unrealistic charge of \(charge)"
                print(err)
            }
        }
        
        if rad != 0 {
            atom.setSpinMultiplicity(rad)
        }
        
        atom.setAtomicNum(element)
        atom.setIsotope(UInt(isotope))
        
        if arom {
            atom.setAromatic()
        }
        
        if _rxnrole > 1 {
            let pi = MKPairData<Int>()
            pi.setAttribute("rxnrole")
            pi.setValue(_rxnrole)
            atom.setData(pi)
        }
                
        if _prev != 0 { // need to add bond
            guard let prevatom = mol.getAtom(_prev) else { return false }
            if arom && prevatom.isAromatic() && _order == 0 {
                mol.addBond(_prev, mol.numAtoms(), 1, Int(OB_AROMATIC_BOND)) // this will be kekulized later
            } else {
                mol.addBond(_prev, mol.numAtoms(), _order == 0 ? 1 : _order)
            }
            // store up/down
            if _updown == BondUpChar || _updown == BondDownChar {
                _upDownMap[mol.getBond(_prev, mol.numAtoms())!] = _updown
            }
            if chiralWatch { // if tetrahedral atom, set previous as from atom
                _tetrahedralMap[atom]!.from_or_towrds = .from(mol.getAtom(_prev)!.getId())
                if canHaveLonePair(element) { // Handle chiral lone pair as in X[S@@](Y)Z
                    _chiralLonePair[UInt(mol.numAtoms())] = "1" // First of the refs
                }
            }
            if squarePlanarWatch { // if squareplanar atom, set previous atom as first ref
                _squarePlanarMap[atom]!.refs[0] = mol.getAtom(_prev)!.getId()
            }
            insertTetrahedralRef(mol, id: atom.getId())
            insertSquarePlanarRef(mol, id: atom.getId())
        } else { // Handle chiral lone pair as in [S@@](X)(Y)Z
            if chiralWatch && canHaveLonePair(element) { // Handle chiral lone pair (only S at the moment)
                _chiralLonePair[UInt(mol.numAtoms())] = "0" // 'from' atom
            }
        }
        
        //set values
        _prev = mol.numAtoms()
        _order = 0
        _updown = " "
        
        if hcount > 0 {
            if chiralWatch {
                insertTetrahedralRef(mol, id: .ImplicitRef)
            }
            if squarePlanarWatch {
                insertSquarePlanarRef(mol, id: .ImplicitRef)
            }
        }
        
        _hcount.append(hcount)
        chiralWatch = false
        squarePlanarWatch = false
        
        return true
    }
    
    func capExternalBonds(_ mol: MKMol) -> Bool {
        if _extbond.isEmpty {
            return true
        }
        var atom: MKAtom
        for bond in _extbond {
            atom = mol.newAtom()
            atom.setAtomicNum(0)
            mol.addBond(bond.prev, atom.getIdx(), bond.order)
            if bond.updown == BondUpChar || bond.updown == BondDownChar {
                _upDownMap[mol.getBond(bond.prev, atom.getIdx())!] = bond.updown
            }
            guard let refbond = atom.getBond(mol.getAtom(bond.prev)!) else { return false }
            var xbd: MKExternalBondData
            if mol.hasData(MKGenericDataType.ExternalBondData) {
                xbd = mol.getData(MKGenericDataType.ExternalBondData) as! MKExternalBondData
            } else {
                xbd = MKExternalBondData()
                xbd.setOrigin(.fileformatInput)
                mol.setData(xbd)
            }
            xbd.setData(atom, refbond, bond.digit)
            //this data gets cleaned up in mol.Clear.
        }
        return true
    }
    
    func parseExternalBond(_ mol: MKMol) -> Bool {
        var digit: Int = 0
        var str = ""
        _ptr.inc()
        // check for bond order indicators CC&=1.C&1
        
        switch _ptr.cur() {
        case "-":
            _order = 1
            _ptr.inc()
        case "=":
            _order = 2
            _ptr.inc()
        case  "#":
            _order = 3
            _ptr.inc()
        case  "$":
            _order = 4
            _ptr.inc()
        case  ";":
            _order = 5
            _ptr.inc()
        case  "/":
            _order = 1
            _updown = BondDownChar
            _ptr.inc()
        case  "\\":
            _order = 1
            _updown = BondUpChar
            _ptr.inc()
        default:
            break
        }
        
        if _ptr.cur() == "%" { // external bond indicator > 10
            _ptr.inc()
            str.append(_ptr.cur())
            _ptr.inc()
            str.append(_ptr.cur())
        } else { // simple single digit external bond indicator
            str.append(_ptr.cur())
        }
        digit = Int(str) ?? 0 // convert indicator to digit
        // check for dot disconnect closures
        var upDown: Character = " "
        var bondOrder: Int = 0
        for (i, bond) in _extbond.enumerated() {
            if bond.digit == digit {
                upDown = (_updown > bond.updown) ? _updown : bond.updown
                bondOrder = (_order > bond.order) ? _order : bond.order
                mol.addBond(bond.prev, _prev, bondOrder)
                // store up/down
                if upDown == BondUpChar || upDown == BondDownChar {
                    _upDownMap[mol.getBond(bond.prev, _prev)!] = upDown
                }
                // after adding a bond to atom "_prev"
                // search to see if atom is bonded to a chiral atom
                insertTetrahedralRef(mol, id: .Ref(bond.prev - 1))
                insertSquarePlanarRef(mol, id: .Ref(bond.prev - 1))
                _extbond.remove(at: i)
                _updown = " "
                _order = 0
                return true
            }
        }
        // since no closures save another ext bond
        let extBond = ExternalBond(digit: digit, prev: _prev, order: _order, updown: _updown)
        _extbond.append(extBond)
        _order = 0
        _updown = " "
        return true
    }
    
    func parseRingBond(_ mol: MKMol) -> Bool {
        // The ring closure must be associated with a 'prev' atom
        guard mol.getAtom(_prev) != nil else {
            fatalError("Number not parsed correctly as a ring bond")
            // TODO: probably make this non fatal error and return false
        }
        
        var digit: Int = 0
        if _ptr.cur() == "%" {
            _ptr.inc()
            if _ptr.cur() == "(" { // %(NNN) extension to OpenSMILES
                _ptr.inc()
                let start = _ptr.cur()
                while _ptr.cur().isNumber {
                    digit *= 10
                    digit += _ptr.cur().wholeNumberValue!
                    _ptr.inc()
                    if _ptr.cur().wholeNumberValue! - start.wholeNumberValue! > 5 {
                        fatalError("Ring closure numbers with more than 5 digits are not supported")
                    }
                }
                if _ptr.cur() != ")" {
                    fatalError("Matching close parenthesis not found for ring closure number")
                }
            } else { // % followed by two-digit ring closure
                if !_ptr.cur().isNumber || !(_ptr[1]!.isNumber) {
                    fatalError("Two digits expected after %")
                }
                digit = (_ptr.cur().wholeNumberValue!) * 10 + _ptr[1]!.wholeNumberValue! + 1
                _ptr.inc()
            }
        } else {
            digit = _ptr.cur().wholeNumberValue!
        }
        
        var upDown: Character = " "
        var bondOrder: Int = 0
        for (i, bond) in _rclose.enumerated() {
            if bond.digit == digit {
                // Check for self-bonding, e.g. C11
                if bond.prev == _prev {
                    // TODO: Enter Error Handling here
                    print("ERROR: Invalid SMILES: Ring closures imply atom bonded to itself.")
                    return false
                }
                upDown = (_updown > bond.updown) ? _updown : bond.updown
                bondOrder = (_order > bond.order) ? _order : bond.order
                // Check if this ring closure bond may be aromatic and set order accordingly
                var aromatic_bond: Bool = false
                if bondOrder == 0 {
                    guard let atom1 = mol.getAtom(bond.prev) else { return false }
                    guard let atom2 = mol.getAtom(_prev) else { return false }
                    if atom1.isAromatic() && atom2.isAromatic() {
                        aromatic_bond = true
                    }
                }
                mol.addBond(bond.prev, _prev, bondOrder == 0 ? 1 : bondOrder, aromatic_bond ? Int(OB_AROMATIC_BOND): 0, insertpos: bond.numConnections)
                // store up/down
                if upDown == BondUpChar || upDown == BondDownChar {
                    _upDownMap[mol.getBond(bond.prev, _prev)!] = upDown
                }
                
                // For assigning cis/trans in the presence of bond closures, we need to
                // remember info on all bond closure bonds.
                var sb: StereoRingBond = StereoRingBond()
                sb.updown.append(_updown)
                sb.atoms.append(mol.getAtom(_prev)!)
                sb.updown.append(bond.updown)
                sb.atoms.append(mol.getAtom(bond.prev)!)
                _stereorbond[mol.getBond(bond.prev, _prev)!] = sb // Store for later
                // after adding a bond to atom "_prev"
                // search to see if atom is bonded to a chiral atom
                // need to check both _prev and bond->prev as closure is direction independent
                insertTetrahedralRef(mol, id: .Ref(_prev - 1))
                insertSquarePlanarRef(mol, id: .Ref(_prev - 1))
                // FIXME: needed for squareplanar too??
                if let ChiralSearchIndex = _tetrahedralMap.keys.firstIndex(where: { $0 == mol.getAtom(bond.prev) }) {
                    if _tetrahedralMap[_tetrahedralMap.keys[ChiralSearchIndex]] != nil {
                        let insertpos = bond.numConnections - 1
                        switch insertpos {
                        case -1:
                            if _tetrahedralMap[_tetrahedralMap.keys[ChiralSearchIndex]]!.from_or_towrds.from != .NoRef {
                                print("Warning: Overwriting previous from reference id.")
                            }
                            _tetrahedralMap[_tetrahedralMap.keys[ChiralSearchIndex]]!.from_or_towrds = .from(mol.getAtom(_prev)!.getId())
                        case 0, 1, 2:
                            if _tetrahedralMap[_tetrahedralMap.keys[ChiralSearchIndex]]!.refs[insertpos] != .NoRef {
                                print("Warning: Overwriting previously set reference id.")
                            }
                            _tetrahedralMap[_tetrahedralMap.keys[ChiralSearchIndex]]!.refs[insertpos] = mol.getAtom(_prev)!.getId()
                        default:
                            print("Warning: Tetrahedral stereo specified for atom with more than 4 connections.")
                        }
                    }
                }
                //CM ensure neither atoms in ring closure is a radical center
                
                guard let patom = mol.getAtom(_prev) else { return false }
                patom.setSpinMultiplicity(0)
                guard let patom = mol.getAtom(bond.prev) else { return false }
                patom.setSpinMultiplicity(0)
                _rclose.remove(at: i)
                _updown = " "
                _order = 0
                return true
            }
        }
        
        //since no closures save another rclose bond
        
        guard let atom = mol.getAtom(_prev) else { return false }
        let ringClosure: RingClosureBond = RingClosureBond(digit: digit, prev: _prev, order: _order, updown: _updown, numConnections: numConnections(atom))
        //store position to insert closure bond
        _rclose.append(ringClosure)
        _order = 0
        _updown = " "
        
        return true
    }
}


// MARK:  Canical SMILES Helpers

/*----------------------------------------------------------------------
 * CLASS: OBBondClosureInfo: For recording bond-closure digits as
 * work progresses on canonical SMILES.
 ----------------------------------------------------------------------*/

class MKBondClosureInfo {
    var toatom: MKAtom      // second atom in SMILES order
    var fromatom: MKAtom    // first atom in SMILES order
    var bond: MKBond
    var ringdigit: Int
    var is_open: Int        // TRUE if SMILES processing hasn't reached 'toatom' yet
    // TODO consider making this a computed property to make it conform to Bool
    
    init(toatom: MKAtom, fromatom: MKAtom, bond: MKBond, ringdigit: Int, is_open: Int) {
        self.toatom = toatom
        self.fromatom = fromatom
        self.bond = bond
        self.ringdigit = ringdigit
        self.is_open = is_open
    }
}

/*----------------------------------------------------------------------
 * CLASS: OBCanSmiNode: A Tree structure, each node of which is an atom in
 * the tree being built to write out the SMILES.
 ----------------------------------------------------------------------*/

class MKCanSmiNode {
    var _atom: MKAtom
    var _parent: MKAtom?
    var _child_nodes: [MKCanSmiNode] = []
    var _child_bonds: [MKBond] = []
    
    init(atom: MKAtom) {
        _atom = atom
        _parent = nil
    }
    
    func size() -> Int {
        return _child_nodes.count
    }
    
    func setParent(_ parent: MKAtom) {
        _parent = parent
    }
    
    func getAtom() -> MKAtom {
        return _atom
    }
    
    func getParent() -> MKAtom? {
        return _parent
    }
    
    func getChildAtom(_ i: Int) -> MKAtom {
        return _child_nodes[i].getAtom()
    }
    
    func getChildBond(_ i: Int) -> MKBond {
        return _child_bonds[i]
    }
    
    func getChildNode(_ i: Int) -> MKCanSmiNode {
        return _child_nodes[i]
    }
    
    func addChildNode(_ node: MKCanSmiNode, _ bond: MKBond) {
        _child_nodes.append(node)
        _child_bonds.append(bond)
    }
}

struct OutOptions {
    var isomeric: Bool
    var kekulesmi: Bool
    var showatomclass: Bool
    var showexplicitH: Bool
    var smarts: Bool
    var ordering: String = "" // This is a pointer to the string in the original map
}

/*----------------------------------------------------------------------
 * CLASS OBMol2Cansmi - Declarations
 ----------------------------------------------------------------------*/
class MKMol2Cansmi {
    var _atmorder: [Int] = []
    var _uatoms: Bitset = Bitset()
    var _ubonds: Bitset = Bitset()
    var _vopen: [MKBondClosureInfo] = []
    var _bcdigit: Int = 0 // Unused unless option "R" is specified
    var _cistrans: [MKCisTransStereo] = []
    var _unvisited_cistrans: [MKCisTransStereo] = []
    var _isup: [MKBond: Bool] = [:]
    
    var _canonicalOutput: Bool = false // regular or canonical SMILES
    
    var _pmol: MKMol
    var _stereoFacade: MKStereoFacade
    var _pconv: MKConversion
    
    var _endatom: MKAtom?
    var _startatom: MKAtom?
    
    var options: OutOptions
    
    init(_ options: OutOptions, _ pmol: MKMol, _ canonical: Bool, _ pconv: MKConversion) {
        _pmol = pmol
        _canonicalOutput = canonical
        _pconv = pconv
        _stereoFacade = MKStereoFacade(pmol)
        _endatom = nil
        _startatom = nil
        self.options = options
    }
    
    func createCisTrans(_ mol: MKMol) {
        guard let vdata: [MKGenericData] = mol.getAllData(.StereoData) else {
            print("ERROR: no stereodata found for mol")
            return
        }
        for data in vdata {
            if (data as? MKStereoBase)?.getType() != .CisTrans { continue }
            if let ct = data as? MKCisTransStereo {
                if ct.getConfig().specified {
                    let config = ct.getConfig()
                    if let dbl_bond = mol.getBond(mol.getAtomById(config.begin)!, mol.getAtomById(config.end)!) {
                        // Do not output cis/trans bond symbols for double bonds in rings of size IMPLICIT_CIS_RING_SIZE or less
                        let boundringsize = MKBondGetSmallestRingSize(dbl_bond, IMPLICIT_CIS_RING_SIZE)
                        if boundringsize == 0 { // either not in ring at all, or not in small ring
                            _cistrans.append(ct)
                        }
                    }
                }
            }
        }
        _unvisited_cistrans = _cistrans // Make a copy of _cistrans
    }
    
    func getCisTransBondSymbol(_ bond: MKBond, _ node: MKCanSmiNode) -> String {
        // Given a cis/trans bond and the node in the SMILES tree, figures out
        // whether to write a '/' or '\' symbol.
        // See the comments smilesformat.cpp: FixCisTransBonds().
        //
        // The OBCanSmiNode is the most-recently-written atom in the SMILES string
        // we're creating.  If it is the double-bonded atom, then the substituent
        // follows, so that "up" means '/' and "down" means '\'.  If the OBCanSmiNode
        // atom is the single-bonded atom then the double-bonded atom comes next,
        // in which case "up" means '\' and "down" means '/'.
        //
        // Note that the story is not so simple for conjugated systems where
        // we need to take into account what symbol was already used.
        
        let atom = node.getAtom()
        let nbr_atom = bond.getNbrAtom(atom)
        guard let mol = atom.getParent() else { return "" }
        
        // If this bond is in two different obcistransstereos (e.g. a conjugated system)
        // choose the one where the dbl bond atom is *atom (i.e. the one which comes first)
        
        var dbl_bond_first: Bool = false
        if atom.hasDoubleBond() {
            if nbr_atom.hasDoubleBond() {
                // Check whether the atom is a center in any CisTransStereo. If so,#
                // then this CisTransStereo takes precedence over any other
                for ChiralSearch in _cistrans {
                    let cfg = ChiralSearch.getConfig()
                    if atom.getId() == cfg.begin || atom.getId() == cfg.end {
                        // **I don't think I need to check whether it has a bond with nbr_atom**
                        dbl_bond_first = true
                        break
                    }
                }
            } else {
                dbl_bond_first = true
            }
        }
        
        // Has the symbol for this bond already been set?
        if !_isup.contains(where: { $0.key === bond }) { // No it hasn't
            var endatom: Ref
            var centeratom: Ref
            if dbl_bond_first {
                if atom.isAromatic() {
                    for bond in atom.getBondIterator()! {
                        if bond.isAromatic() && bond.getBondOrder() == 2 {
                            return ""
                        }
                    }
                }
                endatom = nbr_atom.getId()
                centeratom = atom.getId()
            } else {
                if nbr_atom.isAromatic() {
                    for bond in nbr_atom.getBondIterator()! {
                        if bond.isAromatic() && bond.getBondOrder() == 2 {
                            return ""
                        }
                    }
                }
                endatom = atom.getId()
                centeratom = nbr_atom.getId()
            }
            
            
            for ChiralSearch in _unvisited_cistrans {
                let cfg = ChiralSearch.getConfig(.ShapeU)
                if cfg.refs.contains(endatom) && (cfg.begin == centeratom || cfg.end == centeratom) { // Atoms endatom and centeratom are in this OBCisTransStereo
                    var refbonds: [MKBond?] = [nil, nil, nil, nil]
                    refbonds[0] = mol.getBond(mol.getAtomById(cfg.refs[0])!, mol.getAtomById(cfg.begin)!)
                    
                    if cfg.refs[1] != .ImplicitRef { // Could be a hydrogen
                        refbonds[1] = mol.getBond(mol.getAtomById(cfg.refs[1])!, mol.getAtomById(cfg.begin)!)
                    }
                    
                    if cfg.refs[2] != .ImplicitRef { // Could be a hydrogen
                        refbonds[2] = mol.getBond(mol.getAtomById(cfg.refs[2])!, mol.getAtomById(cfg.end)!)
                    }
                    
                    if cfg.refs[3] != .ImplicitRef { // Could be a hydrogen
                        refbonds[3] = mol.getBond(mol.getAtomById(cfg.refs[3])!, mol.getAtomById(cfg.end)!)
                    }
                    
                    // What symbol would the four refs use if before the dbl bond?
                    let config: [Bool] = [true, false, false, true]
                    var use_same_config: Bool = true
                    // (The actual config used will be config ^ use_same_config)
                    // Make sure that the symbol for this bond is true. This ensures
                    // a canonical string, so that it's always C/C=C/C and not C\C=C\C.
                    for i in 0..<4 {
                        if refbonds[i] == bond {
                            if !config[i] {
                                use_same_config = false
                                break
                            }
                        }
                    }
                    
                    for i in 0..<4 {
                        if _isup.contains(where: { $0.key == refbonds[i] }) { // We have already set this one (conjugated bond)
                            if _isup[refbonds[i]!] == Bool(Int(config[i]) ^ Int(use_same_config)) {
                                use_same_config = !use_same_config
                                break
                            }
                        }
                    }
                    for i in 0..<4 {
                        if refbonds[i] != nil {
                            _isup[refbonds[i]!] = Bool(Int(config[i]) ^ Int(use_same_config))
                        }
                    }
                    if let eraseIndex = _unvisited_cistrans.firstIndex(where: { $0 === ChiralSearch }) {
                        _unvisited_cistrans.remove(at: eraseIndex)
                    } // this one should be a memory match since we are in a loop-closure
                    break
                }
            }
        }
        
        // If ChiralSearch didn't find the bond, we can't set this symbol
        if (!_isup.contains(where: { $0.key == bond })) {
            return ""
        }
        
        if (dbl_bond_first) { // double-bonded atom is first in the SMILES
            if (_isup[bond]!) {
                return "/"
            } else {
                return "\\"
            }
        } else { // double-bonded atom is second in the SMILES
            if (_isup[bond]!) {
                return "\\"
            } else {
                return "/"
            }
        }
    }
    
    /***************************************************************************
     * FUNCTION: AtomIsChiral
     *
     * DESCRIPTION:
     *       Returns TRUE if the atom is genuinely chiral, that is, it meets
     *       the criteria from OBAtom::IsChiral, and additionally it actually
     *       has a connected hash or wedge bond.
     *
     *       We arbitrarily reject chiral nitrogen because for our purposes there's
     *       no need to consider it.
     *
     *       NOTE: This is a simplistic test.  When the full SMILES canonicalization
     *       includes chiral markings, this should check the symmetry classes
     *       of the neighbors, not the hash/wedge bonds.
     ***************************************************************************/
    func atomIsChiral(_ atom: MKAtom) -> Bool {
        let atomid = atom.getId().intValue
        return _stereoFacade.hasTetrahedralStereo(atomid) || _stereoFacade.hasSquarePlanarStereo(atomid)
    }
    
    /***************************************************************************
     * FUNCTION: BuildCanonTree
     *
     * DESCRIPTION:
     *       Builds the SMILES tree, in canonical order, for the specified
     *       molecular fragment.
     ***************************************************************************/
    @discardableResult
    func buildCanonTree(_ mol: MKMol, _ frag_atoms: inout Bitset, _ canonical_order: [Int], _ node: MKCanSmiNode) -> Bool {
        let atom = node.getAtom()
        
        // Create a vector of neighbors sorted by canonical order, but favor
        // double and triple bonds over single and aromatic.  This causes
        // ring-closure digits to avoid double and triple bonds.
        //
        // Since there are typically just one to three neighbors, we just do a
        // ordered insertion rather than sorting.
        var favor_multiple: Bool = true  // Visit 'multiple' bonds first
        if options.ordering != "" {
            favor_multiple = false // Visit in strict canonical order (if using user-specified order)
        }
        var sort_nbrs: [MKAtom] = []
        var set_nbr: Bool = false
        for nbr in atom.getNbrAtomIterator()! {
            let idx = atom.getIdx()
            if _uatoms[idx] || !frag_atoms.contains(idx) { continue }
            guard let nbr_bond = atom.getBond(nbr) else { fatalError("ERROR: Could not find bond") }
            _ = nbr_bond.getBondOrder()
            let new_needs_bsymbol = needsBondSymbol(nbr_bond)
            
            for (i, ai) in sort_nbrs.enumerated() {
                guard let bond = atom.getBond(ai) else { break }
                _ = bond.getBondOrder()
                let sorted_needs_bsymbol = needsBondSymbol(bond)
                if favor_multiple && new_needs_bsymbol && !sorted_needs_bsymbol {
                    sort_nbrs.insert(nbr, at: i)
                    break
                }
                if (!favor_multiple || new_needs_bsymbol == sorted_needs_bsymbol) && canonical_order[idx-1] < canonical_order[ai.getIdx()-1] {
                    sort_nbrs.insert(nbr, at: i)
                    set_nbr = true
                    break
                }
            }
            if !set_nbr {
                sort_nbrs.append(nbr)
            }
            set_nbr = false
        }
        
        _uatoms.add(atom.getIdx()) // mark the aotm as visited
        
        if _endatom != nil {
            if !_uatoms.contains(_endatom!.getIdx()) && sort_nbrs.count > 1 {
                // If you have specified an _endatom, the following section rearranges
                // sort_nbrs as follows:
                //   - if a branch does not lead to the end atom, move it to the front
                //     (i.e. visit it first)
                //   - otherwise move it to the end
                // This section is skipped if sort_nbrs has only a single member, or if
                // we have already visited _endatom.
                var children: [MKAtom] = []
                myFindChildren(mol, &children, _uatoms, _endatom!)
                
                var front: [MKAtom] = []
                var end: [MKAtom] = []
                for nbr in sort_nbrs {
                    if !children.contains(nbr) && nbr != _endatom {
                        front.append(nbr)
                    } else {
                        end.append(nbr)
                    }
                }
                sort_nbrs = front
                sort_nbrs.append(contentsOf: end)
            }
        }
        for nbr in sort_nbrs {
            let idx = nbr.getIdx()
            if _uatoms[idx] { continue }
            guard let bond = atom.getBond(nbr) else { fatalError("ERROR: Could not find bond") }
            _ubonds.add(Int(bond.getIdx()))
            let next = MKCanSmiNode(atom: nbr)
            next.setParent(atom)
            node.addChildNode(next, bond)
            buildCanonTree(mol, &frag_atoms, canonical_order, next)
        }
        
        return true
    }
    
    /***************************************************************************
    * FUNCTION: CreateFragCansmiString
    *
    * DESCRIPTION:
    *       Selects the "root" atom, which will be first in the SMILES, then
    *       builds a tree in canonical order, and finally generates the SMILES.
    *       If there are then atoms that haven't been visited (i.e. a molecule
    *       with disconnected parts), selects a new root from the remaining
    *       atoms and repeats the process.
    ***************************************************************************/
    func createFragCansmiString(_ mol: MKMol, _ frag_atoms: inout Bitset, _ buffer: inout String) {
        
        var symmetry_classes: [Ref] = []
        var canonical_order: [Ref] = []
        
        symmetry_classes.reserveCapacity(mol.numAtoms())
        canonical_order.reserveCapacity(mol.numAtoms())
        
        // Remember the desired endatom, if specified
        var pp: String? = _pconv.isOption("l")
        var atom_idx: Int = pp != nil ? Int(pp!)! : 0
        if atom_idx >= 1 && atom_idx <= mol.numAtoms() {
            guard let _endatom = mol.getAtom(atom_idx) else { fatalError("ERROR: Could not find atom") }
        }
        // Was a start atom specified?
        pp = _pconv.isOption("f")
        atom_idx  = pp != nil ? Int(pp!)! : 0
        if atom_idx >= 1 && atom_idx <= mol.numAtoms() {
            guard let _startatom = mol.getAtom(atom_idx) else { fatalError("ERROR: Could not find atom") }
        }
        
        // Was an atom ordering specified?
        var ppo = Optional(options.ordering)
        var s_atom_order: [String] = []
        var atom_order: [Int] = []
        if ppo != nil {
            s_atom_order = ppo!.components(separatedBy: "-()")
            if s_atom_order.count != mol.numHeavyAtoms() {
                ppo = nil
            } else {
                for s_atom in s_atom_order {
                    atom_order.append(Int(s_atom)!)
                }
                atom_idx = atom_order[0]
                if atom_idx >= 1 && atom_idx <= mol.numAtoms() {
                    guard let _startatom = mol.getAtom(atom_idx) else { fatalError("ERROR: Could not find atom") }
                }
            }
        }
        // Was Universal SMILES requested?
        var universal_smiles: Bool = _pconv.isOption("U")
        if universal_smiles {
            let parsedOkay = parseInChI(mol, &atom_order)
            if !parsedOkay {
                universal_smiles = false
            }
        }
        if _canonicalOutput {
            // Find the (dis)connected fragments.
            var visited = Bitset()
            var fragments: [Bitset] = []
            for i in 0..<mol.numAtoms() {
                if !frag_atoms.contains(i+1) || visited.contains(i+1) { continue }
                if let atom = mol.getAtom(i+1) {
                    fragments.append(getFragment(atom, frag_atoms))
                    visited |= fragments.last!
                } else {
                    fatalError("ERROR: Could not find atom")
                }
            }
            // Determine symmetry classes for each disconnected fragment separately
            symmetry_classes = Array(repeating: 0, count: mol.numAtoms())
            for i in 0..<fragments.count {
                let gs = MKGraphSym(mol, &fragments[i])
                var tmp: [Ref] = []
                gs.getSymmetry(&tmp)
                for j in 0..<mol.numAtoms() {
                    if fragments[i].contains(j+1) {
                        symmetry_classes[j] = tmp[j]
                    }
                }
            }
            // Was a canonicalization timeout given?
            var maxSeconds: Int = 5
            let timeoutString: String? = _pconv.isOption("T")
            if timeoutString != nil {
                if let seconds = Int(timeoutString!) {
                    maxSeconds = seconds
                } else {
                    print("ERROR: Canonicalization timeout should be a number")
                    maxSeconds = 5
                }
            }
            var symclasses = symmetry_classes.map { UInt($0.intValue) }
            var canorder = canonical_order.map { UInt($0.intValue) }
            canonicalLabels(mol, &symclasses, &canorder, frag_atoms, maxSeconds)
        } else {
            if _pconv.isOption("C") { // "C" == "anti-canonical form"
                randomLabels(mol, frag_atoms, &symmetry_classes, &canonical_order)
            } else if ppo != nil || universal_smiles { // user-specified or InChI canonical labels
                canonical_order = Array(repeating: 0, count: mol.numAtoms())
                symmetry_classes = Array(repeating: 0, count: mol.numAtoms())
                var idx = 3 // Start the labels at 3 (to leave space for special values 0, 1 and 2)
                for i in 0..<atom_order.count {
                    if canonical_order[atom_order[i] - 1] == 0 { // Ignore ring closures (for "U")
                        canonical_order[atom_order[i] - 1] = Ref(integerLiteral: idx)
                        symmetry_classes[atom_order[i] - 1] = Ref(integerLiteral: idx)
                        idx += 1
                    }
                }
                for i in 0..<canonical_order.count {
                    if canonical_order[i] == 0 { // Explicit hydrogens
                        if let atom = mol.getAtom(i+1) {
                            if atom.getAtomicNum() == MKElements.Hydrogen.atomicNum && atom.getIsotope() != 0 { // [2H] or [3H]
                                canonical_order[i] = .Ref(Int(atom.getIsotope()) - 1) // i.e. 1 or 2
                                symmetry_classes[i] = canonical_order[i]
                            }
                        } else {
                            fatalError("ERROR: Could not find atom")
                        }
                    }
                }
            } else {
                standardLabels(mol, frag_atoms, &symmetry_classes, &canonical_order)
            }
        }
        
        // OUTER LOOP: Handles dot-disconnected structures and reactions.  Finds the
        // lowest unmarked canorder atom in the current reaction role, and starts there
        // to generate a SMILES.
        // Repeats until no atoms remain unmarked.

        var new_rxn_role = false // flag to indicate whether we have started a new reaction role
        var isrxn: Bool = mol.isReaction()
        let rxn = MKReactionFacade(mol)
        var rxnrole: UInt = 1 
        while true {
            if _pconv.isOption("R") {
                _bcdigit = 0 // Reset the bond closure index for each disconnected component 
            }

            // It happens that the lowest canonically-numbered atom is usually a place to start the 
            // canonical SMILES 
            var root_atom: MKAtom? = nil 
            var lowest_canorder = 999999
            // If we specified a startatom_idx & it's in this fragment, use it to start the fragment
            if (_startatom != nil) {
                if !_uatoms[_startatom!.getIdx()] && frag_atoms.contains(_startatom!.getIdx()) &&
                    (!isrxn || rxn.getRole(atom: _startatom!).rawValue == rxnrole) {
                    root_atom = _startatom
                }
            }
            if root_atom == nil {
                for atom in mol.getAtomIterator() {
                    let idx = atom.getIdx()
                    if ( //atom.getAtomicNum() != MKElements.Hydrogen       // don't start with a hydrogen
                        !_uatoms[idx]          // skip atoms already used (for fragments)
                        && frag_atoms.contains(idx)// skip atoms not in this fragment
                        && (!isrxn || rxn.getRole(atom: atom).rawValue == rxnrole) // skip atoms not in this rxn role
                        //&& !atom->IsChiral()    // don't use chiral atoms as root node
                        && canonical_order[idx-1].intValue < lowest_canorder) {
                        root_atom = atom
                        lowest_canorder = canonical_order[idx-1].intValue
                    }
                }
                // For Inchified or Universal SMILES, if the start atom is an [O-] attached to atom X, choose any =O attached to X instead.
                //          Ditto for [S-] and =S.
                if (_pconv.isOption("I") || universal_smiles)
                    && (root_atom != nil) && root_atom!.getFormalCharge() == -1  && root_atom!.getExplicitDegree() == 1
                    && root_atom!.hasSingleBond() && (root_atom!.getAtomicNum() == MKElements.Oxygen.atomicNum || root_atom!.getAtomicNum() == MKElements.Sulfur.atomicNum) {
                    if let central_atom = root_atom?.getNbrAtomIterator()?.next() {
                        for nbr in central_atom.getNbrAtomIterator()! {
                            if (root_atom == nbr) { continue }
                            if (nbr.getAtomicNum() == root_atom!.getAtomicNum() && nbr.getExplicitDegree() == 1 && nbr.hasDoubleBond()) {
                                root_atom = nbr
                                break
                            }
                        }
                    }
                }
            }
            
            // No atom found?  If it's not a reaction, then we've done all fragments.
            // If it is, then increment the rxn role and try again.
            if root_atom == nil {
                if mol.isReaction() {
                    rxnrole += 1
                    if rxnrole == 4 {
                        break
                    }
                    buffer += ">"
                    new_rxn_role = true
                    continue
                } else {
                    break
                }
            }
            // Clear out closures in case structure is dot disconnected
            //      _atmorder.clear();
            _vopen.removeAll()
            // Dot disconnected structure or new rxn role?
            if new_rxn_role {
                new_rxn_role = false
            } else if !buffer.isEmpty {
                buffer += "."
            }
            // root = new OBCanSmiNode (root_atom);
            let root = MKCanSmiNode(atom: root_atom!)
            let symclasses = symmetry_classes.map { $0.intValue }
            let canorder = canonical_order.map { $0.intValue }
            buildCanonTree(mol, &frag_atoms, canorder, root)
            toCansmilesString(root, &buffer, &frag_atoms, symclasses, canorder)

        }
    }
    
    
    /***************************************************************************
     * FUNCTION: GetTetrahedralStereo
     *
     * DESCRIPTION:
     *       If the atom is chiral, return either "@" or "@@". Otherwise 0.
     ***************************************************************************/
    func getTetrahedralStereo(_ node: MKCanSmiNode, _ chiral_neighbors: [MKAtom?], _ symmetry_classes: [Int]) -> String {
        // If not enough chiral neighbors were passed in, we're done
        if chiral_neighbors.count < 4 { return "" }
        // If atom is not a tetrahedral center, we're done
        let atom = node.getAtom()
        guard let ts = _stereoFacade.getTetrahedralStereo(atom.getId().intValue) else { return "" }
        // get the Config struct defining the stereochemistry
        let atomConfig = ts.getConfig()
        // Unspecified or unknown stereochemistry
        if !atomConfig.specified || (atomConfig.specified && atomConfig.winding == .UnknownWinding) { return "" }
        // create a Config struct with the chiral_neighbors in canonical output order
        var canonRefs: Refs = []
        
        for atomit in chiral_neighbors[1...] {
            if atomit != nil {
                canonRefs.append(atomit!.getId())
            } else {
                canonRefs.append(.ImplicitRef)
            }
        }
        
        var canConfig: MKTetrahedralStereo.Config = MKTetrahedralStereo.Config()
        canConfig.center = atom.getId()
        
        if chiral_neighbors[0] != nil {
            canConfig.from_or_towrds = .from(chiral_neighbors[0]!.getId())
        } else { // Handle a chiral lone pair, represented by a NULL OBAtom* in chiral_neighbors
            canConfig.from_or_towrds = .from(.ImplicitRef)
        }
        
        canConfig.refs = canonRefs
        // config is clockwise
        if atomConfig == canConfig {
            return "@@"
        } else {
            return "@"
        }
    }
    
    func getSquarePlanarStereo(_ node: MKCanSmiNode, _ chiral_neighbors: [MKAtom?], _ symmetry_classes: [Int]) -> String {
        // If not enough chiral neighbors were passed in, we're done
        if chiral_neighbors.count < 4 { return "" }
        if chiral_neighbors.contains(where: { $0 == nil }) { return "" }
        let atom = node.getAtom()
        
        // If atom is not a square-planar center, we're done
        guard let sp = _stereoFacade.getSquarePlanarStereo(atom.getId().intValue) else {
            return ""
        }
        
        // get the Config struct defining the stereochemistry
        let atomConfig = sp.getConfig()
        if !atomConfig.specified {
            return ""
        }
        // create a Config struct with the chiral_neighbors in canonical output order
        let canonRefs: Refs = MKStereo.makeRefs(chiral_neighbors[0]!.getId(), chiral_neighbors[1]!.getId(), chiral_neighbors[2]!.getId(), chiral_neighbors[3]!.getId())
        var canConfig: MKSquarePlanarStereo.Config = MKSquarePlanarStereo.Config()
        canConfig.center = atom.getId()
        canConfig.refs = canonRefs
        
        // canConfig is U shape
        if atomConfig == canConfig {
            return "@SP1"
        }
        canConfig.shape = .Shape4
        if atomConfig == canConfig {
            return "@SP2"
        }
        canConfig.shape = .ShapeZ
        if atomConfig == canConfig {
            return "@SP3"
        }
        return ""
    }
    
    /***************************************************************************
     * FUNCTION: GetSmilesElement
     *
     * DESCRIPTION:
     *       Writes the symbol for an atom, e.g. "C" or "[NH2]" or "[C@@H]".
     *
     * RETURNS: true (always)
     ***************************************************************************/
    @discardableResult
    func getSmilesElement(_ node: MKCanSmiNode, _ chiral_neighbors: [MKAtom?], _ symmetry_classes: [Int], _ buffer: inout String) -> Bool {
        var symbol: String = ""
        
        var bracketElement: Bool = false
        //        var normalValence: Bool = true
        //        var writeExplicitHydrogens: Bool = false
        
        let atom = node.getAtom()
        let element = atom.getAtomicNum()
        
        // Handle SMILES Valence model, and explicit and implicit hydrogens
        if isOutsideOrganicSubset(element) {
            bracketElement = true
        }
        
        var numExplicitHsToSuppress: Int = 0
        // Don't suppress any explicit Hs attached if the atom is an H itself (e.g. [H][H]) or -xh was specified
        if atom.getAtomicNum() != MKElements.Hydrogen.atomicNum && !options.showexplicitH {
            for nbr in atom.getNbrAtomIterator()! {
                if nbr.getAtomicNum() == MKElements.Hydrogen.atomicNum && (!options.isomeric || nbr.getIsotope() == 0) && nbr.getExplicitDegree() == 1 && nbr.getFormalCharge() == 0 && (!options.showatomclass || nbr.getData("Atom Class") == nil) {
                    numExplicitHsToSuppress += 1
                }
            }
        }
        
        var numImplicitHs = 0
        if options.smarts {
            if numExplicitHsToSuppress > 0 {
                bracketElement = true
                numImplicitHs = numExplicitHsToSuppress
            }
        } else {
            numImplicitHs = Int(atom.getImplicitHCount()) + numExplicitHsToSuppress
            if !bracketElement {
                if element == 0 { // asterisk is always hypervalent but we don't bracket it unless has Hs
                    if numImplicitHs > 0 {
                        bracketElement = true
                    }
                } else {
                    let bosum = Int(atom.getExplicitValence()) - numExplicitHsToSuppress
                    let implicitValence = smilesValence(element, bosum, false)
                    let defaultNumImplicitHs = implicitValence - bosum
                    if implicitValence == 0 || // hypervalent
                        numImplicitHs != defaultNumImplicitHs // undervalent
                        || (!options.kekulesmi && element != 6 && atom.isAromatic() && numImplicitHs != 0) // aromatic nitrogen/phosphorus
                    {
                        bracketElement = true
                    }
                }
            }
        }
        
        if atom.getFormalCharge() != 0 ||  // charged elements
            (options.isomeric && atom.getIsotope() != 0) ||  // isotopes
            (options.showatomclass && atom.getData("Atom Class") != nil) { // If the molecule has Atom Class data and -xa option set and atom has data
            bracketElement = true
        }
        
        var stereo: String? = nil
        if getSmilesValence(atom) > 2 && options.isomeric {
            stereo = getTetrahedralStereo(node, chiral_neighbors, symmetry_classes)
            if stereo == nil {
                stereo = getSquarePlanarStereo(node, chiral_neighbors, symmetry_classes)
            }
        }
        
        if stereo != nil {
            bracketElement = true
        }
        
        if !bracketElement {
            // ordinary non-bracketed element
            if element != 0 {
                symbol = MKElements.getSymbol(atom.getAtomicNum())
                if (!options.kekulesmi && atom.isAromatic()) || // aromatic atom
                    (atom.getSpinMultiplicity() != 0 && _pconv.isOption("r"))
                { // radical centers lowercase if r option is set
                    buffer += symbol[0] + " "
                    if symbol[1] != " " || symbol[1] != "" {
                        buffer += symbol[1]
                    }
                } else {
                    buffer += symbol
                }
            } else {
                // Atomic number zero - either '*' or an external atom
                var external: Bool = false
                if let externalBonds: MKExternalBondData = atom.getParent()!.getData("extBonds") as? MKExternalBondData {
                    for externalBond in externalBonds.vexbonds {
                        if externalBond.atom == atom {
                            external = true
                            buffer += "&"
                            if let bond = externalBond.bond {
                                if bond.getBondOrder() == 2 && !bond.isAromatic() {
                                    buffer += "="
                                }
                                if bond.getBondOrder() == 2 && bond.isAromatic() {
                                    buffer += ":"
                                }
                                if bond.getBondOrder() == 3 {
                                    buffer += "#"
                                }
                                if bond.getBondOrder() == 4 {
                                    buffer += "$"
                                }
                                buffer += String(externalBond.idx)
                                break
                            }
                        }
                    }
                    if !external {
                        buffer += "*"
                    }
                }
            }
            return true
        }
        
        // Bracket atoms, e.g. [Pb], [OH-], [C@]
        buffer += "["
        let iso = atom.getIsotope()
        if options.isomeric && iso != 0 {
            if iso >= 10000 { // max 4 characaters
                print("Warning: Isotope value \(iso) is too large for SMILES")
            } else {
                buffer += String(iso)
            }
        }
        if atom.getAtomicNum() == 0 {
            buffer += "*"
        } else {
            if atom.getAtomicNum() == MKElements.Hydrogen.atomicNum && options.smarts {
                buffer += "#1"
            } else {
                let elem = atom.getAtomicNum()
                let symbol = MKElements.getSymbol(elem)
                if symbol == "" {
                    buffer += "#\(elem)"
                } else if (!options.kekulesmi && atom.isAromatic()) { // aromatic atom
                    buffer += symbol[0] + " "
                    if symbol[1] != " " || symbol[1] != "" {
                        buffer += symbol[1]
                    }
                } else {
                    buffer += symbol
                }
            }
        }
        
        // If chiral, append '@' or '@@'...unless we're creating a SMARTS ("s") and it's @H or @@H
        if stereo != nil && !(options.smarts && atom.getImplicitHCount() > 0){
            buffer += stereo!
        }
        // Add extra hydrogens.
        
        var hcount = numImplicitHs
        if hcount > 0 && (atom == _endatom || (atom == _startatom && (options.ordering == ""))) { // Leave a free valence for attachment
            hcount -= 1
        }
        if hcount > 0 {
            if options.smarts && stereo == nil {
                for i in 0..<hcount {
                    buffer += "!H\(i)"
                }
            } else {
                buffer += "H"
                if hcount > 1 {
                    buffer += String(hcount)
                }
            }
        }
        
        // Append charge to the end
        
        let charge = atom.getFormalCharge()
        if charge != 0 {
            if charge > 0 {
                buffer += "+"
            } else {
                buffer += "-"
            }
            if abs(charge) > 1 {
                buffer += String(abs(charge))
            }
        }
        
        //atom class e.g. [C:2]
        
        if options.showatomclass {
            if let data = atom.getData("Atom Class") as? MKPairData<Int> {
                if let ac = data.getValue() {
                    if ac >= 0 {
                        buffer += ":"
                        buffer += String(ac)
                    }
                }
            }
        }
        
        buffer += "]"
        
        return true;
    }
    
    /***************************************************************************
     * FUNCTION: GetSmilesValence
     *
     * DESCRIPTION:
     *       This is like GetHvyDegree(), but it returns the "valence" of an
     *       atom as it appears in the SMILES string.  In particular, hydrogens
     *       count if they will appear explicitly -- see IsSuppressedHydrogen()
     *       above.
     ***************************************************************************/
    func getSmilesValence(_ atom: MKAtom) -> Int {
        if atom.getAtomicNum() == MKElements.Hydrogen.atomicNum {
            return atom.getExplicitDegree()
        }
        if options.showexplicitH {
            return atom.getExplicitDegree()
        }
        return atom.getNbrAtomIterator()!.map({($0.getAtomicNum() != MKElements.Hydrogen.atomicNum ||
                                                $0.getIsotope() != 0 ||
                                                $0.getExplicitDegree() != 1) ? 1: 0 }).reduce(0, +)
    }
    
    /***************************************************************************
     * FUNCTION: GetUnusedIndex
     *
     * DESCRIPTION:
     *       Returns the next available bond-closure index for a SMILES.
     *
     *       You could just do this sequentially, not reusing bond-closure
     *       digits, thus (chosen by Option("R")):
     *
     *               c1cc2ccccc2cc1          napthalene
     *               c1ccccc1c2ccccc2        biphenyl
     *
     *       But molecules with more than ten rings, this requires the use of
     *       two-digit ring closures (like c1ccccc1C...c%11ccccc%11).  To help
     *       avoid digit reuse, this finds the lowest digit that's not currently
     *       "open", thus
     *
     *               c1cc2ccccc2cc1          napthalene (same)
     *               c1ccccc1c1ccccc1        biphenyl (reuses "1")
     *
     ***************************************************************************/
    
    func getUnusedIndex() -> Int {
        if _pconv.isOption("R") {
            // Keep incrementing the bond closure digits (for each connected component)
            _bcdigit += 1
            return _bcdigit
        }
        // TODO: make sure this works as intended. The original code did an inefficient search,
        // but maybe that was because the numbers can appear out of order
        var idx = _vopen.map { $0.ringdigit }.max() ?? 0 // if _vopen is empty, use 0
        idx += 1
        return idx
    }
    
    /***************************************************************************
     * FUNCTION: GetCanonClosureDigits
     *
     * DESCRIPTION:
     *       Given an atom, returns the ring-closure digits for that atom, in
     *       the form of a vector of digit/OBBond* pair.  Some of the digits may
     *       be for newly-opened rings (the matching digit occurs later in the
     *       SMILES string), and some may be for closing rings (the matching
     *       digit occurred earlier in the string).
     *
     *       Canonicalization requires that atoms with more than one digit
     *       have the digits assigned in a canonical fashion.  For example,
     *       the SMILES  "CC12(NCCC2)CCC1" and "CC12(NCCC1)CCC2" are the
     *       same molecule; we need to assign the digits of the first "C12"
     *       such that it always comes out one way or the other.
     *
     *       This needs to find closing bonds (ring bonds already assigned a
     *       digit) and opening bonds (ring bonds not encountered yet).
     *
     *    Closing Bonds:
     *       This is easy: open bonds are already stored in the _vopen vector,
     *       in canonical order.  Just find open bonds to this atom and copy
     *       them from _vopen to our return vector.
     *
     *    Opening Bonds:
     *       This function looks through the bonds for this atoms and finds
     *       any that aren't on the _ubonds "used" list, (and also are non-H
     *       and are in this fragment).  Any such bonds must be ring-closure
     *       bonds.  If there is more than one, they are ordered by the
     *       canonical order of the bonds' neighbor atoms; that is, the bond
     *       to the lowest canonical-ordered neighbor is assigned the first
     *       available number, and upwards in neighbor-atom canonical order.
     ***************************************************************************/
    func getCanonClosureDigits(_ atom: MKAtom, _ frag_atoms: inout Bitset, _ canonical_order: [Int]) -> [MKBondClosureInfo] {
        var vp_closures: [MKBondClosureInfo] = []
        var vbonds: [MKBond] = []
        
        // Find new ring-closure bonds for this atom
        var setNewBond: Bool = false
        for bond1 in atom.getBondIterator()! {
            // Is this a ring-closure neighbor?
            if _ubonds.contains(Int(bond1.getIdx())) { continue }
            let nbr1 = bond1.getNbrAtom(atom)
            // Skip hydrogens before checking canonical_order
            // PR#1999348
            if ((nbr1.getAtomicNum() == MKElements.Hydrogen.atomicNum && isSuppressedHydrogen(nbr1)) ||
                (!frag_atoms.contains(nbr1.getIdx()))) { continue }
            let nbr1_canorder = canonical_order[nbr1.getIdx() - 1]
            
            // Insert into the bond-vector in canonical order (by neightbor atom order)
            for (i, bi) in vbonds.enumerated() {
                let nbr2 = bi.getNbrAtom(atom)
                let nbr2_canorder = canonical_order[nbr2.getIdx() - 1]
                if nbr1_canorder < nbr2_canorder {
                    vbonds.insert(bond1, at: i)
                    setNewBond = true
                    break
                }
            }
            if !setNewBond {
                vbonds.append(bond1)      // highest one (or first one) - append to end
            }
            setNewBond = false
        }
        
        // If we found new open bonds, assign a bond-closure digits to each one,
        // add it to _vopen, and add it to the return vector.
        for bond1 in vbonds {
            _ubonds.add(Int(bond1.getIdx()))
            let digit = getUnusedIndex()
            // let bo = bond1.isAromatic() ? 1 : bond1.getBondOrder() // CJ: why was this line added?  bo is never used?
            let newAdd = MKBondClosureInfo(toatom: bond1.getNbrAtom(atom), fromatom: atom, bond: bond1, ringdigit: digit, is_open: 1) // interpret 1 as True
            _vopen.append(newAdd)
            vp_closures.append(newAdd)
        }
        
        // Now look through the list of open closure-bonds and find any to this
        // atom (but watch out for the ones we just added).  For each one found,
        // add it to the return vector, and erase it from _vopen.
        
        if !_vopen.isEmpty {
            var j = 0
            while j < _vopen.count {
                if _vopen[j].toatom == atom {
                    let bci = _vopen[j]
                    _vopen.remove(at: j)                         // take bond off "open" list
                    bci.is_open = 0                              // mark it "closed" false ===0 in this case
                    vp_closures.append(bci)                      // and add it to this atom's list
                    j = 0                                        // reset iterator
                }
                else {
                    j += 1
                }
            }
        }
        
        return vp_closures
    }
    
    /***************************************************************************
     * FUNCTION: IsSuppressedHydrogen
     *
     * DESCRIPTION:
     *       For a hydrogen atom, returns TRUE if the atom is not [2H] or [3H], only
     *       has one bond, and is not bonded to another hydrogen.
     *
     *       NOTE: Return value is nonsensical if you pass it a non-hydrogen
     *       atom.  Presumably, you're calling this because you've found a
     *       hydrogen and want to know if it goes in the SMILES.
     ***************************************************************************/
    func isSuppressedHydrogen(_ atom: MKAtom) -> Bool {
        if atom.getIsotope() != 0 { return false } // Deuterium or Tritium
        if atom.getExplicitDegree() != 1 { return false } // Not exactly one bond
        for nbr in atom.getNbrAtomIterator()! {
            if nbr.getAtomicNum() == 1 { return false } // Bonded to another hydrogen
        }
        return true
    }
    
    
    /***************************************************************************
    * FUNCTION: ToCansmilesString
    *
    * DESCRIPTION:
    *       Recursively writes the canonical SMILES string to a buffer.  Writes
    *       this node, then selects each of the child nodes (in canonical
    *       order) and writes them.
    *
    *       Chirality is the tricky bit here.  Before we can write out a chiral
    *       atom, we have to "look ahead" to determine the order in which the
    *       neighbor atoms are/will be written.
    *
    *       The SMILES language defines the order-of-appearance of a ring-closure
    *       bond as the position of the digit, in the SMILES, not the actual atom.
    *       For example, the fragments N[C@H](C)Br, and N[C@H]1(Br)CCCC1 have
    *       the same chiral center, because the "1" in the second one is a "stand
    *       in" for the "C" in the first, even though the actual carbon atom appears
    *       after the Bromine atom in the second string.
    ***************************************************************************/
    func toCansmilesString(_ node: MKCanSmiNode, _ buffer: inout String, _ frag_atoms: inout Bitset, _ symmetry_classes: [Int], _ canonical_order: [Int]) {
        let atom = node.getAtom()
        var chiral_neighbors: [MKAtom?] = [] 

        // Get the ring-closure digits in canonical order.  We'll use these in
        // two places: First, for figuring out chirality, then later for writing
        // the actual ring-closure digits to the string.

        let vclose_bonds = getCanonClosureDigits(atom, &frag_atoms, canonical_order)

        // First thing: Figure out chirality.  We start by creating a vector of the neighbors
        // in the order in which they'll appear in the canonical SMILES string.  This is more
        // complex than you'd guess because of implicit/explicit H and ring-closure digits.

        // Don't include chiral symbol on _endatom or _startatom.
        // Otherwise, we end up with C[C@@H](Br)(Cl), where the C has 4 neighbours already
        // and we cannot concatenate another SMILES string without creating a 5-valent C.

        let is_chiral: Bool = atomIsChiral(atom)
        if is_chiral && atom != _endatom && atom != _startatom {
            // If there's a parent node, it's the first atom in the ordered neighbor-vector
            // used for chirality.
            if let parent = node.getParent() {
                chiral_neighbors.append(parent)
            }
            // Next for chirality order will be hydrogen -- since it occurs
            // inside the atom's [] brackets, it's always before other neighbors.
            //
            // Note that we check the regular neighbor list, NOT the canonical
            // SMILES tree, because hydrogens normally aren't part of the canonical
            // SMILES, but we still need them to figure out chirality.
            //
            // There are two cases: it's explicit in the OBMol object but should be
            // written inside the brackets, i.e. "[C@H]", or it is explicit and
            // must be outside the brackets, such as for deuterium.  (A hydrogen
            // that will appear explicitly in the SMILES as a separate atom is
            // treated like any other atom when calculating the chirality.)

            if !options.showexplicitH {
                for nbr in atom.getNbrAtomIterator()! {
                    if nbr.getAtomicNum() == MKElements.Hydrogen.atomicNum && isSuppressedHydrogen(nbr) {
                        chiral_neighbors.append(nbr)
                        break // quit loop: only be one H if atom is chiral
                    }
                }
            }
            // Handle implict H by adding a NULL OBAtom*
            if atom.getImplicitHCount() == 1 {
                chiral_neighbors.append(nil)
            }

            // Ok, done with H. Now we need to consider the case where there is a chiral
            // lone pair. If it exists (and we won't know for sure until we've counted up
            // all the neighbours) it will go in here
            let lonepair_location = chiral_neighbors.count
            // Ok, done with all that. Next in the SMILES will be the ring-closure characters.
            // So we need to find the corresponding atoms and add them to the list.
            // (We got the canonical ring-closure list earlier.)
            if !vclose_bonds.isEmpty {
                for i in 0..<vclose_bonds.count {
                    let bond = vclose_bonds[i].bond
                    let nbr = bond.getNbrAtom(atom)
                    chiral_neighbors.append(nbr)
                }
            }

            // Finally, add the "regular" neighbors, the "child" nodes in the
            // canonical-SMILES tree, to the chiral-neighbors list.
            for child in 0..<node.size() {
                let nbr = node.getChildAtom(child)
                chiral_neighbors.append(nbr)
            }
            // Handle a chiral lone-pair on a sulfur, by inserting a NULL OBAtom* at the
            // appropriate location
            if chiral_neighbors.count == 3 && atom.getAtomicNum() == MKElements.Sulfur.atomicNum { // Handle sulfur
                chiral_neighbors.insert(nil, at: lonepair_location)
            }
        }

        // Write the current atom to the string
        // GetSmilesElement(node, chiral_neighbors, symmetry_classes, buffer);
        getSmilesElement(node, chiral_neighbors, symmetry_classes, &buffer)

        _atmorder.append(atom.getIdx()) //store the atom ordering

        // Write ring-closure digits
        if !vclose_bonds.isEmpty {
            for bci in vclose_bonds {
                if bci.is_open == 0 {
                    // Ring closure
                    var bs = ""
                    // Only get symbol for ring closures on the dbl bond
                    if hasStereoDblBond(bci.bond, node.getAtom()) {
                        bs = getCisTransBondSymbol(bci.bond, node)
                    }
                    if !bs.isEmpty {
                        // append "/" or "\"
                        buffer += bs
                    } else {
                        switch bci.bond.getBondOrder() {
                        case 1:
                            if !bci.bond.isAromatic() && bci.bond.isInRing() && bci.bond.getBeginAtom().isAromatic() && bci.bond.getEndAtom().isAromatic() {
                                buffer += "-"
                            }
                        case 2:
                            if options.kekulesmi || !bci.bond.isAromatic() {
                                buffer += "="
                            }
                        case 3:
                            buffer += "#"
                        case 4:
                            buffer += "$"
                        default:
                            break
                        }
                    }
                } else {
                    // Ring opening
                    var bs = ""
                    // Only get symbol for ring openings on the dbl bond
                    if !hasStereoDblBond(bci.bond, bci.bond.getNbrAtom(node.getAtom())) {
                        bs = getCisTransBondSymbol(bci.bond, node)
                    }
                    if !bs.isEmpty {
                        // append "/" or "\"
                        buffer += bs
                    }
                }
                if bci.ringdigit > 9 {
                    buffer += "%"
                    if bci.ringdigit > 99 {
                        buffer += "("
                    }
                    buffer += String(bci.ringdigit)
                    if bci.ringdigit > 99 {
                        buffer += ")"
                    }
                } else {
                    buffer += String(bci.ringdigit)
                }
            }
        }

        // Write child bonds, then recursively follow paths to child nodes
        // to print the SMILES for each child branch.
        for i in 0..<node.size() {
            let bond = node.getChildBond(i)
            if i+1 < node.size() || node.getAtom() == _endatom {
                buffer += "("
            }
            switch bond.getBondOrder() {
            case 1:
                let cc = getCisTransBondSymbol(bond, node)
                if !cc.isEmpty {
                    buffer += cc
                } else {
                    // Write a single bond symbol if not aromatic but end atoms are both aromatic
                    // This will speed up reading as it will avoid ring perception around line 563 (bond->IsInRing())
                    // TODO: Consider making the test for IsInRing() an option
                    if !bond.isAromatic() && bond.isInRing() && bond.getBeginAtom().isAromatic() && bond.getEndAtom().isAromatic() {
                        buffer += "-"
                    }
                }
            case 2:
                if options.kekulesmi || !bond.isAromatic() {
                    buffer += "="
                }
            case 3:
                buffer += "#"
            case 4:
                buffer += "$"
            default:
                break
            }
            
            toCansmilesString(node.getChildNode(i), &buffer, &frag_atoms, symmetry_classes, canonical_order)
            
            if i+1 < node.size() || node.getAtom() == _endatom {
                buffer += ")"
            }
        }
    }
    
    func hasStereoDblBond(_ bond: MKBond?, _ atom: MKAtom?) -> Bool {
        // This is a helper function for determining whether to
        // consider writing a cis/trans bond symbol for bond closures.
        // Returns TRUE only if the atom is connected to the cis/trans
        // double bond. To handle the case of conjugated bonds, one must
        // remember that the ring opening preceded the closure, so if the
        // ring opening bond was on a stereocenter, it got the symbol already.
        
        guard bond != nil else { return false }
        guard atom != nil else { return false }
        
        let nbr_atom = bond!.getNbrAtom(atom!)
        var stereo_dbl = false
        if atom!.hasDoubleBond() {
            stereo_dbl = true
            if nbr_atom.hasDoubleBond() {
                // Check whether the nbr_atom is a begin or end in any CisTransStereo. If so,
                // then the ring opening already had the symbol.
                for ct in _cistrans {
                    let cfg = ct.getConfig()
                    if nbr_atom.getId() == cfg.begin || nbr_atom.getId() == cfg.end {
                        // I don't think I need to check whether it has a bond with atom
                        stereo_dbl = false
                        break
                    }
                }
            }
        }
        
        return stereo_dbl
    }
    
    //! Adaptation of OBMol::FindChildren to allow a vector of OBAtoms to be passed in
    //  MOVE THIS TO OBMOL FOR 2.4
    //         TODO: substitute this with the moleucule classes implementation :: mol.findChildren
    func myFindChildren(_ mol: MKMol, _ children: inout [MKAtom], _ seen: Bitset, _ end: MKAtom) {
        var curr = Bitset()
        var next = Bitset()
        var used = Bitset(seen)
        used.add(end.getIdx())
        curr.add(end.getIdx())
        children.removeAll()
        while true {
            next.removeAll()
            for i in curr {
                guard let atom = mol.getAtom(i) else { fatalError("Cannot find atom in molecule : myFindChildren") }
                for nbr in atom.getNbrAtomIterator()! {
                    if !used[nbr.getIdx()] {
                        children.append(nbr)
                        next.add(nbr.getIdx())
                        used.add(nbr.getIdx())
                    }
                }
            }
            if next.isEmpty() {
                break
            }
            curr = next
        }
    }
    
    func getOutputOrder(_ outorder: inout String) {
        // std::vector<int>::iterator it = _atmorder.begin();
        // if (it != _atmorder.end()) {
        // char tmp[15];
        // snprintf(tmp, 15, "%d", *it);
        // outorder += tmp;
        // ++it;
        // for (; it != _atmorder.end(); ++it) {
        //     snprintf(tmp, 15, "%d", *it);
        //     outorder += ' ';
        //     outorder += tmp;
        // }
        // }
        if !_atmorder.isEmpty {
            outorder += String(_atmorder[0])
            for i in 1..<_atmorder.count {
                outorder += " "
                outorder += String(_atmorder[i])
            }
        }
    }

    // Returns canonical label order
    func parseInChI(_ mol: MKMol, _ atom_order: inout [Int]) -> Bool {
        let MolConv = MKConversion()
        MolConv.setOutFormat("InChI")
        MolConv.setAuxConv(nil) //temporary until a proper OBConversion copy constructor written
        let newstream = OutputStringStream()
        MolConv.setOutStream(newstream)
        // I'm sure there's a better way of preventing InChI warning output
        MolConv.addOption("w", .OUTOPTIONS)
        MolConv.addOption("a", .OUTOPTIONS)
        MolConv.addOption("X", .OUTOPTIONS, "RecMet FixedH")
        MolConv.write(mol)


        let splitlines = newstream.string.components(separatedBy: .whitespacesAndNewlines)

        var split: [String]
        var split_aux: [String]
        var aux_part: String
        
        let rm_start = splitlines[0].firstIndex(of: Character("/r"))

        if rm_start == nil {
            split = splitlines[0].components(separatedBy: "/")
            aux_part = splitlines[1] // Use the normal labels
        } else {
            let tmp = splitlines[0].substring(fromIndex: rm_start!)
            split = tmp.components(separatedBy: "/")
            split.insert("", at: 0)
            let rm_start_b = splitlines[1].firstIndex(of: Character("/R:"))
            aux_part = splitlines[1].substring(fromIndex: rm_start_b!) // Use the reconnected metal labels
        }

        split_aux = aux_part.components(separatedBy: "/")

        // Parse the canonical labels

        var canonical_labels: [[Int]] = []
        var s_components: [String] = []
        var s_atoms: [String] = []

        var tmp = split_aux[2].substring(fromIndex: 2)
        s_components = tmp.components(separatedBy: ";")
        for s_component in s_components {
            s_atoms = s_component.components(separatedBy: ",")
            var atoms: [Int] = []
            for s_atom in s_atoms {
                atoms.append(Int(s_atom)!)
            }
            canonical_labels.append(atoms)
        }

        // Adjust the canonical labels if necessary using a /F section
        if let f_start = aux_part.firstIndex(of: Character("/F:")) {
            tmp = aux_part.substring(fromIndex: f_start.utf16Offset(in: aux_part) + 3)
            split_aux = tmp.components(separatedBy: "/")
            s_components = split_aux[0].components(separatedBy: ";")
            var new_canonical_labels: [[Int]] = []
            var total = 0
            for it in s_components {
                // e.g. "1,2,3;2m" means replace the first component by "1,2,3"
                //                       but keep the next two unchanged
                if it.last == "m" {
                    var mult: Int
                    if it.count == 1 {
                        mult = 1
                    } else {
                        mult = Int(it.substring(toIndex: it.length - 1))!
                    }
                    // TODO: make sure this doesn't need to be inclusive of total+multi
                    new_canonical_labels.append(contentsOf: canonical_labels[total..<total+mult])
                    total += mult
                } else {
                    s_atoms = it.components(separatedBy: ",")
                    var atoms: [Int] = []
                    for itb in s_atoms {
                        atoms.append(Int(itb)!)
                    }
                    new_canonical_labels.append(atoms)
                    total += 1
                }
            }
            canonical_labels = new_canonical_labels
        }

        // Flatten the canonical_labels
        for it in canonical_labels {
            atom_order.append(contentsOf: it)
        }
        return true
    }
    
    
}

// Do we need to write out a bond symbol for this bond?
// No - if it's aromatic
// Otherwise, yes if the bond order is not 1
// If the bond order *is* 1, then only if the bond is in a ring and between aromatic atoms
func needsBondSymbol(_ bond: MKBond) -> Bool {
    if bond.isAromatic() { return false }
    switch bond.getBondOrder() {
    case 1:
        if bond.isInRing() && bond.getBeginAtom().isAromatic() && bond.getEndAtom().isAromatic() {
            return true
        }
        return false
    default: // bond orders != 1
        return true
    }
}


/****************************************************************************
 * FUNCTION: StandardLabels
 *
 * DESCRIPTION:
 *        Creates a set of non-canonical labels for the fragment atoms
 * ***************************************************************************/

func standardLabels(_ mol: MKMol, _ frag_atoms: Bitset, _ symmetry_classes: inout [Ref], _ labels: inout [Ref]) {
    for atom in mol.getAtomIterator() {
        if frag_atoms[atom.getIdx()] {
            labels.append(.Ref(atom.getIdx() - 1))
            symmetry_classes.append(.Ref(atom.getIdx() - 1))
        }
        else {
            labels.append(.ImplicitRef) //to match situation when canonical ordering. Just a big number?
            symmetry_classes.append(.ImplicitRef)
        }
    }
}

/***************************************************************************
 * FUNCTION: RandomLabels
 *
 * DESCRIPTION:
 *    Creates a set of random labels for the fragment atoms.  Primarily
 *    for testing: you can create a bunch of random SMILES for the same
 *    molecule, and use those to test the canonicalizer.
 ***************************************************************************/

func randomLabels(_ mol: MKMol, _ frag_atoms: Bitset, _ symmetry_classes: inout [Ref], _ labels: inout [Ref]) {
    let natoms = mol.numAtoms()
    let used = Bitset()
    
    for atom in mol.getAtomIterator() {
        if frag_atoms[atom.getIdx()] {
            var r = Int.random(in: 0..<natoms)
            while used[r] {
                r = (r + 1) % natoms // find an unused number
            }
            used.add(r)
            labels.append(.Ref(r))
            symmetry_classes.append(.Ref(r))
        }
        else {
            labels.append(.ImplicitRef) //to match situation when canonical ordering. Just a big number?
            symmetry_classes.append(.ImplicitRef)
        }
    }
}

/**
* Helper function for getFragment below.
*/
func addNbrs(_ fragment: inout Bitset, _ atom: MKAtom, _ mask: Bitset) {
    for nbr in atom.getNbrAtomIterator()! {
        if !mask[nbr.getIdx()] { continue } // skip atoms not in mask
        if fragment[nbr.getIdx()] { continue }  // skip visited atoms
        fragment.add(nbr.getIdx()) // add the neighbor atom to the fragment
        addNbrs(&fragment, nbr, mask) // recurse...
    }
}

/**
* Create an OBBitVec objects with bets set for the fragment consisting of all
* atoms for which there is a path to atom without going through skip. These
* fragment bitvecs are indexed by atom idx (i.e. OBAtom::GetIdx()).
*/

func getFragment(_ atom: MKAtom, _ mask: Bitset) -> Bitset {
    var fragment = Bitset()
    fragment.add(atom.getIdx())
    addNbrs(&fragment, atom, mask)
    return fragment
}

/***************************************************************************
* FUNCTION: CreateCansmiString
*
* DESCRIPTION:
*       Writes the canonical SMILES for a molecule or molecular fragment
*       to the given buffer.
*
*       frag_atoms represents atoms in a fragment of the molecule; the
*       SMILES will contain those atoms only.
*
*       (Note: This is an ordinary public C++ function, not a member
*       of any class.)
*
***************************************************************************/
func createCansmiString(_ mol: inout MKMol, _ buffer: inout String, _ frag_atoms: inout Bitset, _ pConv: MKConversion) {
    let canonical: Bool = pConv.isOption("c")

    let options = OutOptions(isomeric: !pConv.isOption("i"), kekulesmi: pConv.isOption("k"),
                             showatomclass: pConv.isOption("a"),
                             showexplicitH: pConv.isOption("h"), smarts: pConv.isOption("s"),
                             ordering: pConv.isOption("o")!)
    
    let m2s = MKMol2Cansmi(options, mol, canonical, pConv)

    if options.isomeric {
        perceiveStereo(&mol)
        m2s.createCisTrans(mol) // No need for this if not iso
    } else {
        // Not isomeric - be sure there are no Z coordinates, clear
        // all stereo-center and cis/trans information.
        for bond in mol.getBondIterator() {
            bond.setHash(false)
            bond.setWedge(false)
        }
    }

    if !options.showexplicitH {
        // If the fragment includes explicit hydrogens, exclude them.
        // They won't appear in the SMILES anyway (unless they're attached to
        // a chiral center, or it's something like [H][H]).
        for atom in mol.getAtomIterator() {
            if frag_atoms[atom.getIdx()] && atom.getAtomicNum() == MKElements.Hydrogen.atomicNum
                && (!options.isomeric || m2s.isSuppressedHydrogen(atom)) {
                frag_atoms.remove(atom.getIdx())
            }
        }
    }
    m2s.createFragCansmiString(mol, &frag_atoms, &buffer)
    if pConv.isOption("O") {
        // This atom order data is useful not just for canonical SMILES
        // Could also save canonical bond order if anyone desires
        var canData: MKPairData<String>
        if mol.hasData("SMILES Atom Order") {
            // Create new OBPairData
            canData = MKPairData()
            canData.setAttribute("SMILES Atom Order")
            canData.setOrigin(.local)
            mol.setData(canData)
        }
        else {
            // Recanonicalizing - update existing new OBPairData
            canData = mol.getData("SMILES Atom Order") as! MKPairData
        }
        var atmorder = ""
        m2s.getOutputOrder(&atmorder)
        canData.setValue(atmorder)
    }

}


class FIXFormat: MKMoleculeFormat {
    
    override init() {
        super.init() 
        MKConversion.registerFormat("fix", self)
    }
    
    required init(_ id: String, _ isDefault: Bool) {
        fatalError("init(_:_:) has not been implemented")
    }
    
    override func description() -> String? {
        return "SMILES FIX format\n  No comments yet\n" 
    }

    override func specificationURL() -> String {
        return "" 
    }

    override func flags() -> Int {
        return NOTREADABLE
    }

    override func writeMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        guard var mol = pOb as? MKMol else { return false }
        //Define some references so we can use the old parameter names
        let ofs = pConv.getOutStream()

        var buffer = ""

        let options = OutOptions(isomeric: !pConv.isOption("i"), kekulesmi: pConv.isOption("k"),
                                 showatomclass: pConv.isOption("a"),
                                 showexplicitH: pConv.isOption("h"), smarts: pConv.isOption("s"),
                                 ordering: pConv.isOption("o")!)
        let m2s = MKMol2Cansmi(options, mol, true, pConv)
        
        var allbits = Bitset()
        
        for a in mol.getAtomIterator() {
            allbits.add(a.getIdx())
        }
        
        if mol.numAtoms() > 0 {
            createCansmiString(&mol, &buffer, &allbits, pConv)
        }
        // add newline to end of buffer 
        buffer += "\n"
        do {
            try ofs?.write(data: buffer.data(using: .utf8)!)
        } catch {
            print("Error writing to output stream")
        }
        var orderString: String = ""
        m2s.getOutputOrder(&orderString)
        let canonical_order = orderString.components(separatedBy: .whitespacesAndNewlines)

        for j in 0..<mol.numConformers() {
            mol.setConformer(j)
            for index in 0..<canonical_order.count {
                let atomIdx = Int(canonical_order[index])!
                guard let atom = mol.getAtom(atomIdx) else {
                    fatalError("Could not get atom \(atomIdx)")
                }
                var coords = String(format: "%9.3f %9.3f %9.3f", atom.getX(), atom.getY(), atom.getZ())
                coords += "\n"
                do {
                    try ofs?.write(data: coords.data(using: .utf8)!)
                } catch {
                    print("Error writing coords to output stream")
                }
            }
        }
        return true 
    }

}
