//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/24/23.
//

import Foundation
import Collections


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
class SMIBaseFormat: MKMoleculeFormat {

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
    override func readMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        var pmol = pOb.castAndClear(true) as! MKMol
        var ifs = pConv.getInStream()
        var ln: String
        var title: String
        var smiles: String
        var pos: String.Index
        
        //Ignore lines that start with #
        while(ifs != nil && ifs!.peek() == "#") {
            if(!getline(ifs, ln))
                return false;
        }

        //Get title
        if(getline(ifs, ln)) {
            pos = ln.find_first_of(" \t");
            if(pos!=string::npos)
            {
                smiles = ln.substr(0,pos);
                title = ln.substr(pos+1);
                Trim(title);
                pmol->SetTitle(title.c_str());
            } else {
                smiles = ln
            }
        }

        pmol.setDimension(0)
        var sp = MKSmilesParser(preserve_aromaticity: (pConv.isOption("a", .INOPTIONS) != nil))
        if !(pConv.isOption("S", .INOPTIONS) == nil){
            pmol.setChiralityPerceived()
        }

        return sp.smiToMol(pmol, smiles) //normal return
    }

    override func writeMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        fatalError()
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
        var ifs = pConv.getInStream()
        fatalError()
    }
    
    private func getInchifiedSMILESMolecules(_ mol: MKMol, _ useFixedHRecMet: Bool) { }
    
}

class SMIFormat: SMIBaseFormat {
    
    override init() {
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
        var atoms: Array<MKAtom> 
        var updown: Array<Character> 
    }
    
    var _updown: Character
    var _order: Int
    var _prev: Int
    var _rxnrole: Int
    // const char *_ptr; ??
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
    var chiralWatch: Bool  // set when a tetrahedral atom is read
    var _tetrahedralMap: OrderedDictionary<MKAtom, MKTetrahedralStereo.Config>  // map of tetrahedral atoms and their data
    var _upDownMap: OrderedDictionary<MKBond, Character>  // store the '/' & '\' as they occurred in smiles
    var _chiralLonePair: OrderedDictionary<UInt, Character> // for atoms with potential chiral lone pairs, remember when the l.p. was encountered
    var squarePlanarWatch: Bool  // set when a square planar atom is read
    var _squarePlanarMap: OrderedDictionary<MKAtom, MKSquarePlanarStereo.Config>
        
    init(preserve_aromaticity: Bool = false) {
        self._preserve_aromaticity = preserve_aromaticity
        self._rxnrole = 1
    }

    
    func smiToMol(_ mol: MKMol, _ name: String) -> Bool {
        _vprev.clear();
        _rclose.clear();
        _prev=0;
        chiralWatch=false;
        squarePlanarWatch = false;

        // We allow the empty reaction (">>") but not the empty molecule ("")
        if (!ParseSmiles(mol, s) || (!mol.IsReaction() && mol.NumAtoms() == 0))
        {
            mol.Clear();
            return(false);
        }

        // TODO: Is the following a memory leak? - there are return statements above
        map<OBAtom*, OBTetrahedralStereo::Config*>::iterator i;
        for (i = _tetrahedralMap.begin(); i != _tetrahedralMap.end(); ++i)
        delete i->second;
        _tetrahedralMap.clear();

        map<OBAtom*, OBSquarePlanarStereo::Config*>::iterator j;
        for (j = _squarePlanarMap.begin(); j != _squarePlanarMap.end(); ++j)
        delete j->second;
        _squarePlanarMap.clear();

        mol.SetAutomaticFormalCharge(false);

        return(true);
    } 

    func parseSmiles(_ mol: MKMol, _ smiles: String) -> Bool {
        mol.SetAromaticPerceived(); // Turn off perception until the end of this function
        mol.BeginModify();

        for (_ptr=smiles.c_str();*_ptr;_ptr++)
        {
        switch(*_ptr)
        {
        case '\r':
            if (*(_ptr+1) == '\0') // may have a terminating '\r' due to Windows line-endings
            break;
            return false;
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9':
        case '%':  //ring open/close
            if (_prev == 0)
            return false;
            if (!ParseRingBond(mol))
            return false;
            break;
        case '&': //external bond
            if (_prev == 0)
            return false;
            if (!ParseExternalBond(mol))
            return false;
            break;
        case '.':
            _prev=0;
            break;
        case '>':
            _prev = 0;
            _rxnrole++;
            if (_rxnrole == 2) {
            mol.SetIsReaction();
            // Handle all the reactant atoms
            // - the remaining atoms will be handled on-the-fly
            FOR_ATOMS_OF_MOL(atom, mol) {
                OBPairInteger *pi = new OBPairInteger();
                pi->SetAttribute("rxnrole");
                pi->SetValue(1);
                atom->SetData(pi);
            }
            }
            else if (_rxnrole == 4) {
            stringstream errorMsg;
            errorMsg << "Too many greater-than signs in SMILES string";
            std::string title = mol.GetTitle();
            if (!title.empty())
                errorMsg << " (title is " << title << ")";
            errorMsg << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
            return false;
            }
            break;
        case '(':
            _vprev.push_back(_prev);
            break;
        case ')':
            if(_vprev.empty()) //CM
            return false;
            _prev = _vprev.back();
            _vprev.pop_back();
            break;
        case '[':
            if (!ParseComplex(mol))
            {
            mol.EndModify();
            mol.Clear();
            return false;
            }
            break;
        case '-':
            if (_prev == 0)
            return false;
            _order = 1;
            break;
        case '=':
            if (_prev == 0)
            return false;
            _order = 2;
            break;
        case '#':
            if (_prev == 0)
            return false;
            _order = 3;
            break;
        case '$':
            if (_prev == 0)
            return false;
            _order = 4;
            break;
        case ':':
            if (_prev == 0)
            return false;
            _order = 0; // no-op
            break;
        case '/':
            if (_prev == 0)
            return false;
            _order = 1;
            _updown = BondDownChar;
            break;
        case '\\':
            if (_prev == 0)
            return false;
            _order = 1;
            _updown = BondUpChar;
            break;
        default:
            if (!ParseSimple(mol))
            {
            mol.EndModify();
            mol.Clear();
            return false;
            }
        } // end switch
        } // end for _ptr

        // place dummy atoms for each unfilled external bond
        if(!_extbond.empty())
        CapExternalBonds(mol);

        // Check to see if we've balanced out all ring closures
        // They are removed from _rclose when matched
        if (!_rclose.empty()) {
        mol.EndModify();
        mol.Clear();

        stringstream errorMsg;
        errorMsg << "Invalid SMILES string: " << _rclose.size() << " unmatched "
                << "ring bonds." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false; // invalid SMILES since rings aren't properly closed
        }

        // Check to see if we've the right number of '>' for reactions
        if (_rxnrole > 1 && _rxnrole !=3) {
        mol.EndModify();
        stringstream errorMsg;
        errorMsg << "Invalid reaction SMILES string: only a single '>' sign found (two required to be valid).";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false; // invalid SMILES since rings aren't properly closed
        }
        if (mol.IsReaction()) {
        OBReactionFacade facade(&mol);
        facade.AssignComponentIds();
        }

        // Apply the SMILES valence model
        FOR_ATOMS_OF_MOL(atom, mol) {
        unsigned int idx = atom->GetIdx();
        int hcount = _hcount[idx - 1];
        if (hcount == -1) { // Apply SMILES implicit valence model
            unsigned int bosum = 0;
            FOR_BONDS_OF_ATOM(bond, &(*atom)) {
            bosum += bond->GetBondOrder();
            }
            unsigned int impval = SmilesValence(atom->GetAtomicNum(), bosum);
            unsigned int imph = impval - bosum;
            if (imph > 0 && atom->IsAromatic())
            imph--;
            atom->SetImplicitHCount(imph);
        }
        else // valence is explicit e.g. [CH3]
            atom->SetImplicitHCount(hcount);
        }

        mol.EndModify(false);

        // Unset any aromatic bonds that *are not* in rings where the two aromatic atoms *are* in a ring
        // This is rather subtle, but it's correct and reduces the burden of kekulization
        FOR_BONDS_OF_MOL(bond, mol) {
        if (bond->IsAromatic() && !bond->IsInRing()) {
            if (bond->GetBeginAtom()->IsInRing() && bond->GetEndAtom()->IsInRing())
            bond->SetAromatic(false);
        }
        }

        // TODO: Only Kekulize if the molecule has a lower case atom
        bool ok = OBKekulize(&mol);
        if (!ok) {
        stringstream errorMsg;
        errorMsg << "Failed to kekulize aromatic SMILES";
        std::string title = mol.GetTitle();
        if (!title.empty())
            errorMsg << " (title is " << title << ")";
        errorMsg << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        // return false; // Should we return false for a kekulization failure?
        }

        // Add the data stored inside the _tetrahedralMap to the atoms now after end
        // modify so they don't get lost.
        if(!_tetrahedralMap.empty()) {
        OBAtom* atom;
        map<OBAtom*, OBTetrahedralStereo::Config*>::iterator ChiralSearch;
        for(ChiralSearch = _tetrahedralMap.begin(); ChiralSearch != _tetrahedralMap.end(); ++ChiralSearch) {
            atom = ChiralSearch->first;
            OBTetrahedralStereo::Config *ts = ChiralSearch->second;
            if (!ts)
            continue;
            if (ts->refs.size() != 3)
            continue;
            if (ts->refs[2] == OBStereo::NoRef) {
            // This happens where there is chiral lone pair or where there simply aren't enough connections
            // around a chiral atom. We handle the case where there is a S with a chiral lone pair.
            // All other cases are ignored, and raise a warning. (Note that S can be chiral even without
            // a lone pair, think of C[S@](=X)(=Y)Cl.

            // We have remembered where to insert the lone pair in the _chiralLonePair map
            map<unsigned int, char>::iterator m_it = _chiralLonePair.find(atom->GetIdx());
            if (CanHaveLonePair(atom->GetAtomicNum()) && m_it != _chiralLonePair.end()) {
                ts->refs[2] = ts->refs[1]; ts->refs[1] = ts->refs[0];
                if (m_it->second == 0) { // Insert in the 'from' position
                ts->refs[0] = ts->from;
                ts->from = OBStereo::ImplicitRef;
                }
                else // Insert in the refs[0] position
                ts->refs[0] = OBStereo::ImplicitRef;
            }
            else { // Ignored by Open Babel
                stringstream ss;
                ss << "Ignoring stereochemistry. Not enough connections to this atom. " << mol.GetTitle();
                obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
                continue;
            }
            }

            // cout << "*ts = " << *ts << endl;
            OBTetrahedralStereo *obts = new OBTetrahedralStereo(&mol);
            obts->SetConfig(*ts);
            mol.SetData(obts);
        }
        }

        // Add the data stored inside the _squarePlanarMap to the atoms now after end
        // modify so they don't get lost.
        if(!_squarePlanarMap.empty()) {
        OBAtom* atom;
        map<OBAtom*, OBSquarePlanarStereo::Config*>::iterator ChiralSearch;
        for(ChiralSearch = _squarePlanarMap.begin(); ChiralSearch != _squarePlanarMap.end(); ++ChiralSearch) {
            atom = ChiralSearch->first;
            OBSquarePlanarStereo::Config *sp = ChiralSearch->second;
            if (!sp)
            continue;
            if (sp->refs.size() != 4)
            continue;

            // cout << "*ts = " << *ts << endl;
            OBSquarePlanarStereo *obsp = new OBSquarePlanarStereo(&mol);
            obsp->SetConfig(*sp);
            mol.SetData(obsp);
        }
        }

        if (!_preserve_aromaticity)
        mol.SetAromaticPerceived(false);

        createCisTrans(mol);

        return(true);
    

     bool ParseSimple(OBMol&);
     bool ParseComplex(OBMol&);
     bool ParseRingBond(OBMol&);
     bool ParseExternalBond(OBMol&);
     bool CapExternalBonds(OBMol &mol);
     int  NumConnections(OBAtom *, bool isImplicitRef=false);
            func createCisTrans(_ mol: inout MKMol) {}
     char SetRingClosureStereo(StereoRingBond rcstereo, OBBond* dbl_bond);
     void InsertTetrahedralRef(OBMol &mol, unsigned long id);
     void InsertSquarePlanarRef(OBMol &mol, unsigned long id);
    
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
    
}



