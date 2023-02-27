


import Foundation
import Collections

class MKMoleculeFormat: MKFormat {
    
    typealias NameIndexType = Dictionary<String, Int>
    
    static var optionsRegistered: Bool = false
    static var iMols: OrderedDictionary<String, MKMol> = [:]
    static var _jmol: MKMol? = nil           //!< Accumulates molecules with the -j option
    static var molArray: [MKMol]? = nil      //!< Used in --separate option
    static var storedMolsReady: Bool = false //!< Used in --separate option
    static var _pDesc: MKDescriptor?
    
    override init() {
        super.init()
        
        if !MKMoleculeFormat.optionsRegistered {
            MKMoleculeFormat.optionsRegistered = true
            
            MKConversion.registerOptionParam("b",          self, 0, .INOPTIONS)
            MKConversion.registerOptionParam("s",          self, 0, .INOPTIONS)
            MKConversion.registerOptionParam("title",      self, 1, .GENOPTIONS)
            MKConversion.registerOptionParam("addtotitle", self, 1, .GENOPTIONS)
            MKConversion.registerOptionParam("property",   self, 2, .GENOPTIONS)
            MKConversion.registerOptionParam("C",          self, 0, .GENOPTIONS)
            MKConversion.registerOptionParam("j",          self, 0, .GENOPTIONS)
            MKConversion.registerOptionParam("join",       self, 0, .GENOPTIONS)
            MKConversion.registerOptionParam("separate",   self, 0, .GENOPTIONS)
            
            //The follow are OBMol options, which should not be in OBConversion.
            //But here isn't entirely appropriate either, since one could have
            //OBMol formats loaded but which don't derived from this class.
            //However, this possibility is remote.
            
            MKConversion.registerOptionParam("s",      nil, 1, .GENOPTIONS)
            MKConversion.registerOptionParam("v",      nil, 1, .GENOPTIONS)
            MKConversion.registerOptionParam("h",      nil, 0, .GENOPTIONS)
            MKConversion.registerOptionParam("d",      nil, 0, .GENOPTIONS)
            MKConversion.registerOptionParam("b",      nil, 0, .GENOPTIONS)
            MKConversion.registerOptionParam("c",      nil, 0, .GENOPTIONS)
            MKConversion.registerOptionParam("p",      nil, 1, .GENOPTIONS)
            MKConversion.registerOptionParam("t",      nil, 0, .GENOPTIONS)
            MKConversion.registerOptionParam("k",      nil, 0, .GENOPTIONS)
            MKConversion.registerOptionParam("filter", nil, 1, .GENOPTIONS)
            MKConversion.registerOptionParam("add",    nil, 1, .GENOPTIONS)
            MKConversion.registerOptionParam("delete", nil, 1, .GENOPTIONS)
            MKConversion.registerOptionParam("append", nil, 1, .GENOPTIONS)
            
        }
        
        
    }
    
    //! Static routine,  which can be called from elsewhere
    static func readChemObjectImpl(_ pConv: MKConversion, _ format : MKFormat) -> Bool {
        
    }
    
    static func writeChemObjectImpl(_ pConv: MKConversion, _ format : MKFormat) -> Bool {
        
    }
    
    override func readChemObject(_ pConv: MKConversion) -> Bool {
        return MKMoleculeFormat.readChemObjectImpl(pConv, self)
    }
    
    override func writeChemObject(_ pConv: MKConversion) -> Bool {
        return MKMoleculeFormat.writeChemObjectImpl(pConv, self)
    }

    ///Applies output options to molecule. Returns false to terminate output.
    static func doOutputOptions(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        if pConv.isOption("addoutindex", .GENOPTIONS) != nil {
            let newS = "\(pOb.getTitle()) \(pConv.getOutputIndex())"
            pOb.setTitle(newS)
        }
        
        if var pmol = pOb as? MKMol {
            if pConv.isOption("writeconformers", .GENOPTIONS) != nil {
                //The last conformer is written in the calling function
                var c = 0
                while c < pmol.numConformers() - 1 {
                    pmol.setConformer(c)
                    if !pConv.getOutFormat().writeMolecule(pOb, pConv) {
                        break
                    }
                    c += 1
                }
                pmol.setConformer(c)
            }
        }
        return true
    }

    /// \name Routines to handle the -C option for combining data from several OBMols
    //@{
    //! Defer output of a molecule until later, so it can be combined with others
    //! \return Success, or false if no molecule was read.
    static func  deferMolOutput(_ pmol: MKMol, _ pConv: MKConversion, _ pF: MKFormat) -> Bool {
        
    }
    //! Write out all molecules queued with DeferMolOutput
    static func  outputDeferredMols(_ pConv: MKConversion) -> Bool {
        
    }
    //! Delete the list of queued molecules from DeferMolOutput
    static func  deleteDeferredMols() -> Bool {
        //Empties IMols, deteting the OBMol objects whose pointers are stored there, for c++ that is
        iMols.removeAll()
        return false
    }
    //! \return the OBMol which combines @p pFirst and @p pSecond (i.e.)
    static func  makeCombinedMolecule(_ pFirst: MKMol, _ pSecond: MKMol) -> MKMol {
        
    }

    //!When sent an OBReaction object, output all the constituent molecules
    static func outputMolsFromReaction(_ pReact: MKReaction, _ pConv: MKConversion, _ pFormat: MKFormat) -> Bool {}
    
    //////////////////////////////////////////////////////////////////
    /** Attempts to read the index file datafilename.obindx successively
        from the following directories:
        - the current directory
        - that in the environment variable BABEL_DATADIR or in the macro BABEL_DATADIR
        if the environment variable is not set
        - in a subdirectory of the BABEL_DATADIR directory with the version of OpenBabel as its name
        An index of type NameIndexType is then constructed. NameIndexType is defined
        in obmolecformat.h as std::unordered_map. It is searched by
        @code
        NameIndexType::iterator itr = index.find(molecule_name);
        if(itr!=index.end())
        unsigned pos_in_datafile = itr->second;
        @endcode
        pos_in_datafile is used as a parameter in seekg() to read from the datafile

        If no index is found, it is constructed from the datafile by reading all of
        it using the format pInFormat, and written to the directory containing the datafile.
        This means that this function can be used without worrying whether there is an index.
        It will be slow to execute the first time, but subsequent uses get the speed benefit
        of indexed access to the datafile.

        The serialization and de-serialization of the NameIndexType is entirely in
        this routine and could possibly be improved. Currently re-hashing is done
        every time the index is read.
    **/
    static func readNameIndex(_ index: NameIndexType, _ datafilename: String, _ pInFormat: MKFormat) -> Bool {
        
    }

    //! \return the type of data converted by this format (here, OBMol)
    override func getType() -> String {
        return MKMol.self.className()
    }

}

