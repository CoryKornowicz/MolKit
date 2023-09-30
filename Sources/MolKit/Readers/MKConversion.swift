//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/24/23.
//

import Foundation
import Collections

///@brief Three types of options set on the the command line by -a? , -x? , or -?
enum Option_Type {
    case INOPTIONS
    case OUTOPTIONS
    case GENOPTIONS
    case ALL
    
    static var allCases: [Option_Type] {
        return [.GENOPTIONS, .INOPTIONS, .OUTOPTIONS]
    }
}

//*************************************************
/// @brief Class to convert from one format to another.

public class MKConversion {

    fileprivate typealias MKAMapType = OrderedDictionary<String, Int>
    
    // Utilize base datatype as String, could be replaced with a class later that handles streaming to a live Stream and uses InputFileHandler and OutputFileHandler and conforms to StreamDelegate
//    MARK: StreamState
    private struct StreamState {
        
        var pStream: FileHandlerProtocol?
        var owndedStreams: [FileHandlerProtocol?]
        
        init() {
            self.pStream = nil
            self.owndedStreams = []
        }
        
        //save the current input state to this streamstate and clear conv
        mutating func pushInput(_ conv: MKConversion) {
            precondition(owndedStreams.count == 0, "ownedStreams be empty")
            
            pStream = conv.pInput as? any FileHandlerProtocol
            
            conv.ownedInStreams.forEach { ifp in
                owndedStreams.append(ifp as? FileHandlerProtocol)
            }
            
            conv.pInput = nil
            conv.ownedInStreams = []
        }
        
        //restore state, blowing away whatever is in conv
        mutating func popInout(_ conv: MKConversion) {
            conv.setInStream(nil)
            conv.pInput = pStream as? InputFileHandler
            
            precondition(conv.ownedInStreams.count == 0, "ownedStreams be empty")
            
            owndedStreams.forEach { fH in
                conv.ownedInStreams.append(fH as! InputFileHandler)
            }
            
            pStream = nil
            owndedStreams = []
        }
        
        mutating func pushOutput(_ conv: MKConversion) {
            precondition(owndedStreams.count == 0, "ownedStreams be empty")
            
            pStream = conv.pOutput as? FileHandler
            
            conv.ownedOutStreams.forEach { ofp in
                owndedStreams.append(ofp as! OutputFileHandler)
            }
            
            
            conv.pOutput = nil
            conv.ownedOutStreams = []
        }
        
        mutating func popOutput(_ conv: MKConversion) {
            conv.setOutStream(nil)
            conv.pOutput = pStream as? OutputFileHandler
            
            precondition(conv.ownedOutStreams.count == 0, "ownedStreams be empty")
            
            owndedStreams.forEach { fH in
                conv.ownedOutStreams.append(fH as! OutputFileHandler)
            }
            
            pStream = nil
            owndedStreams = []
        }
        
        func isSet() -> Bool { return self.pStream != nil || pStream?.streamStatus == .none }
    }
    
//    MARK: Class Variables
    
    var inFilename: String = ""
    var outFilename: String = ""
    
    var pInput: InputFileHandlerProtocol? //input stream, may be filtered
    var ownedInStreams: [InputFileHandlerProtocol] = []
    
    var pOutput: OutputFileHandlerProtocol?  //output stream, may have filters applied
    var ownedOutStreams: [OutputFileHandlerProtocol] = []
    
    static var pDefaultFormat: MKFormat?
    var pInFormat: MKFormat?
    var pOutFormat: MKFormat?
    
    var OptionsArray: OrderedDictionary<Option_Type, OrderedDictionary<String, String>> = [.GENOPTIONS: [:], .INOPTIONS: [:], .OUTOPTIONS: [:]]
    
    static var opa : OrderedDictionary<Option_Type, OrderedDictionary<String, String>> = [.GENOPTIONS: [:], .INOPTIONS: [:], .OUTOPTIONS: [:]]
    
    var Index: Int = 0
    var StartNumber: Int = 1
    var EndNumber: Int = 0
    var Count: Int = -1
    var m_IsFirstInput: Bool = false
    var m_IsLast: Bool = false
    var MoreFilesToCome: Bool = false
    var OneObjectOnly: Bool = false
    var ReadyToInput: Bool = false
    var SkippedMolecules: Bool = false     /// skip molecules using -f and -l

    //unlike the z and zin options, these are not sticky - setting formats will reset them
    var inFormatGzip: Bool = false
    var outFormatGzip: Bool = false

    var pOb: MKBase?
    var wInpos: Int?  ///<position in the input stream of the object being written
    var rInpos: Int?  ///<position in the input stream of the object being read
    var wInlen: Int?  ///<length in the input stream of the object being written
    var rInlen: Int?  ///<length in the input stream of the object being read
    
    var pAuxConv: MKConversion? ///<Way to extend OBConversion?

    var SupportedInputFormat:  [String] = [] ///< list of supported input format
    var SupportedOutputFormat: [String] = []///< list of supported output format
    
    private func defaultInit() {
        self.pInput = nil
        self.pOutput = nil
        self.pInFormat = nil
        self.pOutFormat = nil
        self.m_IsFirstInput = true
        self.m_IsLast = true
        self.MoreFilesToCome = false
        self.OneObjectOnly = false
        self.SkippedMolecules = false
        self.inFormatGzip = false
        self.outFormatGzip = false
        self.pOb = nil
        self.wInpos = 0
        self.wInlen = 0
        self.pAuxConv = nil
        self.inFilename = ""
        self.outFilename = ""
        
        self.OptionsArray[.INOPTIONS] = [:]
        self.OptionsArray[.OUTOPTIONS] = [:]
        self.OptionsArray[.GENOPTIONS] = [:]
        
        MKConversion.registerOptionParam("f", nil, 1, .GENOPTIONS)
        MKConversion.registerOptionParam("l", nil, 1, .GENOPTIONS)

    }
    
    public init(_ inStream: InputFileHandler? = nil, _ outStream: OutputFileHandler? = nil) {
        defaultInit()
        setInStream(inStream)
        setOutStream(outStream)
    }
    
    public init(_ inFilename: String, _ outFilename: String) {
        defaultInit()
        openInAndOutFiles(inFilename, outFilename)
    }
    
    
//    MARK: Functions
    
    // Collection of Formats
    /// @brief Called once by each format class
    /// Class information on formats is collected by making an instance of the class
    /// derived from OBFormat(only one is usually required). RegisterFormat() is called
    /// from its constructor.
    ///
    /// If the compiled format is stored separately, like in a DLL or shared library,
    /// the initialization code makes an instance of the imported OBFormat class.
    @discardableResult
    static func registerFormat(_ ID: String, _ pFormat: MKFormat, _ MIME: String? = nil) -> Int {
        return pFormat.registerFormat(ID, MIME)
    }

    ///@brief Searches registered formats
    static func findFormat(_ ID: String) -> MKFormat? {
        return MKFormat.findType(ID)
    }
    /// @brief Searches registered formats for an ID the same as the file extension
    static func formatFromExt(_ filename: String) -> MKFormat? {
        var isgzip = false
        return MKConversion.formatFromExt(filename, &isgzip)
    }
    
    static func formatFromExt(_ filename: String, _ isgzip: inout Bool) -> MKFormat? {
//        find last occurence of "."
        isgzip = false
        
        var fileExt: String = ""
        
        if let dotIndex = filename.lastIndex(of: ".") {
//            assert there is no "/" after this index
            if filename[dotIndex...].lastIndex(of: "/") == nil {
                fileExt = String(filename[dotIndex...]).trimmingCharacters(in: .init(charactersIn: "."))
                if fileExt != "gz" { // gzip file
                    return findFormat(fileExt)
                } else {
                    fatalError("Gzip support is not implemented yet")
                }
            }
        }
        // Check the filename if no extension (e.g. VASP does not use extensions):
        if let lastBreak = filename.lastIndex(of: "/") {
            return findFormat(String(filename[lastBreak...]))
        }
        
        return findFormat(fileExt)
    }
   /// @brief Searches registered formats for a MIME the same as the chemical MIME type passed
    static func formatFromMIME(_ MIME: String) -> MKFormat? {
        return MKFormat.formatFromMIME(MIME)
    }
    
    static func getNextFormat(_ itr: any IteratorProtocol<MKPlugin>, _ str: String, _ pFormat: MKFormat) -> Bool {
        fatalError()
    }
    
    /// @name Information
      //@{
    static func description() -> String {
        return """
Conversion options
-f <#> Start import at molecule # specified
-l <#> End import at molecule # specified
-e Continue with next object after error, if possible
"""
// -k Attempt to translate keywords TODO: Will need to implement with locale support
    } //generic conversion options
    //@}
    
    /// These return a filtered stream for reading/writing (possible filters include compression, decompression, and newline transformation)
    /// @name Parameter get and set
    //@{
    public func getInStream() -> (any InputFileHandlerProtocol)? { return pInput ?? nil }
    public func getOutStream() -> (any OutputFileHandlerProtocol)? { return pOutput ?? nil }
    
    /// @brief Set input stream.  If takeOwnership is true, will deallocate when done.
    /// If isGzipped is true, will treat as a gzipped stream regardless of option settings,
    //  if false, then will be treated as gzipped stream only if z/zin is set.
    /// Set input stream, removing/deallocating previous stream if necessary.
    /// If takeOwnership is true, takes responsibility for freeing pIn
    public func setInStream(_ pIn: InputFileHandlerProtocol?, _ takeOwnership: Bool = false) {
        // clear and deallocate any existing streams
        ownedInStreams.removeAll()
        pInput = nil
        
        if pIn != nil {
            if takeOwnership {
                ownedInStreams.append(pIn!)
            }
            
            self.pInput = pIn
            //        LIBZ support to come in the future
            
            //always transform newlines if input isn't binary/xml
            // TODO: Will need to fix for other file types ??
            if((pInFormat != nil) && !(((pInFormat?.flags() ?? 0) & (READBINARY | READXML)) == 0) && (pIn! as! FileHandler != FileHandle.standardInput)) //avoid filtering stdin as well
            {
                //            LEInStream *leIn = new LEInStream(*pInput);
                ownedInStreams.append(pInput!)
                //            pInput = leIn;
            }
        }
    }
    
    /// Set output stream, removing/deallocating previous stream if necessary.
    /// If takeOwnership is true, takes responsibility for freeing pOut
    /// Be aware that if the output stream is gzipped format, then this outstream
    /// either needs to be replaced (e.g., SetOutStream(NULL)) or the OBConversion
    /// destroyed before the underlying OutputFileHandler is deallocated.
    public func setOutStream(_ pOut: (any OutputFileHandlerProtocol)?, _ takeOwnership: Bool = false) {
        ownedOutStreams.removeAll()
        pOutput = nil
        if pOut != nil {
            if takeOwnership {
                ownedOutStreams.append(pOut!)
            }
        }
        pOutput = pOut
        // libz support to come in the future
    }
    
    /// Sets the formats from their ids, e g CML
    /// //////////////////////////////////////////////////////
    /// If inID is NULL, the input format is left unchanged. Similarly for outID
    /// Returns true if both formats have been successfully set at sometime
    @discardableResult
    public func setInAndOutFormats(_ inID: String, _ outID: String, _ ingzip: Bool = false, _ outgzip: Bool = false) -> Bool {
        return setInFormat(inID, isgzip: ingzip) && setOutFormat(outID, isgzip: outgzip)
    }
    
    @discardableResult
    public func setInAndOutFormats(_ pIn: MKFormat, _ pOut: MKFormat, _ ingzip: Bool = false, _ outgzip: Bool = false) -> Bool {
        return setInFormat(pIn, isgzip: ingzip) && setOutFormat(pOut, isgzip: outgzip)
    }
    
    /// Sets the input format from an id e.g. CML
    @discardableResult
    public func setInFormat(_ inID: String, isgzip: Bool = false) -> Bool {
        inFormatGzip = isgzip
        if let format = MKConversion.findFormat(inID) {
            pInFormat = format
            return ((pInFormat!.flags() & NOTREADABLE) == 0)
        } else {
            return false
        }
    }
    
    @discardableResult
    public func setInFormat(_ pIn: MKFormat, isgzip: Bool = false) -> Bool {
        inFormatGzip = isgzip
        pInFormat = pIn
        return ((pInFormat!.flags() & NOTREADABLE) == 0)
    }

    /// Sets the output format from an id e.g. CML
    @discardableResult
    public func setOutFormat(_ outID: String, isgzip: Bool = false) -> Bool {
        outFormatGzip = isgzip
        if let format = MKConversion.findFormat(outID) {
            pOutFormat = format
            return ((pInFormat!.flags() & NOTREADABLE) == 0)
        } else {
            return false
        }
    }
    
    public func setOutFormat(_ pOut: MKFormat, isgzip: Bool = false) -> Bool {
        outFormatGzip = isgzip
        pOutFormat = pOut
        return ((pOutFormat!.flags() & NOTREADABLE) == 0)
    }

    func getInFormat() -> MKFormat? { return pInFormat }
    func getOutFormat() -> MKFormat? { return pOutFormat }

    func getInGzipped() -> Bool { return inFormatGzip }
    func getOutGzipped() -> Bool { return outFormatGzip }

    func getInFilename() -> String { return inFilename }
    func getOutFilename() -> String { return outFilename }

    ///Get the position in the input stream of the object being read
    func getInPos() -> Int? { return wInpos }

    ///Get the length in the input stream of the object being read
    func getInLen() -> Int? { return wInlen }

    /// \return a default title which is the filename
    func getTitle() -> String { return inFilename }

    ///@brief Extension method: deleted in ~OBConversion()
    func getAuxConv() -> MKConversion? { return pAuxConv }
    func setAuxConv(_ pConv: MKConversion?) { pAuxConv = pConv }
    
    //@}
      /** @name Option handling
       Three types of Option provide information and control instructions to the
       conversion process, INOPTIONS, OUTOPTIONS, GENOPTIONS, and are stored in each
       OBConversion object in separate maps. Each option has an id and an optional
       text string. They are set individually by AddOption() or (rarely) collectively
       in SetOptions(). Options cannot be altered but can be replaced with AddOption()
       and deleted with RemoveOption(), which, however, should be used in an op derived
       from OBOp (because of iterator invalidation).

       If the "Convert" interface is used, the GENOPTIONS are acted upon in the
       OBBase::DoTransformations() functions (currently only OBMol has one). This
       happens after the object has been input but before it has been output.
       All the options are available to input and output formats, etc. via the IsOption()
       function, and the interpretation of any text string needs to be done subsequently.

       In the commandline interface, options with single character ids are are indicated
       like -s, and those with multiple character ids like --gen3D. An option may have
       one or more parameters which appear, space separated, in the option's text string.
       With babel, unless the option is at the end of the command, it is necessary for
       the number of its parameters to be exactly that specified in RegisterOptionParam().
       The default is 0, but if it is more, and babel is likely to be used, this function
       should be called in the constructor of a format or op.
       With obabel (or the GUI), it is not necessary to call RegisterOptionParam().

       New GENOPTIONS can be defined (as plugins) using the class OBOp.

       It is customary for a format or op to document any INOPTIONS or OUTPTIONS it
       uses in its Description() function. As well as providing documentation during
       use, this is also parsed by the GUI to construct its checkboxes,etc., so it is
       advisable to give new Descriptions the same form as existing ones.

       Some conversion options, such as -f, -l, -m, are unlikely to be used in
       programming, but are listed in OBConversion::Description().  The built-in
       GENOPTIONS for OBMol objects are listed in OBMol::ClassDescription() which
       is in transform.cpp and also in this documentation under AddOption().
       */
    //@{

    ///@brief Determine whether an option is set. \return NULL if option not and a pointer to the associated text if it is
    func isOption(_ opt: String, _ opttype: Option_Type = .OUTOPTIONS) -> String? {
        //Returns NULL if option not found or a pointer to the text if it is
        if let pos = OptionsArray[opttype]![opt] {
            return pos
        } else {
            return nil
        }
    }
    
    ///@brief Determine whether an option is set. \return NULL if option not and a pointer to the associated text if it is
    func isOption(_ opt: String, _ opttype: Option_Type = .OUTOPTIONS) -> Bool {
        //Returns NULL if option not found or a pointer to the text if it is
        if OptionsArray[opttype]![opt] != nil {
            return true
        } else {
            return false
        }
    }

    ///@brief Access the map with option name as key and any associated text as value
    func getOptions(_ opttyp: Option_Type) -> OrderedDictionary<String, String>? { return OptionsArray[opttyp] }

    ///@brief Set an option of specified type, with optional text
    func addOption(_ opt: String, _ opttype: Option_Type = .OUTOPTIONS, _ text: String? = nil) {
        if text == nil {
            OptionsArray[opttype]![opt] = ""
        } else {
            OptionsArray[opttype]![opt] = text!
        }
    }

    func removeOption(_ opt: String, _ opttype: Option_Type) -> Bool {
        OptionsArray[opttype]?.removeValue(forKey: opt) != nil
    }

    ///@brief Set several single character options of specified type from string like ab"btext"c"ctext"
    func setOptions(_ opt: String, _ opttype: Option_Type) {
        if opt == "" { // clear all
            OptionsArray[opttype]?.removeAll()
        }
        let optionIterator = Iterator<Character>(Array(opt))
        
        while !optionIterator.isEmpty() {
            guard let ch = optionIterator.next() else { break }
            if optionIterator.peek() == "\"" {
                optionIterator.ignore()
                var txt = String(optionIterator.constructToEnd())
                if let pos = txt.firstIndex(of: "\"") {
                    txt.remove(at: pos)
                    OptionsArray[opttype]?[String(ch)] = txt
                    optionIterator.ignore(by: pos.utf16Offset(in: txt) + 2)
                } else {
                    //options is illformed
                    return
                }
            } else {
                OptionsArray[opttype]![String(ch)] = ""
            }
        }
    }

    ///@brief For example -h takes 0 parameters; -f takes 1. Call in a format constructor.
    static func registerOptionParam(_ name: String, _ pFormat: MKFormat?, _ numberParams: Int = 0, _ typ: Option_Type = .OUTOPTIONS) {
//        MARK: YIKESSSS this is a rough closure
        MKConversion.optionParamArray(typ) { optionParamsArray in
            if let pos = optionParamsArray.first(where: { (key: String, value: Int) in
                key == name
            }) {
                if pos.value != numberParams {
                    var description: String = "API"
                    if pFormat != nil {
                        description = pFormat!.description() ?? "nil description for \(pos.key)"
                    }
    //                Throw error
                    print("""
                        The number of parameters needed by option \(pos.key) in
                        \(description.substring(toIndex: (description.firstIndex(of: "\n")?.utf16Offset(in: description))!))
                        differs from an earlier registration.
                        """)
                }
            } else {
                optionParamsArray.updateValue(numberParams, forKey: name)
            }
        }
    }

    /// \return the number of parameters registered for the option, or 0 if not found
    static func getOptionParams(_ name: String, _ typ: Option_Type) -> Int {
        return optionParamArray(typ) { optionParamsArray in
            if let pos = optionParamsArray.first(where: { (key: String, value: Int) in
                key == name
            })?.value {
                return pos
            } else {
                return 0
            }
        }
    }

    ///@brief Copies the options (by default of all types) from one OBConversion Object to another.
    func copyOptions(_ pSourceConv: MKConversion, _ opttype: Option_Type = .ALL) {
        if opttype == .ALL {
            for i in Option_Type.allCases {
                OptionsArray[i] = pSourceConv.OptionsArray[i]
            }
        } else {
            OptionsArray[opttype] = pSourceConv.OptionsArray[opttype]
        }
    }

    /// @name Supported file format
    //@{
    // @brief Set and return the list of supported input format
    public func getSupportedInputFormat() -> [String] {
        var vlist: [String] = []
        var param: String? = "in"
        MKPlugin.ListAsVector("formats", &param, &vlist)
        return vlist
    }
    // @brief Set and return the list of supported output format
    public func getSupportedOutputFormat() -> [String] {
        var vlist: [String] = []
        var param: String? = "out"
        MKPlugin.ListAsVector("formats", &param, &vlist)
        return vlist
    }
    //@}

    // MARK: Conversion
    //@{
    /// @brief Conversion for single input and output stream
    /// //////////////////////////////////////////////////////
    /// Convert molecules from is into os.  If either is null, uses existing streams.
    /// If streams are specified, they do _not_ replace any existing streams.
    func convert(_ iss: InputFileHandler, _ oss: OutputFileHandler) -> Int {
        var savedIn = StreamState()
        var savedOut = StreamState()
        
        if iss.streamStatus == .open {
            savedIn.pushInput(self)
            setInStream(iss, false)
        }
        
        if oss.streamStatus == .open {
            savedOut.pushOutput(self)
            setOutStream(oss)
        }
        
        let count = convert()
        
        if savedIn.isSet() { savedIn.popInout(self) }
        if savedOut.isSet() { savedOut.popOutput(self) }
        return count
    }

    /// @brief Conversion with existing streams
    /// ////////////////////////////////////////////////////
    /// Actions the "convert" interface.
    ///    Calls the OBFormat class's ReadMolecule() which
    ///     - makes a new chemical object of its chosen type (e.g. OBMol)
    ///     - reads an object from the input file
    ///     - subjects the chemical object to 'transformations' as specified by the Options
    ///     - calls AddChemObject to add it to a buffer. The previous object is first output
    ///       via the output Format's WriteMolecule(). During the output process calling
    /// IsFirst() and GetIndex() (the number of objects including the current one already output.
    /// allows more control, for instance writing \<cml\> and \</cml\> tags for multiple molecule outputs only.
    ///
    ///    AddChemObject does not save the object passed to it if it is NULL (as a result of a DoTransformation())
    ///    or if the number of the object is outside the range defined by
    ///    StartNumber and EndNumber.This means the start and end counts apply to all chemical objects
    ///    found whether or not they    are output.
    ///
    ///    If ReadMolecule returns false the input conversion loop is exited.
    ///
    func convert() -> Int {
        
        guard pInput != nil else {
            print("InputFileHandler is nil")
            return 0
        }
        
        guard pInFormat != nil else {
            print("Input format is not set")
            return 0
        }
        
        Count = 0
        
        if !setStartAndEnd() {
            return 0
        }
        
        ReadyToInput = true
        m_IsLast = false
        pOb = nil
        wInlen = 0
        
        if ((pInFormat!.flags() & READONEONLY) != 0) {
            OneObjectOnly=true
        }
        
        //Input loop
        while (ReadyToInput && pInput!.streamStatus != .error) { //Possible to omit? && pInStream->peek() != EOF
        
            if((pInput! as! FileHandler) == FileHandle.standardInput) {
                if(pInput!.peek() == -1) { //Cntl Z Was \n but interfered with piping
                    break
                }
            } else {
                do {
                    rInpos = try pInput!.tellg()
                } catch {
                    fatalError("Error in retrieving pInput reading position!!")
                }
            }
            var ret: Bool = false
            do {
                ret = try pInFormat!.readChemObject(self)
                setFirstInput(false)
            } catch {
                if((isOption("e", .GENOPTIONS) == nil) && !OneObjectOnly) {
//                    TODO: maybe this should be a fatal error
                    print("Convert failed with an exception")
                    return Index // the number we've actually output so far
                }
            }
            
            if(!ret) {
                //error or termination request: terminate unless
                // -e option requested and successfully can skip past current object
                if((isOption("e", .GENOPTIONS) == nil) || pInFormat!.skipObjects(0, self) != 1) {
                    break
                }
            }
            if(OneObjectOnly) {
                break
            }
            // Objects supplied to AddChemObject() which may output them after a delay
            //ReadyToInput may be made false in AddChemObject()
            // by WriteMolecule() returning false  or by Count==EndNumber
        }
        
        //Output last object
        m_IsLast = !MoreFilesToCome
        
        //Output is always occurs at the end with the --OutputAtEnd option
        let oae: Bool = isOption("OutputAtEnd", .GENOPTIONS) != nil
        if(pOutFormat != nil && (!oae || m_IsLast)) {
            if((oae || (pOb != nil)) && !pOutFormat!.writeChemObject(self)) {
                Index -= 1
            }
        }
        //Put AddChemObject() into non-queue mode
        Count = -1
        EndNumber = 0
        StartNumber = 0
        pOb = nil //leave tidy
        MoreFilesToCome=false
        OneObjectOnly=false
                
        return Index
    }

    /// @brief Conversion with multiple input/output files:
    /**
       Makes input and output streams, and carries out normal,
       batch, aggregation, and splitting conversion.

       Normal
       Done if FileList contains a single file name and OutputFileName
       does not contain a *.

       Aggregation
       Done if FileList has more than one file name and OutputFileName does
       not contain * . All the chemical objects are converted and sent
       to the single output file.

       Splitting
       Done if FileList contains a single file name and OutputFileName
       contains a * . Each chemical object in the input file is converted
       and sent to a separate file whose name is OutputFileName with the
       * replaced by 1, 2, 3, etc.  OutputFileName must have at least one
       character other than the * before the extension.
       For example, if OutputFileName is NEW*.smi then the output files are
       NEW1.smi, NEW2.smi, etc.

       Batch Conversion
       Done if FileList has more than one file name and contains a * .
       Each input file is converted to an output file whose name is
       OutputFileName with the * replaced by the inputfile name without its
       path and extension.
       So if the input files were inpath/First.cml, inpath/Second.cml
       and OutputFileName was NEW*.mol, the output files would be
       NEWFirst.mol, NEWSecond.mol.

       If FileList is empty, the input stream that has already been set
       (usually in the constructor) is used. If OutputFileName is empty,
       the output stream already set is used.

       On exit, OutputFileList contains the names of the output files.

       Returns the number of Chemical objects converted.
    */
    func fullConvert(_ fileList: [String], _ outputFileName: String, _ outputFileList: [String]) -> Int {
        
        var pIs: InputStringStream?
        var pOs: OutputStringStream?

        var ssOut: InputStringStream?
        var ssIn: OutputStringStream?

        var pInStream: InputStream?
        var pOutStream: OutputStream?
        
        var hasMultipleOutputFile: Bool = false
        var count: Int = 0 
        setFirstInput()
        var commonInFormat: Bool = pInFormat != nil ? true : false // whether set in calling routine
        
        // TODO: Fill in once the file handler protocols are in place
        
        //OUTPUT

        if outputFileName.isEmpty {
            
        }

        return Count
    }
    //@}

     /// @name Conversion loop control
    //@{
    ///////////////////////////////////////////////////////
    ///    Called by ReadMolecule() to deliver an object it has read from an input stream.
    /// Used in two modes:
    ///  - When Count is negative it is left negative and the routine is just a store
    ///    for an OBBase object.  The negative value returned tells the calling
    ///    routine that no more objects are required.
    ///  - When count is >=0, probably set by Convert(), it acts as a queue of 2:
    ///    writing the currently stored value before accepting the supplied one. This delay
    ///    allows output routines to respond differently when the written object is the last.
    ///    Count is incremented with each call, even if pOb=NULL.
    ///    Objects are not added to the queue if the count is outside the range
    ///    StartNumber to EndNumber. There is no upper limit if EndNumber is zero.
    ///    The return value is Count ((>0) or 0 if WriteChemObject returned false.
    func addChemObject(_ pOb1: MKBase) -> Int {
        if Count < 0 {
            pOb = pOb1
            return Count
        }
        Count += 1
        if Count >= StartNumber {
            if Count == EndNumber {
                ReadyToInput = false
            }
            do {
                rInlen = pInput != nil ? try pInput!.tellg() - rInpos! : 0
            } catch {
                print("error in calling tellg")
                return 0
            }
            if pOb != nil && pOutFormat != nil {
                if !pOutFormat!.writeChemObject(self) {
                    Index -= 1
                    pOb = nil
                    return 0
                }
                if ((pOutFormat!.flags() & WRITEONEONLY) != 0) {
                    let errmsg: String = """
                    WARNING: You are attempting to convert a file
                    with multiple molecule entries into a format
                    which can only store one molecule. The current
                    output will only contain the first molecule.
                    
                    To convert this input into multiple separate
                    output files, with one molecule per file, try:
                    obabel [input] [output] -m
                    
                    To pick one particular molecule
                    (e.g., molecule 4), try:
                    obabel -f 4 -l 4 [input] [output]
                    """
                    print(errmsg)
                    ReadyToInput = false
                    pOb = nil
                    return Count
                }
            }
            pOb = pOb1
            wInpos = rInpos
            wInlen = rInlen
        }
        return Count
    }
    ///< @brief Retrieve from internal array during output
    func getChemObject() -> MKBase? {
        if pOb != nil {
            self.Index += 1
            return pOb!
        } else { return nil }
    }
    ///< @brief True if no more objects to be output
    func isLast() -> Bool {
        return m_IsLast
    }
      ///< @brief True if the first input object is being processed
    func isFirstInput() -> Bool {
        return m_IsFirstInput
    }
      ///< @brief Setwhether or not is the first input
    func setFirstInput(_ b: Bool = true) {
        m_IsFirstInput = b
    }
    
    ///< @brief Retrieves number of ChemObjects that have been actually output
    ///  //////////////////////////////////////////////////////
    ///Returns the number of objects which have been output or are currently being output.
    ///The outputindex is incremented when an object for output is fetched by GetChemObject().
    ///So the function will return 1 if called from WriteMolecule() during output of the first object.
    func getOutputIndex() -> Int {
        return Index
    }
      ///< @brief Sets output index (maybe to control whether seen as first object)
    func setOutputIndex(_ indx: Int) {
        Index = indx
    }
      ///<@brief Used with multiple input files. Off by default.
    func setMoreFilesToCome() {
        MoreFilesToCome = true
    }
      ///< @brief Used with multiple input files. Off by default.
    func setOneObjectOnly(_ b: Bool = true) {
        OneObjectOnly = b
        m_IsLast = true
    }
      ///< @brief Synonym for SetOneObjectOnly()
    func setLast(_ b: Bool = true) { setOneObjectOnly(b) }
      ///< @brief True if no more files to be read
    func isLastFile() -> Bool { return !MoreFilesToCome }

    /// @brief Number of objects read and processed
    /// Incremented after options are processed, so 0 for first object.  Returns -1 if Convert interface not used. 
    func getCount() -> Int { return Count}
    //@}

    /// @name Convenience functions
    ///The default format is set in a single OBFormat class (generally it is OBMol)
    static func getDefaultFormat() -> MKFormat? {
        if let res = MKFormat.findType(nil) {
            if type(of: res) == MKFormat.self {
                return res
            }
        }
        return nil
    }
    
    /// @brief Outputs an object of a class derived from OBBase.
    /// Part of "API" interface.
    /// The output stream can be specified and the change is retained in the OBConversion instance
    @discardableResult
    func write<T: MKBase>(_ pOb: T, _ pout: OutputFileHandler? = nil) -> Bool {
        
        if(pout != nil) { setOutStream(pout, false) }

        guard pOutFormat != nil else { return false }
        guard pOutput != nil else { return false }

        // Set the locale for number parsing to avoid locale issues: PR#1785463
//        obLocale.SetLocale();
//        // Also set the C++ stream locale
//        locale originalLocale = pOutput->getloc(); // save the original
//        locale cNumericLocale(originalLocale, "C", locale::numeric);
//        pOutput->imbue(cNumericLocale);

        // Increment the output counter.
        // This is done *before* the WriteMolecule because some of
        // the format plugins initialized when GetOutputIndex() == 1.
        // This matches the original Convert(), which increments
        // the count for GetChemObject() before the write.
        Index += 1

        // The actual work is done here
        let success = pOutFormat!.writeMolecule(pOb, self)

        // return the C locale to the original one
//        obLocale.RestoreLocale();
//        // Restore the C++ stream locale too
//        pOutput->imbue(originalLocale);

        return success
    }

    /// @brief Outputs an object of a class derived from OBBase as a string
      /// Part of "API" interface.
      /// The output stream is temporarily changed to the string and then restored
      /// This method is primarily intended for scripting languages without "stream" classes
      /// The optional "trimWhitespace" parameter allows trailing whitespace to be removed
      /// (e.g., in a SMILES string or InChI, etc.)
    func writeString<T: MKBase>(_ pOb: T, _ trimWhitespace: Bool = true) -> String { 
        let newStream: OutputStringStream = OutputStringStream()
        var temp: String = "" 
        if pOutFormat != nil {
            var savedOut = StreamState()
            savedOut.pushOutput(self)
            // The StreamState doesn't save all of the properties so
            // do it manually here.
            // Set/reset the Index to 0 so that any initialization
            // code in the formatters will be executed.
            let oldIndex = Index
            Index = 0
            // We'll only send one object, so save those properties too.
            let oldOneObjectOnly = OneObjectOnly
            let oldm_IsLast = m_IsLast
            setOneObjectOnly(true)
            setOutStream(newStream, false)
            write(pOb)
            savedOut.popOutput(self)
            // Restore the other stream properties
            m_IsLast = oldm_IsLast
            OneObjectOnly = oldOneObjectOnly
            Index = oldIndex
        } 

        temp = newStream.string
        if (trimWhitespace) { // trim the trailing whitespace
            temp = temp.trimmingCharacters(in: .whitespaces)
        }
        return temp
    }
    
    /// @brief Outputs an object of a class derived from OBBase as a file (with the supplied path)
      /// Part of "API" interface.
      /// The output stream is changed to the supplied file and the change is retained in the
      /// OBConversion instance.
      /// This method is primarily intended for scripting languages without "stream" classes
    func writeFile<T: MKBase>(_ pOb: T, _ filePath: String) -> Bool {
        if pOutFormat == nil {
            //attempt to autodetect format
            pOutFormat = MKConversion.formatFromExt(filePath, &outFormatGzip)
        }
        var ofs: OutputFileHandler
        do {
            ofs = try OutputFileHandler(path: URL(string: filePath)!, mode: "w")
        } catch {
            print("Cannot write to \(filePath)")
            return false
        }
        setOutStream(ofs, true)
        // Set/reset the Index so that any initialization code
        // in the formatters will be executed.
        Index = 0
        // We can't touch the Last property because only the caller
        // knows if the first molecule is also the last molecule.
        return write(pOb)
    }
    
    
    /// @brief Manually closes and deletes the output stream
    /// The file is closed anyway when in the OBConversion destructor or when WriteFile
    /// is called again.
    /// \since version 2.1
    func closeOutFile() {
        setOutStream(nil)
    }
    
    /// @brief Reads an object of a class derived from OBBase into pOb.
    
    /// Part of "API" interface.
    /// The input stream can be specified and the change is retained in the OBConversion instance
    /// \return false and pOb=NULL on error
    open func read<T: MKBase>(_ pOb: inout T, _ pin: InputFileHandler? = nil) -> Bool {
//        if(pin != nil) {
          //for backwards compatibility, attempt to detect a gzip file
//    #ifdef HAVE_LIBZ
//            if(!inFormatGzip && pInFormat && zlib_stream::isGZip(*pin))
//          {
//            inFormatGzip = true;
//          }
//    #endif
//          SetInStream(pin, false);
//        }
//        TODO: Gzip input
        guard pInFormat != nil else { return false }
        guard pInput != nil else { return false }

        //mysterious line to ensure backwards compatibility
        //previously, even an open istream would have the gzip check applied
        //this meant that a stream at the eof position would end up in an error state
        //code has come to depend on this behavior
//        if(pInput-) pInput->get();

        // Set the locale for number parsing to avoid locale issues: PR#1785463
//        obLocale.SetLocale();

        // Also set the C++ stream locale
//        locale originalLocale = pInput->getloc(); // save the original
//        locale cNumericLocale(originalLocale, "C", locale::numeric);
//        pInput->imbue(cNumericLocale);

        // skip molecules if -f or -l option is set
        if (!SkippedMolecules) {
            Count = 0 // make sure it's 0
            if(!setStartAndEnd()) {
               return false
            }
            SkippedMolecules = true
        }

        // catch last molecule acording to -l
        Count += 1
        var success: Bool = false
        if (EndNumber==0 || Count <= EndNumber) {
            success = pInFormat!.readMolecule(pOb, self)
        }
        
        // return the C locale to the original one
        //        obLocale.RestoreLocale();
        // Restore the original C++ locale as well
        //        pInput->imbue(originalLocale);
        
        // If we failed to read, plus the stream is over, then check if this is a stream from ReadFile
        if (!success && !(pInput!.streamStatus == .error) && ownedInStreams.count > 0) {
            if let inFstream = ownedInStreams[0] as? InputFileHandler {
                do {
                    try inFstream.close()
                } catch {
                    print("unable to close inFstream")
                }
            }
            // We will free the stream later, but close the file now
        }
        
        return success
    }
    
    
    /// @brief Reads an object of a class derived from OBBase into pOb from the supplied string
    
    /// Part of "API" interface.
    /// \return false and pOb=NULL on error
    /// This method is primarily intended for scripting languages without "stream" classes
    /// Any existing input stream will be replaced by stringstream.
    @discardableResult
    open func readString<T: MKBase>(_ pOb: inout T, _ input: String) -> Bool {
        let inStream = InputStringStream(_wrappedBuffer: input)
        setInStream(inStream, true)
        return read(&pOb)
    }
    
    /// @brief Reads an object of a class derived from OBBase into pOb from the file specified
    
    /// Part of "API" interface.
    /// The output stream is changed to the supplied file and the change is retained in the
    /// OBConversion instance. For multi-molecule files, the remaining molecules
    /// can be read by repeatedly calling the Read() method.
    /// \return false and pOb=NULL on error
    /// This method is primarily intended for scripting languages without "stream" classes
    @discardableResult
    open func readFile<T: MKBase>(_ pOb: inout T, _ filePath: String) -> Bool {
        if pInFormat == nil {
            //attempt to autodetect format
            pInFormat = MKConversion.formatFromExt(filePath, &inFormatGzip)
        }
        // save the filename
        inFilename = filePath
        var ifs: InputFileHandler
        do {
            ifs = try InputFileHandler(path: URL(string: inFilename)!, mode: "r")

        } catch {
            print("Cannot read from \(inFilename)")
            return false
        }
        // libz support coming in the future 
        setInStream(ifs, true)
        return read(&pOb)
    }

    /// Part of the "Convert" interface.
    /// Open the files and update the streams in the OBConversion object.
    /// This method is primarily intended for scripting languages without "stream" classes
    /// and will usually followed by a call to Convert().
    /// Will set format from file extension if format has not already been set.
    /// Files will be opened even if format cannot be determined, but not if file path is empty.
    /// \return false if unsuccessful.
    @discardableResult
    func openInAndOutFiles(_ infilepath: String, _ outfilepath: String) -> Bool {
        if pInFormat == nil {
            //attempt to autodetect format
            pInFormat = MKConversion.formatFromExt(infilepath, &inFormatGzip)
        }
        // TODO: these might need to be rb instead of r
        var ifstream: InputFileHandler
        do {
            ifstream = try InputFileHandler(path: URL(string: infilepath)!, mode: "r")
        } catch {
            print("Cannot read from \(infilepath)")
            return false
        }
        setInStream(ifstream, true)
        inFilename = infilepath

        if outfilepath.isEmpty { // don't open an outfile with an emptry name
            return true
        }

        if pOutFormat == nil {
            //attempt to autodetect format
            pOutFormat = MKConversion.formatFromExt(outfilepath, &outFormatGzip)
        }

        var ofstream: OutputFileHandler
        do {
            ofstream = try OutputFileHandler(path: URL(string: outfilepath)!, mode: "w")
        } catch {
            print("Cannot write to \(outfilepath)")
            return false
        }
        setOutStream(ofstream, true)
        outFilename = outfilepath

        return true
    }
    
    /// @brief Sends a message like "2 molecules converted" to clog
    /// The type of object is taken from the TargetClassDescription
    /// of the specified class (or the output format if not specified)and
    /// is appropriately singular or plural.
    func reportNumberConverted(_ count: Int, _ pFormat: inout MKFormat?) {
        //Send info message to clog. This constructed from the TargetClassDescription
        //of the specified class (or the output format if not specified).
        //Get the last word on the first line of the description which should
        //be "molecules", "reactions", etc and remove the s if only one object converted
        if pFormat == nil {
            pFormat = pOutFormat
        }
        var objectname = pFormat!.targetClassDescription()
        let pos = objectname.firstIndex(of: Character("\n"))
        if pos == nil {
            objectname = String(objectname[..<objectname.endIndex])
        } else {
            objectname = String(objectname[..<pos!])
        }
        if count == 1 {
            objectname = String(objectname[..<objectname.index(before: objectname.endIndex)])
        }
        let pos2 = objectname.lastIndex(of: Character(" "))
        if pos2 == nil {
            objectname = String(objectname[objectname.startIndex...])
        } else {
            objectname = String(objectname[objectname.index(after: pos2!)...])
        }
        print("\(count) \(objectname) converted")
    }
    
    /// \return the number of objects in the InputFileHandler,
    /// or -1 if error or if SkipObjects for the input format is not implemented
    /// Adjusts for the value of -f and -l options (first and last objects).
    func numInputObjects() -> Int {
        
        guard let ifs = getInStream() else {
            return -1
        }
        let pos = try! ifs.tellg()

        //check that the input format supports SkipObjects()
        if getInFormat()!.skipObjects(0, self) == 0 {
            print("Input format does not have a SkipObjects function.")
            return -1
        }
        //counts objects only between the values of -f and -l options
        var nfirst = 1
        var nlast = Int.max
        if let p = isOption("f", .GENOPTIONS) { // extra parens to indicate truth value
            nfirst = Int(p)!
        }
        if let p = isOption("l", .GENOPTIONS) { // extra parens to indicate truth value
            nlast = Int(p)!
        }
        try! ifs.seekg(to: 0) //rewind
        //Compressed files currently show an error here.***TAKE CHANCE: RESET ifs****
        // ifs.clear();

        guard let pFormat = getInFormat() else {
            return -1
        }
        var count: Int = 0
        //skip each object but stop after nlast objects
        // while(ifs && pFormat->SkipObjects(1, this)>0  && count<nlast)
        // ++count;
        while !ifs.isEOF && pFormat.skipObjects(1, self) > 0 && count < nlast {
            count += 1
        }
        // ifs.clear(); //clear eof
        try! ifs.seekg(to: pos) //restore old position

        count -= nfirst-1
        return count
    }

//    MARK: Protected methods
    ///Replaces * in BaseName by InFile without extension and path
    static func batchFileName(_ baseName: String, _ inFile: String) -> String {
        //Replaces * in BaseName by InFile without extension and path
        var ofname = baseName
        if let pos = ofname.firstIndex(of: Character("*")) {
            //Replace * by input filename
            let posdot: String.Index = inFile.lastIndex(of: Character(".")) ?? inFile.endIndex
            // } else  {
            //     // if libz support is added, this will need to be changed
            // }
            // ofname.replace(pos,1, InFile, posname+1, posdot-posname-1);
            if let posname = inFile.lastIndex(where: { $0 == Character("/") || $0 == Character("\\") }) {
                ofname.replaceSubrange(pos...pos, with: inFile[posname...posdot])
            } else {
                ofname.replaceSubrange(pos...pos, with: inFile[...posdot])
            }
        }
        return ofname
    }

    ///Replaces * in BaseName by Count
    static func incrementedFilename(_ baseName: String, _ Count: Int) -> String {
        var ofname = baseName
        if let pos = ofname.firstIndex(of: Character("*")) {
            let num = String(Count)
            ofname.replaceSubrange(pos...pos, with: num) // TODO: this might need to be an insert instead
        }
        return ofname
    }

    ///Checks for misunderstandings when using the -m option
    static func checkForUnintendedBatch(_ inFile: String, _ outFile: String) -> Bool {
        let infile1 = inFile.substring(toIndex: (inFile.lastIndex(of: ".")?.utf16Offset(in: inFile))!)
        let infile2 = outFile.substring(toIndex: (outFile.lastIndex(of: ".")?.utf16Offset(in: outFile))!)
        
        if infile1 == infile2 {
            print("ERROR: This was a batch operation. For splitting, use non-empty base name for the output files")
        }
        
        if inFile == outFile { return false }
        return true
    }
        
    func setStartAndEnd() -> Bool {
        var TempStartNumber = 0
        var p: String? = isOption("f", .GENOPTIONS)
        if p != nil {
            do {
                StartNumber = try Int(p!, format: .number)
            } catch {
                print("StartNumber could not be parsed from \(p!)", FileHandle.standardError)
            }
            if StartNumber > 1 {
                TempStartNumber = StartNumber
//                 Try to skip objects now
                let ret = pInFormat!.skipObjects(StartNumber-1, self)
                if ret == -1 { return false } // error?
                if ret == 1 {
                    Count = StartNumber - 1
                    StartNumber = 0
                }
            }
        }
        
        p = isOption("l", .GENOPTIONS)
        
        if p != nil {
            do {
                EndNumber = try Int(p!, format: .number)
            } catch {
                print("StartNumber could not be parsed from \(p!)", FileHandle.standardError)
            }
            if (TempStartNumber != 0) && EndNumber<TempStartNumber {
                EndNumber = TempStartNumber
            }
        }
        
        return true
    }

    fileprivate static func optionParamArray<T>(_ typ: Option_Type, completion: @escaping (inout MKAMapType) -> (T)) -> T {
        do {
            var vals = try MKConversion.opa[typ]!.mapValues({ try Int($0, format: .number)})
            return completion(&vals)
        } catch {
            fatalError("Unable to convert params. Need to upgrade parser or error handler")
        }
    }
    
    func openAndSetFormat(_ setFormat: Bool, _ iss: InputFileHandler, _ ss: InputFileHandler? = nil) -> Bool {
        // opens file using inFilename and sets pInFormat if requested
        if ss != nil && inFilename[0] == "-" {
            // inFilename is actually -:SMILES
            inFilename = inFilename.substring(fromIndex: 2) // cut out "-:"
            if setFormat || setInFormat("smi") {
                return true
            }
        } else if !setFormat {
            pInFormat = MKConversion.formatFromExt(inFilename, &inFormatGzip)
            if pInFormat == nil {
                if let punctuationPoint = inFilename.lastIndex(of: ".") {
                    let ext = inFilename.substring(fromIndex: punctuationPoint)
                    var errorMsg = "Could not parse input format \(ext) for file \(inFilename)"
                    MKLogger.throwError(#function, errorMsg: errorMsg)
                }
                return false
            }
        }
        
        do {
            try iss.open(inFilename, mode: "r")
        } catch {
            print(error)
            return false
        }
        
        return true
    }
    
}


// Comments retained from OpenBabel
/**OBConversion maintains a list of the available formats,
provides information on them, and controls the conversion process.

A conversion is carried out by the calling routine, usually in a
user interface or an application program, making an instance of
OBConversion. It is loaded with the in and out formats, any options
and (usually) the default streams for input and output. Then either
the Convert() function is called, which allows a single input file
to be converted, or the extended functionality of FullConvert()
is used. This allows multiple input and output files, allowing:
- aggregation      - the contents of many input files converted
and sent to one output file;
- splitting        - the molecules from one input file sent to
separate output files;
- batch conversion - each input file converted to an output file.

These procedures constitute the "Convert" interface. OBConversion
and the user interface or application program do not need to be
aware of any other part of OpenBabel - mol.h is not \#included. This
allows any chemical object derived from OBBase to be converted;
the type of object is decided by the input format class.
However,currently, almost all the conversions are for molecules of
class OBMol.
///
OBConversion can also be used with an "API" interface
called from programs which manipulate chemical objects. Input/output is
done with the Read() and Write() functions which work with any
chemical object, but need to have its type specified. (The
ReadMolecule() and WriteMolecule() functions of the format classes
can also be used directly.)


Example code using OBConversion

<b>To read in a molecule, manipulate it and write it out.</b>

Set up an istream and an ostream, to and from files or elsewhere.
(cin and cout are used in the example). Specify the file formats.

@code
OBConversion conv(&cin,&cout);
if(conv.SetInAndOutFormats("SMI","MOL"))
{
    OBMol mol;
    if(conv.Read(&mol))
    // ...manipulate molecule

    conv->Write(&mol);
}
@endcode

A two stage construction is used to allow error handling
if the format ID is not recognized. This is necessary now that the
formats are dynamic and errors are not caught at compile time.
OBConversion::Read() uses a pointer to OBBase, so that, in addition
to OBMol, other kinds of objects, such as reactions, can also be handled
if the format routines are written appropriately.

<b>To make a molecule from a SMILES string.</b>
@code
std::string SmilesString;
OBMol mol;
stringstream ss(SmilesString)
OBConversion conv(&ss);
if(conv.SetInFormat("smi") && conv.Read(&mol))
    // ...
@endcode

An alternative way is more convenient if using bindings from another language:
@code
std::string SmilesString;
OBMol mol;
OBConversion conv;
if(conv.SetInFormat("smi") && conv.ReadString(&mol, SmilesString))
    // ...
@endcode

<b>To do a file conversion without manipulating the molecule.</b>

@code
#include <openbabel/obconversion.h> //mol.h is not needed
...set up an istream is and an ostream os
OBConversion conv(&is,&os);
if(conv.SetInAndOutFormats("SMI","MOL"))
{
    conv.AddOption("h",OBConversion::GENOPTIONS); //Optional; (h adds expicit hydrogens)
    conv.Convert();
}
@endcode

<b>To read a multi-molecule file if using bindings from another language</b>

The first molecule should be read using ReadFile, and subsequent molecules using Read,
as follows:
@code
#include <openbabel/obconversion.h> //mol.h is not needed
OBConversion conv;
OBMol mol;
bool success = conv.SetInFormat("sdf");
if(success)
{
    bool notatend = conv.ReadFile(&mol, "myfile.sdf");
    // Do something with mol
while(notatend)
{
        notatend = conv.Read(&mol);
    // Do something with mol
}
}
@endcode

<b>To add automatic format conversion to an existing program.</b>

The existing program inputs from the file identified by the
const char* filename into the istream is. The file is assumed to have
a format ORIG, but other formats, identified by their file extensions,
can now be used.

@code
ifstream ifs(filename); //Original code

OBConversion conv;
OBFormat* inFormat = conv.FormatFromExt(filename);
OBFormat* outFormat = conv.GetFormat("ORIG");
istream* pIn = &ifs;
stringstream newstream;
if(inFormat && outFormat)
{
    conv.SetInAndOutFormats(inFormat,outFormat);
    conv.Convert(pIn,&newstream);
    pIn=&newstream;
}
//else error; new features not available; fallback to original functionality

...Carry on with original code using pIn
@endcode
*/
