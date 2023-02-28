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
}

//*************************************************
/// @brief Class to convert from one format to another.

class MKConversion {

    
    fileprivate typealias MKAMapType = OrderedDictionary<String, Int>
    
    // Utilize base datatype as String, could be replaced with a class later that handles streaming to a live Stream and uses InputStream and OutputStream and conforms to StreamDelegate
//    MARK: StreamState
    private struct StreamState {
        
        var pStream: String?
        var owndedStreams: [String] = []
        
        init() {
            self.pStream = nil
        }
        
        func pushInput(_ conv: MKConversion) {}
        func popInout(_ conv: MKConversion) {}
        
        func pushOutput(_ conv: MKConversion) {}
        func popOuput(_ conv: MKConversion) {}
        
        
        func isSet() -> Bool { return self.pStream != nil }
    }
    
//    MARK: Class Variables
    
    var inFilename: String = ""
    var outFilename: String = ""
    
    var pInput: String?  //input stream, may be filtered
    var ownedInStreams: [String] = []
    
    var pOutput: String?  //output stream, may have filters applied
    var ownedOutStreams: [String] = []
    
    static var pDefaultFormat: MKFormat?
    var pInFormat: MKFormat?
    var pOutFormat: MKFormat?
    
    var OptionsArray: OrderedDictionary<Option_Type, Dictionary<String, String?>> = .init(minimumCapacity: 3)
    
    var Index: Int?
    var StartsNumber: Int?
    var EndNumber: Int?
    var Count: Int?
    var m_IsFirstInput: Bool 
    var m_IsLast: Bool 
    var MoreFilesToCome: Bool 
    var OneObjectOnly: Bool 
    var ReadyToInput: Bool 
    var SkippedMolecules: Bool     /// skip molecules using -f and -l

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
    
//
//    init(_ inFilename: String, _ outFilename: String) {
//        self.pInput = nil
//        self.pOutput = nil
//        self.pInFormat = nil
//        self.pOutFormat = nil
//        self.Index = 0
//        self.StartsNumber = 1
//        self.EndNumber = 0
//        self.Count = -1
//        self.m_IsFirstInput = true
//        self.m_IsLast = true
//        self.MoreFilesToCome = false
//        self.OneObjectOnly = false
//        self.SkippedMolecules = false
//        self.inFormatGzip = false
//        self.outFormatGzip = false
//        self.pOb = nil
//        self.wInpos = 0
//        self.wInlen = 0
//        self.pAuxConv = nil
//        self.inFilename = inFilename
//        self.outFilename = outFilename
//
//        MKConversion.registerOptionParam("f", nil, 1, .GENOPTIONS)
//        MKConversion.registerOptionParam("l", nil, 1, .GENOPTIONS)
//
//        setInStream(inFilename)
//        setOutStream(outFilename)
//
//        openInAndOutFiles(inFilename, outFilename)
//    }
    internal init(inFilename: String = "", outFilename: String = "", pInput: String? = nil, ownedInStreams: [String] = [], pOutput: String? = nil, ownedOutStreams: [String] = [], pInFormat: MKFormat? = nil, pOutFormat: MKFormat? = nil, OptionsArray: OrderedDictionary<Option_Type, Dictionary<String, String?>> = .init(minimumCapacity: 3), Index: Int? = nil, StartsNumber: Int? = nil, EndNumber: Int? = nil, Count: Int? = nil, m_IsFirstInput: Bool, m_IsLast: Bool, MoreFilesToCome: Bool, OneObjectOnly: Bool, ReadyToInput: Bool, SkippedMolecules: Bool, inFormatGzip: Bool = false, outFormatGzip: Bool = false, pOb: MKBase? = nil, wInpos: Int? = nil, rInpos: Int? = nil, wInlen: Int? = nil, rInlen: Int? = nil, pAuxConv: MKConversion? = nil, SupportedInputFormat: [String] = [], SupportedOutputFormat: [String] = []) {
        self.inFilename = inFilename
        self.outFilename = outFilename
        self.ownedInStreams = ownedInStreams
        self.ownedOutStreams = ownedOutStreams
        self.OptionsArray = OptionsArray
        self.ReadyToInput = ReadyToInput
        self.rInpos = rInpos
        self.rInlen = rInlen
        self.SupportedInputFormat = SupportedInputFormat
        self.SupportedOutputFormat = SupportedOutputFormat
        
        self.pInput = nil
        self.pOutput = nil
        self.pInFormat = nil
        self.pOutFormat = nil
        self.Index = 0
        self.StartsNumber = 1
        self.EndNumber = 0
        self.Count = -1
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

        
        MKConversion.registerOptionParam("f", nil, 1, .GENOPTIONS)
        MKConversion.registerOptionParam("l", nil, 1, .GENOPTIONS)

        setInStream(inFilename)
        setOutStream(outFilename)
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
    static func registerFormat(_ ID: String, _ pFormat: MKFormat, _ MIME: String? = nil) -> Int {
        return pFormat.registerFormat(ID, MIME)
    }

    
    ///@brief Searches registered formats
    static func findFormat(_ ID: String) -> MKFormat {
        fatalError()
    }
    /// @brief Searches registered formats for an ID the same as the file extension
    static func formatFromExt(_ filename: String) -> MKFormat {
        fatalError()
    }
    static func formatFromExt(_ filename: String, _ isgzip: Bool) -> MKFormat {
        fatalError()
    }
   /// @brief Searches registered formats for a MIME the same as the chemical MIME type passed
    static func formatFromMIME(_ MIME: String) -> MKFormat {
        fatalError()
    }
    
    static func getNextFormat(_ itr: any IteratorProtocol<MKPlugin>, _ str: String, _ pFormat: MKFormat) -> Bool {
        fatalError()
    }
    
    /// @name Information
      //@{
    static func description() -> String {
        fatalError()
    } //generic conversion options
    //@}
    
    /// These return a filtered stream for reading/writing (possible filters include compression, decompression, and newline transformation)
    /// @name Parameter get and set
    //@{
    func getInStream() -> String { return pInput ?? "" }
    func getOutStream() -> String { return pOutput ?? "" }
    
    /// @brief Set input stream.  If takeOwnership is true, will deallocate when done.
    /// If isGzipped is true, will treat as a gzipped stream regardless of option settings,
    //  if false, then will be treated as gzipped stream only if z/zin is set.
    func setInStream(_ pIn: String, _ takeOwnership: Bool = false) {
        // clear and deallocate any existing streams
        self.pInput = pIn
    }
    
    func setOutStream(_ pOut: String, _ takeOwnership: Bool = false) {
        self.pOutput = pOut
    }
    
    /// Sets the formats from their ids, e g CML
    func setInAndOutFormats(_ inID: String, _ outID: String, ingzip: Bool = false, outgzip: Bool = false) -> Bool { return false }
    func setInAndOutFormats(_ pIn: MKFormat, _ pOut: MKFormat, ingzip: Bool = false, outgzip: Bool = false) -> Bool { return false }
    
    /// Sets the input format from an id e.g. CML
    func setInFormat(_ inID: String, isgzip: Bool = false) -> Bool {
        fatalError()
    }
    func setInFormat(_ pIn: MKFormat, isgzip: Bool = false) -> Bool {
        fatalError()
    }

    /// Sets the output format from an id e.g. CML
    func setOutFormat(_ outID: String, isgzip: Bool = false) -> Bool {
        fatalError()
    }
    func setOutFormat(_ pOut: MKFormat, isgzip: Bool = false) -> Bool {
        fatalError()
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
    func getTitle() -> String {
        fatalError()
    }

    ///@brief Extension method: deleted in ~OBConversion()
    func getAuxConv() -> MKConversion? { return pAuxConv }
    func setAuxConv(_ pConv: MKConversion) { pAuxConv = pConv }
    
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
        if var pos = OptionsArray[opttype]![opt] {
            return pos
        } else {
            return nil
        }
    }

    ///@brief Access the map with option name as key and any associated text as value
    func getOptions(_ opttyp: Option_Type) -> [String: String?] { return OptionsArray[opttyp]! }

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
        fatalError()
    }

    ///@brief For example -h takes 0 parameters; -f takes 1. Call in a format constructor.
    static func registerOptionParam(_ name: String, _ pFormat: MKFormat?, _ numParams: Int, _ typ: Option_Type = .OUTOPTIONS) {
        fatalError()
    }

    /// \return the number of parameters registered for the option, or 0 if not found
    static func getOptionParams(_ name: String, _ typ: Option_Type) -> Int {
        fatalError()
    }

    ///@brief Copies the options (by default of all types) from one OBConversion Object to another.
    func copyOptions(_ pConv: MKConversion, _ opttype: Option_Type = .ALL) {
        fatalError()
    }

    /// @name Supported file format
    //@{
    // @brief Set and return the list of supported input format
    func getSupportedInputFormat() -> [String] {
        fatalError()
    }
    // @brief Set and return the list of supported output format
    func getSupportedOutputFormat() -> [String] {
        fatalError()
    }
    //@}

    // MARK: Conversion
    //@{
    /// @brief Conversion for single input and output stream
    func convert(_ iss: String, _ oss: String) -> Int {
        fatalError()
    }

    /// @brief Conversion with existing streams
    func convert() -> Int {
        fatalError()
    }

    /// @brief Conversion with multiple input/output files:
    /// makes input and output streams, and carries out normal, batch, aggregation, and splitting conversion.
    func fullConvert(_ fileList: [String], _ outputFileName: String, _ outputFileList: [String]) -> Int {
        fatalError()
    }
    //@}

     /// @name Conversion loop control
    //@{
    ///< @brief Adds to internal array during input
    func addChemObject(_ pOb: MKBase) -> Int {
        fatalError()
    }
    ///< @brief Retrieve from internal array during output
    func getChemObject() -> MKBase? {
        fatalError()
    }
    ///< @brief True if no more objects to be output
    func isLast() -> Bool {
        fatalError()
    }
      ///< @brief True if the first input object is being processed
    func isFirstInput() -> Bool {
        fatalError()
    }
      ///< @brief Setwhether or not is the first input
    func setFirstInput(_ b: Bool = true) {
        fatalError()
    }
      ///< @brief Retrieves number of ChemObjects that have been actually output
    func getOutputIndex() -> Int {
        fatalError()
    }
      ///< @brief Sets output index (maybe to control whether seen as first object)
    func setOutputIndex(_ indx: Int) {
        fatalError()
    }
      ///<@brief Used with multiple input files. Off by default.
    func setMoreFilesToCome() {
        fatalError()
    }
      ///< @brief Used with multiple input files. Off by default.
    func setOneObjectOnly(_ b: Bool = true) {
        fatalError()
    }
      ///< @brief Synonym for SetOneObjectOnly()
    func setLast(_ b: Bool = true) { setOneObjectOnly(b) }
      ///< @brief True if no more files to be read
    func isLastFile() -> Bool { return !MoreFilesToCome }

    /// @brief Number of objects read and processed
    /// Incremented after options are processed, so 0 for first object.  Returns -1 if Convert interface not used. 
    func getCount() -> Int { return Count ?? 0}
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
    func write(_ pOb: MKBase, _ pout: String?) -> Bool {
        fatalError()
    }

    /// @brief Outputs an object of a class derived from OBBase as a string
      /// Part of "API" interface.
      /// The output stream is temporarily changed to the string and then restored
      /// This method is primarily intended for scripting languages without "stream" classes
      /// The optional "trimWhitespace" parameter allows trailing whitespace to be removed
      /// (e.g., in a SMILES string or InChI, etc.)
    func writeString(_ pOb: MKBase, _ trimWhitespace: Bool = true) -> String {
        fatalError()
    }
    
    /// @brief Outputs an object of a class derived from OBBase as a file (with the supplied path)
      /// Part of "API" interface.
      /// The output stream is changed to the supplied file and the change is retained in the
      /// OBConversion instance.
      /// This method is primarily intended for scripting languages without "stream" classes
    func writeFile(_ pOb: MKBase, _ filePath: String) -> Bool {
        fatalError()
    }
    
    
    /// @brief Manually closes and deletes the output stream
    /// The file is closed anyway when in the OBConversion destructor or when WriteFile
    /// is called again.
    /// \since version 2.1
    func closeOutFile() {
        fatalError()
    }
    
    /// @brief Reads an object of a class derived from OBBase into pOb.
    
    /// Part of "API" interface.
    /// The input stream can be specified and the change is retained in the OBConversion instance
    /// \return false and pOb=NULL on error
    func read(_ pOb: MKBase, _ pin: String? = nil) -> Bool {
        fatalError()
    }
    
    
    /// @brief Reads an object of a class derived from OBBase into pOb from the supplied string
    
    /// Part of "API" interface.
    /// \return false and pOb=NULL on error
    /// This method is primarily intended for scripting languages without "stream" classes
    /// Any existing input stream will be replaced by stringstream.
    func readString(_ pOb: MKBase, _ input: String) -> Bool {
        fatalError()
    }
    
    /// @brief Reads an object of a class derived from OBBase into pOb from the file specified
    
    /// Part of "API" interface.
    /// The output stream is changed to the supplied file and the change is retained in the
    /// OBConversion instance. For multi-molecule files, the remaining molecules
    /// can be read by repeatedly calling the Read() method.
    /// \return false and pOb=NULL on error
    /// This method is primarily intended for scripting languages without "stream" classes
    func readFile(_ pOb: MKBase, _ filePath: String) -> Bool {
        fatalError()
    }

    /// Part of the "Convert" interface.
    /// Open the files and update the streams in the OBConversion object.
    /// This method is primarily intended for scripting languages without "stream" classes
    /// and will usually followed by a call to Convert().
    /// Will set format from file extension if format has not already been set.
    /// Files will be opened even if format cannot be determined, but not if file path is empty.
    /// \return false if unsuccessful.
    func openInAndOutFiles(_ infilepath: String, _ outfilepath: String) -> Bool {
        fatalError()
    }
    
    /// @brief Sends a message like "2 molecules converted" to clog
    /// The type of object is taken from the TargetClassDescription
    /// of the specified class (or the output format if not specified)and
    /// is appropriately singular or plural.
    func reportNumberConverted(_ count: Int, _ pFormat: MKFormat? = nil) {
        fatalError()
    }
    
    /// \return the number of objects in the inputstream,
    /// or -1 if error or if SkipObjects for the input format is not implemented
    /// Adjusts for the value of -f and -l options (first and last objects).
    func numInputObjects() -> Int {
        fatalError()
    }

//    MARK: Protected methods
    ///Replaces * in BaseName by InFile without extension and path
    static func batchFileName(_ baseName: String, _ inFile: String) -> String {
        fatalError()
    }
    ///Replaces * in BaseName by Count
    static func incrementedFilename(_ baseName: String, _ Count: Int) -> String {
        fatalError()
    }
    ///Checks for misunderstandings when using the -m option
    static func checkForUnintendedBatch(_ inFile: String, _ outFile: String) -> Bool {
        fatalError()
    }
    
    func clearInStreams() {}
    
    func setStartAndEnd() -> Bool {
        fatalError()
    }

    fileprivate static func optionParamArray(_ typ: Option_Type) -> MKAMapType {
        fatalError()
    }
    
    func openAndSetFormat(_ setFormat: Bool, _ iss: String, _ ss: String? = nil) -> Bool {
        fatalError()
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
