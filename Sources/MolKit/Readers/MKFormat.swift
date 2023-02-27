

class MKFormat: MKPlugin, MKPluginProtocol {
    
        
    private let NOTREADABLE     = 0x01
    private let READONEONLY     = 0x02
    private let READBINARY      = 0x04
    private let ZEROATOMSOK     = 0x08
    private let NOTWRITABLE     = 0x10
    private let WRITEONEONLY    = 0x20
    private let WRITEBINARY     = 0x40
    private let READXML         = 0x80
    private let DEPICTION2D     = 0x100
    private let DEFAULTFORMAT   = 0x4000
    
        
    var pMime: String = ""
    
    static var Default: MKFormat?
    static var map: PluginMapType<MKFormat> = PluginMapType<MKFormat>()
    static var formatsMIMEMap: PluginMapType = PluginMapType<MKFormat>()
    
//    override init() {
////        if isDefault || MKFormat.map.isEmpty {
////            MKFormat.Default = self
////        }
////        if MKFormat.map.map({ $0.0 == _id ? 1 : 0}).reduce(0, +) == 0 {
////            MKFormat.map[_id] = self
////            MKPlugin.pluginMap[typeID()] = self
////        }
//
//    }
    
    override init() {
        super.init()
    }
    
    required init(_ id: String, _ isDefault: Bool) {
        super.init()
//        _id = "formats"
    }
    
    override func typeID() -> String {
        return "formats"
    }
    
    func getMap() -> PluginMapType<MKFormat> {
        return MKFormat.map
    }
    
    func registerFormat(_ ID: String, _ MIME: String? = nil) -> Int {
        var mp = getMap()
        mp.updateValue(self, forKey: ID)
        if MIME != nil {
            MKFormat.formatsMIMEMap[MIME!] = self
        }
        if ((flags() & DEFAULTFORMAT) != 0) {
            MKFormat.Default = self
        }
        
        //ensure "formats" is registered as a plugin
        MKPlugin.pluginMap[typeID()] = self
        self._id = ID
        return getMap().count
    }
    
    /// @brief The "API" interface Read function.
    /// Reads a single object.
    /// Does not make a new object on the heap;
    /// can be used with a pointer to an chem object on the heap or the stack.
    func readMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        print("HIER")
        print("Not a valid input format")
        return false
    }
    
    /// @brief The "Convert" interface Read function.
    /// Possibly reads multiple new objects on the heap and subjects them
    /// to its DoTransformations() function, which may delete them again.
    /// Sends result to pConv->AddChemObject()
    func readChemObject(_ pConv: MKConversion) -> Bool {
        print("Not a valid input format")
        return false
    }
    
    /// @brief The "API" interface Write function.
    /// Writes a single object
    /// Does not delete the object;
    /// can be used with a pointer to an chem object on the heap or the stack.
    /// \return false on error.
    func writeMolecule(_ pOb: MKBase, _ pConv: MKConversion) -> Bool {
        print("Not a valid output format")
        return false
    }
    
    /// @brief The "Convert" interface Write function.
    /// Writes a single object
    /// Deletes the object after writing
    /// \return false on error
    func writeChemObject(_ pConv: MKConversion) -> Bool {
        print("Not a valid output format")
        return false
    }
    
    /// @brief Information on this format. Printed out in response to -Hxxx option where xxx id the id of the format.
    
    /// Must be provided by each format class.
    /// Can include a list of command line Options. These may be used to construction
    /// check boxes, radio buttons etc for GUI interface.
    override func description() -> String? {
        return ""
    }
    
    /// @brief A decription of the chemical object converted by this format.
    /// If not provided, the object type used by the default format is used (usually MKMol).
    func targetClassDescription() -> String {
        if let res = T.findType(nil) {
            if type(of: res) == MKFormat.self {
                return res.targetClassDescription()
            }
        }
        return ""
    }
    /// \return the type of chemical object used by the format.
    /// Defaults to that used by the default format. Useful for checking
    /// that a format can handle a particular object.
    func getType() -> String { return "" } // ??
    
    /// @brief Web address where the format is defined.
    func specificationURL() -> String { return "" }
    
    /// @brief Chemical MIME type associated with this file type (if any)
    func getMIMEType() -> String { return pMime }
    
    
    override func makeInstance(_ v: [String]) -> MKFormat? {
        return nil
    }
    
    override func display(_ txt: inout String, _ param: inout String, _ ID: String?) -> Bool {
        return false
    }
    
    static func findType(_ ID: String?) -> MKFormat? {
        if ID == nil {
            return MKFormat.Default
        }
        return MKPlugin.baseFindType(MKFormat.map, ID!) as? MKFormat
    }
    
    /// @brief Decribes the capabilities of the format (Read only etc.)
    /// Currently, can be a bitwise OR of any of the following
    /// NOTREADABLE READONEONLY NOTWRITABLE WRITEONEONLY DEFAULTFORMAT
    /// READBINARY WRITEBINARY READXML
    func flags() -> Int { return 0 }
    
    /// @brief Skip past first n objects in input stream (or current one with n=0)
    /// \return 1 on success, -1 on error and 0 if not implemented
    func skipObjects(_ n: Int, _ pConv: MKConversion) -> Int { return 0 }
    
    static func formatFromMIME(_ MIME: String) -> MKFormat? {
        return MKFormat.formatsMIMEMap[MIME]
    }
    
    /* Functions provided by the superclass

      ///Constructor that registers the ID of the format
      Not currently used for formats
      MKFormat(const char* ID, bool IsDefault=false);

      ///Returns the sub-type associated with the ID, or the default subtype if ID NULL or empty.
      static MKFormat* FindType(const char* ID);
     */
     // We will need to see if this class acts appropriately. Else, we will need to override those methods, which is a pain
    
}
