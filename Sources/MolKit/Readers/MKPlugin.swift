//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/24/23.
//

import Foundation
import Collections

//Case insensitive string comparison for PluginMapType key **
//Maps of thistype are used to store
// (a)a list of the plugin types in OBPlugin, and
// (b)a list of the sub-types in each type class derived from OBPlugin.
public typealias PluginMapType<T: MKPlugin> = OrderedDictionary<String, T>

protocol MKPluginProtocol {
        
    associatedtype T: MKPlugin
    
    static var Default: T? { get set }
    static var map: PluginMapType<T> { get set }
        
    static func findType(_ ID: String?) -> T?
//    static func findType(_ ID: String?) -> MKPlugin? {
//        if ID == nil {
//            return T.Default
//        }
//        return MKPlugin.baseFindType(getMap(), ID!)
//    }

        
//    ///Get a pointer to a plugin from its type and ID. Return NULL if not found.
//    ///If Type is NULL, search all types. Not cast to Type*
//    static func getPlugin(_ type: String?, _ ID: String) -> MKPluginProtocol?
    
    init(_ id: String, _ isDefault: Bool)
//        self._id = id
//        if isDefault || T.map.isEmpty {
//            T.Default = self
//        }
//        if T.map.map({ $0.0 == _id ? 1 : 0}).reduce(0, +) == 0 {
//            T.map[_id] = self
//            MKPlugin.pluginMap[typeID()] = self
//        }
//    }
    
    ///Returns the map of the subtypes
    func getMap() -> PluginMapType<T>
    
    ///\brief Returns a reference to the map of the plugin types.
    /// Is a function rather than a static member variable to avoid initialization problems.
}

public class MKPlugin: Equatable {
        
    var _id: String = ""
    
    static public var allPluginsLoaded: Int = 0
    static public var pluginMap: PluginMapType<MKPlugin> = PluginMapType<MKPlugin>()
            
    func getID() -> String { return _id }
    
    func description() -> String? { return nil }
    
    func typeID() -> String {
        return "plugins"
    }
    
    func display(_ txt: inout String, _ param: inout String?, _ ID: String?) -> Bool {
        if ID != nil {
            txt = ID!
        } else {
            txt = getID()
        }
        txt += "    "
//        MARK: needs locale updating here
        if param != nil && param == "verbose" {
            guard let description = description() else { return false }
            txt += description
            txt += "\n"
        } else {
            guard let description = description() else { return false }
            txt += MKPlugin.firstLine(description)
        }
        return true
    }
    
    func makeInstance(_ v: [String]) -> MKPlugin? {
        fatalError("not implemented in base class")
    }
    
//    init(_ id: String, _ isDefault: Bool) {
//        self._id = id
//        if isDefault || MKPlugin.map.isEmpty {
//            MKPlugin.Default = self
//        }
//        if MKPlugin.pluginMap.map({ $0.0 == _id ? 1 : 0}).reduce(0, +) == 0 {
//            T.map[_id] = self
//            MKPlugin.pluginMap[typeID()] = self
//        }
//    }
    
    static func getPlugin(_ Type: String?, _ ID: String) -> MKPlugin? {
        if Type != nil || Type != "" {
            return MKPlugin.baseFindType(MKPlugin.getTypeMap(Type!), ID)
        }
        
        if MKPlugin.allPluginsLoaded == 0 {
            MKPlugin.loadAllPlugins()
        }
        
        //When Type==NULL, search all types for matching ID and stop when found
        for plug in MKPlugin.pluginMap {
            let mapp = plug.value.getMap()
            if (type(of: mapp) == PluginMapType<MKPlugin>.self) {
                let res = MKPlugin.baseFindType(mapp!, ID)
                if res != nil { return res }
            }
        }
        return nil
    }
    
    func getMap<T: MKPlugin>() -> PluginMapType<T>? { return MKPlugin.pluginMap as! PluginMapType<T>? }

    ///Load all plugins (formats, fingerprints, forcefields etc.)
    static func loadAllPlugins() {
        
//        MARK: Add more plugins here or elsewhere?!
        
        MKPlugin.allPluginsLoaded = 1
        let pdef = MKPlugin.getPlugin("loaders", "define")
        if pdef != nil {
            let vec: [String] = ["", "define", "plugindefines.txt"]
            pdef?.makeInstance(vec)
        }
    }
    
    ///Returns the map of a particular plugin type, e.g. GetMapType("fingerprints")
    static func getTypeMap(_ plugID: String) -> PluginMapType<MKPlugin> {
        if MKPlugin.allPluginsLoaded == 0 {
            MKPlugin.loadAllPlugins()
        }
        if let plug = MKPlugin.pluginMap.first(where: { (key: String, value: MKPlugin) in
            key == plugID
        }) {
            return plug.value.getMap()!
        } else {
            return MKPlugin.pluginMap //error: type not found; return plugins map
        }
    }
    
    static func getPluginIterator(_ pluginID: String) -> PluginMapType<MKPlugin>? {
        if pluginID != "plugins" || MKPlugin.getTypeMap(pluginID) != MKPlugin.pluginMap {
            let typemap = MKPlugin.getTypeMap(pluginID)
            return typemap
        } else {
            return MKPlugin.pluginMap
        }
    }
    
    ///\brief Returns the type with the specified ID, or NULL if not found.
    ///Needs to be cast to the appropriate class in the calling routine.
    static func baseFindType<T>(_ map: PluginMapType<T>, _ ID: String) -> MKPlugin? where T: MKPlugin {
        // Make sure the plugins are loaded
        if MKPlugin.allPluginsLoaded == 0 {
            MKPlugin.loadAllPlugins()
        }
        
        if let res = map.first(where: { (key: String, value: MKPlugin) in
            key == ID
        })?.value {
            return res
        } else {
            return nil
        }
    }
    
    func List(_ pluginID: String, _ param: inout String?, _ os: inout any TextOutputStream) { // Needs to be a real output stream
        var vlist: [String] = []
        if !MKPlugin.ListAsVector(pluginID, &param, &vlist) {
            os.write("\(pluginID) is not a recognized plugin type. Those with instances of sub-types loaded are:\n")
        }
//        for each string in vlist, append it to os with "\n" as the separator
        for val in vlist {
            os.write(val)
            os.write("\n")
        }
    }
    
    @discardableResult
    static func ListAsVector(_ pluginID: String?, _ param: inout String?,_ vlist: inout [String]) -> Bool {
        var ret = true
        
        // Make sure the plugins are loaded
        if (MKPlugin.allPluginsLoaded == 0) {
            MKPlugin.loadAllPlugins()
        }
        
        if(pluginID != nil) {
            if(pluginID!.first != "0" && pluginID != "plugins") {
                //List the sub classes of the specified type
                guard let itr = MKPlugin.pluginMap[pluginID!] else { return false }
                let onlyIDs: Bool = param == "ids"
                    //Get map of plugin type (like OBFingerprint) and output its contents
                guard let Map = itr.getMap() else { return false }
                for itrM in Map {
                    if (itrM.key == "_" ) { continue } //no listing when ID starts with '_'
                    if (onlyIDs) {
                        vlist.append(itrM.key)
                    } else {
                        var txt: String = ""
                        if(itrM.value.display(&txt, &param, itrM.key)) {
                            vlist.append(txt)
                        }
                    }
                }
                return true
            }
        } else {
            ret = false //asked for a type not available; provide plugin types instead
            //List the plugin types
            for itr in MKPlugin.pluginMap {
                vlist.append(itr.key)
            }
        }
        
        return ret
    }
    
    func ListAsString(_ pluginID: String,_ param: inout String?) -> String {
        var ss: TextOutputStream = ""
        List(pluginID, &param, &ss)
        return ss as! String
    }
    
    
    static func firstLine(_ text: String) -> String {
//        Return the first set of characters up until the newline character
        let newline = text.firstIndex(of: "\n") ?? text.endIndex
        return String(text[..<newline])
    }
    
    
    public static func == (_ lhs: MKPlugin, _ rhs: MKPlugin) -> Bool {
        return lhs.typeID() == rhs.typeID() && lhs._id == rhs._id
    }
    
}
