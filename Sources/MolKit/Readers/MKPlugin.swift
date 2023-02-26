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
typealias PluginMapType<T> = Dictionary<String,T>


protocol MKPluginProtocol {
    
    associatedtype T = AnyClass
    
    static var Default: T? { get set }
    static var map: PluginMapType<T> { get set }
    var _id: String { get set }
    func getID() -> String
    
    static func findType(_ ID: String?) -> (any MKPluginProtocol)?
    
    ///Required description of a sub-type
    func description() -> String?
    
    ///Redefined by each plugin type: "formats", "fingerprints", etc.
    func typeID() -> String
    
    ///Write information on a plugin class to the string txt.
      ///Return false if not written.
      ///The default implementation outputs:
      /// the ID, a tab character, and the first line of the Description.
      ///The param string can be used in derived types to provide different outputs.
    static func display(_ txt: inout String, _ param: inout String, _ ID: String?) -> Bool
    
    ///Make a new instance of the class.
      ///See OpTransform, OBGroupContrib, SmartsDescriptor classes for derived versions.
      ///Usually, the first parameter is the classname, the next three are
      ///parameters(ID, filename, description) for a constructor, and the rest data.
    func makeInstance(_ v: [String]) -> (any MKPluginProtocol)?
    
//    ///Get a pointer to a plugin from its type and ID. Return NULL if not found.
//    ///If Type is NULL, search all types. Not cast to Type*
//    static func getPlugin(_ type: String?, _ ID: String) -> MKPluginProtocol?
    
    ///Returns the map of the subtypes
    static func getMap() -> PluginMapType<T>
    
    ///\brief Returns a reference to the map of the plugin types.
    /// Is a function rather than a static member variable to avoid initialization problems.
    
}

class MKPlugin: MKPluginProtocol {
   
    var _id: String = ""
    
    static var allPluginsLoaded: Int = 0
    static var pluginMap: PluginMapType<any MKPluginProtocol> = PluginMapType<any MKPluginProtocol>()
    
    static var map: PluginMapType<MKPlugin> = PluginMapType<MKPlugin>()
    static var Default: MKPlugin?
    
    init() {}
    
    init(_ id: String, _ isDefault: Bool) {
        self._id = id
        if isDefault || MKPlugin.map.isEmpty {
            MKPlugin.Default = self
        }
        if MKPlugin.map.map({ $0.0 == _id ? 1 : 0}).reduce(0, +) == 0 {
            MKPlugin.map[_id] = self
            MKPlugin.pluginMap[typeID()] = self
        }
    }
    
    func getID() -> String { return _id }
    
    
    func description() -> String? {
        return nil
    }
    
    func typeID() -> String {
        return "plugins"
    }
    
    static func display(_ txt: inout String, _ param: inout String, _ ID: String?) -> Bool {
//        MARK IMPL and throw proper error
        return false
    }
    
    func makeInstance(_ v: [String]) -> (any MKPluginProtocol)? {
        return nil
    }
    
    static func getPlugin(_ Type: String?, _ ID: String) -> (any MKPluginProtocol)? {
        if Type != nil || Type != "" {
            return MKPlugin.baseFindType(MKPlugin.getTypeMap(Type!), ID)
        }
        
        if MKPlugin.allPluginsLoaded == 0 {
            MKPlugin.loadAllPlugins()
        }
        
        //When Type==NULL, search all types for matching ID and stop when found
        for plug in MKPlugin.pluginMap {
            var mapp = (type(of: plug.value)).getMap()
            if (type(of: mapp) == PluginMapType<MKPlugin>.self) {
                let res = MKPlugin.baseFindType((mapp as! PluginMapType<MKPlugin>), ID)
                if res != nil { return res }
            }
        }
        return nil
    }
//
    
    static func getMap() -> PluginMapType<MKPlugin> {
        return MKPlugin.map
    }

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
    static func getTypeMap(_ plugID: String) -> PluginMapType<any MKPluginProtocol> {
        if MKPlugin.allPluginsLoaded == 0 {
            MKPlugin.loadAllPlugins()
        }
        if MKPlugin.pluginMap.first(where: { (key: String, value: any MKPluginProtocol) in
            key == plugID && type(of: value) == MKPlugin.self
        }) != nil {
            return (MKPlugin.getMap())
        } else {
            return MKPlugin.pluginMap //error: type not found; return plugins map
        }
    }
    
    static func findType(_ ID: String?) -> (any MKPluginProtocol)? {
        if ID == nil {
            return MKPlugin.Default
        }
        return MKPlugin.baseFindType(getMap(), ID!)
    }
    
    ///\brief Returns the type with the specified ID, or NULL if not found.
    ///Needs to be cast to the appropriate class in the calling routine.
    static func baseFindType(_ map: PluginMapType<any MKPluginProtocol>, _ ID: String) -> (any MKPluginProtocol)? {
        // Make sure the plugins are loaded
        if MKPlugin.allPluginsLoaded == 0 {
            MKPlugin.loadAllPlugins()
        }
        
        if let res = map.first(where: { (key: String, value: any MKPluginProtocol) in
            key == ID
        })?.value {
            return res
        } else {
            return nil
        }
    }
    
}
