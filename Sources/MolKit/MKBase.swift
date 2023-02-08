//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 1/31/23.
//

import Foundation

// Generate UUID in UInt format
func generateUUID() -> Int {
    let uuid = UUID().uuidString
    let uuidData = uuid.data(using: .utf8)!
    let uuidUInt = uuidData.withUnsafeBytes { $0.load(as: Int.self) }
    return uuidUInt
}

public enum MKBaseID {
    case NoId
    case _id(Int)

    var rawValue: Int {
        switch self {
        case .NoId:
            return -1
        case ._id(let id):
            return id
        }
    }
}

/// Concrete Iterators implement various traversal algorithms. These classes
/// store the current traversal position at all times.
public class MKIterator<T>: IteratorProtocol, Sequence {

    private let collection: [T]
    private var index = 0

    init(_ collection: [T]) {
        self.collection = collection
    }

    public func next() -> T? {
        defer { index += 1 }
        return index < collection.count ? collection[index] : nil
    }
}


public class MKBase: NSObject {
    
    var _vdata: [MKGenericData]? = nil
    var _flags: UInt = 0

    public override init() {
        super.init()
        self._vdata = []
        self._flags = 0
    }
    
    public init(_ base: MKBase) {
        super.init()
        self.clone(base)
    }
    
    public func clone(_ base: MKBase) {
        self._vdata = base._vdata
        self._flags = base._flags
    }

    func clear() { 
        self._vdata = []
        self._flags = 0
    }

    func getFlag() -> UInt {
        return self._flags
    }
    
    func setFlag(_ flag: UInt) {
        self._flags |= flag
    }
    
    func unsetFlag(_ flag: UInt) {
        self._flags &= ~flag
    }
    
    func hasFlag(_ flag: UInt) -> Bool {
        return (self._flags & flag) != 0
    }

    //! \return whether the generic attribute/value pair exists
    func hasData(_ attr: String) -> Bool {
        if let vdata = self._vdata {
            for data in vdata {
                if data.attr == attr {
                    return true
                }
            }
        }
        return false
    }
    //! \return whether the generic attribute/value pair exists, for a given MKGenericDataType
    func hasData(_ type: MKGenericDataType) -> Bool {
        if let vdata = self._vdata {
            for data in vdata {
                if data.type == type {
                    return true
                }
            }
        }
        return false
    }

    func dataSize() -> Int {
        if let vdata = self._vdata {
            return vdata.count
        }
        return 0
    }

    //! \return the first matching data for a given type from OBGenericDataType
    //!    or NULL if nothing matches
    func getData(_ type: MKGenericDataType) -> MKGenericData? {
        if let vdata = self._vdata {
            for data in vdata {
                if data.type == type {
                    return data
                }
            }
        }
        return nil
    }

    //! \return any data matching the given attribute name or NULL if nothing matches
    func getData(_ attr: String) -> MKGenericData? {
        if let vdata = self._vdata {
            for data in vdata {
                if data.attr == attr {
                    return data
                }
            }
        }
        return nil
    }

    //! \return the all matching data for a given type from OBGenericDataType
    //!    or an empty vector if nothing matches
    func getDataVector(_ type: MKGenericDataType) -> [MKGenericData]? {
        var tempvdata: [MKGenericData] = []
        if let vdata = self._vdata {
            for data in vdata {
                if data.type == type {
                    tempvdata.append(data)
                }
            }
        }
        return tempvdata
    }

    //! \return all data, suitable for iterating
    func getDataVector() -> [MKGenericData]? {
        return self._vdata
    }

    //! \return all data with a specific origin, suitable for iterating
    func getDataVector(_ origin: DataOrigin) -> [MKGenericData]? {
        var tempvdata: [MKGenericData] = []
        if let vdata = self._vdata {
            for data in vdata {
                if data.source == origin {
                    tempvdata.append(data)
                }
            }
        }
        return tempvdata
    }

    //   //! \return An iterator pointing to the beginning of the data
    //   OBDataIterator  BeginData()
    //     { return(_vdata.begin());        }
    //   //! \return An iterator pointing to the end of the data
    //   OBDataIterator  EndData()
    //     { return(_vdata.end());          }

    //! Adds a data object; does nothing if d==NULL
    func setData(_ d: MKGenericData) {
        if self._vdata == nil {
            self._vdata = [MKGenericData]()
        }
        self._vdata!.append(d)
    }
    
    //! Delete any data matching the given OBGenericDataType
    func deleteData(_ type: MKGenericDataType) {
        if let vdata = self._vdata {
            for (idx, data) in vdata.enumerated() {
                if data.type == type {
                    self._vdata?.remove(at: idx)
                }
            }
        }
    }

    //! Delete the given generic data from this object
    //   void                              DeleteData(OBGenericData*);
    func deleteData(_ d: MKGenericData) {
        if let vdata = self._vdata {
            for (idx, data) in vdata.enumerated() {
                if data == d {
                    self._vdata?.remove(at: idx)
                }
            }
        }
    }

    //! Delete all of the given generic data from this object
    func deleteData(_ v: [MKGenericData]) {
        if let vdata = self._vdata {
            for d in v {
                for (idx, data) in vdata.enumerated().reversed() {
                    if data == d {
                        self._vdata?.remove(at: idx)
                        break
                    }
                }
            }
        }
    }

    //! Deletes the generic data with the specified attribute, returning false if not found
    func deleteData(_ attr: String) -> Bool {
        if let vdata = self._vdata {
            for (idx, data) in vdata.enumerated() {
                if data.attr == attr {
                    self._vdata?.remove(at: idx)
                    return true
                }
            }
        }
        return false
    }   

    // func doTransformations(_ map: Dictionary<String,String>, _ conversion: OBConversion)

    static func classDescription() -> String {
        return ""
    }

    deinit {
        _vdata = nil
    }

}
