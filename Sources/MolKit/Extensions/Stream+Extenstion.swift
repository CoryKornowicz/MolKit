//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/28/23.
//

import Foundation


public struct PosixError : Error, CustomStringConvertible {
    public var code: Int32
    
    public init(code: Int32) {
        self.code = code
    }
    
    public var description: String {
        let str = String(cString: strerror(code), encoding: .utf8) ?? ""
        return "\(str)(\(code))"
    }
    
    public static var current: PosixError {
        return PosixError(code: Darwin.errno)
    }
}

public protocol InputFileHandlerProtocol {
    func read(maxSize: Int) throws -> Data
}

public protocol OutputFileHandlerProtocol {
    func write(data: Data) throws
}

public protocol FileHandlerProtocol {
    var isEOF: Bool { get }
    func tellg() throws -> Int
    func seekg(to position: Int) throws
}

public class FileHandler : FileHandlerProtocol {
    
    public enum Closer {
        case close
        case none
    }
    
    var filepath: URL
    var streamStatus: Stream.Status = .notOpen
    
    public convenience init(path: URL, mode: String) throws {
        guard let handle = Darwin.fopen(path.path, mode) else {
            throw PosixError.current
        }
        self.init(handle: handle, closer: .close, path)
        self.filepath = path
        streamStatus = .open
    }
    
    public init(handle: UnsafeMutablePointer<FILE>,
                closer: Closer, _ path: URL)
    {
        self.handle = handle
        self.closer = closer
        self.filepath = path
    }
    
    public let handle: UnsafeMutablePointer<FILE>
    public let closer: Closer
    
    deinit {
        switch closer {
        case .close:
            do {
                try close()
                self.streamStatus = .closed
            } catch {
                print("FileHandle.close() failed")
                self.streamStatus = .error
            }
        case .none:
            self.streamStatus = .notOpen
        }
    }
    
    public func close() throws {
        guard Darwin.fclose(handle) == 0 else {
            self.streamStatus = .error
            throw PosixError.current
        }
    }
    
    public var isEOF: Bool {
        return Darwin.feof(handle) != 0
    }
    
    public func tellg() throws -> Int {
        let position = Darwin.ftell(handle)
        if position == -1 {
            self.streamStatus = .error
            throw PosixError.current
        }
        return position
    }

    public func seekg(to position: Int) throws {
        let status = Darwin.fseek(handle, position, SEEK_SET)
        guard status == 0 else {
            self.streamStatus = .error
            throw PosixError.current
        }
    }
        
    static func == (_ lhs: FileHandler, _ rhs: FileHandle) -> Bool {
        return fileno(lhs.handle) == rhs.fileDescriptor
    }
    
    static func != (_ lhs: FileHandler, _ rhs: FileHandle) -> Bool {
        return fileno(lhs.handle) != rhs.fileDescriptor
    }
    
    
}

class InputFileHandler: FileHandler {}

extension InputFileHandler: InputFileHandlerProtocol {
    public func read(maxSize: Int) throws -> Data {
        self.streamStatus = .reading
        var data = Data(count: maxSize)
        let readSize = data.withUnsafeMutableBytes { (bytes) in
            Darwin.fread(bytes, 1, maxSize, handle)
        }
        if readSize < maxSize, !isEOF {
            throw PosixError(code: Darwin.ferror(handle))
        }
        data.count = readSize
        self.streamStatus = .open
        return data
    }

    public func readIntoString() throws -> String {
        self.streamStatus = .reading
        // Max string length is 2^31-1
        var data = Data(count: Int.max)
        let readSize = data.withUnsafeMutableBytes { (bytes) in
            Darwin.fread(bytes, 1, Int.max, handle)
        }
        data.count = readSize
        self.streamStatus = .open
        guard let stringChar = String.init(data: data, encoding: String.Encoding.utf8) else { return "\0" }
        return stringChar
    }
    
    public func peek() -> Int {
//        returns a copy of the next byte
        self.streamStatus = .reading
        var data = Data(count: 1)
        let readSize = data.withUnsafeMutableBytes { (bytes) in
            Darwin.fread(bytes, 1, 1, handle)
        }
        Darwin.fseek(handle, -1, SEEK_CUR)
        data.count = readSize
        self.streamStatus = .open
        
        guard let stringInt = String.init(data: data, encoding: String.Encoding.utf8) else { return 0 }

        return Int.init(stringInt) ?? 0
    }
    
    public func peek() -> String {
//        returns a copy of the next byte
        self.streamStatus = .reading
        var data = Data(count: 1)
        let readSize = data.withUnsafeMutableBytes { (bytes) in
            Darwin.fread(bytes, 1, 1, handle)
        }
        Darwin.fseek(handle, -1, SEEK_CUR)
        data.count = readSize
        self.streamStatus = .open
        
        guard let stringChar = String.init(data: data, encoding: String.Encoding.utf8) else { return "\0" }

        return stringChar
    }


}


class OutputFileHandler: FileHandler{}

extension OutputFileHandler: OutputFileHandlerProtocol {
    public func write(data: Data) throws {
        self.streamStatus = .writing
        let writtenSize = data.withUnsafeBytes { (bytes) in
            Darwin.fwrite(bytes, 1, data.count, handle)
        }
        if writtenSize != data.count {
            self.streamStatus = .error
            throw PosixError(code: Darwin.ferror(handle))
        }
        self.streamStatus = .atEnd
    }

}


class StringStream: FileHandlerProtocol {
    var _wrappedBuffer: String = ""
    
    var isEOF: Bool {
        return false
    }
    
    func tellg() throws -> Int {
        return self._wrappedBuffer.length
    }
    
    func seekg(to position: Int) throws {
        return
    }
    
    init() {}
    
    init(_wrappedBuffer: String) {
        self._wrappedBuffer = _wrappedBuffer
    }
}

class InputStringStream: StringStream, InputFileHandlerProtocol {
    func read(maxSize: Int) throws -> Data {
        if maxSize > self._wrappedBuffer.length {
            return self._wrappedBuffer.data(using: .utf8)!
        } else {
            return self._wrappedBuffer.substring(toIndex: maxSize).data(using: .utf8)!
        }
    }
    
    func read(maxLength: Int) -> String {
        if maxLength < self._wrappedBuffer.length {
            return self._wrappedBuffer.substring(toIndex: maxLength)
        } else {
            return self._wrappedBuffer
        }
    }
    
}

class OutputStringStream: StringStream, OutputFileHandlerProtocol {

    var string: String {
        return self._wrappedBuffer
    }

    func write(data: Data) throws {
        // convert data to string
        if let dataStr: String = String.init(data: data, encoding: .utf8) {
            self._wrappedBuffer += dataStr
        }
    }
    
    func write(data: String) {
        self._wrappedBuffer += data
    }
    
}




//class QFile {
//
//    init(fileURL: URL) {
//        self.fileURL = fileURL
//    }
//
//    deinit {
//        // You must close before releasing the last reference.
//        precondition(self.file == nil)
//    }
//
//    let fileURL: URL
//
//    private var file: UnsafeMutablePointer<FILE>? = nil
//
//    func open() throws {
//        guard let f = fopen(fileURL.path, "r") else {
//            throw NSError(domain: NSPOSIXErrorDomain, code: Int(errno), userInfo: nil)
//        }
//        self.file = f
//    }
//
//    func close() {
//        if let f = self.file {
//            self.file = nil
//            let success = fclose(f) == 0
//            assert(success)
//        }
//    }
//
//    func readLine(maxLength: Int = 1024) throws -> String? {
//        guard let f = self.file else {
//            throw NSError(domain: NSPOSIXErrorDomain, code: Int(EBADF), userInfo: nil)
//        }
//        var buffer = [CChar](repeating: 0, count: maxLength)
//        guard fgets(&buffer, Int32(maxLength), f) != nil else {
//            if feof(f) != 0 {
//                return nil
//            } else {
//                throw NSError(domain: NSPOSIXErrorDomain, code: Int(errno), userInfo: nil)
//            }
//        }
//        return String(cString: buffer)
//    }
//}



extension Character {
    // mutating func to edit character by offsetting its unicode scalar value
    func offset(by offset: Int) -> Character {
        guard let unicode = unicodeScalars.first else { return self }
        let value = unicode.value
        let offsetValue = value + UInt32(offset)
        guard let scalar = UnicodeScalar(offsetValue) else { return self }
        return Character(scalar)
    }
}
