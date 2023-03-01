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
        case .none: break
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

