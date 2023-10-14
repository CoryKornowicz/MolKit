//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/4/23.
//

import Foundation

extension String {
    
    func splitFileType() -> (String, String) {
        
        let components = self.components(separatedBy: ".")
        if components.count == 2 {
            return (components[0], components[1])
        } else {
            return (String(components[0..<components.endIndex-1].joined()), components[-1])
        }
    }
    
    func removeWhiteSpaceAndUnderscore() -> String {
        var outString: String = ""
        for char in self {
            if char != " " && char != "_" {
                outString += String(char)
            }
        }
        return outString
    }

    var length: Int {
        return count
    }

    subscript (i: Int) -> String {
        return self[i ..< i + 1]
    }

    func substring(fromIndex: Int) -> String {
        return self[min(fromIndex, length) ..< length]
    }
    
    func substring(fromIndex: Self.Index) -> String {
        return self[min(fromIndex.utf16Offset(in: self), length) ..< length]
    }

    func substring(toIndex: Int) -> String {
        return self[0 ..< max(0, toIndex)]
    }
    
    func substring(toIndex: Self.Index) -> String {
        return self[0 ..< max(0, toIndex.utf16Offset(in: self))]
    }

    subscript (r: Range<Int>) -> String {
        let range = Range(uncheckedBounds: (lower: max(0, min(length, r.lowerBound)),
                                            upper: min(length, max(0, r.upperBound))))
        let start = index(startIndex, offsetBy: range.lowerBound)
        let end = index(start, offsetBy: range.upperBound - range.lowerBound)
        return String(self[start ..< end])
    }
    
    func toScalar() -> (any Scalar)? {
        return NumberFormatter().number(from: self)?.doubleValue as? (any Scalar)
    }
    
    func toDouble() -> Double? {
        return NumberFormatter().number(from: self)?.doubleValue
    }
    
    func toInt() -> Int? {
        return NumberFormatter().number(from: self)?.intValue
    }
    
    var isNumber: Bool {
        let digitsCharacters = CharacterSet(charactersIn: "0123456789")
        return CharacterSet(charactersIn: self).isSubset(of: digitsCharacters)
    }

    func index<S: StringProtocol>(of string: S, options: String.CompareOptions = []) -> Index? {
        range(of: string, options: options)?.lowerBound
    }
    func endIndex<S: StringProtocol>(of string: S, options: String.CompareOptions = []) -> Index? {
        range(of: string, options: options)?.upperBound
    }
    func indices<S: StringProtocol>(of string: S, options: String.CompareOptions = []) -> [Index] {
        ranges(of: string, options: options).map(\.lowerBound)
    }
    func ranges<S: StringProtocol>(of string: S, options: String.CompareOptions = []) -> [Range<Index>] {
        var result: [Range<Index>] = []
        var startIndex = self.startIndex
        while startIndex < endIndex,
            let range = self[startIndex...]
                .range(of: string, options: options) {
                result.append(range)
                startIndex = range.lowerBound < range.upperBound ? range.upperBound :
                    index(range.lowerBound, offsetBy: 1, limitedBy: endIndex) ?? endIndex
        }
        return result
    }
    
}


extension URL {
    func foreachRow(_ rowParcer:((String, Int) -> Void) )
    {
        //Here we should use path not the absoluteString (wich contains file://)
        let path = self.path
        let m = "r"
        guard let cfilePath = (path as NSString).utf8String else {return}
        
        //Open file with specific mode (just use "r")
        guard let file = fopen(cfilePath, m)
        else {
            print("fopen can't open file: \"\(path)\", mode: \"\(m)\"")
            return
        }
        
        //Row capacity for getline()
        var cap = 0
        
        var row_index = 0
        
        //Row container for getline()
        var cline:UnsafeMutablePointer<CChar>? = nil
        
        //Free memory and close file at the end
        defer{free(cline); fclose(file)}
                    
        while getline(&cline, &cap, file) > 0
        {
            if let crow = cline,
               // the output line may contain '\n' that's why we filtered it
//               let s = String(utf8String: crow)?.filter({($0.asciiValue ?? 0) >= 32})
                let s = String(utf8String: crow)?.trimmingCharacters(in: .newlines)
            {
                rowParcer(s, row_index)
            }
            
            row_index += 1
        }
    }
    
    func foreachRow(_ rowParcer: ((String, Int) throws -> Void)) rethrows {
        //Here we should use path not the absoluteString (wich contains file://)
        let path = self.path
        let m = "r"
        guard let cfilePath = (path as NSString).utf8String else {return}
        
        //Open file with specific mode (just use "r")
        guard let file = fopen(cfilePath, m)
        else {
            print("fopen can't open file: \"\(path)\", mode: \"\(m)\"")
            return
        }
        
        //Row capacity for getline()
        var cap = 0
        
        var row_index = 0
        
        //Row container for getline()
        var cline:UnsafeMutablePointer<CChar>? = nil
        
        //Free memory and close file at the end
        defer{free(cline); fclose(file)}
                    
        while getline(&cline, &cap, file) > 0
        {
            if let crow = cline,
               // the output line may contain '\n' that's why we filtered it
//               let s = String(utf8String: crow)?.filter({($0.asciiValue ?? 0) >= 32})
                let s = String(utf8String: crow)?.trimmingCharacters(in: .newlines)
            {
                try rowParcer(s, row_index)
            }
            
            row_index += 1
        }
    }
    
    func foreachRow(offset: off_t?, _ rowParcer:((String, Int) -> Void) )
    {
        //Here we should use path not the absoluteString (wich contains file://)
        let path = self.path
        let m = "r"
        guard let cfilePath = (path as NSString).utf8String else {return}
        
        //Open file with specific mode (just use "r")
        guard let file = fopen(cfilePath, m)
        else {
            print("fopen can't open file: \"\(path)\", mode: \"\(m)\"")
            return
        }
        
        //Row capacity for getline()
        var cap = 0
        
        var row_index = 0
        
        //Row container for getline()
        var cline:UnsafeMutablePointer<CChar>? = nil
        
        //Free memory and close file at the end
        defer{free(cline); fclose(file)}
        
        if offset != nil {
            fseeko(file, offset!, SEEK_SET)
        }
                    
        while getline(&cline, &cap, file) > 0
        {
            if let crow = cline,
               // the output line may contain '\n' that's why we filtered it
//               let s = String(utf8String: crow)?.filter({($0.asciiValue ?? 0) >= 32})
                let s = String(utf8String: crow)?.trimmingCharacters(in: .newlines)
            {
                rowParcer(s, row_index)
            }
            
            row_index += 1
        }
    }
}

