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

    func substring(toIndex: Int) -> String {
        return self[0 ..< max(0, toIndex)]
    }

    subscript (r: Range<Int>) -> String {
        let range = Range(uncheckedBounds: (lower: max(0, min(length, r.lowerBound)),
                                            upper: min(length, max(0, r.upperBound))))
        let start = index(startIndex, offsetBy: range.lowerBound)
        let end = index(start, offsetBy: range.upperBound - range.lowerBound)
        return String(self[start ..< end])
    }
    
    typealias Scalar = BinaryFloatingPoint & FloatingPoint
    
    func toScalar() -> (any Scalar)? {
        return NumberFormatter().number(from: self)?.doubleValue as? (any Scalar)
    }

}

extension URL

{
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
               let s = String(utf8String: crow)?.filter({($0.asciiValue ?? 0) >= 32})
            {
                rowParcer(s, row_index)
            }
            
            row_index += 1
        }
    }
}
