//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 9/25/23.
//

import Foundation

enum MKLoggingLevel {
    case mkError
    case mkWarning
    case mkInfo
    case mkAuditMsg
    case mkDebug
}

class MKLogger {
    
    
    static func throwError(_ method: String = #function, errorMsg: String) {
        let msg = MKLogger.decorateWithFunctionName(method) + errorMsg
        //write to stderror
        msg.data(using: .utf8).map(FileHandle.standardError.write)
    }
    
}

extension MKLogger {
    
    private static func decorateWithFunctionName(_ funcName: String) -> String {
        return "### originating function: \(funcName) ###"
    }
    
}
