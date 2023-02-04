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
    
}
