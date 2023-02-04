//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/4/23.
//

import Foundation

protocol MKGlobalDataBase {
    
    var _filename: String { get }
    var _subdir: String { get }
        
    func getSize() -> Int
        
    func parseLine(_ line: String)
}


