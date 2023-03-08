//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/4/23.
//

import Foundation

public class MKGlobalDataBase {
    
    var _filename: String
    var _subdir: String
    
    init(fileName: String, subDir: String) {
        self._filename = fileName
        self._subdir = subDir
        //TODO: should set locale here globally
    }
    
    func getSize() -> Int {
        return 0
    }
        
//    Need to override implementation
    func readFile() { return }
    
}


