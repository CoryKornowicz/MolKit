//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 3/7/23.
//

import Foundation

protocol Copying {
    init(instance: Self)
}

extension Copying {
    func mutcopy() -> Self {
        return Self.init(instance: self)
    }
}
