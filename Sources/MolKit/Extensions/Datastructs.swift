//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 3/7/23.
//

import Foundation

// LIFO
struct Stack<T> {
    
    
    var items: [T] = []
    
    init() {}
    
    init(_ items: [T]) {
        self.items = items
    }
    
    mutating func push(_ item: T) {
        self.items.insert(item, at: self.items.endIndex)
    }
    
    @discardableResult
    mutating func pop() -> T? {
        if items.isEmpty { return nil }
        return self.items.removeLast()
    }
    
    func peek() -> T? {
        return self.items.last
    }
    
    func isEmpty() -> Bool {
        return self.items.isEmpty
    }
    
}
