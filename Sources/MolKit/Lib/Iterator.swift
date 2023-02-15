//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/15/23.
//

import Foundation

/// Concrete Iterators implement various traversal algorithms. These classes
/// store the current traversal position at all times.
public class MKIterator<T>: IteratorProtocol, Sequence {
    
    private var collection: [T]
    private var index = 0

    init(_ collection: [T]) {
        self.collection = collection
    }

    public func next() -> T? {
        defer { index += 1 }
        return index < collection.count ? collection[index] : nil
    }
    
    public func append(_ newElement: T) {
        self.collection.append(T.self as! T)
    }
    
    
}
