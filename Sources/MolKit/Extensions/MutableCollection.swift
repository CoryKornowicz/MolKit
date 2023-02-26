//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/7/23.
//

import Foundation
import Collections


extension MutableCollection {
  mutating func updateEach(_ update: (inout Element) -> Void) {
    for i in indices {
      update(&self[i])
    }
  }
}

enum SubscriptError: Error {
    case outOfBounds
    case greaterThanZero
    case lessThanLastIndex
}

extension Collection where Indices.Iterator.Element == Index {
    
    public subscript(safelyAccess ind: Index) -> Iterator.Element {
        get { return (ind as! Int) < 0 ? self[index(endIndex, offsetBy: (ind as! Int) - 1)] : self[ind]}
        
        
    }

}

extension Sequence where Iterator.Element: Hashable {
    func unique() -> [Iterator.Element] {
        var seen: Set<Iterator.Element> = []
        return filter { seen.insert($0).inserted }
    }
}


// MARK: Simple Iterator class
public class Iterator<T>: IteratorProtocol, Sequence where T: Equatable {
    
    private var collection: [T]
    private var index = 0

    var first: T? {
        self.collection.first
    }
    
    var last: T? {
        self.collection.last
    }
    
    
    init(_ collection: [T]) {
        self.collection = collection
    }

    public func tellg() -> Int {
        return self.index
    }

    public func seekg(_ pos: Int) {
        self.index = pos
    }

    public func next() -> T? {
        defer { index += 1 }
        return index < collection.count ? collection[index] : nil
    }
    
    public func ignore() {
        index += 1
    }
    
    public func ignore(by: Int) {
        index += by
    }

    public func ignore(until: Character) {
        while self.peek() != until && !self.isEmpty(){
            self.ignore()
        }
    }
    
    public func peek(by: Int = 0) -> T {
        if self.index + by < self.collection.count {
            return self.collection[self.index + by]
        } else {
            let newIdx = (self.index + by) % self.collection.count
            self.index = newIdx
            return self.collection[newIdx]
        }
    }

    public func unget() {
        if self.index > 0 {
            self.index -= 1
        }
    }

    public func parseDouble() -> Double? {
        var numBuffer: String = ""
        var char = self.peek()
        while char.isNumber || char == "." {
            numBuffer.append(char)
            self.ignore()
            char = self.peek()
        }
        return numBuffer.toDouble()
    }
    
    public func nextElement(_ elem: T, updateIndex: Bool = false) -> T? {
        if let currElemIndex = self.collection.firstIndex(where: {$0 == elem}) {
//            peek at the next index
            if currElemIndex < self.collection.endIndex - 1 {
                if updateIndex {
                    self.index = currElemIndex + 1
                    return self.next()
                } else {
                    return self.collection[currElemIndex+1]
                }
            } else {
                return nil
            }
        }
        return nil
    }
    
    public func nextUntil(_ elem: T) -> ArraySlice<T> {
        let elemIndx = self.collection.firstIndex(where: {$0 == elem})
        // return substring of index to elemIndx
        if let elemIndx = elemIndx {
            let subArr = self.collection[self.index..<elemIndx]
            self.index = elemIndx
            return subArr.last
        }
    }  
    
    public func append(_ newElement: T) {
        self.collection.append(T.self as! T)
    }
    
    subscript (_ idx: Int) -> T? {
        if idx < self.collection.count { return self.collection[idx] }
        else { return self.collection.last }
    }

    func isEmpty() -> Bool {
        return self.collection.isEmpty || self.index >= self.collection.count
    }

    func setEmpty() {
        self.index = self.collection.count
    }
    
}
