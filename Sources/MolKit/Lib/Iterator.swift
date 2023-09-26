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

    var curr: T? {
        get {
            index < collection.count ? collection[index] : nil
        }
        set {
            guard newValue != nil else { return }
            collection[index] = newValue!
        }
    }
    
    public func next() -> T? {
        defer { index += 1 }
        return index < collection.count ? collection[index] : nil
    }
    
    public func append(_ newElement: T) {
        self.collection.append(T.self as! T)
    }
    
    static func += (_ lhs: MKIterator<T>, _ rhs: Int) {
        lhs.index += rhs
    }
    
    static func -= (_ lhs: MKIterator<T>, _ rhs: Int) {
        lhs.index -= rhs
    }
    
}


/** \class OBMolAtomDFSIter obiter.h <openbabel/obiter.h>

      \since version 2.1

      To facilitate iteration through all atoms in a molecule, without resorting
      to atom indexes (which <strong>will</strong> change in the future), a
      variety of iterator methods are provided.

      This class provides a depth-first search ordering of atoms. When one
      connected component is exhausted, the iterator will start at another until
      all atoms are visited. No guarantee is made as to the ordering of
      iteration through connected components.

      The iterator maintains an internal stack and list of visited
      atoms. As such it may not be appropriate for memory-constrained
      situations when iterating through large molecules.

      Use of this iterator has been made significantly easier by a series
      of macros in the obiter.h header file:

      \code
      \#define FOR_DFS_OF_MOL(a,m)     for( OBMolAtomDFSIter     a(m); a; ++a )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      FOR_DFS_OF_MOL(a, mol)
      {
         // The variable a behaves like OBAtom* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*a

      }
      \endcode
  **/

postfix operator ++

public class MKAtomDFSIterator: IteratorProtocol {
    
    var _parent: MKMol 
    var _ptr: MKAtom?
    var _notVisited: MKBitVec = MKBitVec()
    var _stack: Stack<MKAtom> = Stack<MKAtom>()

    init(_ mol: MKMol, _ startIndex: Int = 1) {
        self._parent = mol
        if let atom = mol.getAtom(startIndex) {
            self._ptr = atom
        } else {
            fatalError("Cannot find atom in mol")
        }
        
        self._notVisited.resize(UInt32(_parent.numAtoms()))
        self._notVisited.setRangeOn(0, UInt32(_parent.numAtoms() - 1))
        
        _notVisited.setBitOff(_ptr!.getIdx() - 1)
        
        if let nbrs = _ptr!.getNbrAtomIterator() {
            for a in nbrs {
                _stack.push(a)
                _notVisited.setBitOff(a.getIdx() - 1)
            }
        }
    }
    
    func current() -> MKAtom? {
        return self._ptr
    }

    @discardableResult
    public func next() -> MKAtom? {
        if self._stack.isEmpty() {
            return nil
        } else {
            self._ptr = self._stack.pop()
        }
        return self._ptr
    }
    
    func getIdx() -> Int? {
        return self._ptr?.getIdx()
    }
    
    func isEmpty() -> Bool {
        return self._ptr == nil
    }
    
    static postfix func ++ (_ iter: inout MKAtomDFSIterator) {
        
        if !iter._stack.isEmpty() {
            iter._ptr = iter._stack.pop()
        } else { // are there any disconnected subgraphs?
            let next = iter._notVisited.firstBit()
            if next != iter._notVisited.endBit() {
                iter._ptr = iter._parent.getAtom(next + 1)
                iter._notVisited.setBitOff(next)
            } else {
                iter._ptr = nil
            }
        }
        if iter._ptr != nil {
            if let nbrs = iter._ptr?.getNbrAtomIterator() {
                for a in nbrs {
                    if iter._notVisited[a.getIdx() - 1] {
                        iter._stack.push(a)
                        iter._notVisited.setBitOff(a.getIdx() - 1)
                    }
                }
            }
        }
    }
    
    
}

//public class MKBFSIterator<T>: IteratorProtocol, Sequence {
//
//    private var collection: [T]
//    private var index = 0
//
//    init(_ collection: [T]) {
//        self.collection = collection
//    }
//
//    public func next() -> T? {
//        defer { index += 1 }
//        return index < collection.count ? collection[index] : nil
//    }
//
//    public func append(_ newElement: T) {
//        self.collection.append(T.self as! T)
//    }
//
//}


