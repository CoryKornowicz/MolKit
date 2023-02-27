//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/26/23.
//

import Foundation

//!\brief Used to store chemical reactions (i.e., reactants -> products)
//!
//! Reactants and products stored as smart pointers to molecules stored elsewhere.
//!
//! For performing actual reaction transformations (i.e., deleting atoms,
//! changing bonds, etc.) use the OBChemTsfm class.

class MKReaction: MKBase {
    
    var _reactants: [MKMol] = []
    var _products: [MKMol] = []
    var _agents: [MKMol] = []
    var _ts: MKMol? = nil
    var _title: String = ""
    var _comment: String = ""
    var _reversible: Bool = false
    
    func numReactants() -> Int {
        return _reactants.count
    }

    func numProducts() -> Int {
        return _products.count
    }

    func numAgents() -> Int {
        return _agents.count
    }

    func addReactant(_ sp: MKMol) {
        _reactants.append(sp)
    }

    func addProduct(_ sp: MKMol) {
        _products.append(sp)
    }

    func setTransitionState(_ sp: MKMol) {
        _ts = sp
    }

    func addAgent(_ sp: MKMol) {
        _agents.append(sp)
    }

    func getReactant(_ i: Int) -> MKMol? {
        if i < 0 || i >= _reactants.count {
            return nil
        }
        return _reactants[i]
    }

    func getProduct(_ i: Int) -> MKMol? {
        if i < 0 || i >= _products.count {
            return nil
        }
        return _products[i]
    }

    func getAgent(_ i: Int) -> MKMol? {
        if i < 0 || i >= _agents.count {
            return nil
        }
        return _agents[i]
    }

    func getTransitionState() -> MKMol? {
        return _ts
    }
    
    override func getTitle() -> String {
        return _title
    }
    
    func getComment() -> String {
        return _comment
    }
    
    func setComment(_ comment: String) {
        _comment = comment
    }
    
    override func setTitle(_ title: String) {
        _title = title
    }
    
    func isReversible() -> Bool {
        return _reversible
    }
    
    func setReversible(_ b: Bool = true) {
        _reversible = b
    }
    
    override func classDescription() -> String {
        return "reactions"
    }
    
    func clear() -> Bool {
        _reactants.removeAll()
        _products.removeAll()
        _agents.removeAll()
        _ts = nil
        _title = ""
        _comment = ""
        _reversible = false
        return true
    }
    
}
