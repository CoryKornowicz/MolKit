//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation

  /*\brief Assigns atom types, hybridization, and formal charges
    The MKAtomTyper class is designed to read in a list of atom typing
    rules and apply them to molecules. The code that performs atom
    typing is not usually used directly as atom typing, hybridization
    assignment, and charge are all done
    automatically when their corresponding values are requested of
    atoms:
    \code
    atom->GetType();
    atom->GetFormalCharge();
    atom->GetHyb();
    \endcode
  */

class MKAtomTyper {
    
    static let sharedInstance = MKAtomTyper()
    
    private init() {}
    
    func assignTypes(_ mol: MKMol) {
        
    }
    
    func assignHyb(_ mol: MKMol) {
        
    }
    
}
committment