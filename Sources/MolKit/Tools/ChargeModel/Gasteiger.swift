//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 3/19/23.
//

import Foundation

class GasteigerCharges : MKChargeModel {
    
    required init(_ id: String, _ isDefault: Bool) {
        super.init(id, false)
    }
    
    override class func dipoleScalingFactor() -> Double {
        return 3.4927
    }
    
    override func description() -> String? {
        return "Assign Gasteiger-Marsili sigma partial charges"
    }
    
    override func computeCharges(_ mol: MKMol) -> Bool {
        mol.setPartialChargesPerceived()
        let gc: MKGastChrg = MKGastChrg()
        let returnValue = gc.assignPartialCharges(mol)
        fillChargeVectors(mol)
        return returnValue
    }
}
