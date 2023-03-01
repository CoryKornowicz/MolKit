

// Define Globals 

import Foundation
import Surge

public let ONE_OVER_SQRT3: Double = 0.57735026918962576451
public let SQRT_TWO_THIRDS: Double = 0.81649658092772603272

public class MolKit {

    static let _SpaceGroups = MKSpaceGroups()
    static let _AtomicHeatOfFormationTable = MKAtomicHeatOfFormationTable()
    static let _TypeTable = MKTypeTable()
    static let _RingTyper = MKRingTyper()
    static let _BondTyper = MKBondTyper()
    static let _AtomTyper = MKAtomTyper()
    static let _AromTyper = MKAromaticTyper()
    static let _PhModel = MKPhModel()
    
    static let _plugin_ids: [String] = []
    
    static let theSMIFormat: SMIFormat = SMIFormat()
    static let theCANSMIFormat: CANSMIFormat = CANSMIFormat()
    
    init() {
        
//        Enable Plugins??
//        MolKit._plugin_ids.append([MolKit.theSMIFormat,
//                                   MolKit.theCANSMIFormat])
    }
    
}
