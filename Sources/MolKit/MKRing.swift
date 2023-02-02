

import Foundation 

class MKRing: MKBase {


    func size() -> Int {
        return 0
    }

    //! \return Whether @p i as an atom index is in this ring
    func isInRing(_ i: Int) -> Bool {
    //   return(_pathset.BitIsSet(i));
        return false
    }

    func pathSize() -> Int {
        return 0
    }
    
    func isMember(_ bond: MKBond) -> Bool { return false }

} 
