

import Foundation 

class MKRing: MKBase {


    func size() -> UInt {
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

} 