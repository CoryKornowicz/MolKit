

import Foundation 

enum MKAtomFlag: UInt {
    case Aromatic = 2    // 1 << 1
    case Ring = 16       // 1 << 4
    case Closure = 1024  // 1 << 10
}

enum MKAtomStereoFlag: UInt {
    case Wedge = 4           // 1 << 2
    case Hash = 8            // 1 << 3
    case WedgeOrHash =  2048 // 1 << 11
    case CisOrTrans = 4096   // 1 << 12
}

class MKBond: MKBase {


    // TODO: Need to write implementations

    func getBondOrder() -> UInt { return 0}

    func getBeginAtom() -> MKAtom? {
        return nil
    }

    func getEndAtom() -> MKAtom? {
        return nil
    }

    func getNbrAtom(_ atom: MKAtom) -> MKAtom {
        return MKAtom()
    }

    func isInRing() -> Bool {
        return false
    }

    func isEster() -> Bool {
        return false
    }

}
