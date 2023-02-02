

import Foundation 
import simd

class MKBuilder {

    static let sharedInstance = MKBuilder()

    private init() {}

    func getNewBondVector(_ atom: MKAtom, _ length: Double) -> SIMD3<Double> {
        return SIMD3<Double>(0, 0, 0)
    }

}