

import Foundation
import Surge
import simd

class MKBuilder {

    static let sharedInstance = MKBuilder()

    private init() {}

    func getNewBondVector(_ atom: MKAtom, _ length: Double) -> Vector<Double> {
        return Vector<Double>([0, 0, 0])
    }

}
