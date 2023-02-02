

import Foundation 

class MKStereoFacade {

    let mol: MKMol

    init(_ mol: MKMol) {
        self.mol = mol
    }

    func hasTetrahedralStereo(_ id: Int) -> Bool {
        return false
    }

}