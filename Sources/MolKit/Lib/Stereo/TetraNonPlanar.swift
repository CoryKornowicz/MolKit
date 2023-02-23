

import Foundation 

class MKTetraNonPlanarStereo: MKStereoBase {
    override init(_ mol: MKMol) {
        super.init(mol)
    }

    
    func toConfig<T: ConfigNonPlanar>(_ cfg: T) -> T {
        return cfg
    }


}

