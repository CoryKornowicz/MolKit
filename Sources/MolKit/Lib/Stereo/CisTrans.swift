

import Foundation 

class MKCisTransStereo: MKTetraPlanarStereo {
    
    
    class Config: ConfigPlanar {
        var begin: Ref
        var end: Ref
        var shape: MKStereo.Shape
        var refs: Refs
        var specified: Bool
        
        init(begin: Ref, end: Ref, shape: MKStereo.Shape, refs: Refs, specified: Bool = true) {
            self.begin = begin
            self.end = end
            self.shape = shape
            self.refs = refs
            self.specified = true
        }
        
        init() {
            self.begin = .NoRef
            self.end = .NoRef
            self.refs = []
            self.shape = .ShapeU
            self.specified = true
        }
    }
    
    
    override func getType() -> MKStereo.TType {
        return .CisTrans
    }
    
    
}
