

import Foundation 

public class MKTetraNonPlanarStereo: MKStereoBase {
    
    override init(_ mol: MKMol) {
        super.init(mol)
    }
    
    
    /**
     * Convert a @p ConfigType struct from any View/Winding to the
     * desired representation.
     *
     * This is a template method which works on ConfigType structs from
     * OBTetraNonPlanar subclasses. The subclasses can decide what data
     * member are needed to store the stereogenic unit (i.e. 1 atom for
     * tetrahedral, 3 for allene like, ...) and still use this generic
     * method to handle the real stereochemistry.
     *
     * A ConfigType struct should at least have the following data members:
     * @code
     * class SomeNonPlanarStereo : public TetraNonPlanarStereo
     * {
     *   public:
     *     struct Config
     *     {
     *       // constructor(s) are recommended!
     *
     *       // subclass specific stereogenic unit
     *       ...
     *
     *       union {
     *         unsigned long from;
     *         unsigned long towards;
     *       };
     *       OBStereo::Refs refs;
     *       OBStereo::Winding winding;
     *       OBStereo::View view;
     *     };
     * };
     * @endcode
     *
     * @param cfg A ConfigType struct from a OBTetraNonPlanar subclass.
     * @param from_or_towards The desired from/towards reference id (see @p view)
     * @param winding The desired winding.
     * @param view The desired viewing direction.
     *
     * @return The ConfigType struct with desired from/towards, winding and view.
     */
    static func toConfig<T: ConfigNonPlanar, U: ConfigNonPlanar>(_ cfg: T, _ fromTowards: from_or_towrds, _ winding: MKStereo.Winding = .Clockwise, _ view: MKStereo.View = .ViewFrom) -> U {
//        Swift ensures cases with defined values are uniquely treated, thus the need for reduncant code.
        switch cfg.from_or_towrds {
        case .from(let ref):
            if ref == .NoRef {
                print("MKTetraNonPlanarStereo.toConfig : Invalid from in ConfigType struct.")
                return U()
            }
        case .towards(let ref):
            if ref == .NoRef {
                print("MKTetraNonPlanarStereo.toConfig : Invalid from in ConfigType struct.")
                return U()
            }
        }
        
        if cfg.refs.count != 3 {
            print("MKTetraNonPlanarStereo.toConfig : Invalid refs size.")
            return U()
        }
        
        // copy the internal refs
        var result: U = cfg as! U
        
        result.center = cfg.center
        result.refs = cfg.refs
        result.from_or_towrds = fromTowards
        result.winding = winding
        result.view = view
        
        // keep track of the permuations by using the oddness
        var odd: Bool = false
        
        // find id
        if (cfg.from_or_towrds != fromTowards) {
            // move id to front and remove it = 1 permutation
            for i in 0..<3 {
                if (cfg.refs[i] == fromTowards) {
                    result.refs[i] = cfg.from_or_towrds.refValue
                    break
                }
            }
            // 1 permutation perfromed --> odd = true
            odd = true
        }
        
        // clockwise <-> anti-clockwise : odd = true
        if (winding == cfg.winding) {
            odd = !odd
        }
        // ViewFrom <-> ViewTowards : odd = true
        if (view == cfg.view) {
            odd = !odd
        }
        // make sure we actually found id
        if (result.refs.count == 3) {
            if (odd) {
                MKStereo.permutate(&result.refs, 1, 2)
            }
            return result
        }
        
        print("MKTetraNonPlanarStereo.toConfig : Parameter id not found in internal refs.")
        return result
    }
    /**
     * Change the winding of the ConfigType struct while maintaining the stereochemistry.
     */
    static func changeWinding<T: ConfigNonPlanar>(_ cfg: inout T) {
        cfg.winding = (cfg.winding == MKStereo.Winding.Clockwise) ? MKStereo.Winding.AntiClockwise : MKStereo.Winding.Clockwise
        MKStereo.permutate(&cfg.refs, 1, 2)
    }
    /**
     * Change the view of the ConfigType struct while maintaining the stereochemistry.
     */
    
    static func changeView<T: ConfigNonPlanar>(_ cfg: inout T) {
        cfg.view = (cfg.view == MKStereo.View.ViewFrom) ? MKStereo.View.ViewTowards : MKStereo.View.ViewFrom
        MKStereo.permutate(&cfg.refs, 1, 2)
    }
    /**
     * Invert the stereochemistry of the ConfigType struct.
     */
    static func invert<T: ConfigNonPlanar>(_ cfg: inout T) {
        MKStereo.permutate(&cfg.refs, 1, 2)
    }
    
}

