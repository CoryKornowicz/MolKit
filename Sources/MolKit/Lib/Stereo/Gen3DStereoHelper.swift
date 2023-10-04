

/**
* @brief Helper class for 3D coordinate generation.
*
* This class can be used to check stereochemistry when generating 3D
* coordinates. This class also keeps track of unspecified stereochemistry
* to ensure this information is not lost when calling StereoFrom3D().
*/

class MKGen3DStereoHelper {


    var m_inputSmiles: String = ""
    var m_unspecifiedTetrahedral: [Ref] = []
    var m_unspecifiedCisTrans: [Ref] = []

    init() {}
    /**
    * @brief Store stereochemical information for later comparison.
    */
    func setup(_ mol: MKMol) {
        m_unspecifiedTetrahedral = []
        m_unspecifiedCisTrans = []

        // Store canonical SMILES of original molecules 
        let conv = MKConversion()
        conv.setOutFormat("can")
        m_inputSmiles = conv.writeString(mol, true)
        
        // Keep track of unspecified stereochemistry
        let facade = MKStereoFacade(mol)

        let tetrahedral: [MKTetrahedralStereo] = facade.getAllTetrahedralStereo()
        for i in 0..<tetrahedral.count {
            let cfg = tetrahedral[i].getConfig()
            if !cfg.specified {
                m_unspecifiedTetrahedral.append(cfg.center)
            }
        }
        
        let cistrans = facade.getAllCisTransStereo()
        for i in 0..<cistrans.count {
            let cfg = cistrans[i].getConfig()
            let begin = mol.getAtomById(cfg.begin)
            let end = mol.getAtomById(cfg.end)
            if begin == nil || end == nil {
                continue
            }
            let bond = mol.getBond(begin!, end!)
            if bond == nil {
                continue
            }
            if !cfg.specified {
                m_unspecifiedCisTrans.append(bond!.getId())
            }
        }
    }

    /**
    * @brief Check the stereochemistry.
    *
    * This function will perceive stereochemistry from 3D and compare this
    * with the stereochemistry that was stored when Setup() was called.
    *
    * @return True if the stereochemistry is correct.
    */
    func check(_ mol: inout MKMol) -> Bool {
        // Perceive stereo from 3D coords
        stereoFrom3D(&mol, true) // true  = force
        // // Make sure to respect previously unspecifed stereochemistry
        let facade = MKStereoFacade(mol)
        for i in 0..<m_unspecifiedTetrahedral.count {
            let ts = facade.getTetrahedralStereo(m_unspecifiedTetrahedral[i].intValue)
            if ts == nil {
                continue
            }
            var cfg = ts!.getConfig()
            cfg.specified = false
            ts!.setConfig(cfg)
        }
        for i in 0..<m_unspecifiedCisTrans.count {
            let ct = facade.getCisTransStereo(m_unspecifiedCisTrans[i].intValue)
            if ct == nil {
                continue
            }
            let cfg = ct!.getConfig()
            cfg.specified = false
            ct!.setConfig(cfg)
        }
        // // Generate canonical SMILES with stereochemistry perceived from 3D coords.
        let conv: MKConversion = MKConversion()
        conv.setOutFormat("can")
        let predictedSmiles: String = conv.writeString(mol, true)
        return m_inputSmiles == predictedSmiles
    }

}
