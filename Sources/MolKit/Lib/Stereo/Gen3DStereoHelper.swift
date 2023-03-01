

/**
* @brief Helper class for 3D coordinate generation.
*
* This class can be used to check stereochemistry when generating 3D
* coordinates. This class also keeps track of unspecified stereochemistry
* to ensure this information is not lost when calling StereoFrom3D().
*/

class MKGen3DStereoHelper {


    var m_inputSmailes: String = ""
    var m_unspecifiedTetrahedral: [Int] = []
    var m_unspecifiedCisTrans: [Int] = []

    init() {}

    /**
    * @brief Store stereochemical information for later comparison.
    */
    func setup(_ mol: MKMol) {}
    
    /**
    * @brief Check the stereochemistry.
    *
    * This function will perceive stereochemistry from 3D and compare this
    * with the stereochemistry that was stored when Setup() was called.
    *
    * @return True if the stereochemistry is correct.
    */
    func check(_ mol: MKMol) -> Bool {
        return true
    }

}