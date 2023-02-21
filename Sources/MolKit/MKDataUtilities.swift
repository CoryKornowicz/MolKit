
import Foundation

public let HARTEE_TO_KCALPERMOL = 627.509469
public let HARTREE_TO_KJPERMOL = 2625.49962
public let KJPERMOL_TO_KCALPERMOL = 1.0/4.184
public let RYDBERG_TO_KCALPERMOL = 313.755026
public let ELECTRONVOLT_TO_KCALPERMOL = 23.060538

/*! \brief
 * Convenience function to extract thermochemistry from a molecule structure
 *
 * \param[in] mol          The molecule structure
 * \param[in] bVerbose     If true will print information 
 * \param[inout] Nsymm     If not zero and differing from the rotational symmetry
 *                         in the input molecule, corrections to the entropy and
 *                         free energy will be applied. If zero will hold the symmetry
 *                         number from the input molecule on return.
 * \param[out] temperature The temperature
 * \param[out] DeltaHf0    Enthalpy of formation at T = 0
 * \param[out] DeltaHfT    Enthalpy of formation at T
 * \param[out] DeltaGfT    Gibbs energy of formation at T
 * \param[out] DeltaSfT    Entropy of formation at T
 * \param[out] S0T         Standard entropy at T
 * \param[out] CVT         Heat capacity at T and constant Volume
 * \param[out] Scomponents Translational, Rotational and Vibrational components of S0
 * \return true if all values were found, false otherwise.
 */

func extract_thermochemistry(mol: MKMol,
				    bVerbose:Bool,
				    Nsymm: Int,
				    Nrotbonds: Int,
				    dbdt: Double,
				    temperature: Double,
				    DeltaHf0: Double,
				    DeltaHfT: Double,
				    DeltaGfT: Double,
				    DeltaSfT: Double,
				    S0T: Double,
				    CVT: Double,
				    CPT: Double,
				    Scomponents: [Double],
				    ZPVE: Double) -> Bool {
    return false
}