//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/18/23.
//

import Foundation


/**
* Convert any reference to @p atomId in a stereo object to an OBStereo::ImplicitRef.
* This function is called from OBMol::DeleteHydrogens()
* (via OBMol::DeleteHydrogen()) to remove any explicit references to a
* hydrogen atom that has been deleted. However, the code is not specific
* to hydrogen atoms and could be used for other atoms.
*
* @param mol The molecule
* @param atomId The Id of the atom to be converted to an OBStereo::ImplicitRef
* @since version 2.3
*/
func stereoRefToImplicit(_ mol: MKMol, _ atomId: Ref) {

}
/**
* Convert any reference to an OBStereo::ImplicitRef attached to @p centerId
* in a stereo object to an explicit reference to @p newId.
* This function is called from OBMol::AddHydrogens() and
* OBMol::AddHydrogen() to convert any implicit references to a
* hydrogen atom that has just been added. However, the code is not specific
* to hydrogen atoms and could be used for other atoms.
*
* @param mol The molecule
* @param centerId The Id of the atom to which the new explicit atom is attached
* @param newId The Id of the atom which was previously an OBStereo::ImplicitRef
* @since version 2.4
*/
func implicitRefToStereo(_ mol: MKMol, _ centerId: Ref, _ newId: Ref) {

}
