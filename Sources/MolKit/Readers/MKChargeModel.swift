//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 3/19/23.
//

import Foundation
import Surge 
/** \class OBChargeModel chargemodel.h <openbabel/chargemodel.h>
 \brief Atomic partial charge models
 \since version 2.3
 
 Classes derived from OBChargeModel implement different atomic partial
 charge models. It is intended to allow assinging partial charges
 beyond the traditional Gasteiger-Marsili sigma charges previously used
 in Open Babel. A --partialcharge method is provided for the obabel
 command-line, allowing you to override the Gasteiger charge assignment
 and use other charge models.
 
 The advantage of plugin classes is that no existing code has to be modified
 when a new class is added. You can list those that are present by
 obabel -L charges
 or from a menu item in the GUI.
 
 Any OBChargeModel derived class works like other plugins and needs to
 to have a constructor, a function returning a short description, and a
 ComputeCharges() function which does the work. A single global
 instance of the class needs to be instantiated to define the ID, by
 which the class is subsequently accessed.
 
 Once ComputeCharges() has been called, the atoms of the molecule can
 be queried for partial or formal charges using
 OBAtom::GetPartialCharge() or in vector form from the model itself:
 
 \code
 OBMol inputMolecule;
 OBChargeModel *mmffCharges = OBChargeModel::FindType("mmff94");
 const std::vector<double> partialCharges;
 if (mmffCharges && mmffCharges->ComputeCharges(inputMolecule)) {
 partialCharges = mmffCharges->GetPartialCharges();
 }
 \endcode
 
 Note: Formal charges are also returned as floating point values, since
 some charge models consider delocalized charges (e.g., 0.5 for an O in
 a carboxylate CO2- group).
 
 \code
 OBChargeModel *gasteiger = OBChargeModel::FindType("gasteiger");
 if (gasteiger) {
 cout << " gasteiger: " << dipoleMagnitude(gasteiger->GetDipoleMoment(mol));
 }
 \endcode
 
 By default, Open Babel 2.3 includes Gasteiger and MMFF94 partial
 charges. If the Eigen matrix library is found when compiling, the QEq
 and QTPIE methods will be added. Future releases will likely add
 additional charge models, including the EEM method.
 
 */

class MKChargeModel: MKPlugin, MKPluginProtocol {
    
    static var Default: MKChargeModel?
    
    static var map: PluginMapType<MKChargeModel> = PluginMapType<MKChargeModel>()
    
    var m_partialCharges: [Double] = []
    var m_formalCharges: [Double] = [] 
        
    static func findType(_ ID: String?) -> MKChargeModel? {
        if ID == nil {
            return MKChargeModel.Default
        }
        return MKPlugin.baseFindType(MKChargeModel.map, ID!) as? MKChargeModel
    }
    
    required init(_ id: String, _ isDefault: Bool) {
        super.init()
        self._id = _id
        if isDefault || MKChargeModel.map.isEmpty {
            MKChargeModel.Default = self
        }
        if MKChargeModel.map.map({ $0.0 == _id ? 1 : 0}).reduce(0, +) == 0 {
            MKChargeModel.map[_id] = self
            MKPlugin.pluginMap[typeID()] = self
        }
    }
    
    func getMap() -> PluginMapType<MKChargeModel> {
        MKChargeModel.map
    }
    
    override func typeID() -> String {
        return "charges"
    }

    func fillChargeVectors(_ mol: MKMol) {
        m_partialCharges = []
        m_partialCharges.reserveCapacity(mol.numAtoms())
        m_formalCharges = []
        m_formalCharges.reserveCapacity(mol.numAtoms())
        
        for atom in mol.getAtomIterator() {
            m_partialCharges.append(atom.getPartialCharge())
            m_formalCharges.append(Double(atom.getFormalCharge()))
        }
    }

    /// Provide a scaling factor for the dipole moment -- ideally calibrated from many molecules
    class func dipoleScalingFactor() -> Double {
        return 1.0
    }

    /// \return whether partial charges were successfully assigned to this molecule
    /// \note The method should fill m_partialCharges and m_formalCharges as well
    func computeCharges(_ mol: MKMol) -> Bool {
        return false
    }

    func computeCharges(_ mol: MKMol, _ args: [String]) -> Bool {
        return self.computeCharges(mol)
    }

    /// \return a vector of the formal charges on each atom, indexed from 0
    /// This method returns floating point formal charges since some
    /// charge models consider fractional charges (e.g., 0.5 for an
    /// oxygen in a carboxylate CO2- group).
    /// \note If ComputeCharges() has not been called, this will return an empty vector
    func getFormalCharges() -> [Double] {
        return m_formalCharges
    }

    /// \return a vector of the partial charges on each atom, indexed from 0
    /// \note If ComputeCharges() has not been called, this will return an empty vector
    func getPartialCharges() -> [Double] {
        return m_partialCharges
    }

    /// \return a vector of the dipole moment from this molecule
    func getDipoleMoment(_ mol: MKMol) -> Vector<Double> {
        var dipoleMoment = Vector<Double>.init(dimensions: 3, repeatedValue: 0.0)
        
        if self.computeCharges(mol) {
            for atom in mol.getAtomIterator() {
                dipoleMoment += atom.getVector() * atom.getPartialCharge()
            }
        }
        
        return dipoleMoment * MKChargeModel.dipoleScalingFactor()
    }
    

}

