//
//  MKGastChrg.swift
//  
//
//  Created by Cory Kornowicz on 2/1/23.
//

import Foundation

public let OB_GASTEIGER_DENOM: Double = 20.02
public let OB_GASTEIGER_DAMP: Double =  0.5
public let OB_GASTEIGER_ITERS: Int = 6

struct MKGasteigerState {
    var a: Double, b: Double, c: Double, denom: Double, chi: Double, q: Double = 0.0
    
    init() {
        self.a = 0.0
        self.b = 0.0
        self.c = 0.0
        self.denom = 0.0
        self.chi = 0.0
        self.q = 0.0
    }

    init(a: Double, b: Double, c: Double, denom: Double, chi: Double, q: Double) {
        self.a = a
        self.b = b
        self.c = c
        self.denom = denom
        self.chi = chi
        self.q = q
    }

    mutating func setValues(_ a: Double, _ b: Double, _ c: Double, _ q: Double) {
        self.a = a
        self.b = b
        self.c = c
        self.denom = a+b+c
        self.q = q
    }
}


/*! \class OBGastChrg molchrg.h <openbabel/molchrg.h>
    \brief Assigns Gasteiger partial charges

    The OBGastChrg class is responsible for the assignment of partial
    charges to a molecule according to the Gasteiger charge model (sigma). When
    the partial charge of an atom is requested and partial charges do not
    yet exist for the molecule the OBGastChrg class will automatically be
    called to assign charges for the molecule. If charges have been read
    in for a molecule the following code shows how to force the
    recalculation of partial charges:
    \code
    OBMol mol;

    mol.UnsetPartialChargesPerceived();
    FOR_ATOMS_IN_MOL(atom, mol)
    {
       cout << "atom number = " << atom->GetIdx();
       cout << " charge = " << atom->GetPartialCharge() << endl;
    }
    \endcode
    Formal charges are used as seed values of the initial charge of atoms,
    and the partial charge is propagated to neighboring atoms. For
    example, quaternary amines would have a +1 charge, and the effect of
    the positive charge would be felt by neighboring atoms according to
    the Gasteiger model (sigma).

    For more information, see:
    J. Gasteiger & M. Marsili, "A New Model for Calculating Atomic Charges in Molecules" Tetrahedron Lett., (1978) 3181-3184.
  */

class MKGastChrg {
    
    var _gsv: Array<MKGasteigerState> = []
    
    init() {}

    //! Set initial partial charges in @p mol
    //! Carbonyl O => -0.5
    //! Phosphate O => -0.666
    //! Sulfate O => -0.5
    //! All other atoms are set to have their initial charge from their formal charge
    private func initialPartialCharges(_ mol: MKMol) {
        for atom in mol.getAtomIterator() {
            if atom.isCarboxylOxygen() {
                atom.setPartialCharge(-0.500)
            } else if atom.isPhosphateOxygen() && atom.getHeavyDegree() == 1 {
                atom.setPartialCharge(-0.666)
            } else if atom.isSulfateOxygen() {
                atom.setPartialCharge(-0.500)
            } else {
                atom.setPartialCharge(Double(atom.getFormalCharge()))
            }
        }
    }
    
    private func gasteigerSigmaChi(_ atom: MKAtom, _ a: inout Double, _ b: inout Double, _ c: inout Double) -> Bool {
        
        var count: Int = 0
        var val: Array<Double> = [0.0, 0.0, 0.0]
        switch atom.getAtomicNum() {
        case 1:
            val[0] = 0.37
            val[1] = 7.17
            val[2] = 12.85
        case 6:
            if atom.getHyb() == 3 {
                val[0] = 0.68
                val[1] = 7.98
                val[2] = 19.04
            }
            if atom.getHyb() == 2 {
                val[0] = 0.98
                val[1] = 8.79
                val[2] = 19.62
            }
            if atom.getHyb() == 1 {
                val[0] = 1.67
                val[1] = 10.39
                val[2] = 20.57
            }
        case 7:
            if atom.getHyb() == 3 {
                if atom.getExplicitDegree() == 4 || atom.getFormalCharge() != 0 {
                    val[0] = 0.0
                    val[1] = 0.0
                    val[2] = 23.72
                } else {
                    val[0] = 2.08
                    val[1] = 11.54
                    val[2] = 23.72
                }
            }
            
            if atom.getHyb() == 2 {
                if atom.getType() == "Npl" || atom.getType() == "Nam" {
                    val[0] = 2.46
                    val[1] = 12.32
                    val[2] = 24.86
                } else {
                    val[0] = 2.57
                    val[1] = 12.87
                    val[2] = 24.87
                }
            }
            
            if atom.getHyb() == 1 {
                val[0] = 3.71
                val[1] = 15.68
                val[2] = 27.11
            }
        case 8:
            if atom.getHyb() == 3 {
                val[0] = 2.65
                val[1] = 14.18
                val[2] = 28.49
            }
            if atom.getHyb() == 2 {
                val[0] = 3.75
                val[1] = 17.07
                val[2] = 31.33
            }
        case 9:
            val[0] = 3.12
            val[1] = 14.66
            val[2] = 30.82
        case 15: //P
            val[0] = 1.62
            val[1] = 8.90
            val[2] = 18.10
        case 16: //S
            count = Int(atom.countFreeOxygens())
            if count == 0 || count == 1 {
                val[0] = 2.39
                val[1] = 10.14
                val[2] = 20.65
            }
            if count > 1 {
                val[0] = 2.39
                val[1] = 12.00
                val[2] = 24.00
            }
        case 17: //Cl
            val[0] = 2.66
            val[1] = 11.00
            val[2] = 22.04
        case 35: //Br
            val[0] = 2.77
            val[1] = 10.08
            val[2] = 19.71
        case 53: //I
            val[0] = 2.90
            val[1] = 9.90
            val[2] = 18.82
        case 13: //Al
            val[0] = 1.06
            val[1] = 5.47
            val[2] = 11.65
        default:
            return false
        }
        if !isNearZero(val[2]) {
            a = val[1]
            b = (val[2] - val[0]) / 2
            c = (val[2] + val[0]) / 2 - val[1]
        } else {
            return false
        }
        return true
    }
    
    func gsvResize(_ int: Int) {
        self._gsv = Array<MKGasteigerState>(repeating: MKGasteigerState(), count: int)
    }
    
    @discardableResult
    func assignPartialCharges(_ mol: MKMol) -> Bool {
        // annotate that partial charges come from Gasteiger
        let dp = MKPairData<String>()
        dp.setAttribute("PartialCharges")
        dp.setValue("Gasteiger")
        dp.setOrigin(.perceived)
        mol.setData(dp)
        
        gsvResize(mol.numAtoms() + 1)
        var a: Double = 0.0
        var b: Double = 0.0
        var c: Double = 0.0
        
        for atom in mol.getAtomIterator() {
            if !gasteigerSigmaChi(atom, &a, &b, &c) {
                return false
            }
            _gsv[atom.getIdx()].setValues(a, b, c, atom.getPartialCharge())
        }
        
        var alpha: Double = 1.0
        var charge: Double
        var denom: Double
        
        for _ in 0..<OB_GASTEIGER_ITERS {
            alpha *= OB_GASTEIGER_DAMP
            
            for j in 1..<_gsv.count {
                charge = _gsv[j].q
                _gsv[j].chi = (_gsv[j].c*charge+_gsv[j].b)*charge+_gsv[j].a
            }
            
            for bond in mol.getBondIterator() {
                let src = bond.getBeginAtom()
                let dst = bond.getEndAtom()
                
                if _gsv[src.getIdx()].chi >= _gsv[dst.getIdx()].chi {
                    if dst.getAtomicNum() == 1 {
                        denom = Double(OB_GASTEIGER_DENOM)
                    } else {
                        denom = _gsv[dst.getIdx()].denom
                    }
                } else {
                    if src.getAtomicNum() == 1 {
                        denom = Double(OB_GASTEIGER_DENOM)
                    } else {
                        denom = _gsv[src.getIdx()].denom
                    }
                }
                
                charge = (_gsv[src.getIdx()].chi - _gsv[dst.getIdx()].chi)/denom
                _gsv[src.getIdx()].q -= alpha*charge
                _gsv[dst.getIdx()].q += alpha*charge
            }
        }
        
        for atom in mol.getAtomIterator() {
            atom.setPartialCharge(_gsv[atom.getIdx()].q)
        }
        
        return true
    }
    
    deinit {
        self._gsv.removeAll()
    }
    
}
