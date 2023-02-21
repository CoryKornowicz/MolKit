//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/20/23.
//

import Foundation


/** \class MKAtomHOF data.h <openbabel/data.h>
      \brief helper class for OBAtomicHeatOfFormationTable

      Stores both theoretical and experimental corrections
      needed to compute the Enthalpy of formation. In order to
      use these you need to perform
      Gaussian G2/G3/G4 or CBS-QB3 calculations.
  **/
struct MKAtomHOF {
    var _element: String
    var _charge: Int
    var _method: String
    var _desc: String
    var _T: Double
    var _value: Double
    var _multiplicity: Int
    var _unit: String
}

/** \class OBAtomicHeatOfFormationTable data.h <openbabel/data.h>
      \brief Atomic Heat of Formation Table

      Contributions of atoms to Enthalpy of Formation calculations performed
      in Gaussian, using the G2/G3/G4 or CBS-QB3 methods respectively.
      The energies produced by Gaussian have to be corrected according to their
      document on Thermochemistry with Gaussian. The data in the file
      BABEL_DATA/atomization_energies.txt supplies this information based on
      single atom calculations with Gaussian and the appropriate method and
      experimental data from Curtiss et al., J. Chem. Phys. 106 (1997) 1063-1079.
  */

class MKAtomicHeatOfFormationTable: MKGlobalDataBase {
    
    var _atomhof: [MKAtomHOF] = []
        
    public init() {
        super.init(fileName: "atomization-energies", subDir: "Data")
        self.readFile()
    }
    
    public override func getSize() -> Int {
        return self._atomhof.count
    }
        
    public override func readFile() {
        //        Try to load contents of file
        guard let filePath = Bundle.module.url(forResource: self._filename, withExtension: "txt", subdirectory: self._subdir) else { return }
        
        filePath.foreachRow { rowContents, lineNum in
            
            if !rowContents.starts(with: "#") {
                let vs = rowContents.split(separator: "|")
                if vs.count >= 8 {
                    let mka = MKAtomHOF(_element: String(vs[0]), _charge: String(vs[1]).toInt()!, _method: String(vs[2]), _desc: String(vs[3]), _T: String(vs[4]).toDouble()!, _value: String(vs[5]).toDouble()!, _multiplicity: String(vs[6]).toInt()!, _unit: String(vs[7]))
                    self._atomhof.append(mka)
                    
                }
            }
        }
    }
    
    static func UnitNameToConversionFactor(_ unit: String) -> Double {
        let p = unit
        switch(p[0]) {
        case "e":
            if (p[1]=="V" && p[2]=="\0") { return ELECTRONVOLT_TO_KCALPERMOL }  // eV
            if (p[1]=="l" && p[2]=="e" && p[3]=="c" && p[4]=="t" && p[5]=="r" && p[6]=="o" && p[7]=="n" &&
                p[8]=="v" && p[9]=="o" && p[10]=="l" && p[11]=="t" && p[12]=="\0") { return ELECTRONVOLT_TO_KCALPERMOL } // electronvolt
        case "k":
            if (p[1]=="J" && p[2]=="/" && p[3]=="m" && p[4]=="o" && p[5]=="l" && p[6]=="\0") { return KJPERMOL_TO_KCALPERMOL }  // kJ/mol
            if (p[1]=="c" && p[2]=="a" && p[3]=="l" && p[4]=="/" && p[5]=="m" && p[6]=="o" && p[7]=="l" && p[8]=="\0") { return 1.0 } // kcal/mol
        case "H":
            if (p[1]=="a" && p[2]=="r" && p[3]=="t" && p[4]=="r" && p[5]=="e" && p[6]=="e" && p[7]=="\0") {return HARTEE_TO_KCALPERMOL } // Hartree
        case "J":
            if (p[1]=="/" && p[2]=="m" && p[3]=="o" && p[4]=="l" && p[5]==" " && p[6]=="K" && p[7]=="\0") { return KJPERMOL_TO_KCALPERMOL } // J/mol K
        case "R":
            if (p[1]=="y" && p[2]=="d" && p[3]=="b" && p[4]=="e" && p[5]=="r" && p[6]=="g" && p[7]=="\0") { return RYDBERG_TO_KCALPERMOL } // Rydberg
        default: break
        }
        print("WARNING: Unknown energy unit in thermochemistry file")
        // Log error
        return 1.0
    }
    
    /** \brief Extract heat of formation and entropy for an atom
           @param elem         The chemical element we're looking for
           @param charge       At this formal charge
           @param method       The method used for computing/measuring
           @param T            The temperature
           @param dhof0        The output energy at 0K
           @param dhof1        The output energy at T
           @param S0T          The entropy at T (it is 0 at 0K)
           \return 1 if the contribution to the Heat of Formation for this atom
           is known at temperature T. If 1 the values
           including all corrections are returned in the dhof variable.
          */
    public func getHeatOfFormation(_ elem: String, _ charge: Int, _ method: String, _ T: Double,
                                   _ dhof0: inout Double, _ dhofT: inout Double, _ S0T: inout Double) -> Int {
        var found = 0
        let Ttol = 0.05 /* Kelvin */
        var Vmodel: Double = 0.0
        var Vdhf: Double = 0.0
        var S0: Double = 0.0
        var HexpT: Double = 0.0
        let desc = "\(method)(0K)"
        for it in _atomhof {
            if ((it._element == elem) && (it._charge == charge)) {
                let eFac = MKAtomicHeatOfFormationTable.UnitNameToConversionFactor(it._unit)
                if (fabs(T - it._T) < Ttol) {
                    if (it._method == "exp") {
                        if (it._desc == "H(0)-H(T)") {
                            HexpT += it._value * eFac
                            found+=1
                        } else if (it._desc == "S0(T)") {
                            S0 += it._value
                            found+=1
                        }
                    }
                } else if ( it._T == 0) {
                    if ((it._method == method) && (it._desc == desc)) {
                        Vmodel += it._value * eFac
                        found+=1
                    }
                    if (it._method == "exp") {
                        if (it._desc == "DHf(T)") {
                            Vdhf += it._value * eFac
                            found+=1
                        }
                    }
                }
            }
        }
        
        if (found == 4) {
            dhof0 = Vdhf-Vmodel
            dhofT = Vdhf-Vmodel-HexpT
            S0T   = -S0/4.184
            return 1
        }

        return 0
    }
}

