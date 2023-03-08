//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 3/7/23.
//

import Foundation
//! \brief Supports a set of rotamer coordinate sets for some number of potentially rotatable bonds
// Further class introduction in rotamer.cpp
class MKRotamerList: MKGenericData {

    //! Number of atoms in the base coordinate set (i.e., OBMol::NumAtoms())
    // unsigned int                         _NBaseCoords;
    // //! Base coordinate sets (i.e., existing conformers to be modified)
    // std::vector<double*>                 _c;
    // //! Individual bond rotors (from an OBRotor object or other)
    // std::vector<std::pair<OBAtom**,std::vector<int> > > _vrotor;
    // //! \brief Index of each rotor's different sampling states ("resolution")
    // //! Usually from OBRotor::GetResolution()
    // std::vector<std::vector<double> >    _vres;
    // //! Individual rotamer states (i.e., the array of rotor settings)
    // std::vector<unsigned char*>          _vrotamer;
    // //! Rotors in rings
    // std::vector<std::vector<int> >       _vrings;
    // //! Dihedral angles of ring bonds
    // std::vector<std::vector<double> >    _vringTors;
    var _NBaseCoords: Int = 0                 //! Number of atoms in the base coordinate set (i.e., OBMol::NumAtoms())
    var _c: [Double]  = []                    //! Base coordinate sets (i.e., existing conformers to be modified)
    var _vrotor: [Pair<[MKAtom], [Int]>] = [] //! Individual bond rotors (from an OBRotor object or other)
    var _vres: [[Double]] = []                //! \brief Index of each rotor's different sampling states ("resolution")
                                              //! Usually from OBRotor::GetResolution()
    var _vrotamer: [[Int]] = []               //! Individual rotamer states (i.e., the array of rotor settings)
    var _vrings: [[Int]] = []                 //! Rotors in rings
    var _vringTors: [[Double]] = []           //! Dihedral angles of ring bonds

    init() {
        super.init("RotamerList", .RotamerList, .any)
        _NBaseCoords = 0 
    }

    





}
