//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/27/23.
//

import Foundation

/** \class OBOp op.h <openbabel/op.h>
      \brief Operations to modify molecules before output
      \since version 2.2

Classes derived from OBOp implement options for the obabel program (for both
its commandline and GUI interfaces). It is intended for options that carry out some
modification on the molecule(or reaction) after it has been input, but before
it is output. An example is the --center option implemented in the OpCenter class
in ops.cpp, which is a duplicate of the built in -c option for centering coordinates.

The advantage of plugin classes is that no existing code has to be modified
when a new class is added. You can list those that are present by
obabel -L ops
or from a menu item in the GUI.

Any OBOp derived class has to have a constructor, a function returning a short description,
and a Do() function which does the work. It also needs a WorksWith() function
which is always the same when operating on OBMol objects. (It is not made a default
to reducecode dependencies.) A single global instance of the class needs to be
instantiated to define the ID, by which the class is subsequently accessed.

OBOp works by two of its static functions being called from code in transform.cpp:
 - OpOptions(OBBase* pOb, OpMap* pOptions) returns a string describing each of the
derivated classes relevant to objects of the class of the OBBase parameter,
for use in the help text and to set checkboxes in the GUI;
 - DoOps(OBBase* pOb) applies each option whose ID is listed in the  Opmap parameter
to the object (ususally an OBMol) in the OBBase parameter.

Options which need parameters are passed these (space delimited) in the text parameter
of the Do() function. They can also access other general options specified on the
command line by examining the the OpMap parameter.

To use an OBOp class from the API it is necessary to use an extra step in case it isn't
present. So to apply the OBOp class with ID gen3D to your mol

\code
OBOp* pOp = OBOp::FindType("gen3D");
if(!pOp)
  ...report error
pOp->Do(mol);
\endcode

  */

typealias OpMap = Dictionary<String, String>

class MKOp: MKPlugin, MKPluginProtocol {
    
    static var Default: MKOp?
    static var map: PluginMapType<MKOp> = PluginMapType<MKOp>()
    
    required init(_ id: String, _ isDefault: Bool) {
        super.init()
        self._id = id
        if isDefault || MKOp.map.isEmpty {
            MKOp.Default = self
        }
        if MKOp.map.map({ $0.0 == _id ? 1 : 0}).reduce(0, +) == 0 {
            MKOp.map[_id] = self
            MKPlugin.pluginMap[typeID()] = self
        }
    }
    
    static func findType(_ ID: String?) -> MKOp? {
        if ID == nil {
            return MKOp.Default
        }
        return MKPlugin.baseFindType(MKOp.map, ID!) as? MKOp
    }
    
    override func description() -> String? {
        return nil
    }
    
    override func typeID() -> String {
        return "ops"
    }
    
    override func display(_ txt: inout String, _ param: inout String?, _ ID: String?) -> Bool {
        return false
    }
    
    override func makeInstance(_ v: [String]) -> MKOp? {
        fatalError()
    }
    
    func getMap() -> PluginMapType<MKOp> {
        return MKOp.map
    }
    
    ///Required function that does the work. Normally return true, unless object is not to be output.
    //NOTE: the parameters were changed in r3532
    open func Do(_ pOb: MKBase, _ pptionText: String? = nil, _ pOptions: OpMap? = nil, _ pConv: MKConversion? = nil) -> Bool {return false}

    /// \return true if this op is designed to work with the class of pOb, e.g. OBMol
    func worksWith(_ pOb: MKBase) -> Bool { return false }

    /// Do something with an array of objects. Used a a callback routine in OpSort, etc.
    func processVec(_ vec: [MKBase]) -> Bool { return false }

    /// \return string describing options, for display with -H and to make checkboxes in GUI
    static func opOptions(_ pOb: MKBase) -> String {
        var s: String = ""
        guard let itrOPs = MKPlugin.getPluginIterator("ops") else { return "" }
        for itr in itrOPs {
            guard let pOp: MKOp = itr.value as? MKOp else { continue }
            //ignore ops with IDs that begin with '_' or have "not displayed in GUI" in their first line of description
            if ((itr.key) == "_") || MKPlugin.firstLine(pOp.description()!).contains("not displayed in GUI") {
                continue
            }
            if(pOp.worksWith(pOb)) {
                s += "--"
                s += itr.key //ID
                s += " "
                s += MKPlugin.firstLine(pOp.description()!) + "\n"
            }
        }
        s += "\n"
        return s
    }

    ///Call Do() of all the OBOps whose ID is a key in the map.
    ///Called from DoTransformations(). The map has general options like -x or --multicharoption
    ///The key is the option name and the value, if any, is text which follows the option name.
    /// In some cases, there may be several parameters, space separated)
    /// \return false indicating object should not be output, if any Do() returns false
    static func doOps(_ pOb: MKBase,_ pOptions: OpMap ,_ pConv: MKConversion) -> Bool {
        for itr in pOptions {
            guard let pOp = findType(itr.key) else { continue }
            if(!pOp.Do(pOb, itr.value, pOptions, pConv)) {
                return false; //Op has decided molecule should not be output
            }
        }
        return true
    }
}

    

