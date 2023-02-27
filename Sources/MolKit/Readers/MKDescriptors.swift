//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/25/23.
//

import Foundation
import OrderedCollections

protocol MKDescriptorProtocol {
    
    func predict(_ pOb: MKBase, _ param: String?) -> Double
    func getStringValue(_ pOb: MKBase, _ svalue: inout String?, _ param: inout String?) -> Double
    func compare(_ pOb: MKBase, _ optionText: Iterator<Character>, _ noEval: Bool, _ param: String?) -> Bool
     
}

/**
OBDescriptor and Filtering

On the command line, using the option --filter filter-string converts only
those molecules which meet the criteria specified in the filter-string. This
is useful to select particular molecules from a set.
It is used like:
babel dataset.sdf outfile.smi --filter "MW>200 SMARTS!=c1ccccc1 PUBCHEM_CACTVS_ROTATABLE_BOND<5"

The identifier , "PUBCHEM_CACTVS_ROTATABLE_BOND" is the name of an attribute
of an OBPairData which has probably been imported from a property in a SDF
or CML file. The identifier names are (currently) case dependent. A comparison
is made with the value in the OBPairData. This is a numeric comparison if both
operands can be converted to numbers (as in the example). If the 5 had been
enclosed in single or double quotes the comparison would have been a string
comparison, which gives a different result in some cases. OBPairData is searched
first to match an identifier.

If there are no OBPair attributes that match, the identifier is taken to be the
ID of an OBDescriptor class object. The class OBDescriptor is the base class
for classes that wrap molecular properties, descriptors or features. In the example
"MW" and "SMARTS" are OBDescriptor IDs and are case independent. They are plugin
classes, like fingerprints, forcefields and formats, so that new molecular features
can be added or old ones removed (to prevent code bloat) without altering old code.
A list of available descriptors is available from the commandline:
babel -L descriptors
or from the functions OBPlugin::List, OBPlugin::ListAsString and OBPlugin::ListAsVector.

The filter-string is interpreted by a static function of OBDescriptor,
FilterCompare(). This identifies the descriptor IDs and then calls a virtual
function, Compare(), of each OBDescriptor class to interpret the rest of relational
expression, for example, ">200", or "=c1ccccc1". The default version of Compare()
is suitable for descriptors, like MW or logP, which return a double from
their Predict() method. Classes like SMARTS which need different semantics
provide their own.

By default, as in the example, OBDescriptor::FilterCompare() would AND each
comparison so that all the comparisons must be true for the test to succeed.
However filter-string could also be a full boolean expression, with &, |, !,
and parenthases allowing any combination of features to be selected.
FilterCompareAs calls itself recursively to give AND precidence over OR and
evaluation is not carried out if not needed.

The aim has been to make interpretation of the filter-string as liberal as
possible, so that AND can be &&, there can be spaces or commas in places
that are reasonable.

The base class, OBDescriptor, uses pointers to OBBase in its functions,
like OBFormat, to improve extendability - reactions could have
features too. It does mean that a dynamic_cast is needed at the start of the
Predict(OBBase* pOb, string*) functions.

To use a particular descriptor, like logP, when programming with the API, use
code like the following:
\code
  OBDescriptor* pDescr = OBDecriptor::FindType("logP");
  if(pDescr)
    double val = pDescr->Predict(mol, param);
\endcode
To add the descriptor ID and the predicted data to OBPairData attached to
the object, use PredictAndSave().

Descriptors can have a string parameter, which each descriptor can interpret
as it wants, maybe, for instance as multiple numeric values. The parameter
is in brackets after the descriptor name, e.g. popcount(FP4). In the above
programming example param is a pointer to a std::string which has a default
value of NULL, meaning no parameter. GetStringValue() and Compare() are similar.

To parse a string for descriptors use GetIdentifier(), which returns both
the ID and the parameter, if there is one.

This facility can be called from the command line.Use the option
--add "descriptor list", which will add the requested descriptors to the
molecule.  They are then visible as properties in SDF and CML formats.
The IDs in the list can be separated by spaces or commas.
All Descriptors will provide an output value as a string through a  virtual
function GetStringValue((OBBase* pOb, string& svalue)) which
assigns the value of a string descriptor(like inchi) to svalue or a string
representation of a numerical property like logP.

The classes MWFilter and TitleFilter illustrate the code that has to be
provided for numerical and non-numerical descriptors.

*/

//! \brief Base class for molecular descriptors
class MKDescriptor: MKPluginProtocol, MKDescriptorProtocol {

    static var Default: MKDescriptor?
    static var map: PluginMapType<MKDescriptor> = PluginMapType<MKDescriptor>()
    
    var _id: String = ""
    
    func getID() -> String {
        return _id
    }
    
    static func findType(_ ID: String?) -> (any MKPluginProtocol)? {
       if ID == nil {
            return MKDescriptor.Default
        }
        return MKPlugin.baseFindType(getMap(), ID!)
    }
    
    func description() -> String? {
        return nil
    }
    
    func typeID() -> String {
        return "descriptors"
    }
    
    ///Write information on a plugin class to the string txt.
    ///If the parameter is a descriptor ID, displays the verbose description for that descriptor only
    /// e.g. babel -L descriptors HBA1
    static func display(_ txt: inout String, _ param: inout String, _ ID: String?) -> Bool {
        //Use the base class version except when the parameter is a descriptor ID.
        //For a parameter which is the matching descriptor set verbose.
        //No display for other descriptors.
        //Allows babel descriptors HBA1
        if !param.isEmpty && (MKDescriptor.findType(param) != nil) {
            if ID == param {
                return false
            }
            param = "verbose"
        }
        return MKPlugin.display(&txt, &param, ID)
    }
    
    func makeInstance(_ v: [String]) -> (any MKPluginProtocol)? {
        return nil 
    }
    
    static func getMap() -> PluginMapType<MKDescriptor> {
        return MKDescriptor.map
    }
    
    func predict(_ pOb: MKBase, _ param: String?) -> Double {
        return Double.nan
    }

    /// \return the value of the descriptor and adds it to the object's OBPairData
    func predictAndSave(_ pOb: MKBase, _ param: inout String?) -> Double {
        var attr = getID()
        var svalue: String? = ""
        var val: Double = getStringValue(pOb, &svalue, &param)
        var dp = pOb.getData(attr) as? MKPairData<String>
        var previouslySet = true 
        if dp == nil {
            previouslySet = false 
            dp = MKPairData<String>()
        }
        dp!.setAttribute(attr)
        dp!.setValue(svalue ?? "")
        dp!.setOrigin(.perceived)
        if !previouslySet {
            pOb.setData(dp!)
        }
        return val
    }

    /// Interprets the --filter option string and returns the combined result of all the comparisons it contains
    /**
        The string has the form:
        PropertyID1 predicate1 [booleanOp] PropertyID2 predicate2 ...
        The propertyIDs are the ID of instances of a OBDescriptor class or
        the Attributes of OBPairData, and contain only letters, numbers and underscores.
        The predicates must start with a punctuation character and are interpreted by
        the Compare function of the OBDescriptor class. The default implementation expects
        a comparison operator and a number, e.g. >=1.3  Whitespace is optional and is ignored.
        Each predicate and this OBBase object (usually OBMol) is passed to
        the Compare function of a OBDescriptor. The result of each comparison
        is combined in a boolean expression (which can include parentheses)
        in the normal way. The AND operator can be & or &&, the OR operator can be
        | or ||, and a unitary NOT is !  The expected operator precedence
        is achieved using recursive calls of the function. If there is no boolean Op, all
        the tests have to return true for the function to return true, i.e. the default is AND.
        If the first operand of an AND is 0, or of an OR is 1, the parsing of the second operand
        continues but no comparisons are done since the result does not matter.
    **/
    static func filterCompare(_ pOb: MKBase, _ optionText: Iterator<Character>, _ noEval: inout Bool) -> Bool {
        while true {
            
        var negate: Bool=false
        var retFromCompare: Bool = false 
        var ret: Bool=false
        var ch: Character? 
        if optionText.isEmpty() { return false }
        
        repeat { 
            ch = optionText.next() 
        } while (ch != nil) ? ch!.isWhitespace : false

        if ch == "!" { 
            negate = true 
            ch = optionText.next()
        }

        if ch == "(" {
            // bracketed expression 
            retFromCompare = filterCompare(pOb, optionText, &noEval) //noEval persists in subsidiary calls
            ch = optionText.next()
            if ch != ")" { 
                print("Missing ')' in filter string")
                return retFromCompare  //missing closing bracket
            } 
        } else { // unbracketed expression 

            if !MKDescriptor.ispunctU(ch!) { //must be start of ID
                optionText.unget()
            } else {
                let mes: String = "Filter string has erroneous character : "
                print(mes + String(ch!))
                optionText.setEmpty()
                return false
            }

            let spair = getIdentifier(optionText)
            var descID = spair.0
            let param = spair.1!
            if descID.isEmpty {
                print("Filter string has no descriptor ID")
                optionText.setEmpty()
                return false // MARK: should show error
            }

            //If there is existing OBPairData use that
            if param.isEmpty && matchPairData(pOb, &descID) {
                var value: String = (pOb.getData(descID) as! MKPairData<String>).getValue()!
                retFromCompare = compareStringWithFilter(optionText, &value, noEval, true)
            } else {
                //if no existing data see if it is an OBDescriptor
                let pDesc = MKDescriptor.findType(descID) as? MKDescriptor
                if pDesc != nil && !noEval {
                    retFromCompare = pDesc!.compare(pOb, optionText, noEval, param)
                } else {
//                     just parse
                    var ch1: Character? = nil 
                    var ch2: Character? = nil
                    var svalue: String = "" 
                    parsePredicate(optionText, &ch1, &ch2, &svalue)
                    //no existing data, not a descriptor result is false meaning "does not exist"
                    retFromCompare = false 
                }
            }
        }

        if negate {
            retFromCompare = !retFromCompare
        }
       
        if !noEval {
            ret = retFromCompare
        }
        
        //Look for boolean operator
        ch = optionText.next()
        if ch == nil {
            return ret // end of filterString
        }

        if ch == ")" {
            optionText.unget()
            return ret // end of bracketed expression
        }

            if !MKDescriptor.ispunctU(ch!) {
            optionText.unget()
        } else {
            if optionText.peek() == ch { //treat && and || as & and |
                optionText.ignore()
            }
        }

        if ch == "|" {
            noEval = ret || noEval
            retFromCompare = filterCompare(pOb, optionText, &noEval)
            return !noEval && (ret || retFromCompare) //always return false if noEval=true;
        } else { //includes & and , and ;
            noEval = !ret //if ret is false keep parsing but don't bother to evaluate
        }
        }//go for next conditional expression
        return false //never come here
    }

    ///Reads list of descriptor IDs and calls PredictAndSave() for each.
    static func addProperties(_ pOb: MKBase, _ DescrList: String) {
        let ss: Iterator<Character> = Iterator([Character](DescrList))
        var pDescr: MKDescriptor? = nil
        while !ss.isEmpty() {
            var spair = getIdentifier(ss)
            pDescr = MKDescriptor.findType(spair.0) as? MKDescriptor 
            if pDescr != nil {
                pDescr!.predictAndSave(pOb, &spair.1)
            } else {
                print("not recognized as a descriptor")
                // TODO: throw error
            }
        }
    }

    ///Deletes all the OBPairDatas whose attribute names are in the list (if they exist).
    static func deleteProperties(_ pOb: MKBase, _ DescrList: String) {
        var vs: [String] = DescrList.components(separatedBy: CharacterSet(charactersIn: "\t\r\n,/-*&;:|%+"))
        for var token in vs {
            if matchPairData(pOb, &token) {
                pOb.deleteData(token)
            }
        }
    }

    //Reads list of descriptor IDs and OBPairData names and returns a list of values
    //each preceded by a space or the first character in the list if it is whitespace or punctuation.
    //Used in OBMol::Transform() to append to title , but that is not done here to avoid
    //having to #include mol.h in this file.
    static func getValues(_ pOb: MKBase, _ DescrList: String) -> String {
        let ss: Iterator<Character> = Iterator<Character>([Character](DescrList))
        var delim: Character = ss.first!
        
        if delim.isWhitespace || ispunctU(delim) {
            ss.ignore()
            if delim == "\\" {
                if ss[1] == "\\" {
                    ss.ignore()
                } else if ss[1] == "t" {
                    delim = "\t"
                    ss.ignore()
                }
            } else {
                delim = " "
            }
        }
        
        var values: String = ""
        var pDescr: MKDescriptor?

        while !ss.isEmpty() {
            var thisvalue: String? = ""
            var spair = getIdentifier(ss)
            //If there is existing OBPairData use that
            if matchPairData(pOb, &spair.0) {
                thisvalue = (pOb.getData(spair.0)! as! MKPairData<String>).getValue()
            } else {
                pDescr = MKDescriptor.findType(spair.0) as? MKDescriptor
                if pDescr != nil {
                    pDescr!.getStringValue(pOb, &thisvalue, &spair.1)
                } else {
                    // Mark throw error
                    print("Descriptor \(spair.0) not recognized as a property or a descriptor")
                    thisvalue = "??"
                }
            }

            values += String(delim) + (thisvalue ?? "??")
        }
        return values
    }

    ///Read an identifier and its parameter from the filter string.
    static func getIdentifier(_ optionText: Iterator<Character>) -> (String, String?) {
        var descID: String = ""
        var param: String? = nil
        var ch: Character? = optionText.next()

        while !optionText.isEmpty() {
            if ch == nil || ch!.isWhitespace || ch == "," {
                break
            }
            if ch == "(" {  // the parameter is in parentheses
                ch = optionText.peek()
                if ch == "\"" || ch == "\'" {
                    // parameter is in parentheses
                    optionText.ignore()
                    param = String(optionText.nextUntil(ch!)!)
                    optionText.ignore(until: ")")
                } else {
                    param = String(optionText.nextUntil(")")!)
                }
                if optionText.isEmpty() {
                    print("Missing ')' in descriptor parameter")
                    descID = ""
                    return (descID, descID)
                }
            } else if MKDescriptor.ispunctU(ch!) {
                optionText.unget()
                break
            } else {
                descID.append(ch!)
            }
            ch = optionText.next()
        }
        optionText.setEmpty()
        return (descID, param)
    }
    
    func getStringValue(_ pOb: MKBase, _ svalue: inout String?, _ param: inout String?) -> Double {
        let val: Double = predict(pOb, param)
        if svalue != nil {
            svalue = String(val)
        }
        return val
    }
    
    func compare(_ pOb: MKBase, _ optionText: Iterator<Character>, _ noEval: Bool, _ param: String?) -> Bool {
        // Scan the optionText until the first occurance of an punctuation character

        let ch1: Character? = optionText.first(where: {MKDescriptor.ispunctU($0)})
        guard let ch1 else { return false }
        let ch2: Character? = optionText.nextElement(ch1, updateIndex: true)
        
//         get number
        var val: Double
        var filterval: Double?
        filterval = optionText.parseDouble()
        guard let filterval else { return false }
        if !optionText.isEmpty() {
            if noEval {
                return false
            }
            val = predict(pOb, param)
            return doComparison(ch1, ch2, val, filterval)
        }
        optionText.setEmpty()
        print("Error in filter string \(optionText) - no value after comparison operator")
        return false
    }
    
    // static double ParsePredicate(std::istream& optionText, char& ch1, char& ch2, std::string& svalue);
    static func parsePredicate(_ optionText: Iterator<Character>, _ ch1: inout Character?, _ ch2: inout Character?, _ svalue: inout String) -> Double {
        var val: Double = Double.nan
        ch2 = nil
        ch1 = nil
        //Get comparison operator
        ch1 = optionText.next()
        if ch1 == nil || (ch1!.isLetter || ch1!.isNumber) || ch1 == "&" || ch1 == "|" || ch1 == ")" {
            //no comparison operator
            optionText.unget()
            optionText.setEmpty() //not an error to reach eof
            ch1 = nil
            return val
        } else {
            if optionText.peek() == "=" {
                ch2 = optionText.next()
            }
        }
        //Try to read a double. Rewind and read as a string
        let spos = optionText.tellg()
        val = optionText.parseDouble() ?? Double.nan
        //only a number when the param has no additional text or only a closing bracket
        if !optionText.isEmpty() && (optionText.peek().isLetter || optionText.peek().isNumber) {
            val = Double.nan
        }

        optionText.setEmpty()
        optionText.seekg(spos)
        readStringFromFilter(optionText, &svalue)
        return val
    }

    /// Reads a string from the filter string  optionally preceded by = or !=
    /** On entry the stringstream position should be just after the ID. On exit it is after the string.
        If there is an error, the stringstream badbit is set.
        Returns false if != found, to indicate negation.
        Can be of any of the following forms:
        mystring  =mystring ==mystring [must be terminated by a space or tab]
        "mystring" 'mystring'  ="mystring" ='mystring' [mystring can contain spaces or tabs]
        !=mystring !="mystring" [Returns false indicating negate]
        There can be spaces or tabs after the operator = == !=
    **/
    static func readStringFromFilter(_ optionText: Iterator<Character>, _ result: inout String) -> Bool {
        var ret = true
        var ch: Character? = optionText.next()
        if ch != nil {
            if ch == "=" || ch == "!" {
                if optionText.peek() != "=" {
                    optionText.unget()
                }
                if ch == "!" {
                    ret = false //to indicate negation
                }
            } else { //no operator
                optionText.unget()
            }
            ch = optionText.next()
            if ch == "\"" || ch == "\'" {
                result = String(optionText.nextUntil(ch!)!) //get quoted text
            } else { // not quoted; get string up to next space or ')'
                optionText.unget()
                result = ""
                ch = optionText.nextWhile({ $0.isWhitespace }) //ignore leading white space
                while ch != nil && !ch!.isWhitespace && ch != ")" {
                    result.append(ch!)
                    ch = optionText.next()
                }
            }
        }
        // TODO: A stream should be used to check for potential errors in parsing 
        return ret
    }

    ///Makes a comparison using the operator and a string read from the filter stream with a provided string.
    /// \return the result of the comparison and true if NoCompOK==true and there is no comparison operator.
    static func compareStringWithFilter(_ optionText: Iterator<Character>, _ sval: inout String, _ noEval: Bool, _ NoCompOK: Bool = false) -> Bool {
        /** Convert to swift
         char ch1=0, ch2=0;
           string sfilterval;
           double filterval = ParsePredicate(optionText, ch1, ch2, sfilterval);
           if(ch1==0 && NoCompOK)
           {
             // there is no comparison operator
             return true; // means that the identifier exists
           }

           stringstream ss(sval);
           double val;
           if((ss >> val) && !IsNan(filterval))
             //Do a numerical comparison if both values are numbers
             return DoComparison(ch1, ch2, val, filterval);
           else
           {
             //Do a string comparison if either the filter or the OBPair value is not a number
             //If sval is quoted remove quotes
             if(sval[0]=='\"' || sval[0]=='\'')
               sval.erase(0,1);
             if(sval[sval.size()-1]=='\"' || sval[sval.size()-1]=='\'')
               sval.erase(sval.size()-1);

             bool leading=false, trailing=false;
             if(sfilterval[0]=='*')
             {
               leading=true;
               sfilterval.erase(0,1);
             }
             if(sfilterval[sfilterval.size()-1]=='*')
             {
               trailing=true;
               sfilterval.erase(sfilterval.size()-1);
             }
             string::size_type pos = sval.find(sfilterval);
             if(pos!=string::npos)
             {
               if(trailing) //needs to be first
                 sval.erase(pos+sfilterval.size()); //erase aftermatch
               if(leading)
                 sval.erase(0, pos); //erase before match
             }
             return DoComparison(ch1, ch2, sval, sfilterval);
           }
         */
        var ch1: Character?
        var ch2: Character?
        var sfilterval: String = ""
        let filterval = parsePredicate(optionText, &ch1, &ch2, &sfilterval)

        if ch1 == nil && NoCompOK {
            // there is no comparison operator
            return true // means that the identifier exists
        }

        var val: Double? = Iterator<Character>([Character](sval)).parseDouble()
        // parse next double value from the string by using Iterator
        if val != nil && !filterval.isNaN {
            // Do a numerical comparison if both values are numbers
            return doComparison(ch1!, ch2, val!, filterval)
        } else {
            // Do a string comparison if either the filter or the OBPair value is not a number
            // If sval is quoted remove quotes
            if sval.first == "\"" || sval.first == "\'" {
                sval.removeFirst()
            }
            if sval.last == "\"" || sval.last == "\'" {
                sval.removeLast()
            }
            
            var leading = false
            var trailing = false

            if sfilterval.first == "*" {
                leading = true
                sfilterval.removeFirst()
            }

            if sfilterval.last == "*" {
                trailing = true
                sfilterval.removeLast()
            }

            if let pos = sval.index(of: sfilterval) {
                if trailing { // needs to be first
                    sval.removeSubrange(pos...sval.endIndex) // erase aftermatch
                }
                if leading { // needs to be second
                    sval.removeSubrange(sval.startIndex...pos) // erase before match
                }
            }
            return doComparison(ch1!, ch2, sval, sfilterval)
        }
    }

    // Treats _ as not a punctuation character and since 2.3.2 also $ # and %
    static func ispunctU(_ ch: Character) -> Bool {
        return ch.isPunctuation && ch != "_" && ch != "$" && ch != "#" && ch != "%"
    }

    /// \return true if s (with or without _ replaced by spaces) is a PairData attribute. On return s is the form which matches.
    static func matchPairData(_ pOb: MKBase, _ s: inout String) -> Bool {
        //If s matches a PairData attribute return true
        //else if s with all '_' replaced by spaces matches return true and s is now the form with spaces
        //else return false.
        if pOb.hasData(s) {
            return true
        }
        // contains _ or space
        if !s.contains("_") && !s.contains(" ") { return false }

        var temp = s
        temp = temp.replacingOccurrences(of: "_", with: " ")
        if pOb.hasData(temp) {
            s = temp
            return true
        }
        return false
    }    
    
}

func doComparison<T: Comparable>(_ ch1: Character, _ ch2: Character?, _ val: T, _ filterval: T) -> Bool {
    switch ch1 {
    case "0", "=": return val == filterval
    case "!": return val != filterval
    case ">":
        guard let ch2 = ch2 else { return false }
        return ch2 == "=" ? val >= filterval : val > filterval
    case "<":
        guard let ch2 = ch2 else { return false }
        return ch2 == "=" ? val <= filterval : val < filterval
    default: return false
    }
}
