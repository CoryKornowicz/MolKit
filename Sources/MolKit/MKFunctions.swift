//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 2/19/23.
//

import Foundation
import Bitset

public func getTypicalValence(_ ele: UInt, _ bosum: UInt, _ charge: Int) -> UInt {
    
    switch ele {
    case 1:
        switch charge {
        case 0:
            if bosum <= 1 { return 1 }
        case 1:
            if bosum == 1 { return 0 }
        default: break
        }
    case 2:
        if charge == 0 {
            if bosum == 0 { return 0 }
        }
    case 3:
        switch charge {
        case 0:
            if bosum <= 1 { return 1 }
        case 1:
            if bosum == 0 { return 0 }
        default: break
        }
    case 4:
        switch charge {
        case 0:
            if bosum <= 2 { return 2 }
        case 1:
            if bosum <= 1 { return 1 }
        case 2:
            if bosum == 0 { return 0 }
        default: break
        }
    case 5:
        switch charge {
        case -2:
            if (bosum <= 3) { return 3 }
        case -1:
            if (bosum <= 4) { return 4 }
        case 0:
            if (bosum <= 3) { return 3 }
        case 1:
            if (bosum <= 2) { return 2 }
        case 2:
            if (bosum <= 1) { return 1 }
        default: break
        }
    case 6:
        switch charge {
        case -2:
            if bosum <= 2 { return 2 }
        case -1:
            if bosum <= 3 { return 3 }
        case 0:
            if bosum <= 4 { return 4 }
        case 1:
            if bosum <= 3 { return 3 }
        case 2:
            if bosum <= 2 { return 2 }
        default: break
        }
    case 7:
        // Note that while N can have valence 5, it doesn't make sense
        // to round up to 5 when adding hydrogens
        switch charge {
        case -2:
            if bosum <= 1 { return 1 }
        case -1:
            if bosum <= 2 { return 2 }
        case 0:
            switch bosum {
            case 0, 1, 2, 3: return 3
            case 4: return 4 // don't round up to 5 for nitrogen
            case 5: return 5
            default: break
            }
        case 1:
            if bosum <= 4 { return 4 }
        case 2:
            if bosum <= 3 { return 3 }
        default: break
        }
    case 8:
        switch (charge) {
        case -2:
            if (bosum == 0) { return 0 }
            break
        case -1:
            if (bosum <= 1) { return 1 }
            break
        case 0:
            if (bosum <= 2) { return 2 }
            break
        case 1:
            switch (bosum) {
            case 0, 1, 2, 3: return 3
            case 4, 5: return 5
            default: break
            }
            break
        default: break
        }
        break
    case 9:
        switch charge {
        case -1:
            if bosum == 0 { return 0 }
        case 0:
            if bosum <= 1 { return 1 }
        case 1:
            if bosum <= 2 { return 2 }
        case 2:
            switch bosum {
            case 0, 1, 2, 3: return 3
            case 4, 5: return 5
            default: break
            }
        default: break
        }
    case 10:
        if (charge == 0) {
            if (bosum == 0) {return 0 }
        }
    case 11:
        switch (charge) {
        case -1:
            if (bosum == 0) { return 0 }
        case 0:
            if (bosum <= 1) { return 1 }
        case 1:
            if (bosum == 0) { return 0 }
        default: break
        }
    case 12:
        switch charge {
        case 0:
            if bosum <= 2 {
                return 2
            }
        case 2:
            if bosum == 0 {
                return 0
            }
        default: break
        }
    case 13:
        switch charge {
        case -2:
            switch bosum {
            case 0, 1, 2, 3: return 3
            case 4, 5: return 5
            default: break
            }
        case -1:
            if bosum <= 4 { return 4}
            break;
        case 0:
            if bosum <= 3 { return 3}
            break;
        case 1:
            if bosum <= 2 { return 2}
            break;
        case 2:
            if bosum <= 1 { return 1}
            break;
        case 3:
            if bosum == 0 { return 0}
            break;
        default: break
        }
        break;
    case 14:
        switch charge {
        case -2:
            if bosum <= 2 { return 2}
            break;
        case -1:
            switch bosum {
            case 0, 1, 2, 3: return 3
            case 4, 5: return 5
            default: break
            }
            break;
        case 0:
            if bosum <= 4 { return 4}
            break;
        case 1:
            if bosum <= 3 { return 3}
            break;
        case 2:
            if bosum <= 2 { return 2}
            break;
        default: break
        }
        break;
    case 15:
        switch charge {
        case -2:
            switch bosum {
            case 0, 1: return 1
            case 2, 3: return 3
            case 4, 5: return 5
            case 6, 7: return 7
            default: break
            }
            break
        case -1:
            switch bosum {
            case 0, 1, 2: return 2
            case 3, 4: return 4
            case 5, 6: return 6
            default: break
            }
            break
        case 0:
            switch bosum {
            case 0, 1, 2, 3: return 3
            case 4, 5: return 5
            default: break
            }
            break
        case 1:
            if bosum <= 4 { return 4 }
            break
        case 2:
            if bosum <= 3 { return 3 }
            break
        default: break
        }
    case 16:
        switch charge {
        case -2:
            if bosum == 0 { return 0 }
            break;
        case -1:
            switch bosum {
            case 0,1: return 1
            case 2,3: return 3
            case 4,5: return 5
            case 6,7: return 7
            default: break
            }
            break;
        case 0:
            switch bosum {
            case 0,1,2: return 2
            case 3,4: return 4
            case 5,6: return 6
            default: break
            }
            break;
        case 1:
            switch bosum {
            case 0,1,2,3: return 3
            case 4,5: return 5
            default: break
            }
            break;
        case 2:
            if bosum <= 4 { return 4 }
            break;
        default: break
        }
        break;
    case 17:
        switch (charge) {
        case -1:
            if (bosum == 0) { return 0 }
            break;
        case 0:
            switch (bosum) {
            case 0, 1: return 1;
            case 2, 3: return 3;
            case 4, 5: return 5;
            case 6, 7: return 7;
            default: break
            }
            break;
        case 1:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case 2:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        default: break
        }
        break;
    case 18:
        if (charge == 0) {
            if (bosum == 0) { return 0 }
        }
        break;
    case 19:
        switch (charge) {
        case -1:
            if (bosum == 0) { return 0 }
            break;
        case 0:
            if (bosum <= 1) { return 1 }
            break;
        case 1:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    case 20:
        switch (charge) {
        case 0:
            if (bosum <= 2) { return 2}
            break;
        case 1:
            if (bosum <= 1) { return 1}
            break;
        case 2:
            if (bosum == 0) { return 0}
            break;
        default: break
        }
    case 31:
        switch (charge) {
        case -2:
            switch (bosum) {
            case 0, 1, 2, 3: return 3
            case 4, 5: return 5
            default: break
            }
            break;
        case -1:
            if (bosum <= 4) { return 4 }
            break;
        case 0:
            if (bosum <= 3) { return 3 }
            break;
        case 1:
            if (bosum == 0) { return 0 }
            break;
        case 2:
            if (bosum <= 1) { return 1 }
            break;
        case 3:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    case 32:
        switch (charge) {
        case -2:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case -1:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        case 0:
            if (bosum <= 4) { return 4 }
            break;
        case 1:
            if (bosum <= 3) { return 3 }
            break;
        case 4:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    case 33:
        switch charge {
        case -3:
            if bosum == 0 {return 0}
        case -2:
            switch bosum {
            case 0, 1: return 1
            case 2, 3: return 3
            case 4, 5: return 5
            case 6, 7: return 7
            default: break
            }
        case -1:
            switch bosum {
            case 0, 1, 2: return 2
            case 3, 4: return 4
            case 5, 6: return 6
            default: break
            }
        case 0:
            switch bosum {
            case 0, 1, 2, 3: return 3
            case 4, 5: return 5
            default: break
            }
        case 1:
            if bosum <= 4 {return 4}
        case 2:
            if bosum <= 3 {return 3}
        default: break
        }
    case 34:
        switch (charge) {
        case -2:
            if (bosum == 0) { return 0 }
            break;
        case -1:
            switch (bosum) {
            case 0, 1:  return 1
            case 2, 3:  return 3
            case 4, 5:  return 5
            case 6, 7:  return 7
            default: break
            }
            break;
        case 0:
            switch (bosum) {
            case 0, 1, 2:  return 2
            case 3, 4:  return 4
            case 5, 6:  return 6
            default: break
            }
            break;
        case 1:
            switch (bosum) {
            case 0, 1, 2, 3:  return 3
            case 4, 5:  return 5
            default: break
            }
            break;
        case 2:
            if (bosum <= 4) { return 4 }
            break;
        default: break
        }
        break;
    case 35:
        switch (charge) {
        case -1:
            if (bosum == 0) { return 0 }
            break;
        case 0:
            switch (bosum) {
            case 0, 1: return 1;
            case 2, 3: return 3;
            case 4, 5: return 5;
            case 6, 7: return 7;
            default: break
            }
            break;
        case 1:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case 2:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        default: break
        }
        break;
    case 36:
        if (charge == 0) {
            switch (bosum) {
            case 0: return 0;
            case 1, 2: return 2;
            default: break
            }
        }
        break;
    case 37:
        switch (charge) {
        case -1:
            if (bosum == 0) { return 0}
            break;
        case 0:
            if (bosum <= 1) { return 1}
            break;
        case 1:
            if (bosum == 0) { return 0}
            break;
        default: break
        }
        break;
    case 38:
        switch (charge) {
        case 0:
            if (bosum <= 2) { return 2}
            break;
        case 1:
            if (bosum <= 1) { return 1}
            break;
        case 2:
            if (bosum == 0) { return 0}
            break;
        default: break
        }
        break;
    case 49:
        switch (charge) {
        case -2:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        case -1:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            default: break
            }
            break;
        case 0:
            if (bosum <= 3) { return 3 }
            break;
        case 1:
            if (bosum == 0) { return 0 }
            break;
        case 2:
            if (bosum <= 1) { return 1 }
            break;
        case 3:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    case 50:
        switch (charge) {
        case -2:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case -1:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        case 0:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            default: break
            }
            break;
        case 1:
            if (bosum <= 3) { return 3 }
            break;
        case 2:
            if (bosum == 0) { return 0 }
            break;
        case 4:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    case 51:
        switch (charge) {
        case -2:
            switch (bosum) {
            case 0, 1: return 1;
            case 2, 3: return 3;
            case 4, 5: return 5;
            case 6, 7: return 7;
            default: break
            }
            break;
        case -1:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case 0:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        case 1:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            default: break
            }
            break;
        case 2:
            if (bosum <= 3) {return 3 }
            break;
        case 3:
            if (bosum == 0) {return 0 }
            break;
        default: break
        }
        break;
    case 52:
        switch (charge) {
        case -2:
            if (bosum == 0) {return 0 }
            break;
        case -1:
            switch (bosum) {
            case 0, 1: return 1;
            case 2, 3: return 3;
            case 4, 5: return 5;
            case 6, 7: return 7;
            default: break
            }
            break;
        case 0:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case 1:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        case 2:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            default: break
            }
            break;
        default: break
        }
        break;
    case 53:
        switch (charge) {
        case -1:
            if (bosum == 0) {return 0 }
            break;
        case 0:
            switch (bosum) {
            case 0, 1: return 1;
            case 2, 3: return 3;
            case 4, 5: return 5;
            case 6, 7: return 7;
            default: break
            }
            break;
        case 1:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case 2:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        default: break
        }
        break;
    case 54:
        if (charge == 0) {
            switch (bosum) {
            case 0: return 0;
            case 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            case 7, 8: return 8;
            default: break
            }
        }
        break;
    case 55:
        switch (charge) {
        case -1:
            if (bosum == 0) {return 0 }
            break;
        case 0:
            if (bosum <= 1) {return 1 }
            break;
        case 1:
            if (bosum == 0) {return 0 }
            break;
        default: break
        }
        break;
    case 56:
        switch (charge) {
        case 0:
            if (bosum <= 2) {return 2 }
            break;
        case 1:
            if (bosum <= 1) {return 1 }
            break;
        case 2:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    case 81:
        if (charge == 0) {
            switch (bosum) {
            case 0, 1: return 1;
            case 2, 3: return 3;
            default: break
            }
        }
        break;
    case 82:
        switch (charge) {
        case -2:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case -1:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        case 0:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            default: break
            }
            break;
        case 1:
            if (bosum <= 3) { return 3 }
            break;
        case 2:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    case 83:
        switch (charge) {
        case -2:
            switch (bosum) {
            case 0, 1: return 1;
            case 2, 3: return 3;
            case 4, 5: return 5;
            case 6, 7: return 7;
            default: break
            }
            break;
        case -1:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case 0:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        case 1:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            default: break
            }
            break;
        case 2:
            if (bosum <= 3) { return 3 }
            break;
        case 3:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    case 84:
        if (charge == 0) {
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
        }
        break;
    case 85:
        switch (charge) {
        case -1:
            if (bosum == 0) { return 0 }
            break;
        case 0:
            switch (bosum) {
            case 0, 1: return 1;
            case 2, 3: return 3;
            case 4, 5: return 5;
            case 6, 7: return 7;
            default: break
            }
            break;
        case 1:
            switch (bosum) {
            case 0, 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            default: break
            }
            break;
        case 2:
            switch (bosum) {
            case 0, 1, 2, 3: return 3;
            case 4, 5: return 5;
            default: break
            }
            break;
        default: break
        }
        break;
    case 86:
        if (charge == 0) {
            switch (bosum) {
            case 0: return 0;
            case 1, 2: return 2;
            case 3, 4: return 4;
            case 5, 6: return 6;
            case 7, 8: return 8;
            default: break
            }
        }
        break;
    case 87:
        switch (charge) {
        case 0:
            if (bosum <= 1) { return 1 }
            break;
        case 1:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    case 88:
        switch (charge) {
        case 0:
            if (bosum <= 2) { return 2 }
            break;
        case 1:
            if (bosum <= 1) { return 1 }
            break;
        case 2:
            if (bosum == 0) { return 0 }
            break;
        default: break
        }
        break;
    default: break
    }
    
    return bosum
}

public func MKAtomAssignTypicalImplicitHydrogens(_ atom: MKAtom) {
    let bosum = atom.getExplicitValence()
    let valence = getTypicalValence(UInt(atom.getAtomicNum()), bosum, atom.getFormalCharge())
    atom.setImplicitHCount(valence - bosum)
}

func MKBondGetSmallestRingSize(_ bond: MKBond, _ bound: Int) -> Int {
    // A bounded BFS from one side of a bond trying to find the other side.
    //   - The required queue is implemented using a std::vector for efficiency
    //   - Note that items are never popped from the queue, I just move a pointer
    //     along to mark the left hand side
    //   - A potential improvement would be to BFS from both sides at the same time
    if !bond.isInRing() { return 0 }

    let start = bond.getBeginAtom()
    let end = bond.getEndAtom()
    var qatoms: [MKAtom] = []
    let numatoms = bond.getParent()!.numAtoms()
    let seen: Bitset = Bitset()
    seen.add(start.getIdx())
    for nbond in start.getBondIterator()! {
        if nbond == bond { continue }
        if !nbond.isInRing() { continue }
        let nbr = nbond.getNbrAtom(start)
        qatoms.append(nbr)
    }
    var depthmarker = qatoms.count // when qstart reaches this, increment the depth
    var qstart = 0
    var depth = 2
    while qatoms.count - qstart > 0 { // While queue not empty
        let curr = qatoms[qstart]
        if qstart == depthmarker {
            depth += 1
            depthmarker = qatoms.count
        }
        qstart += 1
        if seen.contains(curr.getIdx()) {
            continue
        }
        seen.add(curr.getIdx())
        if depth < bound {
            for nbond in curr.getBondIterator()! {
                if !nbond.isInRing() { continue }
                let nbr = nbond.getNbrAtom(curr)
                if nbr == end {
                    return depth + 1
                }
                if !seen.contains(nbr.getIdx()) {
                    qatoms.append(nbr)
                }
            }
        }
    }
    return 0
}
