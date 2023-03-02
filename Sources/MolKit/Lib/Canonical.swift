







func isFerroceneBond(_ bond: MKBond) -> Bool {
    if bond.getBondOrder() != 1 {
        return false
    }

    var Fe: MKAtom? = nil
    var C: MKAtom? = nil

    let begin: MKAtom = bond.getBeginAtom()
    if begin.getAtomicNum() == 26 {
        Fe = begin
    }
    if begin.getAtomicNum() == 6 {
        C = begin
    }

    let end: MKAtom = bond.getEndAtom()
    if end.getAtomicNum() == 26 {
        Fe = end
    }
    if end.getAtomicNum() == 6 {
        C = end
    }

    if Fe == nil || C == nil {
        return false
    }

    if Fe!.getExplicitDegree() < 10 {
        return false
    }

    return C!.hasDoubleBond() && C!.isInRing()
}