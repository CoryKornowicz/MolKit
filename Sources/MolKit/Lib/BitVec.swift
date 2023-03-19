
import Foundation
import Surge

// Use UInt32

let SETWORD: UInt32 = 32
let WORDROLL: UInt32 = 5
let WORDMASK: UInt32 = 31

let STARTWORDS: UInt32 = 10

let bitsoff: [UInt32] = [ 0xFFFFFFFF,0xFFFFFFFE,0xFFFFFFFC,0xFFFFFFF8,0xFFFFFFF0,0xFFFFFFE0,0xFFFFFFC0,
                          0xFFFFFF80,0xFFFFFF00,0xFFFFFE00,0xFFFFFC00,0xFFFFF800,0xFFFFF000,0xFFFFE000,
                          0xFFFFC000,0xFFFF8000,0xFFFF0000,0xFFFE0000,0xFFFC0000,0xFFF80000,0xFFF00000,
                          0xFFE00000,0xFFC00000,0xFF800000,0xFF000000,0xFE000000,0xFC000000,0xF8000000,
                          0xF0000000,0xE0000000,0xC0000000,0x80000000
                        ]

// Used by CountBits
let nibble_bit_count: [UInt16] = [  0, // 0000
                                    1, // 0001
                                    1, // 0010
                                    2, // 0011
                                    1, // 0100
                                    2, // 0101
                                    2, // 0110
                                    3, // 0111
                                    1, // 1000
                                    2, // 1001
                                    2, // 1010
                                    3, // 1011
                                    2, // 1100
                                    3, // 1101
                                    3, // 1110
                                    4  // 1111
                                        ]

func lowBit(_ set: inout UInt32, _ bit: inout Int) {
    if set != 0 {
        bit = 31
        if set != 0x80000000 {
            if (set & 0x0000ffff) != 0 { set &= 0x0000ffff; bit -= 16 }
            if (set & 0x00ff00ff) != 0 { set &= 0x00ff00ff; bit -= 8 }
            if (set & 0x0f0f0f0f) != 0 { set &= 0x0f0f0f0f; bit -= 4 }
            if (set & 0x33333333) != 0 { set &= 0x33333333; bit -= 2 }
            if (set & 0x55555555) != 0 { set &= 0x55555555; bit -= 1 }
        }
    } else { bit = -1 }
}


typealias WordVector = [UInt32]

func WORDSIZE_OF_BITSIZE(_ bit_size: UInt32) -> UInt32 {
    return (( bit_size >> WORDROLL ) + ((( bit_size & WORDMASK ) != 0) ? 1 : 0))
}

public class MKBitVec: Equatable {

    var _size: Int  { _set.count }  // Number of words stored
    var _set: WordVector    // Words

    public init() {
        self._set = WordVector(repeating: 0, count: Int(STARTWORDS))
    }

    public init(_ size_in_bits: UInt32) {
        self._set = WordVector(repeating: 0, count: Int(WORDSIZE_OF_BITSIZE(size_in_bits)))
    }
    
    public init(_ size_in_bits: Int) {
        self._set = WordVector(repeating: 0, count: Int(WORDSIZE_OF_BITSIZE(UInt32(size_in_bits))))
    }
    
    public init(_ vec: MKBitVec) {
        self._set = vec._set
    }

    func setBitOn(_ bit_offset: UInt32) {
        let word_offset = bit_offset >> WORDROLL
        let bit_offset_in_word = bit_offset & WORDMASK
        if word_offset >= _size {
            _ = self.resizeWords(word_offset + 1)
        }
        _set[Int(word_offset)] |= (1 << bit_offset_in_word)
    }
    
    func setBitOff(_ bit_offset: UInt32) {
        let word_offset = bit_offset >> WORDROLL
        let bit_offset_in_word = bit_offset & WORDMASK
        if word_offset < _size {
            _set[Int(word_offset)] &= ~(1 << bit_offset_in_word)
        }
    }
    
    func setBitOn(_ bit_offset: Int) {
        setBitOn(UInt32(bit_offset))
    }
    
    func setBitOff(_ bit_offset: Int) {
        setBitOff(UInt32(bit_offset))
    }
    
    func setBitOn(_ bit_offset: UInt) {
        setBitOn(UInt32(bit_offset))
    }
    

    /** Set the range of bits from \p lo_bit_offset to \p hi_bit_offset to 1
    Increases the size of this bit vector if necessary
    \param[in] lo_bit_offset a zero based offset into the bit vector
    \param[in] hi_bit_offset a zero based offset into the bit vector
    */
    func setRangeOn(_ lo_bit_offset: UInt32, _ hi_bit_offset: UInt32) {
        if (lo_bit_offset > hi_bit_offset) {
            return
        } else if (lo_bit_offset == hi_bit_offset) {
            setBitOn(hi_bit_offset)
            return
        } else {
            let lo_word_offset = lo_bit_offset >> WORDROLL
            let hi_word_offset = hi_bit_offset >> WORDROLL
            let lo_bit_offset_in_word = lo_bit_offset & WORDMASK
            let hi_bit_offset_in_word = hi_bit_offset & WORDMASK

            if hi_word_offset >= _size {
                _ = self.resizeWords(hi_word_offset + 1)
            }

            if lo_word_offset == hi_word_offset {
                for bit in lo_bit_offset_in_word...hi_bit_offset_in_word {
                    _set[Int(lo_word_offset)] |= (1 << bit)
                }
            } else {
                for bit in lo_bit_offset_in_word..<SETWORD {
                    _set[Int(lo_word_offset)] |= (1 << bit)
                }
                for word in lo_word_offset+1..<hi_word_offset {
                    _set[Int(word)] = ~0
                }
                for bit in 0...hi_bit_offset_in_word {
                    _set[Int(hi_word_offset)] |= (1 << bit)
                }
            }
        }
    }

    func setRangeOff(_ lo_bit_offset: UInt32, _ hi_bit_offset: UInt32) {
        if (lo_bit_offset > hi_bit_offset) {
            return
        } else if (lo_bit_offset == hi_bit_offset) {
            setBitOff(hi_bit_offset)
            return
        } else {
            let lo_word_offset = lo_bit_offset >> WORDROLL
            var hi_word_offset = hi_bit_offset >> WORDROLL
            let lo_bit_offset_in_word = lo_bit_offset & WORDMASK
            var hi_bit_offset_in_word = hi_bit_offset & WORDMASK
            
            if lo_word_offset >= _size {
                return
            }

            if hi_word_offset >= _size {
                hi_word_offset = UInt32(_size - 1)
                hi_bit_offset_in_word = SETWORD - 1
            }

            if lo_word_offset == hi_word_offset {
                for bit in lo_bit_offset_in_word...hi_bit_offset_in_word {
                    _set[Int(lo_word_offset)] &= ~(1 << bit)
                }
            } else {
                for bit in lo_bit_offset_in_word..<SETWORD {
                    _set[Int(lo_word_offset)] &= ~(1 << bit)
                }
                for word in lo_word_offset+1..<hi_word_offset {
                    _set[Int(word)] = 0
                }
                for bit in 0...hi_bit_offset_in_word {
                    _set[Int(hi_word_offset)] &= ~(1 << bit)
                }
            }
        }
    }

    func fold(_ new_bit_size: UInt32) {

        let new_word_size = new_bit_size >> WORDROLL
        if _size < new_word_size {
            _ = self.resizeWords(new_word_size)
            return
        }

        var idx: Int = Int(new_word_size)
        var i: Int = 0
        while idx < _size {
            _set[i] |= _set[idx]
            if i+1 < new_word_size {
                i += 1
            } else {
                i = 0
            }
            idx+=1
        }

        _ = self.resizeWords(new_word_size)
    }

    /** Searches the vector for the first true value, starting at the \p last_bit_offset 'th bit
    \param[in] last_bit_offset the bit before the first to consider
	\return the bit offset of the first true bit after \p last_bit_offset, or -1 if there is none
    */
    func nextBit(_ last_bit_offset: Int) -> Int {
        var last_bit_offset = last_bit_offset
        last_bit_offset += 1
        var wrdcnt = last_bit_offset >> WORDROLL
        if wrdcnt >= _size {
            return -1
        }
        if _set[wrdcnt] != 0 {
            var s = _set[wrdcnt] & bitsoff[last_bit_offset & Int(WORDMASK)]
            if s != 0 {
                var bit = 0
                lowBit(&s, &bit)
                if bit != -1 {
                    return Int(bit) + (wrdcnt << WORDROLL)
                }
            }
        }
        wrdcnt += 1
        while wrdcnt < _size {
            if _set[wrdcnt] != 0 {
                var s = _set[safelyAccess: wrdcnt]
                var bit = 0
                lowBit(&s, &bit)
                if bit != -1 {
                    return Int(bit) + (wrdcnt << WORDROLL)
                }
            }
            wrdcnt += 1
        }
        return -1
    }

    func countBits() -> Int {
        var count = 0
        for word in _set {
            var s = word
            while s != 0 {
                count+=Int(nibble_bit_count[Int(s) & 0xf])
                s >>= 4
            }
        }
        return count
    }

//    MARK: Could make this into a map call and reduce call?
    func isEmpty() -> Bool {
        for word in _set {
            if word != 0 {
                return false
            }
        }
        return true
    }
    
    func fromVecInt(_ bit_offsets: [Int]) {
        for bit in bit_offsets {
            self.setBitOn(UInt32(bit))
        }
    }
    
    /** Sets bits on, listed as a string of character-represented integers
     This bit vector is first cleared.
     The format is "[ n0 n1 n2 n3 ... ]".
     The square brackets are optional.
     The whitespace can be SPACE, NEWLINE or HTAB
     For example "[ 1 5 6 9 ]"
     \param[in] line A string containing positive integers
     \param[in] new_bit_size The size that the vector should become
     */
    func fromString(_ line: String, new_bit_size: UInt32) {
        self.clear()
        _ = self.resize(new_bit_size)
        
        let tokens = line.components(separatedBy: .whitespacesAndNewlines)
        for token in tokens {
            if token == "[" {
                continue
            } else if token == "]" {
                break
            }
            
            if let bit = Int(token) {
                if bit >= 0 {
                    self.setBitOn(UInt32(bit))
                } else {
                    print("Negative Bit: \(bit)")
                }
            }
        }
    }
    
    func toVecInt(_ bit_offset: inout [Int]) {
        for i in bit_offset {
            bit_offset[i] = 0
        }
        bit_offset.reserveCapacity(self.countBits())
        var i = -1
        repeat {
            bit_offset.insert(self.nextBit(i), at: bit_offset.count)
            i+=1
        } while self.nextBit(i) != -1
    }

    /// Find the first true bit at or after \p bit_offset
    /** Searches the vector for the first true value, starting at the \p bit_offset 'th bit
        \param[in] bit_offset the first bit to consider
        \return the bit offset of the first true bit at or after \p bit_offset, or -1 if there is none
    */
    func firstBit(_ bit_offset: Int = 0) -> Int {
        return self.bitIsSet(Int(bit_offset)) ? self.nextBit(Int(bit_offset - 1)) : self.nextBit(Int(bit_offset))
    }

    func endBit() -> Int {
        return -1
    }
    
    /// Inverts every bit in the vector
    /** Inverts the entire vector.
        Note that this may give unexpected results, as the vector
        can be considered to end in an arbitrary number of zero bits.
    */
    func negate() {
        _set.updateEach { word in
            word = ~word
        }
    }
    
//    Modifies wordVec in-place
    func getWords(_ vec: inout WordVector) {
        vec.append(contentsOf: self._set)
    }
    
    @discardableResult
    func resize(_ size_in_bits: UInt32) -> Bool {
        return self.resizeWords(WORDSIZE_OF_BITSIZE(size_in_bits))
    }

    func resizeWords(_ size_in_words: UInt32) -> Bool {
        
        if size_in_words <= _size {
            return false
        }
        // Enlarge Vector by filling in with 0's
        self._set.append(contentsOf: WordVector(repeating: 0, count: Int(size_in_words) - _size))
        return true
    }

    /// Asks if the \p bit_offset 'th bit is set
    /** Is the \p bit_offset 'th bit set ?
     \param[in] bit_offset a zero based offset into the bit vector
     \return true if it is set, false otherwise
     */
    func bitIsSet(_ bit_offset: Int) -> Bool {
        var rtn = false
        let word_offset = bit_offset >> WORDROLL
        if word_offset < _size {
            rtn = ((_set[Int(word_offset)] >> bit_offset) & 1) != 0
        }
        return rtn
    }
    
    func bitIsSet(_ bit_offset: MKBaseID) -> Bool {
        var rtn = false
        let word_offset = bit_offset.rawValue >> WORDROLL
        if word_offset < _size {
            rtn = ((_set[Int(word_offset)] >> bit_offset.rawValue) & 1) != 0
        }
        return rtn
    }
    
    func bitIsSet(_ bit_offset: Ref) -> Bool {
        guard bit_offset != .NoRef else { return false }
        guard bit_offset != .ImplicitRef else { return false }
        var rtn = false
        let word_offset = bit_offset.intValue! >> WORDROLL
        if word_offset < _size {
            rtn = ((_set[Int(word_offset)] >> bit_offset.intValue!) & 1) != 0
        }
        return rtn
    }

    func getSize() -> Int { return self._size }

    /** Set all the bits in this vector to zero
        Does not currently change the size of the vector.
    */
    func clear() {
        for i in 0..<_size {
            _set[i] = 0
        }
    }
    
//     Subscript Definition
    public subscript (_ bit_offset: Int) -> Bool {
        return self.bitIsSet(bit_offset)
    }
    
//     Operator Overloads
    
    
    static public func &= (lhs: inout MKBitVec, rhs: MKBitVec) {
        let min = lhs._size < rhs._size ? lhs._size : rhs._size
        for i in 0..<min {
            lhs._set[i] &= rhs._set[i]
        }
    }
    
    static public func |= (lhs: inout MKBitVec, rhs: MKBitVec) {
        if (lhs._size < rhs.getSize()) {
            _ = lhs.resizeWords(UInt32(rhs.getSize()))
        }
        
        for i in 0..<rhs.getSize() {
            lhs._set[i] |= rhs._set[i]
        }
    }
    
    static public func |= (lhs: inout MKBitVec, rhs: Int) {
        if lhs._size <= 0 {
            _ = lhs.resizeWords(1)
        }
        lhs._set[0] |= UInt32(rhs)
    }
    
    static public func ^= (lhs: inout MKBitVec, rhs: MKBitVec) {
        if (lhs._size < rhs.getSize()) {
            _ = lhs.resizeWords(UInt32(rhs.getSize()))
        }
        
        for i in 0..<rhs.getSize() {
            lhs._set[i] ^= rhs._set[i]
        }
    }
    
    public static func -= (lhs: inout MKBitVec, rhs: MKBitVec) {
        if (lhs._size < rhs.getSize()) {
            _ = lhs.resizeWords(UInt32(rhs.getSize()))
        }
        let temp: MKBitVec = lhs ^ rhs
        lhs &= temp
    }

    public static func += (lhs: inout MKBitVec, rhs: MKBitVec) {
        lhs._set.append(contentsOf: rhs._set)
    }

    public static func | (lhs: MKBitVec, rhs: MKBitVec) -> MKBitVec {
        var temp = lhs
        temp |= rhs
        return temp
    }
    
    public static func & (lhs: MKBitVec, rhs: MKBitVec) -> MKBitVec {
        var temp = lhs
        temp &= rhs
        return temp
    }

    public static func ^ (lhs: MKBitVec, rhs: MKBitVec) -> MKBitVec {
        var temp = lhs
        temp ^= rhs
        return temp
    }

    public static func - (lhs: inout MKBitVec, rhs: MKBitVec) -> MKBitVec {
        var temp = lhs ^ rhs
        temp &= lhs
        return temp
    }
    
    public static func == (lhs: MKBitVec, rhs: MKBitVec) -> Bool {
        var i = 0
        if lhs._size < rhs._size {
            for j in 0..<lhs._size {
                if lhs._set[j] != rhs._set[j] { return false }
                i+=1
            }
            for k in i..<rhs._size {
                if rhs._set[k] != 0 { return false }
            }
        } else {
            for j in 0..<rhs._size {
                if lhs._set[j] != rhs._set[j] { return false }
                i += 1
            }
            for k in i..<lhs._size {
                if lhs._set[k] != 0 { return false }
            }
        }
        return true
    }

    /** Return true if \p bv1 i less than \p bv2
     Lexicographical order, with bit vectors written LSB first.
     \param[in] bv1 A bit vector
     \param[in] bv2 Another bit vector
     \return true if equal, false otherwise
     */
    public static func < (lhs: MKBitVec, rhs: MKBitVec) -> Bool {
        var rtn: Bool = false
        var should_continue: Bool = true
        var nextBitLHS = lhs.nextBit(-1)
        var nextBitRHS = rhs.nextBit(-1)
        while should_continue {
            should_continue = false
            
            if nextBitLHS == -1 {
                rtn = (nextBitRHS == -1 ? false : true )
            } else if nextBitRHS == -1 {
                rtn = false
            } else if nextBitRHS < nextBitLHS {
                rtn = true
            } else if nextBitLHS < nextBitRHS {
                return false
            } else {
                nextBitLHS = lhs.nextBit(nextBitLHS)
                nextBitRHS = rhs.nextBit(nextBitRHS)
                should_continue = true
            }
        }
        return rtn
    }
    
    /// The Tanimoto coefficient, which may be regarded as the proportion of the "on-bits" which are shared.
    /// double Tanimoto(const OBBitVec & bv1, const OBBitVec & bv2);
    public static func tanimoto(lhs: MKBitVec, rhs: MKBitVec) -> Double {
        let andBits = Double((lhs & rhs).countBits())
        let orBits = Double((lhs | rhs).countBits())
        return (andBits/orBits)
    }

    public static func / (lhs: MKBitVec, rhs: MKBitVec) -> MKBitVec {
        let min = lhs._size < rhs._size ? lhs._size : rhs._size
        let temp = MKBitVec(UInt32(min))
        for i in 0..<min {
            temp._set[i] = lhs._set[i] & rhs._set[i]
        }
        return temp
    }
    
    public static func % (lhs: MKBitVec, rhs: MKBitVec) -> MKBitVec {
        let min = lhs._size < rhs._size ? lhs._size : rhs._size
        let temp = MKBitVec(UInt32(min))
        for i in 0..<min {
            temp._set[i] = lhs._set[i] % rhs._set[i]
        }
        return temp
    }
    
    public static func %= (lhs: inout MKBitVec, rhs: MKBitVec) {
        let min = lhs._size < rhs._size ? lhs._size : rhs._size
        for i in 0..<min {
            lhs._set[i] %= rhs._set[i]
        }
    }
    
    public static func * (lhs: MKBitVec, rhs: MKBitVec) -> MKBitVec {
        let min = lhs._size < rhs._size ? lhs._size : rhs._size
        let temp = MKBitVec(UInt32(min))
        for i in 0..<min {
            temp._set[i] = lhs._set[i] * rhs._set[i]
        }
        return temp
    }
        
    public static func + (lhs: MKBitVec, rhs: MKBitVec) -> MKBitVec {
        let min = lhs._size < rhs._size ? lhs._size : rhs._size
        let temp = MKBitVec(UInt32(min))
        for i in 0..<min {
            temp._set[i] = lhs._set[i] + rhs._set[i]
        }
        return temp
    }

}
