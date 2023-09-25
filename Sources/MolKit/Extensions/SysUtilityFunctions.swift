import Foundation

func sizeof <T> (_ : T.Type) -> Int
{
    return (MemoryLayout<T>.size)
}

func sizeof <T> (_ : T) -> Int
{
    return (MemoryLayout<T>.size)
}

func sizeof <T> (_ value : [T]) -> Int
{
    return (MemoryLayout<T>.size * value.count)
}

@inline(__always) func boolean_xor(_ lhs: Bool, _ rhs: Bool) -> Bool {
    return lhs != rhs
}

class StandardError: TextOutputStream {
  func write(_ string: String) {
    try! FileHandle.standardError.write(contentsOf: Data(string.utf8))
  }
}
