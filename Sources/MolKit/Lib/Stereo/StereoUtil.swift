

enum NeighborSymmetryClasses: Int {
    // Tetrahedral
    case T1234 = 1234 // 4 different symmetry classes
    case T1123 = 1123 // 3 different symmetry classes, 1 class duplicated (2 times)
    case T1122 = 1122 // 2 different symmetry classes, 1 class duplicated (3 times)
    case T1112 = 1112 // 2 different symmetry classes, each class duplicated (2 times)
    case T1111 = 1111 // 1 symmetry class, duplictaed 4 times
    // CisTrans
    case C12     = 12 // 2 different symmetry classes
    case C11     = 11 // the same symmetry class
  }
  
