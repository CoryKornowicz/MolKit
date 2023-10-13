// swift-tools-version: 5.7
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "MolKit",
    defaultLocalization: "en",
    platforms: [
        .macOS(.v12),
        .iOS(.v16)
    ],
    products: [
        // Products define the executables and libraries a package produces, and make them visible to other packages.
        .library(
            name: "MolKit",
            targets: ["MolKit"])
    ],
    dependencies: [
        // Dependencies declare other packages that this package depends on.
        // .package(url: /* package url */, from: "1.0.0"),
        .package(url: "https://github.com/apple/swift-algorithms", from: "1.0.0"),
        .package(url: "https://github.com/apple/swift-collections", from: "1.0.4"),
        .package(url: "https://github.com/Jounce/Surge.git", .upToNextMajor(from: "2.3.2")),
        .package(url: "https://github.com/CoryKornowicz/SwiftBitset.git", branch: "master")
    ],
    targets: [
        // Targets are the basic building blocks of a package. A target can define a module or a test suite.
        // Targets can depend on other targets in this package, and on products in packages this package depends on.
        .target(
            name: "MolKit",
            dependencies: [.product(name: "Algorithms", package: "swift-algorithms"),
                           "Surge", 
                            .product(name: "Collections", package: "swift-collections"),
                           .product(name: "Bitset", package: "SwiftBitset")],
            resources: [
                    // Apply platform-specific rules.
                    // For example, images might be optimized per specific platform rule.
                    // If path is a directory, the rule is applied recursively.
                    // By default, a file will be copied if no rule applies.
                    // Process file in Sources/MolKit/Data/*
                    .copy("Data")
                  ]),
        .testTarget(
            name: "MolKitTests",
            dependencies: ["MolKit", .product(name: "Algorithms", package: "swift-algorithms"), "Surge", .product(name: "Collections", package: "swift-collections"),
                           .product(name: "Bitset", package: "SwiftBitset")],
            resources: [
                    // Copy Tests/ExampleTests/Resources directories as-is.
                    // Use to retain directory structure.
                    // Will be at top level in bundle.
                    .copy("Data/")
                  ])
    ]
)

