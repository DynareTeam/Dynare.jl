using BinaryProvider

const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

PREPROCESSOR_VERSION = "46a2272ef929a7acbf276fea2c43089bfc7d0e34"
REMOTE_PATH = "https://dynare.adjemian.eu/preprocessor/$PREPROCESSOR_VERSION"

products = Product[
    ExecutableProduct(prefix, "dynare_m", :dynare),
]

download_info = Dict(
    Linux(:i686, :glibc)    => ("$REMOTE_PATH/linux/32/preprocessor.tar.gz", "3e847abfd259e3d21a28a315da322f3b326910e8f1d198efdc0bd200eaea00ec"),
    Linux(:x86_64, :glibc)  => ("$REMOTE_PATH/linux/64/preprocessor.tar.gz", "bd31db8a0b29f43ba63d8df7a96cde5349851bd72c562200dbfa9d37f5606219"),
    MacOS()                 => ("$REMOTE_PATH/osx/64/preprocessor.tar.gz", "001fa0ab957623b1e5b01b8b7563b92f960adda79027f054d1f72df8ad60c30b"),
    Windows(:i686)          => ("$REMOTE_PATH/windows/32/preprocessor.tar.gz", "1be8d8fc8cc798dec811f36525ff5b5fd6bf78e84be659068df30d4cf737eaac"),
    Windows(:x86_64)        => ("$REMOTE_PATH/windows/64/preprocessor.tar.gz", "32485a1dde052e97b7844d731af1f1c0c47c5154f22e29b71789ffe9a1f25e8c"),
)

for p in products
    if platform_key() in keys(download_info)
        url, tarball_hash = download_info[platform_key()]
        install(url, tarball_hash; prefix=prefix, force=true, verbose=true, ignore_platform=true)
    else
        error("Your platform $(Sys.MACHINE) is not supported by Dynare.jl!")
    end
    write_deps_file(joinpath(@__DIR__, "deps.jl"), products; verbose=true)
end
