using BinaryProvider

const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

PREPROCESSOR_VERSION = "70ec52a6cce2ad116b0de68eec505a4092e50a82"
REMOTE_PATH = "https://dynare.adjemian.eu/preprocessor/$PREPROCESSOR_VERSION"

products = Product[
    ExecutableProduct(prefix, "dynare_m", :dynare),
]

download_info = Dict(
    Linux(:i686, :glibc)    => ("$REMOTE_PATH/linux/32/preprocessor.tar.gz", "dcfc2a15c35a6ee916cd0b20cec89a454082ad3aa512ea2dfd62ea2a22e58561"),
    Linux(:x86_64, :glibc)  => ("$REMOTE_PATH/linux/64/preprocessor.tar.gz", "1f9c619c2116542c0810ae2ba8d12525993d03a46d53322e49b3fad35f6c51e7"),
    MacOS()                 => ("$REMOTE_PATH/osx/64/preprocessor.tar.gz", "4d082dbbee46810b7350d9594f1c1d848e37c871f2fa7456be60f2847682d3bf"),
    Windows(:i686)          => ("$REMOTE_PATH/windows/32/preprocessor.tar.gz", "e277c9dc9248ff57476612e42d7fd691ffd72cee1de4aea4926544906cc772e3"),
    Windows(:x86_64)        => ("$REMOTE_PATH/windows/64/preprocessor.tar.gz", "5a8cbe1f58664446c28e074b2e5b79aa23c64ccd84c7b0425d2003a29292af64"),
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
