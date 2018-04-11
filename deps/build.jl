using BinaryProvider

const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

PREPROCESSOR_VERSION = "38fc19f9d33be968f29d8ab73d0633f4ea1f9c33"
REMOTE_PATH = "https://dynare.adjemian.eu/preprocessor/$PREPROCESSOR_VERSION"

products = Product[
    ExecutableProduct(prefix, "dynare_m", :dynare),
]

download_info = Dict(
    Linux(:i686, :glibc)    => ("$REMOTE_PATH/linux/32/preprocessor.tar.gz", "f2256d88b02f3375636e67c36d657b089e9c7d3329d9eb192eaa8475a38d66e4"),
    Linux(:x86_64, :glibc)  => ("$REMOTE_PATH/linux/64/preprocessor.tar.gz", "e9c04c3b1ee693e755431f57801d06cfdae9a24b11cfb02d6b1699eedd51aef4"),
    MacOS()                 => ("$REMOTE_PATH/osx/64/preprocessor.tar.gz", "c4ad5eb94e95fdcd5c8193f099d7cbd53ecda9b044ec5ceec9506dbf21002f98"),
    Windows(:i686)          => ("$REMOTE_PATH/windows/32/preprocessor.tar.gz", "70946d5f2c11a3d2012eb3edb96d65fffd26e1da65519b244a77f4c8a928cf78"),
    Windows(:x86_64)        => ("$REMOTE_PATH/windows/64/preprocessor.tar.gz", "fe3f0643265d895159edb186b4ee74b766b87e544364e9e3a137a51584aa46b5"),
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
