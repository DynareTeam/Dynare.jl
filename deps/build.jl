using BinaryProvider

const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

PREPROCESSOR_VERSION = "c4ae840b207dbc4e61bff315e8eaa28fb742d9ea"
REMOTE_PATH = "https://dynare.adjemian.eu/preprocessor/$PREPROCESSOR_VERSION"

products = Product[
    ExecutableProduct(prefix, "dynare_m", :dynare),
]

download_info = Dict(
    Linux(:i686, :glibc)    => ("$REMOTE_PATH/linux/32/preprocessor.tar.gz", "4971c8443a2a4347b90a002baeab4972857f9dd2863c6b9c646b6d0e399a6f7f"),
    Linux(:x86_64, :glibc)  => ("$REMOTE_PATH/linux/64/preprocessor.tar.gz", "f72b592b910283b141fb3ac0b7435b9dca1133bdaa710bf731150cbd45e77b08"),
    MacOS()                 => ("$REMOTE_PATH/osx/64/preprocessor.tar.gz", "d6b4dea63249996c45dc6249014f796359416e9e16019ab0c01a4c64ec5005a2"),
    Windows(:i686)          => ("$REMOTE_PATH/windows/32/preprocessor.tar.gz", "10c4aa274ace85420d7fbf75fe7239ae8e5accf3d479bbd8f431359650940b56"),
    Windows(:x86_64)        => ("$REMOTE_PATH/windows/64/preprocessor.tar.gz", "66805c5fb097f6972dadaf2f1560555e2c1e51fa4635750b5bfafb23a9ef54d7"),
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
