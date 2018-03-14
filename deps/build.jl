using BinaryProvider

const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

PREPROCESSOR_VERSION = "5727083865753f5abde8bdc0c20eee2b1ed5a501"
REMOTE_PATH = "https://dynare.adjemian.eu/preprocessor/$PREPROCESSOR_VERSION"

products = Product[
    ExecutableProduct(prefix, "dynare_m", :dynare),
]

download_info = Dict(
    Linux(:i686, :glibc)    => ("$REMOTE_PATH/linux/32/preprocessor.tar.gz", "291427358c08e23f7fb7816e5329892695b761c7d7403c9799d13bf6d58fa1a3"),
    Linux(:x86_64, :glibc)  => ("$REMOTE_PATH/linux/64/preprocessor.tar.gz", "b82876a77e819ad89dc7183ab32c0c802a0c2c44b33cc10ec03a647094fb7fd6"),
    MacOS()                 => ("$REMOTE_PATH/osx/64/preprocessor.tar.gz", "e88764e0641fb79bb0f8dbcf06902fc0b698b642a287800e2170a927999321b0"),
    Windows(:i686)          => ("$REMOTE_PATH/windows/32/preprocessor.tar.gz", "8834fbddf862a2546123e188cb1adafab3fa36a7df6950bb326eb743e55afd75"),
    Windows(:x86_64)        => ("$REMOTE_PATH/windows/64/preprocessor.tar.gz", "cfd2eb6bda76c50dd4d2fb789377897859a98925c69158d9adf39dc2deba5cbf"),
)

if any(!satisfied(p; verbose=verbose) for p in products)
    if platform_key() in keys(download_info)
        url, tarball_hash = download_info[platform_key()]
        install(url, tarball_hash; prefix=prefix, force=true, verbose=true, ignore_platform=true)
    else
        error("Your platform $(Sys.MACHINE) is not supported by Dynare.jl!")
    end
    write_deps_file(joinpath(@__DIR__, "deps.jl"), products; verbose=true)
end
