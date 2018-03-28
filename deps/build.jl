using BinaryProvider

const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

PREPROCESSOR_VERSION = "b52ded0eb429b0f24eb7f27cdb1b19835a232bb8"
REMOTE_PATH = "https://dynare.adjemian.eu/preprocessor/$PREPROCESSOR_VERSION"

products = Product[
    ExecutableProduct(prefix, "dynare_m", :dynare),
]

download_info = Dict(
    Linux(:i686, :glibc)    => ("$REMOTE_PATH/linux/32/preprocessor.tar.gz", "861dee9c0c7bbacf30d5553de8110933fa4ac56cd028ff15a7a4f2927912f6c1"),
    Linux(:x86_64, :glibc)  => ("$REMOTE_PATH/linux/64/preprocessor.tar.gz", "e155e7dd8530d9ac80cadf8e41c0ebbfa6297d1b0437bbdc33ddb2867c97348d"),
    MacOS()                 => ("$REMOTE_PATH/osx/64/preprocessor.tar.gz", "8829a829ef5c01a50e5aefd07ca9a0e126cea71b7f15f0e914c3e90bb17acdb3"),
    Windows(:i686)          => ("$REMOTE_PATH/windows/32/preprocessor.tar.gz", "7285f9d95f45eb753dc7f3b614fc85dc9dff7584717973f3c28e967fb3a27d4f"),
    Windows(:x86_64)        => ("$REMOTE_PATH/windows/64/preprocessor.tar.gz", "a05cdfc6388ba4df17a495763fbee6ea15de71de5d01f59a31ba2aa2a2f7c032"),
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
