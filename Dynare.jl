module Dynare
##
 # Copyright (C) 2015 Dynare Team
 #
 # This file is part of Dynare.
 #
 # Dynare is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # Dynare is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
##

export dynare

function dynare(modfile)
    # Process modfile
    println(string("Using ", WORD_SIZE, "-bit preprocessor"))
    preprocessor = string(dirname(@__FILE__()), "/preprocessor", WORD_SIZE, "/dynare_m")
    run(`$preprocessor $modfile`)

    # Temporary: clean up Matlab output
    basename = split(modfile, ".mod", false)
    mfiles = filter(r".*\.m", readdir())
    for file in mfiles
        if isempty(search(file, ".mod"))
            rm(file)
        end
    end
    rm(basename[1], recursive=true)

    # Load module created by preprocessor
    require(basename[1])
end

end
