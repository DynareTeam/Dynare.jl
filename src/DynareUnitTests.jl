module DynareUnitTests

##
 # Copyright (C) 2018 Dynare Team
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

export @trycommand

function clean(basename::String)
    rm(basename; force=true, recursive=true)
    rm("$(basename).jl"; force=true)
    rm("$(basename)Dynamic.jl"; force=true)
    rm("$(basename)Static.jl"; force=true)
    rm("$(basename)SteadyState2.jl"; force=true)
    rm("$(basename)_set_auxiliary_variables.jl"; force=true)
end

end
