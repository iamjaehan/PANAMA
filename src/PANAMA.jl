module PANAMA

greet() = print("Hello World!")

include("..//devel/TOS.jl")
using .TOS
include("..//devel/CTOP.jl")
using .CTOP
include("..//devel/CTB.jl")
using .CTB

end # module PANAMA
