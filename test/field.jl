@testset "field" begin

N = [5,6,7]
ncmp = 4

# Test the size of create_field_array().
F = create_field_array(N, ncmp=ncmp)
@test size(F) == (N..., ncmp)

# Test if field_arr2vec() returns a view.
F .= rand(eltype(F), size(F))
f = field_arr2vec(F, order_cmpfirst=false)
f_cmpfirst = field_arr2vec(F, order_cmpfirst=true)

# Test if order_cmpfirst is working as intended.
@test reshape(transpose(reshape(f_cmpfirst, ncmp, :)), :) == f

end  # @testset "field"
