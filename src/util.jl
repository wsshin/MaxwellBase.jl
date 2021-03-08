# import Base:isapprox, dot

# See http://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array.
# Generate the swap macros at http://pages.ripco.net/~jgamble/nw.html.
# - Number of inputs: length of the input vector
# - Algorithm choices: "Best"
# - Select "Create a set of SWAP macros".
@inline function swap_ind!(ind::AbsArrInteger, v::AbsArr, i::Integer, j::Integer)
    @inbounds if v[ind[j]] < v[ind[i]]
        @inbounds ind[i], ind[j] = ind[j], ind[i]
    end
end

@inline function sort_ind!(ind::MArray{Tuple{2,2,2},<:Integer,3,8}, v::MArray{Tuple{2,2,2},<:Real,3,8})
    # Initialize ind.
    @simd for n = 1:8
        @inbounds ind[n] = n
    end

    # Sort ind.
    swap_ind!(ind,v,1,2); swap_ind!(ind,v,3,4); swap_ind!(ind,v,1,3); swap_ind!(ind,v,2,4);
    swap_ind!(ind,v,2,3); swap_ind!(ind,v,5,6); swap_ind!(ind,v,7,8); swap_ind!(ind,v,5,7);
    swap_ind!(ind,v,6,8); swap_ind!(ind,v,6,7); swap_ind!(ind,v,1,5); swap_ind!(ind,v,2,6);
    swap_ind!(ind,v,2,5); swap_ind!(ind,v,3,7); swap_ind!(ind,v,4,8); swap_ind!(ind,v,4,7);
    swap_ind!(ind,v,3,5); swap_ind!(ind,v,4,6); swap_ind!(ind,v,4,5)

    return ind
end

@inline function sort_ind!(ind::MArray{Tuple{2,2},<:Integer,2,4}, v::MArray{Tuple{2,2},<:Real,2,4})
    # Initialize ind.
    @simd for n = 1:4
        @inbounds ind[n] = n
    end

    # Sort ind.
    swap_ind!(ind,v,1,2); swap_ind!(ind,v,3,4); swap_ind!(ind,v,1,3); swap_ind!(ind,v,2,4); swap_ind!(ind,v,2,3);

    return ind
end

@inline function sort_ind!(ind::MArray{Tuple{2},<:Integer,1,2}, v::MArray{Tuple{2},<:Real,1,2})
    # Initialize ind.
    @simd for n = 1:2
        @inbounds ind[n] = n
    end

    # Sort ind.
    swap_ind!(ind,v,1,2);

    return ind
end

isuniform(mv::MArray{Tuple{2,2,2},<:Any,3,8}) = @inbounds (mv[1]==mv[2]==mv[3]==mv[4]==mv[5]==mv[6]==mv[7]==mv[8])
isuniform(mv::MArray{Tuple{2,2},<:Any,2,4}) = @inbounds (mv[1]==mv[2]==mv[3]==mv[4])
isuniform(mv::MArray{Tuple{2},<:Any,1,2}) = @inbounds (mv[1]==mv[2])
