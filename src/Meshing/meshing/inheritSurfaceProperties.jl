# Updated inheritSurfaceProperties! function with direct inheritance
function inheritSurfaceProperties!(superFace::PolyVolume2D{G}, subFace::PolyVolume2D{G}; from=i, to=j) where {G}
    # the subFace inherits the properties of the superFace
    
    # wall state variables - direct inheritance
    inheritWallProperty!(superFace, subFace, :epsilon, from, to)
    inheritWallProperty!(superFace, subFace, :j_w, from, to)
    inheritWallProperty!(superFace, subFace, :g_a_w, from, to)
    inheritWallProperty!(superFace, subFace, :e_w, from, to)
    inheritWallProperty!(superFace, subFace, :r_w, from, to)
    inheritWallProperty!(superFace, subFace, :g_w, from, to)
    inheritWallProperty!(superFace, subFace, :i_w, from, to)
    inheritWallProperty!(superFace, subFace, :q_in_w, from, to)
    inheritWallProperty!(superFace, subFace, :q_w, from, to)
    inheritWallProperty!(superFace, subFace, :T_in_w, from, to)
    inheritWallProperty!(superFace, subFace, :T_w, from, to)
end