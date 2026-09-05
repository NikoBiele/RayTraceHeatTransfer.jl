"""
    elementFaces(domain)

The subfaces of `domain` in global element order: element `k` is row and column
`k` of the exchange factor matrix, and entry `k` of the solution vectors.

This ordering is the single definition used by the solver write-back and by the
plotting extension. Anything that indexes elements by position should obtain
that position from here rather than re-deriving it.
"""
elementFaces(domain) = [sf for f in domain.facesMesh for sf in f.subFaces]

"""
    elementParents(domain)

Superface index of each element, in the same order as [`elementFaces`](@ref).
"""
elementParents(domain) = [i for (i, f) in enumerate(domain.facesMesh) for _ in f.subFaces]