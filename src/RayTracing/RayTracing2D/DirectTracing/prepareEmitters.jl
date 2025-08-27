function prepare_emitters(rtm::RayTracingMeshOptim{VPF,VVPF,MT,VT,DIII,DII,BOO,FLOA}, gas::GasProperties) where {VPF,VVPF,MT,VT,DIII,DII,BOO,FLOA}
    emitters = Emitter[]
    total_energy = 0.0
    
    for (coarse_index, coarse_face) in enumerate(rtm.coarse_mesh)
        for (fine_index, fine_face) in enumerate(rtm.fine_mesh[coarse_index])
            # Volume emission
            volume_energy = 4 * gas.kappa * STEFAN_BOLTZMANN * fine_face.volume * fine_face.T_in_g^4
            push!(emitters, Emitter(:volume, coarse_index, fine_index, 0, volume_energy))
            total_energy += volume_energy
            
            # Surface emission
            for wall_index in 1:length(fine_face.solidWalls)
                if fine_face.solidWalls[wall_index]
                    surface_energy = fine_face.epsilon[wall_index] * STEFAN_BOLTZMANN * fine_face.area[wall_index] * fine_face.T_in_w[wall_index]^4
                    push!(emitters, Emitter(:surface, coarse_index, fine_index, wall_index, surface_energy))
                    total_energy += surface_energy
                end
            end
        end
    end
    
    # Normalize emitter energies
    for emitter in emitters
        emitter.energy /= total_energy
    end
    
    return emitters, total_energy
end