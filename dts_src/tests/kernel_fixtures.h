#pragma once
#include "check.h"

// MC_Simulation voxelizes at the start of do_Simulation, which tests do not run
static State* make_voxelized_state(const char* input, const char* index = nullptr) {
    State* pState = make_test_state(input, "top.top", 1, index);
    pState->GetVoxelization()->Voxelize(pState->GetMesh()->GetActiveV());
    return pState;
}

// recompute all geometry from positions so the reference does not share the kernel's incremental updates
static double reference_energy(State* pState) {
    pState->GetCurvatureCalculator()->Initialize();
    return pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
}

static std::vector<links*> all_links(State* pState) {
    MESH* mesh = pState->GetMesh();
    std::vector<links*> l = mesh->GetRightL();
    l.insert(l.end(), mesh->GetLeftL().begin(), mesh->GetLeftL().end());
    l.insert(l.end(), mesh->GetEdgeL().begin(), mesh->GetEdgeL().end());
    return l;
}
