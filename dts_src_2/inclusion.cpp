
#include <stdio.h>
#include "inclusion.h"
#include "vertex.h"

inclusion::inclusion(int id, const InclusionType& inctype) : InclusionType(inctype), m_ID(id) {

}
inclusion::~inclusion() {
    
}
void inclusion::Updatevertex(vertex * v) {
    m_pvertex = v;
    return;
}
void inclusion::UpdateLocalDirection(const Vec3D & lo_dir) {
    m_LDirection = lo_dir;
    return;
}
void inclusion::UpdateGlobalDirection(const Vec3D & lg_dir) {
    m_GDirection = lg_dir;
    return;
}
bool inclusion::UpdateGlobalDirectionFromLocal(){
    
    m_GDirection = (m_pvertex->GetL2GTransferMatrix())*m_LDirection;
    
    if(m_GDirection.isbad())
        return false;

    return true;
}
bool inclusion::UpdateLocalDirectionFromGlobal(){
    
    m_LDirection = (m_pvertex->GetG2LTransferMatrix())*m_GDirection;
    if(m_LDirection.isbad())
        return false;

    return true;
}


