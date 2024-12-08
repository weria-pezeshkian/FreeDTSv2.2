#include <stdio.h>
#include "VertexHarmonicBounds.h"
#include "MESH.h"

VertexHarmonicBounds::VertexHarmonicBounds() {
}

VertexHarmonicBounds::~VertexHarmonicBounds() {

}
void VertexHarmonicBounds::AddBondToList(bond *b){
    
    m_VertexBond.push_back(b);
    return;
}
