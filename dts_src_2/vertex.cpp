
#include <stdio.h>
#include "vertex.h"
#include "links.h"
#include "triangle.h"
#include "Nfunction.h"

vertex::vertex(){
    
// default constructor
}
vertex::vertex(int id, double x, double y, double z){

    m_X = x;
    m_Y = y;
    m_Z = z;
    m_ID = id;
    m_SimTimeStep=-1;
    m_Group = 0;
    m_OwnInclusion = false;
    m_GroupName = "system";
    //--- edge vertex
    m_Geodesic_Curvature = 0;
    m_Normal_Curvature = 0;
    m_VertexType = 0;
    m_VLength = 0;                       // length of the vertex
    m_Lambda = 0;                   // line tension
    m_PrincipalCurvature_1 = 0;
    m_PrincipalCurvature_2 = 0;
    
}
vertex::vertex(int id){

    m_ID=id;
    m_X=0;
    m_Y=0;
    m_Z=0;
    m_SimTimeStep=-1;
    m_Group = 0;
    m_OwnInclusion = false;
    m_GroupName = "system";
    
    
    m_Geodesic_Curvature = 0;
    m_Normal_Curvature = 0;
    m_VertexType = 0;
    m_VLength = 0;                       // length of the vertex
    m_Lambda = 0;                   // line tension
    m_PrincipalCurvature_1 = 0;
    m_PrincipalCurvature_2 = 0;
}
vertex::~vertex(){
    
    // destructor
}
void vertex::UpdateGroupName(std::string group_name){
    m_GroupName = group_name;
    
    return;
}
void vertex::UpdateBox(Vec3D *pbox){
    m_pBox=pbox;
    return;
}

void vertex::UpdateDomainID(int domain_id){
    m_DomainID = domain_id;
    return;
}
void vertex::UpdateVID(int i){
    m_ID = i;
    return;
}
void vertex::UpdateVXPos(double x) {
    double boxX = (*m_pBox)(0);

    // Adjust X position if it's outside the box bounds
    if (x >= boxX) {
        m_X = x - boxX;
    } else if (x < 0) {
        m_X = x + boxX;
    } else {
        m_X = x;
    }
}
void vertex::UpdateVYPos(double y) {
    double boxY = (*m_pBox)(1);

    // Adjust Y position if it's outside the box bounds
    if (y >= boxY) {
        m_Y = y - boxY;
    } else if (y < 0) {
        m_Y = y + boxY;
    } else {
        m_Y = y;
    }
}
void vertex::UpdateVZPos(double z) {
    double boxZ = (*m_pBox)(2);

    // Adjust Z position if it's outside the box bounds
    if (z >= boxZ) {
        m_Z = z - boxZ;
    } else if (z < 0) {
        m_Z = z + boxZ;
    } else {
        m_Z = z;
    }
}
void vertex::AddtoLinkList(links* l)
{
    m_VLinkList.push_back(l);
    return;
}
void vertex::UpdateP1Curvature(double p1_curvature){
    
    m_PrincipalCurvature_1 = p1_curvature;
    return;
}
void vertex::UpdateP2Curvature(double p2_curvature){
    
    m_PrincipalCurvature_2 = p2_curvature;
    return;
}
bool vertex::AddtoLinkListCarefully(links* l) {
    // Check if the link already exists in the list
    for (std::vector<links*>::iterator it = m_VLinkList.begin() ; it != m_VLinkList.end(); ++it){
        if (*it == l)
            return false;
    }
    // Link does not exist, add it to the list
    m_VLinkList.push_back(l);
    return true;
}
void vertex::RemoveFromLinkList(links* z) {
    
    m_VLinkList.erase(std::remove(m_VLinkList.begin(), m_VLinkList.end(), z), m_VLinkList.end());
    return;
}
bool vertex::RemoveFromLinkListCarefully(links* l) {
    std::vector<links*>::iterator it = std::find(m_VLinkList.begin(), m_VLinkList.end(), l);
    if (it != m_VLinkList.end()) {
        m_VLinkList.erase(it);
        return true; // Vertex found and removed
    }
    return false; // Vertex not found
}
void vertex::AddtoTraingleList(triangle* z)
{
m_VTraingleList.push_back(z);
}
bool vertex::AddtoTriangleListCarefully(triangle* t) {
    // Check if the link already exists in the list
    for (std::vector<triangle*>::iterator it = m_VTraingleList.begin() ; it != m_VTraingleList.end(); ++it){
        if (*it == t)
            return false;
    }
    // Link does not exist, add it to the list
    m_VTraingleList.push_back(t);
    return true;
}
void vertex::RemoveFromTraingleList(triangle* z)
{
m_VTraingleList.erase(std::remove(m_VTraingleList.begin(), m_VTraingleList.end(), z), m_VTraingleList.end());
}
void vertex::AddtoNeighbourVertex(vertex* z)
{
    m_VNeighbourVertex.push_back(z);
    return;
}
bool vertex::AddtoNeighbourVertexCarefully(vertex* v) {
    // Check if the link already exists in the list
    for (std::vector<vertex*>::iterator it = m_VNeighbourVertex.begin() ; it != m_VNeighbourVertex.end(); ++it){
        if (*it == v)
            return false;
    }
    // Link does not exist, add it to the list
    m_VNeighbourVertex.push_back(v);
    return true;
}
void vertex::RemoveFromNeighbourVertex(vertex* r_v){
    m_VNeighbourVertex.erase(std::remove(m_VNeighbourVertex.begin(), m_VNeighbourVertex.end(), r_v), m_VNeighbourVertex.end());
    return;
}
void vertex::UpdateInclusion(inclusion * inc){
    m_pInclusion=inc;
    return;
}
void vertex::UpdateOwnInclusion(bool own){
    m_OwnInclusion=own;
    return;
}
void vertex::UpdateGroup(int group_id){
    m_Group=group_id;
    return;
}
void vertex::UpdateEnergy(double en){
    m_Energy=en;
    return;
}
//////////// Functions related to Curvature energy and normal
void vertex::UpdateNormal_Area(Vec3D v,double a){
    m_Normal=v;
    m_Area=a;
    return;
}
void vertex::UpdateL2GTransferMatrix(Tensor2 v)
{
    m_T_Local_2_Global = v;
}
void vertex::UpdateG2LTransferMatrix(Tensor2 v)
{
    m_T_Global_2_Local = v;
}
void vertex::UpdateSimTimeStep(int v)
{
    m_SimTimeStep=v;
}
void vertex::UpdateVoxel(Voxel<vertex> * pVoxel){
    
    m_pVoxel = pVoxel;
    return;
}

bool vertex::SetCopy(){
    
    m_OldX = m_X;
    m_OldY = m_Y;
    m_OldZ = m_Z;
    m_OldArea = m_Area;
    m_OldNormal = m_Normal;
    m_OldEnergy = m_Energy;
    m_OldT_Local_2_Global = m_T_Local_2_Global;         //  Local to global transformation matrix
    m_OldT_Global_2_Local = m_T_Global_2_Local;        //  global to local transformation matrix
    m_OldpVoxel = m_pVoxel;
    m_OldGeodesic_Curvature = m_Geodesic_Curvature;          // Edge Vertex Curvature
    m_OldNormal_Curvature = m_Normal_Curvature;          // Edge Vertex Curvature
    m_OldVLength = m_VLength;                       // length of the vertex
    
    m_OldVertexType = m_VertexType;                   // 0 surface vertex; 1 edge vertex;
    m_OldVTraingleList = m_VTraingleList;
    m_OldVLinkList = m_VLinkList;
    m_OldVNeighbourVertex = m_VNeighbourVertex;
    
   // links * m_OldpEdgeLink;
    //links * m_OldpPrecedingEdgeLink;// preceding link at the edge
    //inclusion *m_OldpInclusion;                    // pointer to an inclusion that the vertex hold (could be empty)
    //bool m_OldOwnInclusion;                        // to check if the vertex own any inclusion
    //CNTCell * m_OldCNTCell;                        // a unitcell that the vertex belong to at any point of the simulation, it will be chnage during a simulation
    //m_OldGroup = m_Group;            // Id of a group that the vertex belong too
    //std::string m_OldGroupName;
    // m_OldLambda;                   // line tension

    return true;
}
bool vertex::Reverse2PreviousCopy(){  // reverse the edge to the value set at the time of MakeCopy()
 
    m_X = m_OldX ;
    m_Y = m_OldY ;
    m_Z = m_OldZ ;
    m_Area = m_OldArea ;
    m_Normal = m_OldNormal ;
    m_Energy = m_OldEnergy ;
    m_T_Local_2_Global = m_OldT_Local_2_Global ;         //  Local to global transformation matrix
    m_T_Global_2_Local = m_OldT_Global_2_Local ;        //  global to local transformation matrix
    m_pVoxel = m_OldpVoxel ;
    m_Geodesic_Curvature = m_OldGeodesic_Curvature ;          // Edge Vertex Curvature
    m_Normal_Curvature = m_OldNormal_Curvature ;          // Edge Vertex Curvature
    m_VLength = m_OldVLength ;                       // length of the vertex
    m_VertexType = m_OldVertexType ;                   // 0 surface vertex; 1 edge vertex;
    m_VertexType = m_OldVertexType;                   // 0 surface vertex; 1 edge vertex;
    m_VTraingleList = m_OldVTraingleList;
    m_VLinkList = m_OldVLinkList ;
    m_VNeighbourVertex = m_OldVNeighbourVertex;
    
    return true;
}

bool vertex::CheckVoxel(){
    
//-- obtain the object (vertex)  cell id, with respect to the current cell
        int i = int((m_X)/m_pVoxel->GetXSideVoxel())-m_pVoxel->GetXIndex();
        int j = int((m_Y)/m_pVoxel->GetYSideVoxel())-m_pVoxel->GetYIndex();
        int k = int((m_Z)/m_pVoxel->GetZSideVoxel())-m_pVoxel->GetZIndex();
        //-- check if it has moved too far
        if(i!=0 || j!=0 || k!=0) {
            return false;
        }
    
    return true;
}
bool vertex::UpdateVoxelAfterAVertexMove(){
    
    int i = int(m_X/m_pVoxel->GetXSideVoxel())-m_pVoxel->GetXIndex();
    int j = int(m_Y/m_pVoxel->GetYSideVoxel())-m_pVoxel->GetYIndex();
    int k = int(m_Z/m_pVoxel->GetZSideVoxel())-m_pVoxel->GetZIndex();
    //-- check if it has moved too far
    if(i==0 && j==0 && k==0){
        return true;
    }
    else if(i > 1 || i < -1 || j > 1 || j < -1 || k > 1 || k < -1) {
        std::cout << " ---> error. the object has moved more than one voxel " << std::endl;
        return false;
    }
    //Voxel<vertex>::GetANeighbourCell
    m_pVoxel = m_pVoxel->GetANeighbourCell(i, j, k);

    return true;
}
double vertex::SquareDistanceFromAVertex(vertex* pv2) {
    double dx = pv2->GetVXPos() - m_X;
    double dy = pv2->GetVYPos() - m_Y;
    double dz = pv2->GetVZPos() - m_Z;

    double boxHalfX = (*m_pBox)(0) / 2.0;
    double boxHalfY = (*m_pBox)(1) / 2.0;
    double boxHalfZ = (*m_pBox)(2) / 2.0;

    // Adjust coordinates if outside the periodic boundary
    if (fabs(dx) > boxHalfX) {
        dx = (dx < 0) ? (*m_pBox)(0) + dx : dx - (*m_pBox)(0);
    }
    if (fabs(dy) > boxHalfY) {
        dy = (dy < 0) ? (*m_pBox)(1) + dy : dy - (*m_pBox)(1);
    }
    if (fabs(dz) > boxHalfZ) {
        dz = (dz < 0) ? (*m_pBox)(2) + dz : dz - (*m_pBox)(2);
    }

    // Compute and return squared distance
    return dx * dx + dy * dy + dz * dz;
}
double vertex::SquareDistanceOfAVertexFromAPoint(double X, double Y, double Z, vertex* pv2) {
    double dx = pv2->GetVXPos() - X;
    double dy = pv2->GetVYPos() - Y;
    double dz = pv2->GetVZPos() - Z;

    double boxHalfX = (*m_pBox)(0) / 2.0;
    double boxHalfY = (*m_pBox)(1) / 2.0;
    double boxHalfZ = (*m_pBox)(2) / 2.0;

    // Adjust coordinates if outside the periodic boundary
    if (fabs(dx) > boxHalfX) {
        dx = (dx < 0) ? (*m_pBox)(0) + dx : dx - (*m_pBox)(0);
    }
    if (fabs(dy) > boxHalfY) {
        dy = (dy < 0) ? (*m_pBox)(1) + dy : dy - (*m_pBox)(1);
    }
    if (fabs(dz) > boxHalfZ) {
        dz = (dz < 0) ? (*m_pBox)(2) + dz : dz - (*m_pBox)(2);
    }

    // Compute and return squared distance
    return dx * dx + dy * dy + dz * dz;
}
