

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include "Convert.h"
#include "TSIFile.h"
#include "QFile.h"
#include "WritevtuFiles.h"
#include "VMDOutput.h"

Convert::Convert()
{
}
Convert::Convert(std::vector <std::string> argument)
{

    m_Argument = argument;
    m_Healthy =true;
    m_Seed =36723;
    m_MinFaceAngle = -0.5;
    m_MinVerticesDistanceSquare = 1.0;
    m_MaxLinkLengthSquare = 3.0;
    m_OutputFilename = "out.q";
    m_InputFilename = "in.tsi";
    m_Box(0) = 10;
    m_Box(1) = 10;
    m_Box(2) = 10;
    m_Zoom(0) = 1;
    m_Zoom(1) = 1;
    m_Zoom(2) = 1;
    m_center = false;
    m_CutType = "line";
    m_Translate(0)=0; m_Translate(1)=0; m_Translate(2)=0;
    m_NoVF = 0;
    Nfunction f;
    ExploreArguments();     // read the input data

    std::string ext = m_InputFilename.substr(m_InputFilename.find_last_of(".") + 1);
    QFile Q;
    TSIFile tsi;
    WritevtuFiles vtu;
    VMDOutput  gro;
    if(ext==TSExt)
    {
        m_BluePrint = Q.Read(m_InputFilename);
    }
    else if(ext==TSIExt)
    {
        m_BluePrint = tsi.ReadTSI(m_InputFilename);
    }
    else
    {
        std::cout<<"---> Error: input file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }
    //==============================================================
    
    // make the changes into the mesh
    
        std::vector<Vertex_Map> bvertex = m_BluePrint.bvertex;
        m_Box = m_BluePrint.simbox;
        ExploreArguments();     // read the input data
    
    // first we do the translation
    double min_x = 10000000;
    double max_x = -1000000;
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
    {
        double x = it->y ;

        
        if(x>max_x)
            max_x = x;
        
        if(x<min_x)
            min_x = x;

    }
    double slice = 0;
    std::cout<<" enter the size of the bond slice \n";
    std::cin>>slice;
    
    if(slice<0 || slice>1)
        std::cout<<"error slice size \n";
    
    
    slice = (max_x-min_x)*slice;
    
    std::vector<Vertex_Map> Pole1;
    std::vector<Vertex_Map> Pole2;
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
    {
        double x = it->y ;

        

        if(x>(max_x-slice)){
            Pole1.push_back((*it));
        }
        if(x<(min_x+slice)){
            Pole2.push_back((*it));
        }

    }
    
    std::ofstream ndx;
    ndx.open("ndx.inx");
    int tot_G = Pole1.size() + Pole2.size();
    ndx<<" Group1 "<<tot_G<<"\n";
    for (std::vector<Vertex_Map>::iterator it = Pole1.begin() ; it != Pole1.end(); ++it)
    {

        ndx<<it->id<<"  ";
    }
    for (std::vector<Vertex_Map>::iterator it = Pole2.begin() ; it != Pole2.end(); ++it)
    {

        ndx<<it->id<<"  ";
    }
    std::vector<Inclusion_Map> binclusion = m_BluePrint.binclusion;
    for (std::vector<Inclusion_Map>::iterator it = binclusion.begin() ; it != binclusion.end(); ++it)
    {
        bool isnot = true;
        for (std::vector<Vertex_Map>::iterator it2 = Pole1.begin() ; it2 != Pole1.end(); ++it2)
        {
            if(it2->id == it->vid){
                isnot =false;
            }

        }
        for (std::vector<Vertex_Map>::iterator it2 = Pole2.begin() ; it2 != Pole2.end(); ++it2)
        {
            if(it2->id == it->vid){
                isnot =false;
            }

        }
        if(!isnot){
            it->tid = 2;
        }

    }
    // we do box change
    m_BluePrint.binclusion = binclusion;
    m_BluePrint.simbox = m_Box;
    tsi.WriteTSI(m_OutputFilename, m_BluePrint);

}
Convert::~Convert()
{
    
}
void Convert::ExploreArguments()
{
    Nfunction f;
    for (long i=1;i<m_Argument.size();i=i+2)
    {
        std::string Arg1 = m_Argument.at(i);
        if(Arg1=="-h")
        {
            HelpMessage();
            exit(0);
            break;
        }
        else if(Arg1=="-o")
        {
            m_OutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-in")
        {
            m_InputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-zoom")
        {
            m_Zoom(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Zoom(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Zoom(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;i++;
        }
        else if(Arg1=="-defout")
        {
            m_GeneralOutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-angle")
        {
            m_MinFaceAngle = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-maxDist")
        {
            m_MaxLinkLengthSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-type")
        {
            m_CutType = m_Argument.at(i+1);
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-translate")
        {
            m_Translate(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Translate(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Translate(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;i++;

        }
        else if(Arg1=="-Box")
        {
            m_Box(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Box(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Box(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;i++;
            
        }
        else if(Arg1=="-NoVF")
        {
            m_NoVF = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-center")
        {
            m_center = true;
            i=i-1;
            
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else
        {
            std::cout << "---> Error: "<<Arg1<<" ... bad argument ";
            std::cout<<"\n"<<"For more information and tips run "<< m_Argument.at(0) <<" -h"<<"\n";
            m_Healthy =false;
            exit(0);
            break;
        }
    }

}
void Convert::HelpMessage()
{
    std::cout<<"--------------------copyright: Weria Pezeshkian------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"--------------version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"------------simple example for exacuting  -------------------"<<"\n";
    std::cout<<" ./CNV -in out.q  -o out.tsi"<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  option    type        default            description "<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  -in         string       out.q               input file name, could be tsi or q file formats "<<"\n";
    std::cout<<"  -o          string       out.tsi             output file name, could be tsi/vtu/gro or q file formats "<<"\n";
    std::cout<<"  -type       string       direction           vertex/direction "<<"\n";

    
    
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"------------------ version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
}
