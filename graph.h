#pragma once
#include <list>
//#include <Base\GeoMesh.h>
#include "Heap.h"
#include "Align.h"
// #ifdef NIL
// #define NIL -1
// #endif
#define  NIL -1
#ifndef Inf
#define Inf 2147483647
#endif

#ifndef DInf
#define DInf 1e+16
#endif
//using namespace Base::Geometry;


using namespace std;
// graph 
class adjanode{
public:
	int index;
	double dist;
	adjanode():index(0),dist(0) {}
	friend bool operator>=(const adjanode& le,const adjanode& ri);
	friend bool operator>(const adjanode& le,const adjanode& ri);
	friend bool operator<=(const adjanode& le,const adjanode& ri);
	friend bool operator<(const adjanode& le,const adjanode& ri);
};

class adjalist{
public:
	int vernum;
    list<adjanode>* adjanodelist;
	int* component;
	bool add_undirected_edge(int orin,int dest);
	bool add_directed_edge(int orin,int dest);
	adjalist():vernum(0),adjanodelist(NULL),component(NULL) {}
	adjalist(int _ver);
	adjalist(const adjalist& other);
	adjalist& operator=(const adjalist& other);    
	//traversal
	//DFS
	void DFS(const int v,int tick = 0);
	int DFSComponent();
	//BFS
    void BFS(const int v,int tick = 0);
	int BFSComponent();
	bool clear();
	~adjalist();
};

class GraPath
{
public:
	int vernum;
	int* pathlist;
	double* distlist;
	int source;
	GraPath():vernum(0),pathlist(NULL),source(0),distlist(NULL) {}
    bool GetGraPath(int tar,list<int>& grapath);
	GraPath(int _vernum)
	{
		this->vernum=_vernum;
		pathlist=new int[_vernum];
		memset(pathlist,NIL,_vernum*sizeof(int));
		distlist=new double[_vernum];
		for (int i=0;i<_vernum;i++)
		{
			distlist[i]=DInf;
		}
	}
	GraPath(const GraPath& other);
	GraPath& operator=(const GraPath& other);
	~GraPath();
};
class DijNode
{
public:
	int verindex;
	double dist;
	int times;
	DijNode():verindex(NIL),dist(0),times(0) {}
	DijNode(int _verindex,double _dist, int _times) {verindex=_verindex;dist=_dist;times=_times;}
	friend bool operator<=(const DijNode& le,const DijNode& ri);
	friend bool operator<(const DijNode& le,const DijNode& ri);
	friend bool operator>=(const DijNode& le,const DijNode& ri);
	friend bool operator>(const DijNode& le,const DijNode& ri);
};

int GetHalfEdge(DTriMesh& hemesh,int ver,int nextver);

bool ConvertHETriMeshToAdjalist(DTriMesh* hemesh,adjalist& adlist);

bool Dijkstra(adjalist& meshadja,int source,GraPath& meshpath,int tar=NIL);

struct  verset
{
	int index;
	bool S;
};
bool Dijkstralist(adjalist& meshadja,int source, GraPath& meshpath, vector<verset>& tarlist);
bool Vertex2Edge(DTriMesh& hetrimesh,list<int>& boundaryver,list<int>& boundaryedge);
//get voronoi sample points
//seed is the location in the candidatepoints
bool VoronoiSampling(adjalist& meshadja, int seed, vector<int>& candidatepoints,int samplingsize, vector<int>& samplingpoints);

