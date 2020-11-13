#include "Align.h"
#include "graph.h"
#include <queue>
#include <assert.h>


adjalist::adjalist(const adjalist &other)
{
	 //copy construction function
	this->vernum=other.vernum;
	adjanodelist=new list<adjanode>[other.vernum];
	memcpy(this->adjanodelist,other.adjanodelist,(this->vernum)*sizeof(list<adjanode>));
}

bool operator>=(const adjanode& le,const adjanode& ri)
{
	return (le.dist>=ri.dist);
}
bool operator>(const adjanode& le,const adjanode& ri)
{
	return le.dist>ri.dist;
}
bool operator<=(const adjanode& le,const adjanode& ri)
{
	return le.dist<=ri.dist;
}
bool operator<(const adjanode& le,const adjanode& ri)
{
	return le.dist<ri.dist;
}
adjalist& adjalist::operator =(const adjalist &other)
{
	if (this==&other)
	{
		return (*this);	
	}
	if (adjanodelist!=NULL)
	{
		delete [] adjanodelist;
	}
	this->vernum=other.vernum;
	adjanodelist=new list<adjanode>[other.vernum];
	memcpy(this->adjanodelist,other.adjanodelist,(this->vernum)*sizeof(list<adjanode>));
	return (*this);
}
adjalist::adjalist(int _ver)
{
	vernum = _ver;
	adjanodelist = new list<adjanode> [vernum];
	this->component = NULL;
}
bool adjalist::add_directed_edge(int orin,int dest)
{
	if (orin<0 || orin>=vernum)
	{
		return false;
	}
	if (dest<0 || dest>=vernum)
	{
		return false;
	}
	adjanode _dest;
	_dest.index = dest;
	adjanodelist[orin].push_back(_dest);
	return true;
}

bool adjalist::add_undirected_edge(int orin,int dest)
{
	add_directed_edge(orin,dest);
	return add_directed_edge(dest,orin);
}

void adjalist::DFS(const int v, int tick/* = 0*/)
{
	component[v] = tick;
	for (list<adjanode>::iterator iter = adjanodelist[v].begin();iter!=adjanodelist[v].end();iter++)
	{
		if (component[iter->index]<0)
		{
			DFS(iter->index,tick);
		}
	}
}

int adjalist::DFSComponent()
{
	if (component!=NULL)
	{
		delete[] component;
		component = NULL;
	}
	component = new int[vernum];
	memset(component,-1,vernum*sizeof(int));
	int tick = 0;
	for (int i=0;i<vernum;i++)
	{
		if (component[i]<0)
		{
			this->DFS(i,tick);
			tick++;
		}
	}
	return tick;
}
void adjalist::BFS(const int v,int tick /* = 0 */)
{
	component[v] = tick;
	std::queue<int> q;
	q.push(v);
	while (!q.empty())
	{
		int _v = q.front();
		q.pop();
		for (list<adjanode>::iterator iter=adjanodelist[_v].begin();iter!=adjanodelist[_v].end();iter++)
		{
			if (component[iter->index]<0)
			{
				component[iter->index] = tick;
				q.push(iter->index);
			}
		}
	}
}
int adjalist::BFSComponent()
{
	if (component!=NULL)
	{
		delete[] component;
		component = NULL;
	}
	component = new int[vernum];
	memset(component,-1,vernum*sizeof(int));
	int tick = 0;
	for (int i=0;i<vernum;i++)
	{
		if (component[i]<0)
		{
			this->BFS(i,tick);
			tick++;
		}
	}
	return tick;
}
bool adjalist::clear()
{
	vernum = 0;
	if (adjanodelist!=NULL)
	{
		delete[] adjanodelist;
		adjanodelist = NULL;
	}
	if (component!=NULL)
	{
		delete[] component;
		component = NULL;
	}
	return true;
}
adjalist::~adjalist()
{
	if (this->adjanodelist!=NULL)
	{
		delete[] adjanodelist;
		adjanodelist=NULL;
	}
	if (this->component!=NULL)
	{
		delete[] component;
		component = NULL;
	}
}

int GetHalfEdge(DTriMesh& hemesh, int ver, int nextver)
{
	return 0;
}

bool ConvertHETriMeshToAdjalist(DTriMesh* hemesh,adjalist& adlist)
{
	if (hemesh==NULL)
	{
		return false;
	}
	if (hemesh->n_vertices()==0)
	{
		return false;
	}
	if(adlist.vernum!=0)
	{
		return false;
	}
	adlist.vernum=hemesh->n_vertices();
	adlist.adjanodelist=new list<adjanode>[adlist.vernum];
	for (int i=0;i<adlist.vernum;i++)
	{
		vector<int> adjacentvertices;
		DTriMesh::VertexHandle vi(i);
		for (DTriMesh::VertexVertexIter vj = hemesh->vv_iter(vi); vj.is_valid(); vj++)
		{
			adjacentvertices.push_back(vj->idx());
		}
		//CollectVertexNeighborhoodVertrices(*hemesh,i,adjacentvertices);
		for (vector<int>::iterator iter = adjacentvertices.begin();iter!=adjacentvertices.end();iter++)
		{
			adjanode temp;
			temp.index = (*iter);
			//temp.dist = hemesh->GetVertex(i)->GetCood().Distance(hemesh->GetVertex(*iter)->GetCood());
			temp.dist = (OtoE(hemesh->point(DTriMesh::VertexHandle(i))) - OtoE(hemesh->point(DTriMesh::VertexHandle(*iter)))).norm();
			adlist.adjanodelist[i].push_back(temp);
		}
	}
	return true;
}
 GraPath::GraPath(const GraPath& other)
 {
	vernum=other.vernum;
	pathlist=new int[other.vernum];
	distlist=new double[other.vernum];
	memcpy(pathlist,other.pathlist,(other.vernum)*sizeof(int));
	memcpy(distlist,other.distlist,(other.vernum)*sizeof(double));
 };
GraPath& GraPath::operator=(const GraPath& other)
{
   //检查自赋值
   if (this==&other)
   {
	   return (*this);
   }
   //清空原有内存
   if (pathlist!=NULL)
   {
	   delete[] pathlist;
	   pathlist=NULL;
   }
   if (distlist!=NULL)
   {
	   delete[] distlist;
	   distlist=NULL;
   }
   //赋值
   this->vernum=other.vernum;
   pathlist=new int[other.vernum];
   memcpy(pathlist,other.pathlist,(other.vernum)*sizeof(int));
   //
   distlist=new double[other.vernum];
   memcpy(distlist,other.distlist,(other.vernum)*sizeof(double));
   return *this;
}
GraPath::~GraPath()
{
	if (pathlist!=NULL)
	{
		delete[] pathlist;
		pathlist=NULL;
	}
	if (distlist!=NULL)
	{
		delete[] distlist;
		distlist=NULL;
	}
}
bool operator<=(const DijNode& le,const DijNode& ri)
{
	return le.dist<=ri.dist;
}
bool operator<(const DijNode& le,const DijNode& ri)
{
	return le.dist<ri.dist;
}
bool operator>=(const DijNode& le,const DijNode& ri)
{
	return le.dist>=ri.dist;
}
bool operator>(const DijNode& le,const DijNode& ri)
{
	return le.dist>ri.dist;
}
bool Dijkstra(adjalist& meshadja,int source,GraPath& meshpath,int tar/* =NIL */)
{
	 if (source<0 || source>meshadja.vernum)
	 {
		 return false;
	 }
	 double* dist=new double[meshadja.vernum];
	 int* previous=new int[meshadja.vernum];
	 int* vertimes=new int[meshadja.vernum];
	 //Initialization
	 for (int i=0;i<meshadja.vernum;i++)
	 {
		 dist[i]=DInf;
		 previous[i]=NIL;
		 vertimes[i]=0;
	 }
	 //
	 meshpath.vernum=meshadja.vernum;
	 meshpath.source=source;
	 meshpath.pathlist=new int[meshadja.vernum];
	 memset(meshpath.pathlist,NIL,meshadja.vernum*sizeof(int));
	 //
	 if (meshpath.distlist!=NULL)
	 {
		 delete[] meshpath.distlist;
	 }
	 //initialize
	 meshpath.distlist=new double[meshadja.vernum];
	 for (int i=0;i<meshadja.vernum;i++)
	 {
		 meshpath.distlist[i]=DInf;
	 }
	 //from source to source
	 dist[source]=0;
//	 vertimes[source]=0;
	 Min_Heap<DijNode> Qset;
	 DijNode tempnode(source,0,vertimes[source]);
	 Qset.Insert(tempnode);
	 while (!Qset.IsEmpty())
	 {
		DijNode u=Qset.Pop();
		if (u.times==vertimes[u.verindex])
		{
			  if (u.dist==DInf)
			  {
				  break;
			  }
			  if (u.verindex==tar)
			  {
				  break;
			  }
			  for (list<adjanode>::iterator iter=meshadja.adjanodelist[u.verindex].begin();iter!=meshadja.adjanodelist[u.verindex].end();iter++)
			  {
				   double alt=dist[u.verindex]+iter->dist;
				   if (alt<dist[iter->index])
				   {
					   dist[iter->index]=alt;
					   meshpath.pathlist[iter->index]=u.verindex;
					   vertimes[iter->index]++;
					   DijNode temp;
					   temp.dist=alt;
					   temp.verindex=iter->index;
					   temp.times=vertimes[iter->index];
					   Qset.Insert(temp);
				   }
			  }
		}
	 } 
	 //
	 memcpy(meshpath.distlist,dist,meshadja.vernum*sizeof(double));
	 //
// 	 for (int i=0;i<meshadja.vernum;i++)
// 	 {
// 	  	 meshpath.distlist[i]=dist[i];
// 	 }
	 // 
	 delete[] vertimes;
	 delete[] previous;
	 delete[] dist;
	 return true;
}
bool Dijkstralist(adjalist& meshadja,int source, GraPath& meshpath, vector<verset>& tarlist)
{
	if (source<0 || source>meshadja.vernum)
	{
		return false;
	}
	//
	int* tmp = new int[meshadja.vernum];
	memset(tmp,1,meshadja.vernum*sizeof(int));
	//
	double* dist=new double[meshadja.vernum];
	int* previous=new int[meshadja.vernum];
	int* vertimes=new int[meshadja.vernum];
	//Initialization
	for (int i=0;i<meshadja.vernum;i++)
	{
		dist[i]=DInf;
		previous[i]=NIL;
		vertimes[i]=0;
	}
	//
	meshpath.vernum=meshadja.vernum;
	meshpath.source=source;
	meshpath.pathlist=new int[meshadja.vernum];
	memset(meshpath.pathlist,NIL,meshadja.vernum*sizeof(int));
	if (meshpath.distlist!=NULL)
	{
		delete[] meshpath.distlist;
	}
	//initialize
	meshpath.distlist=new double[meshadja.vernum];
	for (int i=0;i<meshadja.vernum;i++)
	{
		meshpath.distlist[i]=DInf;
	}
	//from source to source
	dist[source]=0;
	tmp[source] = 0;
	//	 vertimes[source]=0;
	Min_Heap<DijNode> Qset;
	DijNode tempnode(source,0,vertimes[source]);
	Qset.Insert(tempnode);
	int tarcount = 0;
	while (!Qset.IsEmpty())
	{
		DijNode u=Qset.Pop();
		if (u.times==vertimes[u.verindex])
		{
			if (u.dist==DInf)
			{
				break;
			}
			tmp[u.verindex] = 0;
// 			if (u.verindex==tar)
// 			{
// 				break;
// 			}
			tarcount=0;
			for (vector<verset>::iterator iter = tarlist.begin();iter!=tarlist.end();iter++)
			{
				if (!(iter->S))
				{
					tarcount += tmp[iter->index];
				}
			}
			if (tarcount == 0)
			{
				break;
			}
			for (list<adjanode>::iterator iter=meshadja.adjanodelist[u.verindex].begin();iter!=meshadja.adjanodelist[u.verindex].end();iter++)
			{
				double alt=dist[u.verindex]+iter->dist;
				if (alt<dist[iter->index])
				{
					dist[iter->index]=alt;
					meshpath.pathlist[iter->index]=u.verindex;
					vertimes[iter->index]++;
					DijNode temp;
					temp.dist=alt;
					temp.verindex=iter->index;
					temp.times=vertimes[iter->index];
					Qset.Insert(temp);
				}
			}
		}
	} 
	//
	memcpy(meshpath.distlist,dist,meshadja.vernum*sizeof(double));
	//
	// 	 for (int i=0;i<meshadja.vernum;i++)
	// 	 {
	// 	  	 meshpath.distlist[i]=dist[i];
	// 	 }
	// 
	delete[] tmp;tmp=NULL;
	delete[] vertimes;
	delete[] previous;
	delete[] dist;
	return true;
}
bool GraPath::GetGraPath(int tar, std::list<int> &grapath)
{
	int u=tar;
	if (u==NIL)
	{
		return false;
	}
	while (pathlist[u]!=NIL)
	{
		grapath.push_front(u);
		u=pathlist[u];
	}
	grapath.push_front(u);
	return true;
}
//int GetHalfEdge(DTriMesh& hemesh, int ver, int nextver)
//{
//	if (ver<0 || ver>hemesh.n_vertices())
//	{
//		return -1;
//	}
//	if (nextver<0 || nextver>hemesh.n_vertices())
//	{
//		return -1;
//	}
//	DTriMesh* temphemesh = &hemesh;
//	//ver---->nextver
//	int outedgeindex = temphemesh->point(DTriMesh::VertexHandle(ver))->GetOneOutGoingHalfedge();
//	int tick = outedgeindex;
//	bool boundary = false;
//	do
//	{
//		int prev = temphemesh->GetHalfedge(tick)->GetPrevHalfedge();
//		tick = temphemesh->GetHalfedge(prev)->GetFlipHalfedge();
//		if (tick < 0)
//		{
//			boundary = true;
//			break;
//		}
//		if (temphemesh->GetHalfedge(tick)->GetDestVertex() == nextver)
//		{
//			return tick;
//		}
//	} while (tick != outedgeindex);
//	if (boundary)
//	{
//		tick = outedgeindex;
//		do
//		{
//			int opp = temphemesh->GetHalfedge(tick)->GetFlipHalfedge();
//			if (opp < 0)
//			{
//				break;
//			}
//			tick = temphemesh->GetHalfedge(opp)->GetNextHalfedge();
//			if (temphemesh->GetHalfedge(tick)->GetDestVertex() == nextver)
//			{
//				return tick;
//			}
//		} while (tick != outedgeindex);
//	}
//	//nextver-->ver
//	outedgeindex = temphemesh->GetHEVertex(nextver)->GetOneOutGoingHalfedge();
//	tick = outedgeindex;
//	boundary = false;
//	do
//	{
//		int prev = temphemesh->GetHalfedge(tick)->GetPrevHalfedge();
//		tick = temphemesh->GetHalfedge(prev)->GetFlipHalfedge();
//		if (tick < 0)
//		{
//			boundary = true;
//			break;
//		}
//		if (temphemesh->GetHalfedge(tick)->GetDestVertex() == ver)
//		{
//			return tick;
//		}
//	} while (tick != outedgeindex);
//	if (boundary)
//	{
//		tick = outedgeindex;
//		do
//		{
//			int opp = temphemesh->GetHalfedge(tick)->GetFlipHalfedge();
//			if (opp < 0)
//			{
//				break;
//			}
//			tick = temphemesh->GetHalfedge(opp)->GetNextHalfedge();
//			if (temphemesh->GetHalfedge(tick)->GetDestVertex() == ver)
//			{
//				return tick;
//			}
//		} while (tick != outedgeindex);
//	}
//	return -1;
//}
bool Vertex2Edge(DTriMesh& hemesh,list<int>& boundaryver,list<int>& boundaryedge)
{
	 if (hemesh.n_vertices()==0)
	 {
		 return false;
	 }
	 if (!boundaryedge.empty())
	 {
		 boundaryedge.clear();
	 }
	 if (boundaryver.empty())
	 {
		 return false;
	 }
	 list<int>::iterator iter=boundaryver.begin();
	 list<int>::iterator jiter=iter;
	 iter++;
	 while(iter!=boundaryver.end())
	 {
		 int heedge=GetHalfEdge(hemesh,(*jiter),(*iter));
		 if (heedge!=(-1))
		 {
			 boundaryedge.push_back(heedge);
		 }
		 jiter=iter;
		 iter++;
	 }
	 return true;
}

bool VoronoiSampling(adjalist& meshadja, int seed, vector<int>& candidatepoints,int samplingsize, vector<int>& samplingpoints)
{
	 //seed 是相对坐标
	 //candidatepoints 是绝对坐标
	if (candidatepoints.empty())
	 {
		 return false;
	 }
	 if (seed<0)
	 {
		 seed = 0;
	 }
	 //create sampling sets
	 //if S[i] is true the i is belong to the sampling set
	 bool* S=new bool[candidatepoints.size()];
	 int* samplingindex=new int[samplingsize];
	 memset(samplingindex,-1,samplingsize*sizeof(int));
	 memset(S,false,candidatepoints.size()*sizeof(bool));
	 double* dis_matrix = NULL;
	 dis_matrix = new double[candidatepoints.size()];
	 vector<verset> canditmp;
	 for (unsigned int i=0;i<candidatepoints.size(); i++)
	 {
		 verset vertmp;
		 vertmp.index = candidatepoints[i];
		 vertmp.S = false;
		 canditmp.push_back(vertmp);
		 dis_matrix[i] = DInf;
	 }
	 S[seed]=true;
	 canditmp[seed].S = true;
	 samplingindex[0] = seed;
	 for (int j=0; j<samplingsize;j++)
	 {
		  GraPath pathtmp;
		  Dijkstralist(meshadja,candidatepoints[samplingindex[j]],pathtmp,canditmp);
		  for (int i=0;i<candidatepoints.size();i++)
		  {
			  if (!S[i])
			  {
				  if (pathtmp.distlist[candidatepoints[i]]<dis_matrix[i])
				  {
					  dis_matrix[i] = pathtmp.distlist[candidatepoints[i]];
				  }
			  }
		  }
		  double furthestdis=0;
		  int furthestindex = -1;
		  for (int i=0;i<candidatepoints.size();i++)
		  {
			  if (!S[i])
			  {
				  if (furthestdis<dis_matrix[i])
				  {
					  furthestdis = dis_matrix[i];
					  furthestindex = i;
				  }
			  }
		  }
		  //
		  assert(furthestindex!=(-1));
		  S[furthestindex] = true;
		  canditmp[furthestindex].S = true;
		  if (j+1<samplingsize)
		  {
			  samplingindex[j+1] = furthestindex;
		  }
		  else
		  {
			  break;
		  }     
	 }
	 if (!samplingpoints.empty())
	 {
		 samplingpoints.clear();
	 }
	 samplingpoints.resize(samplingsize);
	 for (int i=0;i<samplingsize;i++)
	 {
		 samplingpoints[i] = candidatepoints[samplingindex[i]];
	 }

	 delete[] S;S=NULL;
	 delete[] samplingindex;samplingindex=NULL;
	 delete[] dis_matrix;dis_matrix = NULL;
	 return true;
}