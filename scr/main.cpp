#include <string>
#include <direct.h>
//#include <afx.h>

#include "graph.h"
#include "Align.h"
#include <random>
#include <iostream>
using namespace std;

#define DEBUG
//#undef DEBUG

class FeatureSampling
{
public:
	//
	DTriMesh* hemesh;
	vector<int> samplingindex;
	//algorithm properties

	double neighbourradius;
	int numbersamplings;

	vector<Eigen::Vector3d> samplingpoints;
	vector<vector<int>> neighbourindex;
	//construction methods
	FeatureSampling() :neighbourradius(0.0) {}
	FeatureSampling(int _neighbourradius, DTriMesh* _hemesh, int NumberSamplings = 500);
	//FeatureSampling(FeatureSampling& other);
	//FeatureSampling& operator =(FeatureSampling& other);
	//~FeatureSampling();
	//methods
	bool SetHETriMesh(DTriMesh* _hemesh);
	bool SetSamplingIndex(const vector<int>& _samplingindex);
	bool Sampling();
	bool GetNeighbourhood();
	bool clear();
};

double calcarea(std::vector<Eigen::Vector3d> & tranglev)
{
	Eigen::Vector3d edge1 = tranglev[2] - tranglev[0];
	Eigen::Vector3d edge2 = tranglev[1] - tranglev[0];
	return 0.5*(edge1.cross(edge2)).norm();
}

FeatureSampling::FeatureSampling(int _neighbourradius, DTriMesh* _hemesh, int NumberSamplings/* = 500*/)
{
	neighbourradius = _neighbourradius;
	numbersamplings = NumberSamplings;
	//neighourindex.resize(numbersamplings);
	if (_hemesh == NULL)
	{
		cout << "error" << endl;
	}

	this->hemesh = _hemesh;
}

bool FeatureSampling::Sampling()
{
	if (hemesh == NULL)
	{
		return false;
	}
	if (hemesh->n_vertices() <= 0)
	{
		return false;
	}
	double* facearea = new double[hemesh->n_faces()];
	//initialize
	memset(facearea, 0, hemesh->n_faces() * sizeof(double));
	//


	for (int i = 0; i < hemesh->n_faces(); i++)
	{
		std::vector<Eigen::Vector3d> tri3;
		tri3.clear();
		//tri3.resize(3);
		double areatmp = 0;
		Base::Geometry::Vector3D normaltmp;
		DTriMesh::FaceHandle fh(i);
		for (auto it = hemesh->fv_begin(fh); it != hemesh->fv_end(fh); ++it)
		{
			auto vertices = *it;
			Eigen::Vector3d v1 = OtoE(hemesh->point(*it));
			tri3.push_back(v1);
		}
		//Vector3D v01 = hemesh->GetVertex(hemesh->GetFace(i)->GetReference(1))->GetCood()-hemesh->GetVertex(hemesh->GetFace(i)->GetReference(0))->GetCood();
		//Vector3D v02 = hemesh->GetVertex(hemesh->GetFace(i)->GetReference(2))->GetCood()-hemesh->GetVertex(hemesh->GetFace(i)->GetReference(0))->GetCood();
		//normaltmp = v01.cross(v02);
		// areatmp = normaltmp.norm()/2;
		facearea[i] = calcarea(tri3);
		//facearea[i] = areatmp;
		if (i > 0)
		{
			facearea[i] += facearea[i - 1];
		}
	}
	srand(time(NULL));
	std::random_device r;

	// Choose a random mean between 1 and 6
	std::default_random_engine e1(r());
	std::uniform_int_distribution<int> uniform_dist(0, RAND_MAX);
	for (int i = 0; i < numbersamplings; i++)
	{
		
		int randnum = uniform_dist(e1);
		double randratio = ((double)randnum / (double)RAND_MAX);
		double facerand = randratio*facearea[hemesh->n_faces() - 1];
		//find facerand
		int leftind = 0;
		int rightind = hemesh->n_faces() - 1;
		int debugtick = 0;
		while (leftind < rightind)
		{
			int mid = (leftind + rightind) / 2;
			if (facerand <= facearea[mid])
			{
				rightind = mid;
			}
			else
			{
				leftind = mid + 1;
			}
			debugtick++;
			assert(debugtick < hemesh->n_faces());
		}
		assert(leftind == rightind);
		//get the index
		double r1 = ((double)uniform_dist(e1)) / ((double)RAND_MAX);
		double r2 = ((double)uniform_dist(e1)) / ((double)RAND_MAX);
		//
		double para[3];
		para[0] = 1 - sqrt(r1);
		para[1] = (sqrt(r1))*(1 - r2);
		para[2] = (sqrt(r1))*(r2);
		Eigen::Vector3d P1(0.0, 0.0, 0.0);
		DTriMesh::FaceHandle fhh(leftind);
		//for (int k = 0; k<3;k++)
		//{
		int ii = 0;
		for (auto it = hemesh->fv_begin(fhh); it != hemesh->fv_end(fhh); ++it, ii++)
		{
			P1 += para[ii] * OtoE(hemesh->point(*it));
		}
		//P+=para[k]*hemesh->point(hemesh->GetFace(leftind)->GetReference(k))->GetCood();
		//}
		//Vector3D P(P1[0], P1[1], P1[2]);
		samplingpoints.push_back(P1.transpose());
	}
	delete[] facearea;
	return true;
}


void VoronoiSamplingyj(std::string orimesh, std::string voronoipointfile, int samplingsize, int seedid)
{
	//std::string orimesh = "F:\\yangjiee\\yangjie\\tracking\\paper\\22.obj";
	DTriMesh Meshori;
	//int samplingsize = 50;
	vector<int> samplingpoints;
	adjalist adja;
	vector<int> candidatepoints;
	if (!OpenMesh::IO::read_mesh(Meshori, orimesh.c_str()))
	{
		std::cout << "Load reference Mesh Error" << std::endl;
		return;
	}
	for (int i = 0; i < Meshori.n_vertices(); i++)
	{
		candidatepoints.push_back(i);
	}
	ConvertHETriMeshToAdjalist(&Meshori, adja);

	VoronoiSampling(adja, seedid, candidatepoints, samplingsize, samplingpoints);
	std::ofstream voronoi(voronoipointfile.c_str());
	voronoi << samplingsize << endl;

	for (int i = 0 ; i < samplingsize; i++)
	{
		voronoi << "1" << endl << samplingpoints[i] <<endl<< Meshori.point(DTriMesh::VertexHandle(samplingpoints[i]))<<endl;

	}

	voronoi.close();
	cout << "total sampling point size is " << samplingsize << ".\n" << "save file path: " << voronoipointfile.c_str() << endl;
}


void otherSamplingyj(std::string orimesh, std::string voronoipointfile, int samplingsize)
{
	//std::string orimesh = "F:\\yangjiee\\yangjie\\tracking\\paper\\22.obj";
	DTriMesh Meshori;
	//int samplingsize = 50;
	//vector<int> samplingpoints;
	vector<int> candidatepoints;
	if (!OpenMesh::IO::read_mesh(Meshori, orimesh.c_str()))
	{
		std::cout << "Load reference Mesh Error" << std::endl;
		return;
	}
	for (int i = 0; i < Meshori.n_vertices(); i++)
	{
		candidatepoints.push_back(i);
	}
	FeatureSampling samplinghaha(0, &Meshori, samplingsize);
	samplinghaha.Sampling();

	std::ofstream voronoi(voronoipointfile.c_str());
	voronoi << samplingsize << endl;

	for (int i = 0; i < samplingsize; i++)
	{
		voronoi << "1" << endl << i << endl << samplinghaha.samplingpoints[i][0]<<" "<< samplinghaha.samplingpoints[i][1] << " "<< samplinghaha.samplingpoints[i][2] << endl;
	}

	voronoi.close();
	cout << "total sampling point size is " << samplingsize << ".\n" << "save file path: " << voronoipointfile.c_str() << endl;
}


int main(int argc, char *argv[])
{
	//char * input = argv[1];
	//char * output = argv[2];
	//char * pointsize = argv[3];
	//char * seed = argv[4];

	//
	//otherSamplingyj(string(input), string(output), atoi(pointsize));
	//string mesh = "F:\\yangjiee\\yangjie\\tracking\\paper\\2.obj";
	//std::vector<LAs_> las;
	//LAs(mesh, las);
	//calc_las(argv[1]);
	//subdivede(string(input), string(output));
	//std::string inputfile =string(argv[1]);
	//std::string outputfile =string(argv[2]);
	//deleteface(inputfile,outputfile);

	//std::string a=string(argv[1]);
	//std::string b=string(argv[2]);
	//std::string c=string(argv[4]);
	//mkdir(c.c_str());
	//genmodel(a,b,atoi(argv[3]),c);

	//Smoother_1();
	//Smoother_2();
	//Smoother_3();
	//resettime();
	//char * inputfolder = "Y:\\yangjie17\\sig18\\testtvcg\\1";


	//VoronoiSamplingyj(string(input),string(output),atoi(pointsize),atoi(seed));

	//char * inputfolder = argv[1];
	//PreprocessObj2Mat(std::string(inputfolder));

	char * inputfolder = argv[1];
	char * outputfolder = argv[2];

	//testmodele_difference();

	//std::string path = "F:\\yangjiee\\sig17\\model_analysis\\people1";
	//std::string pathout = "F:\\yangjiee\\sig17\\model_analysis\\peoplem";
	//yj_checkmesh(path, pathout);
	//char * inifile = "F:\\yangjiee\\sig17\\asia\\sig170523\\testc3c5\\deform\\0.ini";
	//char * inifile = argv[1];
	//arapdeformm_inter(inifile);
	//arapdeformm_(inifile);
	//system("pause");
	return 0;
}