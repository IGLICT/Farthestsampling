#ifndef new_DALIGN_H
#define new_DALIGN_H

#include <Eigen/Eigen>
#include <vector>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <vector>
#include <set>

#define PI 3.141592653589793
#define Eps 1e-10

//extern int flag;//main 1 .cpp里的int flag=0;这里注释掉，在Align。cpp中改为intflag=0；


struct DTraits : public OpenMesh::DefaultTraits {
	typedef OpenMesh::Vec3d Point; // use double-values points
	typedef OpenMesh::Vec3d Normal; // use double-values points
	VertexAttributes(OpenMesh::Attributes::Normal | OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Normal | OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	typedef OpenMesh::Vec3d MidPoint; // relative with the edges
	//static OpenMesh::Vec3d calc_mid_point(OpenMesh::EdgeHandle & eh)
	//{
	//	OpenMesh::HalfedgeHandle he_it = Meshori.halfedge_handle(e_it, 0);
	//	OpenMesh::HalfedgeHandle oe_it = Meshori.halfedge_handle(e_it, 1);

	//	int num1 = Meshori.from_vertex_handle(he_it).idx();
	//	int num2 = Meshori.to_vertex_handle(he_it).idx();
	//}

};

typedef OpenMesh::TriMesh_ArrayKernelT<DTraits> DTriMesh;

namespace Base
{
	namespace Geometry  // mesh related 
	{
		// definition of triangle mesh structure
		struct MyTraits : OpenMesh::DefaultTraits
		{
			// let point and normal be a vector of doubles
			typedef OpenMesh::Vec3d Point;
			typedef OpenMesh::Vec3d Normal;

			// add face label, e.g. for mesh segmentation 
			FaceTraits
			{
				unsigned int label;
				FaceT() : label(0) {}   // construction
			};

			// add vertex label, e.g. for mesh segmentation
			VertexTraits
			{
				unsigned int label;
				VertexT() : label(0) {}  // construction
			};

			// add vertex normal, color and status
			VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color);

			// add face normal, color and status
			FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color);

			// add edge status
			EdgeAttributes(OpenMesh::Attributes::Status);
		};

		// double precision triangle mesh
		typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyTriMesh;

		// double precision polygonal mesh (mainly for quad mesh)
		typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> MyPolyMesh;
	}
}


namespace Base
{
	namespace Geometry
	{
		//////////////////////////////////////////////////////////////////////////
		// basic definitions related

		// definition of vector based on OpenMesh 
		typedef OpenMesh::Vec3d Vector3D;
		typedef OpenMesh::Vec3f Vector3F;
		typedef OpenMesh::Vec3i Vector3I;
		typedef OpenMesh::Vec3uc Vector3UC;

		typedef OpenMesh::Vec2d Vector2D;
		typedef OpenMesh::Vec2f Vector2F;
		typedef OpenMesh::Vec2i Vector2I;

	}
}

class AffineAlign {
public:
	std::vector<Eigen::Vector3d> p;
	Eigen::Matrix3d AtA;
	AffineAlign(std::vector<Eigen::Vector3d> &v);
	Eigen::Matrix3d calc(const std::vector<Eigen::Vector3d> &v);
	double residual(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v);
	double residualwithoutnormal(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v);
	bool ckt(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v);
};

class RotateAlign {
public:
	std::vector<Eigen::Vector3d> *p, *q;
	double res;

	RotateAlign(std::vector<Eigen::Vector3d>& v1, std::vector<Eigen::Vector3d>& v2);
	Eigen::Matrix3d calc();
	Eigen::Matrix3d svd1();
	double residual(const Eigen::Matrix3d &m);
	bool ckt(Eigen::Matrix3d m);

	static void AlignAtoB(DTriMesh &A, DTriMesh &B);
	static void AlignAtoBCenter(DTriMesh &A, DTriMesh &B);
};

double polarDec(const Eigen::Matrix3d &a, Eigen::Matrix3d &r, Eigen::Matrix3d &s);
double polarDec(const Eigen::Matrix3d &a, Eigen::Matrix3d &r);

OpenMesh::Vec3d EtoO(const Eigen::Vector3d &v);
Eigen::Vector3d OtoE(const OpenMesh::Vec3d &v);

class Rot
{
public:
	Eigen::Vector3d axis;
	double theta;
	double circlek;
	Eigen::Matrix3d logr;
	Eigen::Matrix3d r;
	Rot() { circlek = theta = 0; axis = Eigen::Vector3d::Zero(); }
	void ToLogR();
	double ToAngle();
};

Eigen::Matrix3d exp(Eigen::Matrix3d);
Eigen::Matrix3d log(Eigen::Matrix3d);

Rot logrot(Eigen::Matrix3d, Rot);
Rot logrot(Eigen::Matrix3d);
void logrot(Rot& jrot, Eigen::Matrix3d & jrotation);
#endif
