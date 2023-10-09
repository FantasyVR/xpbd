/*
	PBD实现，只包含距离约束和体积约束，求解器Guass-Seidel。
	作者：Peng Yu
	时间：2023/9/19
*/

#include <igl/opengl/glfw/Viewer.h>
#include <igl/AABB.h>
#include <igl/in_element.h>
#include <igl/signed_distance.h>
#include <igl/barycentric_coordinates.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include <set>
#include <chrono>

std::vector<int> pull_points{811, 812, 813, 814, 815, 809};

// creates small spheres to visualize P on the overlay of the mesh
// Input:
//  P:      #P by 3 coordinates of the centers of spheres
//  N:      #P by 3 normals (the south-north pole direction of the spheres).
//  r: radii of the spheres
//  sphereColors:      #P by 3 - RBG colors per sphere
//  res:    the resolution of the sphere discretization
// extendMesh if to extend the V,T,TC, or to overwrite them
// Output:
//  V:    #V by 3 sphere mesh coordinates
//  T     #T by 3 sphere mesh triangles
//  C:    #T by 3 vertex-based colors
IGL_INLINE bool point_spheres(const Eigen::MatrixXd &P,
							  const Eigen::MatrixXd &normals,
							  const double &r,
							  const Eigen::MatrixXd &sphereColors,
							  const int res,
							  Eigen::MatrixXd &V,
							  Eigen::MatrixXi &T,
							  Eigen::MatrixXd &C)
{
	using namespace Eigen;
	/*V.resize(res*res*P.rows(),3);
	T.resize(2*(res-1)*res*P.rows(),3);
	C.resize(V.rows(),3);*/

	MatrixXd VSphere(res * res, 3);
	MatrixXi TSphere(2 * (res - 1) * res, 3);

	// creating template sphere vertices
	for (int j = 0; j < res; j++)
	{
		double z = r * cos(igl::PI * (double)j / (double(res - 1)));
		for (int k = 0; k < res; k++)
		{
			double x = r * sin(igl::PI * (double)j / (double(res - 1))) * cos(2 * igl::PI * (double)k / (double(res)));
			double y = r * sin(igl::PI * (double)j / (double(res - 1))) * sin(2 * igl::PI * (double)k / (double(res)));
			VSphere.row(j * res + k) << x, y, z;
		}
	}

	for (int j = 0; j < res - 1; j++)
	{
		for (int k = 0; k < res; k++)
		{
			int v1 = j * res + k;
			int v2 = (j + 1) * res + k;
			int v3 = (j + 1) * res + (k + 1) % res;
			int v4 = j * res + (k + 1) % res;
			TSphere.row(2 * (res * j + k)) << v1, v2, v3;
			TSphere.row(2 * (res * j + k) + 1) << v4, v1, v3;
		}
	}

	// std::cout<<"TSphere: "<<TSphere<<std::endl;
	V.resize(VSphere.rows() * P.rows(), 3);
	T.resize(TSphere.rows() * P.rows(), 3);
	C.resize(V.rows(), 3);

	for (int i = 0; i < P.rows(); i++)
	{
		RowVector3d ZAxis = normals.row(i);
		ZAxis.normalize();
		RowVector3d XAxis;
		XAxis << 0.0, -normals(i, 2), normals(i, 1);
		if (XAxis.squaredNorm() < 1e-4)
			XAxis << -normals(i, 2), 0.0, normals(i, 0);
		XAxis.normalize();

		RowVector3d YAxis = ZAxis.cross(XAxis);
		YAxis.rowwise().normalize();

		Matrix3d R;
		R << XAxis, YAxis, ZAxis;
		RowVector3d translation = P.row(i);

		V.block(VSphere.rows() * i, 0, VSphere.rows(), 3) = VSphere * R + translation.replicate(VSphere.rows(), 1);
		T.block(TSphere.rows() * i, 0, TSphere.rows(), 3) = TSphere.array() + VSphere.rows() * i;
		C.block(VSphere.rows() * i, 0, VSphere.rows(), 3) = sphereColors.row(i).replicate(VSphere.rows(), 1);
	}

	return true;
}

Eigen::MatrixXd V_sphere;
Eigen::MatrixXi T_sphere;
Eigen::MatrixXd C_sphere;
// Eigen::MatrixXd P = (Eigen::MatrixXd(1, 3) << 3.9, -18.4, 4.9).finished();
Eigen::MatrixXd P = (Eigen::MatrixXd(1, 3) << 0, 0, 4.9).finished();

void init_spehre(Eigen::MatrixXd &V,
				 Eigen::MatrixXi &T,
				 Eigen::MatrixXd &C)
{
	Eigen::MatrixXd normals = (Eigen::MatrixXd(1, 3) << 0, 1, 0).finished();
	Eigen::MatrixXd sphereColors = (Eigen::MatrixXd(1, 3) << 1, 0, 0).finished();
	const int res = 10;
	bool is_success = point_spheres(P, normals, 1.0, sphereColors, res, V_sphere, T_sphere, C_sphere);
}

using Real = double;
using Vec3d = Eigen::RowVector3d;
using Vec4i = Eigen::RowVector4i;
using Vec3i = Eigen::RowVector3i;
using Vec2i = Eigen::RowVector2i;

class Softbody
{
public:
	Eigen::MatrixXd pos_;
	Eigen::MatrixXd old_pos_;
	Eigen::MatrixXd vel_;
	Eigen::VectorXd w_;
	Eigen::MatrixXi faces_;
	Eigen::MatrixXi edges_;
	Eigen::MatrixXi tets_;
	Eigen::VectorXd lambda_;
	Eigen::VectorXd rest_length_;
	Eigen::VectorXd rest_volume_;

	int num_vertices_{0};
	int num_edges_{0};
	int num_faces_{0};
	int num_tets_{0};
	Real dis_alpha_tilde_ {1.0e-5};
	Real vol_alpha_tilde_ {1.0e-5};
private:
	void read_node_file(const std::string node_file)
	{
		std::ifstream infile(node_file);
		std::string line;
		std::getline(infile, line);
		std::stringstream ss(line);
		int tmp1, tmp2, temp3;
		ss >> num_vertices_ >> tmp1 >> tmp2 >> temp3;
		pos_.resize(num_vertices_, 3);
		int i = 0;
		while (std::getline(infile, line))
		{
			if (line[0] == '#')
				continue;
			std::stringstream ss(line);
			int idx;
			Real a, b, c;
			ss >> idx >> a >> b >> c;
			pos_.row(i++) = Vec3d(a, b, c);
		}
		assert(pos_.rows() == num_vertices_);
		infile.close();
	}
	void read_ele_file(const std::string ele_file)
	{
		std::ifstream infile(ele_file);
		std::string line;
		std::getline(infile, line);
		std::stringstream ss(line);
		int tmp1, tmp2;
		ss >> num_tets_ >> tmp1 >> tmp2;
		tets_.resize(num_tets_, 4);
		int i = 0;
		while (std::getline(infile, line))
		{
			if (line[0] == '#')
				continue;
			std::stringstream ss(line);
			int idx, a, b, c, d;
			ss >> idx >> a >> b >> c >> d;
			tets_.row(i++) = Vec4i(a, b, c, d);
		}
		assert(num_tets_ == tets_.rows());
		infile.close();
	}
	void read_face_file(const std::string face_file)
	{
		std::ifstream infile(face_file);
		std::string line;
		std::getline(infile, line);
		std::stringstream ss(line);
		int tmp1, tmp2;
		ss >> num_faces_ >> tmp1 >> tmp2;
		faces_.resize(num_faces_, 3);
		int i = 0;
		while (std::getline(infile, line))
		{
			if (line[0] == '#')
				continue;
			std::stringstream ss(line);
			int idx, a, b, c;
			ss >> idx >> a >> b >> c;
			faces_.row(i++) = Vec3i(a, b, c);
		}
		assert(num_faces_ == faces_.rows());
		infile.close();
	}
	void extract_edges()
	{
		std::set<std::pair<int, int>> edges;
		for (int i = 0; i < num_tets_; i++)
		{
			int idx0 = tets_(i, 0), idx1 = tets_(i, 1), idx2 = tets_(i, 2), idx3 = tets_(i, 3);
			edges.insert(std::make_pair(std::min(idx0, idx1), std::max(idx0, idx1)));
			edges.insert(std::make_pair(std::min(idx0, idx2), std::max(idx0, idx2)));
			edges.insert(std::make_pair(std::min(idx0, idx3), std::max(idx0, idx3)));
			edges.insert(std::make_pair(std::min(idx1, idx2), std::max(idx1, idx2)));
			edges.insert(std::make_pair(std::min(idx1, idx3), std::max(idx1, idx3)));
			edges.insert(std::make_pair(std::min(idx2, idx3), std::max(idx2, idx3)));
		}
		int eidx = 0;
		edges_.resize(edges.size(), 2);
		for (auto e = edges.begin(); e != edges.end(); e++)
		{
			edges_.row(eidx++) = Vec2i(e->first, e->second);
		}
		num_edges_ = edges_.rows();
	}

public:
	Softbody(const std::string node_file, const std::string ele_file, const std::string face_file)
	{
		read_node_file(node_file);
		read_ele_file(ele_file);
		read_face_file(face_file);
		extract_edges();
	}

	void init_physical_data()
	{
		old_pos_ = pos_;
		vel_ = Eigen::MatrixXd::Zero(num_vertices_, 3);
		w_ = Eigen::VectorXd::Ones(num_vertices_);
	}
	void init_constraints()
	{
		rest_length_.resize(num_edges_);
		for (int i = 0; i < num_edges_; i++)
		{
			int idx0 = edges_(i, 0), idx1 = edges_(i, 1);
			rest_length_[i] = (pos_.row(idx0) - pos_.row(idx1)).norm();
		}
		rest_volume_.resize(num_tets_);
		for (int i = 0; i < num_tets_; i++)
		{
			int idx0 = tets_(i, 0), idx1 = tets_(i, 1), idx2 = tets_(i, 2), idx3 = tets_(i, 3);
			Vec3d p10 = pos_.row(idx1) - pos_.row(idx0);
			Vec3d p20 = pos_.row(idx2) - pos_.row(idx0);
			Vec3d p30 = pos_.row(idx3) - pos_.row(idx0);
			rest_volume_[i] = p30.dot(p10.cross(p20)) / 6.0;
		}
	}

	void init_xpbd(Real dis_alpha, Real vol_alpha, Real h){
		lambda_.resize(num_edges_ + num_tets_);
		lambda_.setZero();
		dis_alpha_tilde_ = dis_alpha / h / h;
		vol_alpha_tilde_ = vol_alpha / h / h / h;
	}

	void init_static_points(const std::vector<int> &static_points)
	{
		for (auto p : static_points)
			w_[p] = 0.0;
	}

	void semi_euler(Real h)
	{
		old_pos_ = pos_;
		for (int i = 0; i < num_vertices_; i++)
		{
			if (w_[i] != 0.0f)
			{
				vel_.row(i) += h * Vec3d(0.0, -9.8, 0.0);
				pos_.row(i) += h * vel_.row(i);
			}
		}
	}
	void solve_distance_constraints()
	{
		for (int i = 0; i < num_edges_; i++)
		{
			int idx0 = edges_(i, 0), idx1 = edges_(i, 1);
			Real w0 = w_[idx0], w1 = w_[idx1];
			Real w_sum = w0 + w1;
			if (w_sum == 0.0)
				continue;

			Vec3d p0p1 = pos_.row(idx0) - pos_.row(idx1);
			Real len = p0p1.norm();
			if (len == 0.0)
				continue;
			Real constraint = len - rest_length_[i];
			Vec3d normal = p0p1.normalized();
			Real delta_lambda = (constraint - lambda_[i] * dis_alpha_tilde_ )/ (w_sum + dis_alpha_tilde_);
			lambda_[i] += delta_lambda;
			Vec3d corr = 1.0 * delta_lambda * normal;
			if (w0 != 0.0)
			{
				pos_.row(idx0) += -w0 * corr;
			}
			if (w1 != 0.0)
			{
				pos_.row(idx1) -= -w1 * corr;
			}
		}
	}
	void solve_volume_constraints()
	{
		for (int i = 0; i < num_tets_; i++)
		{
			int idx0 = tets_(i, 0), idx1 = tets_(i, 1), idx2 = tets_(i, 2), idx3 = tets_(i, 3);
			Vec3d p0 = pos_.row(idx0);
			Vec3d p1 = pos_.row(idx1);
			Vec3d p2 = pos_.row(idx2);
			Vec3d p3 = pos_.row(idx3);
			Real w0 = w_[idx0];
			Real w1 = w_[idx1];
			Real w2 = w_[idx2];
			Real w3 = w_[idx3];

			Real volume = (1.0f / 6.0f) * (p1 - p0).cross(p2 - p0).dot(p3 - p0);

			Vec3d grad0 = (p3 - p1).cross(p2 - p1);
			Vec3d grad1 = (p2 - p0).cross(p3 - p0);
			Vec3d grad2 = (p3 - p0).cross(p1 - p0);
			Vec3d grad3 = (p1 - p0).cross(p2 - p0);

			Real delta_lambda = w0 * grad0.squaredNorm() + w1 * grad1.squaredNorm() + w2 * grad2.squaredNorm() + w2 * grad3.squaredNorm();

			if (fabs(delta_lambda) < 1e-6)
				continue;

			Real constraint = volume - rest_volume_[i];

			delta_lambda = (constraint - lambda_[num_edges_+i] * vol_alpha_tilde_) / (delta_lambda + vol_alpha_tilde_);
			lambda_[num_edges_+i] += delta_lambda;

			Real stiffness = 1.0;
			if (w0 != 0.0f)
				pos_.row(idx0) -= stiffness * grad0 * delta_lambda * w0;
			if (w1 != 0.0f)
				pos_.row(idx1) -= stiffness * grad1 * delta_lambda * w1;
			if (w2 != 0.0f)
				pos_.row(idx2) -= stiffness * grad2 * delta_lambda * w2;
			if (w3 != 0.0f)
				pos_.row(idx3) -= stiffness * grad3 * delta_lambda * w3;
		}
	}

	void collision_response()
	{
		for (int i = 0; i < num_vertices_; i++)
		{
			if (w_[i] != 0.0f)
			{
				if (pos_(i, 1) < -20.0f)
					pos_(i, 1) = -20.0f;
			}
		}
	}

	void update_velocity(Real h)
	{
		for (int i = 0; i < num_vertices_; i++)
		{
			if (w_[i] != 0.0f)
			{
				vel_.row(i) = (pos_.row(i) - old_pos_.row(i)) / h;
			}
		}
	}

	void collision_response(const Eigen::MatrixXd &Q)
	{

		// Choose type of signing to use
		igl::SignedDistanceType sign_type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
		Eigen::VectorXd S;
		Eigen::VectorXi I;
		Eigen::MatrixXd C, N;
		igl::signed_distance(Q, pos_, faces_, sign_type, S, I, C, N);

		for (int i = 0; i < Q.rows(); i++)
		{
			if (S[i] < 0.0)
			{ // collision response
				// std::cout<< "The " << i << "-th point of the sphere is in the soft body!" << std::endl;
				// std::cout<< "The signed distance: " << S[i] << std::endl;
				// std::cout<< "The collision point: "<<Q.row(i)<<std::endl;
				// std::cout<< "The triangle index: " << I[i] << std::endl;
				// std::cout<< "The barycentric corrdinates is "<<C.row(i)<<std::endl;
				// std::cout<< "The collision normal: "<<N.row(i)<<std::endl;

				Vec3i tri_idx = faces_.row(I[i]);
				// Vec3d p0 = pos_.row(tri_idx[0]);
				// Vec3d p1 = pos_.row(tri_idx[1]);
				// Vec3d p2 = pos_.row(tri_idx[2]);
				// Real constraint = S[i]; // signed distance (penetration depth)

				// Vec3d p10 = p1 - p0;
				// Vec3d p20 = p2 - p0;
				// Vec3d normal = p10.cross(p20);
				// Real de = 1.0 / normal.norm();
				// normal.normalize();
				// Eigen::Vector3d n = normal.transpose();
				// Vec3d p = Q.row(i);
				// Eigen::Matrix3d p10x, p20x;
				// p10x << 0.0, -p10[2], p10[1],
				//	p10[2], 0.0, -p10[0],
				//	-p10[1], p10[0], 0.0;
				// p20x << 0.0, -p20[2], p20[1],
				//	p20[2], 0.0, -p20[0],
				//	-p20[1], p20[0], 0.0;
				// Vec3d g1 = (p - p0) * de * (-p20x + n * (n.cross(p20).transpose()));
				// Vec3d g2 = (p - p0) * -de * (-p10x + n * (n.cross(p10).transpose()));
				// Vec3d g0 = -(g1 + g2);

				// Real w0 = w_[tri_idx[0]];
				// Real w1 = w_[tri_idx[1]];
				// Real w2 = w_[tri_idx[2]];

				// Real delta_lambda = w0 * g0.squaredNorm() + w1 * g1.squaredNorm() + w2 * g2.squaredNorm();
				// if(delta_lambda == 0.0) continue;
				// delta_lambda = constraint / delta_lambda;

				// Real stiffness = 0.5f;
				// if (w0 != 0.0f)
				//	pos_.row(tri_idx[0]) += stiffness * delta_lambda * w0 * g0;
				// if (w1 != 0.0f)
				//	pos_.row(tri_idx[1]) += stiffness * delta_lambda * w1 * g1;
				// if (w2 != 0.0f)
				//	pos_.row(tri_idx[2]) += stiffness * delta_lambda * w2 * g2;

				/*--------------first strategy--------------*/
				// Real stiffness = 1.0;
				// Real w0 = w_[tri_idx[0]];
				// Real w1 = w_[tri_idx[1]];
				// Real w2 = w_[tri_idx[2]];
				// if (w0 != 0.0f)
				//	pos_.row(tri_idx[0]) += stiffness * N.row(i) * S[i];
				// if (w1 != 0.0f)
				//	pos_.row(tri_idx[1]) += stiffness * N.row(i) * S[i];
				// if (w2 != 0.0f)
				//	pos_.row(tri_idx[2]) += stiffness * N.row(i) * S[i];

				/*----------------Second strategy-----------------*/
				// Project the point to the surface of the sphere
				Vec3d p0 = pos_.row(tri_idx[0]);
				Vec3d p1 = pos_.row(tri_idx[1]);
				Vec3d p2 = pos_.row(tri_idx[2]);
				Real epsilon = 5.0e-5;
				if ((p0 - P).norm() < 1.0 && w_[tri_idx[0]] != 0.0)
				{
					pos_.row(tri_idx[0]) = P + (1 + epsilon) * (p0 - P).normalized();
				}
				if ((p1 - P).norm() < 1.0 && w_[tri_idx[1]] != 0.0)
				{
					pos_.row(tri_idx[1]) = P + (1 + epsilon) * (p1 - P).normalized();
				}
				if ((p2 - P).norm() < 1.0 && w_[tri_idx[2]] != 0.0)
				{
					pos_.row(tri_idx[2]) = P + (1 + epsilon) * (p2 - P).normalized();
				}
			}
		}

		// Initialize AABB tree
		/*		igl::AABB<Eigen::MatrixXd, 3> tree;
				tree.init(pos_, tets_);
				Eigen::VectorXi I;
				igl::in_element(pos_, tets_, Q, tree, I);

				for (int i = 0; i < Q.rows(); i++) {
					if (I[i] > 0.0) {
						std::cout<<"The " << I[i] <<"-th element is in collision with the sphere!"<<std::endl;
					}
				}*/
	}

	void update(Real h, int maxIte, const Eigen::MatrixXd &Q)
	{
		semi_euler(h);

		lambda_.setZero();

		for (int ite = 0; ite < maxIte; ite++)
		{
			solve_distance_constraints();
			solve_volume_constraints();
			// collision_response(Q);
		}
		collision_response();
		update_velocity(h);
	}

	void write_obj(const std::string obj_file)
	{
		std::ofstream outfile(obj_file);
		for (int i = 0; i < num_vertices_; i++)
		{
			outfile << "v " << pos_(i, 0) << " " << pos_(i, 1) << " " << pos_(i, 2) << std::endl;
		}
		for (int i = 0; i < num_faces_; i++)
		{
			outfile << "f " << faces_(i, 0) + 1 << " " << faces_(i, 1) + 1 << " " << faces_(i, 2) + 1 << std::endl;
		}
		outfile.close();
	}

	int get_num_vertices() { return num_vertices_; }
	int get_num_edges() { return num_edges_; }
	int get_num_faces() { return num_faces_; }
	int get_num_tets() { return num_tets_; }
};

Softbody softbody("../model/liver/liver.node", "../model/liver/liver.ele", "../model/liver/liver.face");
Real h = 0.02;
int maxIte = 5;
bool pause = true;
Real speed = 0.1;
int frame = 0;

void update_sphere(const Vec3d &delta)
{
	P += delta;
	for (int i = 0; i < V_sphere.rows(); i++)
		V_sphere.row(i) += delta;
}

bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
	if (!pause)
	{
		frame++;
		Real height = 0.3 * std::sin(frame * 0.1);
		// update_sphere(Vec3d(0.0, height, 0.0));

		// for (int i = 0; i < pull_points.size(); i++)
		// softbody.pos_.row(pull_points[i]) += Vec3d(0.0, height, 0.0);

		softbody.update(h, maxIte, V_sphere);
	}
	// softbody.write_obj("../model/liver2/f" +std::to_string(frame++) + ".obj");

	viewer.data().clear();
	viewer.data().set_mesh(softbody.pos_, softbody.faces_);
	viewer.core().align_camera_center(softbody.pos_, softbody.faces_);
	viewer.append_mesh();
	viewer.data().set_mesh(V_sphere, T_sphere);
	// viewer.data_list[0].set_colors(Eigen::RowVector3d(1, 0, 0));
	// viewer.data_list[1].set_colors(Eigen::RowVector3d(0, 1, 0));
	return false;
}
bool post_draw(igl::opengl::glfw::Viewer &viewer)
{
	for (auto &data : viewer.data_list)
	{
		data.clear();
	}
	return false;
}

bool key_pressed(igl::opengl::glfw::Viewer &viewer, unsigned int key, int modifiers)
{
	if (key == 'p')
	{
		pause = !pause;
		std::cout << std::boolalpha << pause << std::endl;
	}
	else if (key == 'a')
	{
		update_sphere(Vec3d(-speed, 0, 0));
	}
	else if (key == 'w')
	{
		update_sphere(Vec3d(0.0, speed, 0.0));
	}
	else if (key == 's')
	{
		update_sphere(Vec3d(0.0, -speed, 0.0));
	}
	else if (key == 'd')
	{
		update_sphere(Vec3d(speed, 0.0, 0.0));
	}
	else if (key == 'r')
	{
		update_sphere(Vec3d(0.0, 0.0, speed));
	}
	else if (key == 't')
	{
		update_sphere(Vec3d(0.0, 0.0, -speed));
	}
	else if (key == 'e')
	{
		speed += 0.1;
	}
	else if (key == 'c')
	{
		speed -= 0.1;
	}
	else if (key == 'b')
	{
		std::cout << "sphere poinstion: " << P.row(0) << std::endl;
	}
	return false;
}
int main()
{
	init_spehre(V_sphere, T_sphere, C_sphere);
	softbody.init_physical_data();
	softbody.init_xpbd(1.0e-5, 1.0e-5, h);
	// std::vector<int> static_points;
	// for (int i = 180; i < 200; i++)
	// 	static_points.push_back(i);
	// for (int i = 250; i < 260; i++)
	// 	static_points.push_back(i);
	// for (auto& p : pull_points)
	// 	static_points.push_back(p);

	std::cout << "Number of vertices: " << softbody.get_num_vertices() << std::endl;
	std::cout << "Number of tets: " << softbody.get_num_tets() << std::endl;
	std::cout << "Number of edges: " << softbody.get_num_edges() << std::endl;
	std::cout << "Number of faces: " << softbody.get_num_faces() << std::endl;

	// softbody.init_static_points(static_points);
	softbody.init_constraints();

	igl::opengl::glfw::Viewer viewer;
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_post_draw = &post_draw;
	viewer.callback_key_pressed = &key_pressed;
	viewer.core().is_animating = true;
	viewer.data().set_face_based(false);
	viewer.launch();

	return 0;
}