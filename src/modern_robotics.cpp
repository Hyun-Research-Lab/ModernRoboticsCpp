#include "../include/modern_robotics.h"

/*
 * modernRobotics.cpp
 * Adapted from modern_robotics.py provided by modernrobotics.org
 * Provides useful Jacobian and frame representation functions
 */
#include <Eigen/Dense>
#include <cmath>
#include <vector>

# define M_PI           3.14159265358979323846  /* pi */

namespace mr {

	/* Function: Find if the value is negligible enough to consider 0
	 * Inputs: value to be checked as a double
	 * Returns: Boolean of true-ignore or false-can't ignore
	 */
	bool NearZero(const double val) {
		return (std::abs(val) < .000001);
	}

	/*
	 * Function: Calculate the 6x6 matrix [adV] of the given 6-vector
	 * Input: Eigen::VectorXd (6x1)
	 * Output: Eigen::MatrixXd (6x6)
	 * Note: Can be used to calculate the Lie bracket [V1, V2] = [adV1]V2
	 */
	Matrix<6> ad(const Vector<6> &V) {
		Matrix<3> omgmat = VecToso3(Vector<3>(V(0), V(1), V(2)));

		Matrix<6> result;
		result.topLeftCorner<3, 3>() = omgmat;
		result.topRightCorner<3, 3>() = Matrix<3>::Zero();
		result.bottomLeftCorner<3, 3>() = VecToso3(Vector<3>(V(3), V(4), V(5)));
		result.bottomRightCorner<3, 3>() = omgmat;
		return result;
	}

	/* Function: Returns the skew symmetric matrix representation of an angular velocity vector
	 * Input: Vector<3> 3x1 angular velocity vector
	 * Returns: Eigen::MatrixXd 3x3 skew symmetric matrix
	 */
	Matrix<3> VecToso3(const Vector<3>& omg) {
		Matrix<3> m_ret;
		m_ret << 0, -omg(2), omg(1),
			omg(2), 0, -omg(0),
			-omg(1), omg(0), 0;
		return m_ret;
	}


	/* Function: Returns angular velocity vector represented by the skew symmetric matrix
	 * Inputs: Eigen::MatrixXd 3x3 skew symmetric matrix
	 * Returns: Vector<3> 3x1 angular velocity
	 */
	Vector<3> so3ToVec(const Matrix<3>& so3mat) {
		Vector<3> v_ret;
		v_ret << so3mat(2, 1), so3mat(0, 2), so3mat(1, 0);
		return v_ret;
	}


	/* Function: Translates an exponential rotation into it's individual components
	 * Inputs: Exponential rotation (rotation matrix in terms of a rotation axis
	 *				and the angle of rotation)
	 * Returns: The axis and angle of rotation as [x, y, z, theta]
	 */
	std::tuple<Vector<3>, double> AxisAng3(const Vector<3>& expc3) {
		double theta = expc3.norm();
		return std::make_tuple(expc3 / theta, theta);
	}


	/* Function: Translates an exponential rotation into a rotation matrix
	 * Inputs: exponenential representation of a rotation
	 * Returns: Rotation matrix
	 */
	Matrix<3> MatrixExp3(const Matrix<3>& so3mat) {
		Vector<3> omgtheta = so3ToVec(so3mat);

		Matrix<3> m_ret = Matrix<3>::Identity();
		if (NearZero(so3mat.norm())) {
			return m_ret;
		}
		else {
			double theta = std::get<1>(AxisAng3(omgtheta));
			Matrix<3> omgmat = so3mat * (1 / theta);
			return m_ret + std::sin(theta) * omgmat + ((1 - std::cos(theta)) * (omgmat * omgmat));
		}
	}


	/* Function: Computes the matrix logarithm of a rotation matrix
	 * Inputs: Rotation matrix
	 * Returns: matrix logarithm of a rotation
	 */
	Matrix<3> MatrixLog3(const Matrix<3>& R) {
		double acosinput = (R.trace() - 1) / 2.0;
		Matrix<3> m_ret = Matrix<3>::Zero();
		if (acosinput >= 1)
			return m_ret;
		else if (acosinput <= -1) {
			Vector<3> omg;
			if (!NearZero(1 + R(2, 2)))
				omg = (1.0 / std::sqrt(2 * (1 + R(2, 2))))*Vector<3>(R(0, 2), R(1, 2), 1 + R(2, 2));
			else if (!NearZero(1 + R(1, 1)))
				omg = (1.0 / std::sqrt(2 * (1 + R(1, 1))))*Vector<3>(R(0, 1), 1 + R(1, 1), R(2, 1));
			else
				omg = (1.0 / std::sqrt(2 * (1 + R(0, 0))))*Vector<3>(1 + R(0, 0), R(1, 0), R(2, 0));
			m_ret = VecToso3(M_PI * omg);
			return m_ret;
		}
		else {
			double theta = std::acos(acosinput);
			m_ret = theta / 2.0 / sin(theta)*(R - R.transpose());
			return m_ret;
		}
	}

	/* Function: Combines a rotation matrix and position vector into a single
	 * 				Special Euclidian Group (SE3) homogeneous transformation matrix
	 * Inputs: Rotation Matrix (R), Position Vector (p)
	 * Returns: Matrix of T = [ [R, p],
	 *						    [0, 1] ]
	 */
	Matrix<4> RpToTrans(const Matrix<3>& R, const Vector<3>& p) {
		Matrix<4> m_ret(4, 4);
		m_ret << R, p,
			0, 0, 0, 1;
		return m_ret;
	}


	/* Function: Separates the rotation matrix and position vector from
	 *				the transfomation matrix representation
	 * Inputs: Homogeneous transformation matrix
	 * Returns: std::vector of [rotation matrix, position vector]
	 */
	std::tuple<Matrix<3>, Vector<3>> TransToRp(const Matrix<4>& T) {
		// Get top left 3x3 corner
		Matrix<3> R_out = T.block<3, 3>(0, 0);

		// Get top right 3x1 corner
		Vector<3> p_out = Vector<3>(T(0, 3), T(1, 3), T(2, 3));

		return std::make_tuple(R_out, p_out);
	}


	/* Function: Translates a spatial velocity vector into a transformation matrix
	 * Inputs: Spatial velocity vector [angular velocity, linear velocity]
	 * Returns: Transformation matrix
	 */
	Matrix<4> VecTose3(const Vector<6> &V) {
		// Separate angular (exponential representation) and linear velocities
		Vector<3> exp(V(0), V(1), V(2));
		Vector<3> linear(V(3), V(4), V(5));

		// Fill in values to the appropriate parts of the transformation matrix
		Matrix<4> m_ret;
		m_ret << VecToso3(exp), linear,
			0, 0, 0, 0;

		return m_ret;
	}


	/* Function: Translates a transformation matrix into a spatial velocity vector
	 * Inputs: Transformation matrix
	 * Returns: Spatial velocity vector [angular velocity, linear velocity]
	 */
	Vector<6> se3ToVec(const Matrix<4>& T) {
		Vector<6> m_ret;
		m_ret << T(2, 1), T(0, 2), T(1, 0), T(0, 3), T(1, 3), T(2, 3);
		return m_ret;
	}


	/* Function: Provides the adjoint representation of a transformation matrix
	 *			 Used to change the frame of reference for spatial velocity vectors
	 * Inputs: 4x4 Transformation matrix SE(3)
	 * Returns: 6x6 Adjoint Representation of the matrix
	 */
	Matrix<6> Adjoint(const Matrix<4>& T) {
		
		// T = (R, p)
		Matrix<3> R;
		Vector<3> p;
		std::tie(R, p) = TransToRp(T);

		Matrix<6> ad_ret;
		ad_ret << R, Matrix<3>::Zero(), 
				VecToso3(p) * R, R;
		return ad_ret;
	}


	/* Function: Rotation expanded for screw axis
	 * Inputs: se3 matrix representation of exponential coordinates (transformation matrix)
	 * Returns: 4x4 Matrix representing the rotation
	 */
	Matrix<4> MatrixExp6(const Matrix<4>& se3mat) {
		// Extract the angular velocity vector from the transformation matrix
		Matrix<3> se3mat_cut = se3mat.block<3, 3>(0, 0);
		Vector<3> omgtheta = so3ToVec(se3mat_cut);
		Matrix<4> m_ret;

		// If negligible rotation, m_Ret = [[Identity, angular velocty ]]
		//									[	0	 ,		1		   ]]
		if (NearZero(omgtheta.norm())) {
			// Reuse previous variables that have our required size
			se3mat_cut = Matrix<3>::Identity();
			omgtheta << se3mat(0, 3), se3mat(1, 3), se3mat(2, 3);
			m_ret << se3mat_cut, omgtheta,
				0, 0, 0, 1;
			return m_ret;
		}
		// If not negligible, MR page 105
		else {
			double theta = std::get<1>(AxisAng3(omgtheta));
			Matrix<3> omgmat = se3mat.block<3, 3>(0, 0) / theta;
			Matrix<3> expExpand = Matrix<3>::Identity() * theta + (1 - std::cos(theta)) * omgmat + ((theta - std::sin(theta)) * (omgmat * omgmat));
			Vector<3> linear(se3mat(0, 3), se3mat(1, 3), se3mat(2, 3));
			Vector<3> GThetaV = (expExpand*linear) / theta;
			m_ret << MatrixExp3(se3mat_cut), GThetaV,
				0, 0, 0, 1;
			return m_ret;
		}

	}

	Matrix<4> MatrixLog6(const Matrix<4>& T) {
		Matrix<4> m_ret;

		// Extract the rotation and position vectors
		Matrix<3> R;
		Vector<3> p;
		std::tie(R, p) = TransToRp(T);

		Matrix<3> omgmat = MatrixLog3(R);
		if (NearZero(omgmat.norm())) {
			m_ret << Matrix<3>::Zero(), p,
				0, 0, 0, 0;
		}
		else {
			double theta = std::acos((R.trace() - 1) / 2.0);
			Matrix<3> logExpand1 = Matrix<3>::Identity() - omgmat / 2.0;
			Matrix<3> logExpand2 = (1.0 / theta - 1.0 / std::tan(theta / 2.0) / 2)*omgmat*omgmat / theta;
			Matrix<3> logExpand = logExpand1 + logExpand2;
			m_ret << omgmat, logExpand*p,
				0, 0, 0, 0;
		}
		return m_ret;
	}


	/* Function: Compute end effector frame (used for current spatial position calculation)
	 * Inputs: Home configuration (position and orientation) of end-effector
	 *		   The joint screw axes in the space frame when the manipulator
	 *             is at the home position
	 * 		   A list of joint coordinates.
	 * Returns: Transfomation matrix representing the end-effector frame when the joints are
	 *				at the specified coordinates
	 * Notes: FK means Forward Kinematics
	 */
	Matrix<4> FKinSpace(const Matrix<4>& M, const Eigen::MatrixXd& Slist, const Eigen::VectorXd& thetaList) {
		Matrix<4> T = M;
		for (int i = (thetaList.size() - 1); i > -1; i--) {
			T = MatrixExp6(VecTose3(Slist.col(i)*thetaList(i))) * T;
		}
		return T;
	}

	/*
	 * Function: Compute end effector frame (used for current body position calculation)
	 * Inputs: Home configuration (position and orientation) of end-effector
	 *		   The joint screw axes in the body frame when the manipulator
	 *             is at the home position
	 * 		   A list of joint coordinates.
	 * Returns: Transfomation matrix representing the end-effector frame when the joints are
	 *				at the specified coordinates
	 * Notes: FK means Forward Kinematics
	 */
	Matrix<4> FKinBody(const Matrix<4>& M, const Eigen::MatrixXd& Blist, const Eigen::VectorXd& thetaList) {
		Matrix<4> T = M;
		for (int i = 0; i < thetaList.size(); i++) {
			T = T * MatrixExp6(VecTose3(Blist.col(i)*thetaList(i)));
		}
		return T;
	}


	/* Function: Gives the space Jacobian
	 * Inputs: Screw axis in home position, joint configuration
	 * Returns: 6xn Spatial Jacobian
	 */
	Eigen::MatrixXd JacobianSpace(const Eigen::MatrixXd& Slist, const Eigen::MatrixXd& thetaList) {
		Eigen::MatrixXd Js = Slist;
		Eigen::MatrixXd T = Eigen::MatrixXd::Identity(4, 4);
		Eigen::VectorXd sListTemp(Slist.col(0).size());
		for (int i = 1; i < thetaList.size(); i++) {
			sListTemp << Slist.col(i - 1) * thetaList(i - 1);
			T = T * MatrixExp6(VecTose3(sListTemp));
			// std::cout << "array: " << sListTemp << std::endl;
			Js.col(i) = Adjoint(T) * Slist.col(i);
		}

		return Js;
	}

	/*
	 * Function: Gives the body Jacobian
	 * Inputs: Screw axis in BODY position, joint configuration
	 * Returns: 6xn Bobdy Jacobian
	 */
	Eigen::MatrixXd JacobianBody(const Eigen::MatrixXd& Blist, const Eigen::MatrixXd& thetaList) {
		Eigen::MatrixXd Jb = Blist;
		Eigen::MatrixXd T = Eigen::MatrixXd::Identity(4, 4);
		Eigen::VectorXd bListTemp(Blist.col(0).size());
		for (int i = thetaList.size() - 2; i >= 0; i--) {
			bListTemp << Blist.col(i + 1) * thetaList(i + 1);
			T = T * MatrixExp6(VecTose3(-1 * bListTemp));
			// std::cout << "array: " << sListTemp << std::endl;
			Jb.col(i) = Adjoint(T) * Blist.col(i);
		}
		return Jb;
	}

	Matrix<4> TransInv(const Matrix<4>& T) {
		Matrix<3> R;
		Vector<3> p;
		std::tie(R, p) = TransToRp(T);
		Matrix<3> Rt = R.transpose();
		Vector<3> t = -Rt * p;
		
		Matrix<4> inv = Matrix<4>::Zero();
		inv.block(0, 0, 3, 3) = Rt;
		inv.block(0, 3, 3, 1) = t;
		inv(3, 3) = 1;
		return inv;
	}

	Matrix<3> RotInv(const Matrix<3>& rotMatrix) {
		return rotMatrix.transpose();
	}

	Vector<6> ScrewToAxis(Vector<3> q, Vector<3> s, double h) {
		Vector<6> axis;
		axis.segment(0, 3) = s;
		axis.segment(3, 3) = q.cross(s) + (h * s);
		return axis;
	}

	std::tuple<Vector<6>, double> AxisAng6(const Vector<6>& expc6) {
		double theta = Vector<3>(expc6(0), expc6(1), expc6(2)).norm();
		if (NearZero(theta)) {
			theta = Vector<3>(expc6(3), expc6(4), expc6(5)).norm();
		}
		return std::make_tuple(expc6 / theta, theta);
	}

	Matrix<3> ProjectToSO3(const Matrix<3> & M) {
		Eigen::JacobiSVD<Matrix<3>> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Matrix<3> R = svd.matrixU() * svd.matrixV().transpose();
		if (R.determinant() < 0)
			// In this case the result may be far from M; reverse sign of 3rd column
			R.col(2) *= -1;
		return R;
	}

	Matrix<4> ProjectToSE3(const Matrix<4>& M) {
		Matrix<3> R = M.block<3, 3>(0, 0);
		Vector<3> t = M.block<3, 1>(0, 3);
		Matrix<4> T = RpToTrans(ProjectToSO3(R), t);
		return T;
	}

	double DistanceToSO3(const Matrix<3>& M) {
		if (M.determinant() > 0)
			return (M.transpose() * M - Matrix<3>::Identity()).norm();
		else
			return 1.0e9;
	}

	double DistanceToSE3(const Matrix<4>& T) {
		Matrix<3> matR = T.block<3, 3>(0, 0);
		if (matR.determinant() > 0) {
			Matrix<4> m_ret;
			m_ret << matR.transpose()*matR, Vector<3>::Zero(3),
				T.row(3);
			m_ret = m_ret - Matrix<4>::Identity();
			return m_ret.norm();
		}
		else
			return 1.0e9;
	}

	bool TestIfSO3(const Matrix<3>& M) {
		return std::abs(DistanceToSO3(M)) < 1e-3;
	}

	bool TestIfSE3(const Matrix<4>& T) {
		return std::abs(DistanceToSE3(T)) < 1e-3;
	}
	bool IKinBody(const Eigen::MatrixXd& Blist, const Matrix<4>& M, const Matrix<4>& T,
		Eigen::VectorXd& thetalist, double eomg, double ev) {
		int i = 0;
		int maxiterations = 20;
		Matrix<4> Tfk = FKinBody(M, Blist, thetalist);
		Matrix<4> Tdiff = TransInv(Tfk)*T;
		Eigen::VectorXd Vb = se3ToVec(MatrixLog6(Tdiff));
		Vector<3> angular(Vb(0), Vb(1), Vb(2));
		Vector<3> linear(Vb(3), Vb(4), Vb(5));

		bool err = (angular.norm() > eomg || linear.norm() > ev);
		Eigen::MatrixXd Jb;
		while (err && i < maxiterations) {
			Jb = JacobianBody(Blist, thetalist);
			thetalist += Jb.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Vb);
			i += 1;
			// iterate
			Tfk = FKinBody(M, Blist, thetalist);
			Tdiff = TransInv(Tfk)*T;
			Vb = se3ToVec(MatrixLog6(Tdiff));
			angular = Vector<3>(Vb(0), Vb(1), Vb(2));
			linear = Vector<3>(Vb(3), Vb(4), Vb(5));
			err = (angular.norm() > eomg || linear.norm() > ev);
		}
		return !err;
	}

	bool IKinSpace(const Eigen::MatrixXd& Slist, const Matrix<4>& M, const Matrix<4>& T,
		Eigen::VectorXd& thetalist, double eomg, double ev) {
		int i = 0;
		int maxiterations = 20;
		Matrix<4> Tfk = FKinSpace(M, Slist, thetalist);
		Matrix<4> Tdiff = TransInv(Tfk)*T;
		Vector<6> Vs = Adjoint(Tfk)*se3ToVec(MatrixLog6(Tdiff));
		Vector<3> angular(Vs(0), Vs(1), Vs(2));
		Vector<3> linear(Vs(3), Vs(4), Vs(5));

		bool err = (angular.norm() > eomg || linear.norm() > ev);
		Eigen::MatrixXd Js;
		while (err && i < maxiterations) {
			Js = JacobianSpace(Slist, thetalist);
			thetalist += Js.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Vs);
			i += 1;
			// iterate
			Tfk = FKinSpace(M, Slist, thetalist);
			Tdiff = TransInv(Tfk)*T;
			Vs = Adjoint(Tfk)*se3ToVec(MatrixLog6(Tdiff));
			angular = Vector<3>(Vs(0), Vs(1), Vs(2));
			linear = Vector<3>(Vs(3), Vs(4), Vs(5));
			err = (angular.norm() > eomg || linear.norm() > ev);
		}
		return !err;
	}

	/*
	* Function: This function uses forward-backward Newton-Euler iterations to solve the
	* equation:
	* taulist = Mlist(thetalist) * ddthetalist + c(thetalist, dthetalist) ...
	*           + g(thetalist) + Jtr(thetalist) * Ftip
	* Inputs:
	*  thetalist: n-vector of joint variables
	*  dthetalist: n-vector of joint rates
	*  ddthetalist: n-vector of joint accelerations
	*  g: Gravity vector g
	*  Ftip: Spatial force applied by the end-effector expressed in frame {n+1}
	*  Mlist: List of link frames {i} relative to {i-1} at the home position
	*  Glist: Spatial inertia matrices Gi of the links
	*  Slist: Screw axes Si of the joints in a space frame, in the format
	*         of a matrix with the screw axes as the columns.
	*
	* Outputs:
	*  taulist: The n-vector of required joint forces/torques
	*
	*/
	Eigen::VectorXd InverseDynamics(const Eigen::VectorXd& thetalist, const Eigen::VectorXd& dthetalist, const Eigen::VectorXd& ddthetalist,
									const Eigen::VectorXd& g, const Eigen::VectorXd& Ftip, const std::vector<Eigen::MatrixXd>& Mlist,
									const std::vector<Eigen::MatrixXd>& Glist, const Eigen::MatrixXd& Slist) {
	    // the size of the lists
		int n = thetalist.size();

		Matrix<4> Mi = Matrix<4>::Identity();
		Eigen::MatrixXd Ai = Eigen::MatrixXd::Zero(6,n);
		std::vector<Matrix<6>> AdTi;
		for (int i = 0; i < n+1; i++) {
			AdTi.push_back(Matrix<6>::Zero());
		}
		Eigen::MatrixXd Vi = Eigen::MatrixXd::Zero(6,n+1);    // velocity
		Eigen::MatrixXd Vdi = Eigen::MatrixXd::Zero(6,n+1);   // acceleration

		Vdi.block(3, 0, 3, 1) = - g;
		AdTi[n] = mr::Adjoint(mr::TransInv(Mlist[n]));
		Eigen::VectorXd Fi = Ftip;

		Eigen::VectorXd taulist = Eigen::VectorXd::Zero(n);

		// forward pass
		for (int i = 0; i < n; i++) {
			Mi = Mi * Mlist[i];
			Ai.col(i) = mr::Adjoint(mr::TransInv(Mi))*Slist.col(i);

			AdTi[i] = mr::Adjoint(mr::MatrixExp6(mr::VecTose3(Ai.col(i)*-thetalist(i)))
			          * mr::TransInv(Mlist[i]));

			Vi.col(i+1) = AdTi[i] * Vi.col(i) + Ai.col(i) * dthetalist(i);
			Vdi.col(i+1) = AdTi[i] * Vdi.col(i) + Ai.col(i) * ddthetalist(i)
						   + ad((Vector<6>)Vi.col(i+1)) * Ai.col(i) * dthetalist(i); // this index is different from book!
		}

		// backward pass
		for (int i = n-1; i >= 0; i--) {
			Fi = AdTi[i+1].transpose() * Fi + Glist[i] * Vdi.col(i+1)
			     - ad(Vi.col(i+1)).transpose() * (Glist[i] * Vi.col(i+1));
			taulist(i) = Fi.transpose() * Ai.col(i);
		}
		return taulist;
	}

	/*
	 * Function: This function calls InverseDynamics with Ftip = 0, dthetalist = 0, and
	 *   ddthetalist = 0. The purpose is to calculate one important term in the dynamics equation
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  g: Gravity vector g
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  grav: The 3-vector showing the effect force of gravity to the dynamics
	 *
	 */
	Eigen::VectorXd GravityForces(const Eigen::VectorXd& thetalist, const Eigen::VectorXd& g,
									const std::vector<Eigen::MatrixXd>& Mlist, const std::vector<Eigen::MatrixXd>& Glist, const Eigen::MatrixXd& Slist) {
	    int n = thetalist.size();
		Eigen::VectorXd dummylist = Eigen::VectorXd::Zero(n);
		Eigen::VectorXd dummyForce = Eigen::VectorXd::Zero(6);
		Eigen::VectorXd grav = mr::InverseDynamics(thetalist, dummylist, dummylist, g,
                                                dummyForce, Mlist, Glist, Slist);
		return grav;
	}

	/*
  	 * Function: This function calls InverseDynamics n times, each time passing a
	 * ddthetalist vector with a single element equal to one and all other
	 * inputs set to zero. Each call of InverseDynamics generates a single
	 * column, and these columns are assembled to create the inertia matrix.
	 *
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  M: The numerical inertia matrix M(thetalist) of an n-joint serial
	 *     chain at the given configuration thetalist.
	 */
	Eigen::MatrixXd MassMatrix(const Eigen::VectorXd& thetalist,
                                const std::vector<Eigen::MatrixXd>& Mlist, const std::vector<Eigen::MatrixXd>& Glist, const Eigen::MatrixXd& Slist) {
		int n = thetalist.size();
		Eigen::VectorXd dummylist = Eigen::VectorXd::Zero(n);
		Eigen::VectorXd dummyg = Eigen::VectorXd::Zero(3);
		Eigen::VectorXd dummyforce = Eigen::VectorXd::Zero(6);
		Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n,n);
		for (int i = 0; i < n; i++) {
			Eigen::VectorXd ddthetalist = Eigen::VectorXd::Zero(n);
			ddthetalist(i) = 1;
			M.col(i) = mr::InverseDynamics(thetalist, dummylist, ddthetalist,
                             dummyg, dummyforce, Mlist, Glist, Slist);
		}
		return M;
	}

	/*
  	 * Function: This function calls InverseDynamics with g = 0, Ftip = 0, and
     * ddthetalist = 0.
	 *
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  dthetalist: A list of joint rates
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  c: The vector c(thetalist,dthetalist) of Coriolis and centripetal
	 *     terms for a given thetalist and dthetalist.
	 */
	Eigen::VectorXd VelQuadraticForces(const Eigen::VectorXd& thetalist, const Eigen::VectorXd& dthetalist,
                                const std::vector<Eigen::MatrixXd>& Mlist, const std::vector<Eigen::MatrixXd>& Glist, const Eigen::MatrixXd& Slist) {
		int n = thetalist.size();
		Eigen::VectorXd dummylist = Eigen::VectorXd::Zero(n);
		Eigen::VectorXd dummyg = Eigen::VectorXd::Zero(3);
		Eigen::VectorXd dummyforce = Eigen::VectorXd::Zero(6);
		Eigen::VectorXd c = mr::InverseDynamics(thetalist, dthetalist, dummylist,
                             dummyg, dummyforce, Mlist, Glist, Slist);
		return c;
	}

	/*
  	 * Function: This function calls InverseDynamics with g = 0, dthetalist = 0, and
     * ddthetalist = 0.
	 *
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  Ftip: Spatial force applied by the end-effector expressed in frame {n+1}
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  JTFtip: The joint forces and torques required only to create the
	 *     end-effector force Ftip.
	 */
	Eigen::VectorXd EndEffectorForces(const Eigen::VectorXd& thetalist, const Eigen::VectorXd& Ftip,
								const std::vector<Eigen::MatrixXd>& Mlist, const std::vector<Eigen::MatrixXd>& Glist, const Eigen::MatrixXd& Slist) {
		int n = thetalist.size();
		Eigen::VectorXd dummylist = Eigen::VectorXd::Zero(n);
		Eigen::VectorXd dummyg = Eigen::VectorXd::Zero(3);

		Eigen::VectorXd JTFtip = mr::InverseDynamics(thetalist, dummylist, dummylist,
                             dummyg, Ftip, Mlist, Glist, Slist);
		return JTFtip;
	}

	/*
	 * Function: This function computes ddthetalist by solving:
	 * Mlist(thetalist) * ddthetalist = taulist - c(thetalist,dthetalist)
	 *                                  - g(thetalist) - Jtr(thetalist) * Ftip
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  dthetalist: n-vector of joint rates
	 *  taulist: An n-vector of joint forces/torques
	 *  g: Gravity vector g
	 *  Ftip: Spatial force applied by the end-effector expressed in frame {n+1}
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  ddthetalist: The resulting joint accelerations
	 *
	 */
	Eigen::VectorXd ForwardDynamics(const Eigen::VectorXd& thetalist, const Eigen::VectorXd& dthetalist, const Eigen::VectorXd& taulist,
									const Eigen::VectorXd& g, const Eigen::VectorXd& Ftip, const std::vector<Eigen::MatrixXd>& Mlist,
									const std::vector<Eigen::MatrixXd>& Glist, const Eigen::MatrixXd& Slist) {

		Eigen::VectorXd totalForce = taulist - mr::VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)
                 							 - mr::GravityForces(thetalist, g, Mlist, Glist, Slist)
                                             - mr::EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist);

		Eigen::MatrixXd M = mr::MassMatrix(thetalist, Mlist, Glist, Slist);

		// Use LDLT since M is positive definite
        Eigen::VectorXd ddthetalist = M.ldlt().solve(totalForce);

		return ddthetalist;
	}

	void EulerStep(Eigen::VectorXd& thetalist, Eigen::VectorXd& dthetalist, const Eigen::VectorXd& ddthetalist, double dt) {
		thetalist += dthetalist * dt;
		dthetalist += ddthetalist * dt;
		return;
	}

	Eigen::MatrixXd InverseDynamicsTrajectory(const Eigen::MatrixXd& thetamat, const Eigen::MatrixXd& dthetamat, const Eigen::MatrixXd& ddthetamat,
		const Eigen::VectorXd& g, const Eigen::MatrixXd& Ftipmat, const std::vector<Eigen::MatrixXd>& Mlist, const std::vector<Eigen::MatrixXd>& Glist,
		const Eigen::MatrixXd& Slist) {
		Eigen::MatrixXd thetamatT = thetamat.transpose();
		Eigen::MatrixXd dthetamatT = dthetamat.transpose();
		Eigen::MatrixXd ddthetamatT = ddthetamat.transpose();
		Eigen::MatrixXd FtipmatT = Ftipmat.transpose();

		int N = thetamat.rows();  // trajectory points
		int dof = thetamat.cols();
		Eigen::MatrixXd taumatT = Eigen::MatrixXd::Zero(dof, N);
		for (int i = 0; i < N; ++i) {
			taumatT.col(i) = InverseDynamics(thetamatT.col(i), dthetamatT.col(i), ddthetamatT.col(i), g, FtipmatT.col(i), Mlist, Glist, Slist);
		}
		Eigen::MatrixXd taumat = taumatT.transpose();
		return taumat;
	}

	std::vector<Eigen::MatrixXd> ForwardDynamicsTrajectory(const Eigen::VectorXd& thetalist, const Eigen::VectorXd& dthetalist, const Eigen::MatrixXd& taumat,
		const Eigen::VectorXd& g, const Eigen::MatrixXd& Ftipmat, const std::vector<Eigen::MatrixXd>& Mlist, const std::vector<Eigen::MatrixXd>& Glist,
		const Eigen::MatrixXd& Slist, double dt, int intRes) {
		Eigen::MatrixXd taumatT = taumat.transpose();
		Eigen::MatrixXd FtipmatT = Ftipmat.transpose();
		int N = taumat.rows();  // force/torque points
		int dof = taumat.cols();
		Eigen::MatrixXd thetamatT = Eigen::MatrixXd::Zero(dof, N);
		Eigen::MatrixXd dthetamatT = Eigen::MatrixXd::Zero(dof, N);
		thetamatT.col(0) = thetalist;
		dthetamatT.col(0) = dthetalist;
		Eigen::VectorXd thetacurrent = thetalist;
		Eigen::VectorXd dthetacurrent = dthetalist;
		Eigen::VectorXd ddthetalist;
		for (int i = 0; i < N - 1; ++i) {
			for (int j = 0; j < intRes; ++j) {
				ddthetalist = ForwardDynamics(thetacurrent, dthetacurrent, taumatT.col(i), g, FtipmatT.col(i), Mlist, Glist, Slist);
				EulerStep(thetacurrent, dthetacurrent, ddthetalist, 1.0*dt / intRes);
			}
			thetamatT.col(i + 1) = thetacurrent;
			dthetamatT.col(i + 1) = dthetacurrent;
		}
		std::vector<Eigen::MatrixXd> JointTraj_ret;
		JointTraj_ret.push_back(thetamatT.transpose());
		JointTraj_ret.push_back(dthetamatT.transpose());
		return JointTraj_ret;
	}

	Eigen::VectorXd ComputedTorque(const Eigen::VectorXd& thetalist, const Eigen::VectorXd& dthetalist, const Eigen::VectorXd& eint,
		const Eigen::VectorXd& g, const std::vector<Eigen::MatrixXd>& Mlist, const std::vector<Eigen::MatrixXd>& Glist,
		const Eigen::MatrixXd& Slist, const Eigen::VectorXd& thetalistd, const Eigen::VectorXd& dthetalistd, const Eigen::VectorXd& ddthetalistd,
		double Kp, double Ki, double Kd) {

		Eigen::VectorXd e = thetalistd - thetalist;  // position err
		Eigen::VectorXd tau_feedforward = MassMatrix(thetalist, Mlist, Glist, Slist)*(Kp*e + Ki * (eint + e) + Kd * (dthetalistd - dthetalist));

		Eigen::VectorXd Ftip = Eigen::VectorXd::Zero(6);
		Eigen::VectorXd tau_inversedyn = InverseDynamics(thetalist, dthetalist, ddthetalistd, g, Ftip, Mlist, Glist, Slist);

		Eigen::VectorXd tau_computed = tau_feedforward + tau_inversedyn;
		return tau_computed;
	}

	double CubicTimeScaling(double Tf, double t) {
		double timeratio = 1.0*t / Tf;
		double st = 3 * pow(timeratio, 2) - 2 * pow(timeratio, 3);
		return st;
	}

	double QuinticTimeScaling(double Tf, double t) {
		double timeratio = 1.0*t / Tf;
		double st = 10 * pow(timeratio, 3) - 15 * pow(timeratio, 4) + 6 * pow(timeratio, 5);
		return st;
	}

	Eigen::MatrixXd JointTrajectory(const Eigen::VectorXd& thetastart, const Eigen::VectorXd& thetaend, double Tf, int N, int method) {
		double timegap = Tf / (N - 1);
		Eigen::MatrixXd trajT = Eigen::MatrixXd::Zero(thetastart.size(), N);
		double st;
		for (int i = 0; i < N; ++i) {
			if (method == 3)
				st = CubicTimeScaling(Tf, timegap*i);
			else
				st = QuinticTimeScaling(Tf, timegap*i);
			trajT.col(i) = st * thetaend + (1 - st)*thetastart;
		}
		Eigen::MatrixXd traj = trajT.transpose();
		return traj;
	}
	std::vector<Matrix<4>> ScrewTrajectory(const Matrix<4>& Xstart, const Matrix<4>& Xend, double Tf, int N, int method) {
		double timegap = Tf / (N - 1);
		std::vector<Matrix<4>> traj(N);
		double st;
		for (int i = 0; i < N; ++i) {
			if (method == 3)
				st = CubicTimeScaling(Tf, timegap*i);
			else
				st = QuinticTimeScaling(Tf, timegap*i);
			Matrix<4> Ttemp = MatrixLog6(TransInv(Xstart)*Xend);
			traj.at(i) = Xstart * MatrixExp6(Ttemp*st);
		}
		return traj;
	}

	std::vector<Matrix<4>> CartesianTrajectory(const Matrix<4>& Xstart, const Matrix<4>& Xend, double Tf, int N, int method) {
		double timegap = Tf / (N - 1);
		std::vector<Matrix<4>> traj(N);
		Matrix<3> Rstart, Rend; 
		Vector<3> pstart, pend;
		std::tie(Rstart, pstart) = TransToRp(Xstart);
		std::tie(Rend, pend) = TransToRp(Xend);

		double st;
		for (int i = 0; i < N; ++i) {
			if (method == 3)
				st = CubicTimeScaling(Tf, timegap*i);
			else
				st = QuinticTimeScaling(Tf, timegap*i);
			Matrix<3> Ri = Rstart * MatrixExp3(MatrixLog3(Rstart.transpose() * Rend)*st);
			Vector<3> pi = st*pend + (1 - st)*pstart;
			Matrix<4> traji;
			traji << Ri, pi,
				0, 0, 0, 1;
			traj.at(i) = traji;
		}
		return traj;
	}
	std::vector<Eigen::MatrixXd> SimulateControl(const Eigen::VectorXd& thetalist, const Eigen::VectorXd& dthetalist, const Eigen::VectorXd& g,
		const Eigen::MatrixXd& Ftipmat, const std::vector<Eigen::MatrixXd>& Mlist, const std::vector<Eigen::MatrixXd>& Glist,
		const Eigen::MatrixXd& Slist, const Eigen::MatrixXd& thetamatd, const Eigen::MatrixXd& dthetamatd, const Eigen::MatrixXd& ddthetamatd,
		const Eigen::VectorXd& gtilde, const std::vector<Eigen::MatrixXd>& Mtildelist, const std::vector<Eigen::MatrixXd>& Gtildelist,
		double Kp, double Ki, double Kd, double dt, int intRes) {
		Eigen::MatrixXd FtipmatT = Ftipmat.transpose();
		Eigen::MatrixXd thetamatdT = thetamatd.transpose();
		Eigen::MatrixXd dthetamatdT = dthetamatd.transpose();
		Eigen::MatrixXd ddthetamatdT = ddthetamatd.transpose();
		int m = thetamatdT.rows(); int n = thetamatdT.cols();
		Eigen::VectorXd thetacurrent = thetalist;
		Eigen::VectorXd dthetacurrent = dthetalist;
		Eigen::VectorXd eint = Eigen::VectorXd::Zero(m);
		Eigen::MatrixXd taumatT = Eigen::MatrixXd::Zero(m, n);
		Eigen::MatrixXd thetamatT = Eigen::MatrixXd::Zero(m, n);
		Eigen::VectorXd taulist;
		Eigen::VectorXd ddthetalist;
		for (int i = 0; i < n; ++i) {
			taulist = ComputedTorque(thetacurrent, dthetacurrent, eint, gtilde, Mtildelist, Gtildelist, Slist, thetamatdT.col(i),
				dthetamatdT.col(i), ddthetamatdT.col(i), Kp, Ki, Kd);
			for (int j = 0; j < intRes; ++j) {
				ddthetalist = ForwardDynamics(thetacurrent, dthetacurrent, taulist, g, FtipmatT.col(i), Mlist, Glist, Slist);
				EulerStep(thetacurrent, dthetacurrent, ddthetalist, dt / intRes);
			}
			taumatT.col(i) = taulist;
			thetamatT.col(i) = thetacurrent;
			eint += dt * (thetamatdT.col(i) - thetacurrent);
		}
		std::vector<Eigen::MatrixXd> ControlTauTraj_ret;
		ControlTauTraj_ret.push_back(taumatT.transpose());
		ControlTauTraj_ret.push_back(thetamatT.transpose());
		return ControlTauTraj_ret;
	}
}
