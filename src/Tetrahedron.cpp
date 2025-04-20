#include <SYS/SYS_Types.h>
#include <SOP/SOP_Node.h>
#include <GU/GU_PrimPoly.h>

#include "Tetrahedron.h"
#include <GU/GU_PrimTetrahedron.h>

// -----------------------------------------------------
// 
// -------------------- TETRAHEDRON --------------------
// 
// -----------------------------------------------------
vec3 Tetrahedron::GetCenterOfMass()
{
	if (m_points.size() == 4)
	{
		return 0.25 * (m_points[0] + m_points[1] + m_points[2] + m_points[3]);
	}
	else
	{
		return vec3Zero;
	}
}

// Implement this however
Eigen::Vector3f Tetrahedron::GetVertexDisplacements(int a_myPointIndex) {
	const Eigen::VectorXf* displacements = m_myObj->GetDisplacementVector();

	int index = m_myObj->GetVertexIndex(m_points[a_myPointIndex]);

	Eigen::Vector3f disp;
	disp << (*displacements)(index * 3), 
		(*displacements)((index * 3) + 1), 
		(*displacements)((index * 3) + 2);

	return disp;
}

float Tetrahedron::TetrahedralVolume()
{
	Eigen::Matrix3f A;
	// Calculate the Jacobian
	// the columns of this matrix are the tetrahedrons 1st point coords, subtracted from its 2nd, 3rd, and 4th points coords.
	// the rows correspond to x, y, z
	vec3 diff10 = m_points[1] - m_points[0];
	vec3 diff20 = m_points[2] - m_points[0];
	vec3 diff30 = m_points[3] - m_points[0];

	A << diff10[0], diff20[0], diff30[0],
		diff10[1], diff20[1], diff30[1],
		diff10[2], diff20[2], diff30[2];

	// assign calculated Jacobian to struct member variable J
	J = A;

	// volume is the determinant of the above matrix (whatever the name for that matrix is) divided by 6
	m_V = abs(J.determinant() / 6.0);
	return m_V;
}

void Tetrahedron::Draw(GU_Detail* gdp, int fragmentId, GA_RWHandleS nameHandle)
{
	vec3 p0 = m_points[0];
	vec3 p1 = m_points[1];
	vec3 p2 = m_points[2];
	vec3 p3 = m_points[3];

	GU_PrimTetrahedron* tet = GU_PrimTetrahedron::build(gdp);
	nameHandle.set(tet->getMapOffset(), ("fragment_" + to_string(fragmentId)).c_str());

	GA_Offset off0 = gdp->appendPointOffset();
	GA_Offset off1 = gdp->appendPointOffset();
	GA_Offset off2 = gdp->appendPointOffset();
	GA_Offset off3 = gdp->appendPointOffset();

	gdp->setPos3(off0, UT_Vector3(p0[0], p0[1], p0[2]));
	gdp->setPos3(off1, UT_Vector3(p1[0], p1[1], p1[2]));
	gdp->setPos3(off2, UT_Vector3(p2[0], p2[1], p2[2]));
	gdp->setPos3(off3, UT_Vector3(p3[0], p3[1], p3[2]));

	tet->setPointOffset(0, off0);
	tet->setPointOffset(1, off1);
	tet->setPointOffset(2, off2);
	tet->setPointOffset(3, off3);
}

void Tetrahedron::ComputeCoefficients() {
	a = vector<float>();
	b = vector<float>();
	c = vector<float>();
	d = vector<float>();

	// Use Cramer's rule to solve system of 4 equations of form AX = B
	// equations are the shape function N = a + bx + cy + dz
	// coefficients will be calculated by determinant(Amod) / determinant(A) where A =
	//
	// | 1 x0 y0 z0 |
	// | 1 x1 y1 z1 |
	// | 1 x2 y2 z2 |
	// | 1 x3 y3 z3 |
	//
	// and Amod is A, except one of its columns (which one depends on if we're calculating for coeff a, b, c, or, d)
	// is replaced by...
	//
	// B = [1 0 0 0], [0 1 0 0], [0 0 1 0], or [0 0 0 1] (which one depends on if we're getting coeffs for 1st, 2nd, 3rd, or 4th point)
	//
	Eigen::Matrix4f A;
	A << 1.f, m_points[0][0], m_points[0][1], m_points[0][2],
		1.f, m_points[1][0], m_points[1][1], m_points[1][2],
		1.f, m_points[2][0], m_points[2][1], m_points[2][2],
		1.f, m_points[3][0], m_points[3][1], m_points[3][2];

	// iterate to calculate a, b, c, d, for each point in m_points. i determines which point, inner loop j determines a, b, c, d
	// the location of 1 in vector B corresponds to the first, second, third, fourth points.
	for (int i = 0; i < 4; i++)
	{
		Eigen::Vector4f B = { 0, 0, 0, 0 };
		B[i] = 1;

		// iterate across matrix A, substituting each column with B. first column gets us a, 2nd b, 3rd c, etc.
		for (int j = 0; j < 4; j++)
		{
			Eigen::Matrix4f A_copy = A;
			A_copy.col(j) = B;

			switch (j)
			{
			case 0:
				a.push_back(A_copy.determinant() / A.determinant());
				break;
			case 1:
				b.push_back(A_copy.determinant() / A.determinant());
				break;
			case 2:
				c.push_back(A_copy.determinant() / A.determinant());
				break;
			case 3:
				d.push_back(A_copy.determinant() / A.determinant());
				break;
			}
		}
	}

	// with the coefficients, we can compute the strain-displacement matrix of the tetrahedron
	// it is 12 x 6 and generally of the form
	// | b1 00 00 b2 00 00 b3 00 00 b4 00 00 |
	// | 00 c1 00 00 c2 00 00 c3 00 00 c4 00 |
	// | 00 00 d1 00 00 d2 00 00 d3 00 00 d4 |
	// | c1 b1 00 c2 b2 00 c3 b3 00 c4 b4 00 |
	// | 00 d1 c1 00 d2 c2 00 d3 c3 00 d4 c4 |
	// | d1 00 b1 d2 00 b2 d3 00 b3 d4 00 b4 |
	Eigen::MatrixXf Strain_Disp{
		{ b[0], 0.00, 0.00, b[1], 0.00, 0.00, b[2], 0.00, 0.00, b[3], 0.00, 0.00 },
		{ 0.00, c[0], 0.00, 0.00, c[1], 0.00, 0.00, c[2], 0.00, 0.00, c[3], 0.00 },
		{ 0.00, 0.00, d[0], 0.00, 0.00, d[1], 0.00, 0.00, d[2], 0.00, 0.00, d[3] },
		{ c[0], b[0], 0.00, c[1], b[1], 0.00, c[2], b[2], 0.00, c[3], b[3], 0.00 },
		{ 0.00, d[0], c[0], 0.00, d[1], c[1], 0.00, d[2], c[2], 0.00, d[3], c[3] },
		{ d[0], 0.00, b[0], d[1], 0.00, b[1], d[2], 0.00, b[2], d[3], 0.00, b[3] }
	};

	B = Strain_Disp;

	// Calculate the stiffness matrix for this tetrahedron
	// K_e = V * B^T * D * B
	// V - Volume
	// B - Strain-displacement matrix
	// D - Material stiffness matrix
	K_e = m_V * B.transpose() * (*m_myObj->GetMaterialMatrix()) * B;
}

void Tetrahedron::ComputeStrainTensor() {
	// step 1: compute inverse jacobian of the undeformed tetrahedron
	Eigen::Matrix3f J_inv = J.inverse();

	// step 2: compute shape function gradient matrix. Local gradients of shape functions
	// are constant for a standard tetrahedron
	Eigen::Matrix<float, 3, 4> localGrad;
	localGrad << -1.0, 1.0, 0.0, 0.0,
		-1.0, 0.0, 1.0, 0.0,
		-1.0, 0.0, 0.0, 1.0;

	localGrad = J_inv * localGrad;

	// step 3: compute deformation gradient F for the tetrahedron
	// this is the 3D Identity matrix plus the sum of the product of each vertex displacement 
	// and the sequential columns (transposed to rows) of the local shape function gradient matrix
	Eigen::Matrix3f F = Eigen::Matrix3f::Identity(); // I think F is the Deformation gradient
	for (int i = 0; i < 4; ++i) {
		Eigen::Vector3f delta = GetVertexDisplacements(i);

		// this is adding the product of a 3x1 * 1x3, which gives a 3x3 matrix.
		F += delta * localGrad.col(i).transpose();
	}

	// Step 4: compute Green-Lagrange Strain Tensor
	// Strain Tensor = (1/2) * ((F^T * F) - I)
	m_StrainTensor = 0.5 * ((F.transpose() * F) - Eigen::Matrix3f::Identity());
}

/// <summary>
/// To compute Stress Tensor, we convert Strain Tensor into a 6x1 vector in Voigt notation
/// This is multiplied onto the already created Material Stiffness Matrix to create a 6x1 Stress Tensor in Voigt notation
/// We then convert this to a 3x3 matrix and store it.
/// </summary>
void Tetrahedron::ComputeStrainEnergy() {
	// create strain vector
	Eigen::Matrix<float, 6, 1> strainVec;
	strainVec << m_StrainTensor(0, 0),		// ε_xx
		m_StrainTensor(1, 1),		// ε_yy
		m_StrainTensor(2, 2),		// ε_zz
		2 * m_StrainTensor(1, 2),	// γ_yz (2 * ε_yz)
		2 * m_StrainTensor(0, 2),	// γ_xz (2 * ε_xz)
		2 * m_StrainTensor(0, 1);	// γ_xy (2 * ε_xy)

	// multiply onto Element Stiffness Matrix to get stress vector
	// INSTEAD OF K_e WE MAY NEED TO USE MATERIAL MATRIX FROM PARENT
	Eigen::Matrix<float, 6, 1> stressVec = K_e * strainVec;

	// energy: W = V * (0.5 * dot(strain_voigt, stress_voigt))
	float strainEnergy = 0.5 * strainVec.dot(stressVec);
	m_W = m_V * strainEnergy;
}