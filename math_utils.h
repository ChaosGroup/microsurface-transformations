#pragma once

// A numerical validation for the paper "Microsurface Transformations"
// by Asen Atanasov, Vladimir Koylazov, Rossen Dimov and Alexander Wilkie, EGSR'22

const float PI=3.14159265358979323846f;
inline float rnd(float a=0.0f, float b=1.0f) { return a+(b-a)*rand()/float(RAND_MAX); }

// Follow math utilities copied from V-Ray SDK.

inline float sqr(float x) {
	return x*x;
}

inline float clamp(float x, float a, float b) {
	return std::max(std::min(x, b), a);
}

inline float radicalInverse2(unsigned n) {
	n=(n<<16) | (n>>16);
	n=((n&0x00ff00ff)<<8)	| ((n&0xff00ff00)>>8);
	n=((n&0x0f0f0f0f)<<4)	| ((n&0xf0f0f0f0)>>4);
	n=((n&0x33333333)<<2)	| ((n&0xcccccccc)>>2);
	n=((n&0x55555555)<<1)	| ((n&0xaaaaaaaa)>>1);
	return (float)(double(n)/double(0x100000000));
}

class Vector {
public:
	/// The components of the vector
	union {
		struct { float x, y, z; };
		float f[3];
	}; 

	/// Default constructor - does nothing
	Vector(void) {};

	/// Initializes the vector with the given components
	Vector(float ix, float iy, float iz) { x=float(ix); y=float(iy); z=float(iz); }


	/// Initializes all components of the vector to the given value
	Vector(float f) { x=y=z=float(f); }

	/// Zeroes all components of the vector
	void makeZero(void) { x=y=z=0.0f; }

	/// Sets the components of the vector
	void set(float ix, float iy, float iz) { x=ix; y=iy; z=iz; }

	void set(int i, float v) { f[i]=v; }

	/// Adds the components of the given vector
	void operator +=(const Vector &a) { x+=a.x; y+=a.y; z+=a.z; }

	/// Subtracts the components of the given vector
	void operator -=(const Vector &a) { x-=a.x; y-=a.y; z-=a.z; }

	/// Multiplies all components by the given number
	void operator *=(float mult) { x*=float(mult); y*=float(mult); z*=float(mult); }

	/// Divides all components by the given number
	void operator /=(float f) { x/=float(f); y/=float(f); z/=float(f); }

	/// Reverses the sign of all components
	Vector operator-(void) const { return(Vector(-x, -y, -z)); }

	/// Returns the i-th component (0 for x, 1 for y, 2 for z)
	float& operator [](const int index) { return f[index]; }

	/// Returns the i-th component (0 for x, 1 for y, 2 for z) as a const
	const float& operator [](const int index) const { return f[index]; }

	/// Returns the length of the vector
	float length(void) const { return sqrtf(x*x+y*y+z*z); }

	/// Returns the squared length of the vector
	float lengthSqr(void) const { return x*x+y*y+z*z; }

	/// Makes this vector a unit vector; does not make a check for division by zero.
	/// @note This method uses multiplication by reciprocal length instead of division; if proper division is required, use getUnitVector().
	void makeNormalized(void) {
		float len1=1.0f/length();
		x*=len1;
		y*=len1;
		z*=len1;
	}
};

inline Vector operator *(const Vector &a, float f) {
	Vector res(a);
	res*=f;
	return res;
}

inline Vector operator +(const Vector &a, const Vector &b) {
	Vector res(a);
	res+=b;
	return res;
}

inline Vector operator -(const Vector &a, const Vector &b) {
	Vector res(a);
	res-=b;
	return res;
}

/// \relates Vector
/// Returns a unit vector along the direction of the given one
inline Vector normalize(const Vector &a) { return a*(1.0f/a.length()); }

/// Dot product of two vectors; intermediate calculations are double-precision.
inline double operator *(const Vector &a, const Vector &b) {
	return (double)a.x*(double)b.x+(double)a.y*(double)b.y+(double)a.z*(double)b.z;
}

inline float dotf(const Vector &a, const Vector &b) {
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

/// \relates Vector
/// Returns the mixed product of the three vectors ((a^b)*c)
inline double mixed(const Vector &a, const Vector &b, const Vector &c) {
	return
		((double)a.y*(double)b.z-(double)b.y*(double)a.z)*(double)c.x+
		((double)a.z*(double)b.x-(double)b.z*(double)a.x)*(double)c.y+
		((double)a.x*(double)b.y-(double)b.x*(double)a.y)*(double)c.z;
}

// Used by the Matrix inverse functions to compute the cross-product of two
// Vector's divided by a double-precision number. All calculations are double-precision.
inline Vector crossd(const Vector &a, const Vector &b, double D) {
	return Vector(
		((double)a.y*(double)b.z-(double)b.y*(double)a.z)/D,
		((double)a.z*(double)b.x-(double)b.z*(double)a.x)/D,
		((double)a.x*(double)b.y-(double)b.x*(double)a.y)/D
	);
}

class Matrix {
public:
	Vector f[3]; ///< The three COLUMNS of the matrix.

	/// Return the i-th column of the matrix (i=0,1,2).
	Vector& operator [](const int index) { return f[index]; }

	/// Return the i-th column of the matrix as a const object (i=0,1,2).
	const Vector& operator [](const int index) const { return f[index]; }

	/// Constructor - does not perform any initialization.
	Matrix(void) {}

	/// Constructor to either identity or zero matrix.
	/// @param i If this is 0, the matrix is initialized to the zero matrix; if this is 1, the matrix is
	/// initialized to the identity matrix.
	Matrix(int i) {
		if(i==1) makeIdentity();
		if(i==0) makeZero();
	}

	/// Constructor from the three columns of the matrix.
	/// @param a The first column.
	/// @param b The second column.
	/// @param c The third column.
	Matrix(const Vector &a, const Vector &b, const Vector &c) { f[0]=a; f[1]=b; f[2]=c; }

	/// Set the columns of the matrix.
	/// @param a The first column.
	/// @param b The second column.
	/// @param c The third column.
	void set(const Vector &a, const Vector &b, const Vector &c) { f[0]=a; f[1]=b; f[2]=c; }

	/// Set the i-th column of the matrix.
	/// @param i The index of the column to set (i=0,1,2).
	/// @param a The value for the column.
	void setCol(int i, const Vector &a) { f[i]=a; }

	/// Set the i-th row of the matrix.
	/// @param i The index of the row to set (i=0,1,2).
	/// @param a The value of the row.
	void setRow(int i, const Vector &a) { f[0][i]=a.x; f[1][i]=a.y; f[2][i]=a.z; }

	/// Make the matrix the zero matrix (all elements are zeroes).
	void makeZero(void) { memset(f, 0, sizeof(Vector)*3); }

	/// Make the matrix the identity matrix.
	void makeIdentity(void) {
		f[0].set(1.0f, 0.0f, 0.0f);
		f[1].set(0.0f, 1.0f, 0.0f);
		f[2].set(0.0f, 0.0f, 1.0f);
	}

	/// Make the matrix to be the transpose of its currrent value.
	void makeTranspose(void) {
		float t;
		t=f[0][1]; f[0][1]=f[1][0]; f[1][0]=t;
		t=f[0][2]; f[0][2]=f[2][0]; f[2][0]=t;
		t=f[1][2]; f[1][2]=f[2][1]; f[2][1]=t;
	}

	/// Make the matrix to be the inverse of its current value. Does not check if the
	/// matrix has a zero determinant.
	void makeInverse(void) {
		Matrix r;
		double D=mixed(f[0], f[1], f[2]);
		r.setRow(0, crossd(f[1], f[2], D));
		r.setRow(1, crossd(f[2], f[0], D));
		r.setRow(2, crossd(f[0], f[1], D));
		*this=r;
	}

	/// Multiply all elements in the matrix by the given number.
	void operator *=(float x) { for(int i=0; i<3; i++) f[i]*=x; }

	/// Divide all elements in the matrix by the given number. Does not perform a check for divide by zero.
	void operator /=(float x) { for(int i=0; i<3; i++) f[i]/=x; }
};

/// Vector-matric product which is the same as matrix-vector product with the transposed matrix.
inline Vector operator *(const Vector &a, const Matrix &m) { return Vector(a*m.f[0], a*m.f[1], a*m.f[2]); }

/// Matrix-vector product.
inline Vector operator *(const Matrix &m, const Vector &a) { return m.f[0]*a.x+m.f[1]*a.y+m.f[2]*a.z; }

/// Uniformly map a point from the plane to the unit (hemi)sphere using polar mapping.
/// @param[in] u The first coordinate in the plane. If 0<=u<=1, the result is in the
/// positive hemisphere (z>=0). If -1<=u<0, the result is in the negative hemisphere
/// (z<0). Otherwise the result is undefined.
/// @param[in] v The second coordinate in the plane in the range [0,1].
/// @return A point on the unit (hemi)sphere.
inline Vector getSphereDir(float u, float v) {
	float thetaSin=u;
	float thetaCos=sqrtf(1.0f-thetaSin*thetaSin);
	float phi=2.0f*PI*v;
	return Vector(cosf(phi)*thetaCos, sinf(phi)*thetaCos, thetaSin);
}