#include "aSplineVec3.h"
#include <algorithm>
#include <Eigen\Dense>

#pragma warning(disable:4018)
#pragma warning(disable:4244)


ASplineVec3::ASplineVec3() : mInterpolator(new ABernsteinInterpolatorVec3())
{
}

ASplineVec3::~ASplineVec3()
{
	if (mInterpolator) delete mInterpolator;
}

void ASplineVec3::setFramerate(double fps)
{
	mInterpolator->setFramerate(fps);
}

double ASplineVec3::getFramerate() const
{
	return mInterpolator->getFramerate();
}

void ASplineVec3::setLooping(bool loop)
{
	mLooping = loop;
}

bool ASplineVec3::getLooping() const
{
	return mLooping;
}

void ASplineVec3::setInterpolationType(ASplineVec3::InterpolationType type)
{
	double fps = getFramerate();

	if (mInterpolator) { delete mInterpolator; }
	switch (type)
	{
	case LINEAR: mInterpolator = new ALinearInterpolatorVec3(); break;
	case CUBIC_BERNSTEIN: mInterpolator = new ABernsteinInterpolatorVec3(); break;
	case CUBIC_CASTELJAU: mInterpolator = new ACasteljauInterpolatorVec3(); break;
	case CUBIC_MATRIX: mInterpolator = new AMatrixInterpolatorVec3(); break;
	case CUBIC_HERMITE: mInterpolator = new AHermiteInterpolatorVec3(); break;
	case CUBIC_BSPLINE: mInterpolator = new ABSplineInterpolatorVec3(); break;
	case LINEAR_EULER: mInterpolator = new AEulerLinearInterpolatorVec3(); break;
	case CUBIC_EULER: mInterpolator = new AEulerCubicInterpolatorVec3(); break;
	};

	mInterpolator->setFramerate(fps);
	computeControlPoints();
	cacheCurve();
}

ASplineVec3::InterpolationType ASplineVec3::getInterpolationType() const
{
	return mInterpolator->getType();
}

void ASplineVec3::editKey(int keyID, const vec3& value)
{
	assert(keyID >= 0 && keyID < mKeys.size());
	mKeys[keyID].second = value;
	computeControlPoints();
	cacheCurve();
}

void ASplineVec3::editControlPoint(int ID, const vec3& value)
{
	assert(ID >= 0 && ID < mCtrlPoints.size() + 2);
	if (ID == 0)
	{
		mStartPoint = value;
		computeControlPoints(false);
	}
	else if (ID == mCtrlPoints.size() + 1)
	{
		mEndPoint = value;
		computeControlPoints(false);
	}
	else mCtrlPoints[ID - 1] = value;
	cacheCurve();
}

void ASplineVec3::appendKey(double time, const vec3& value, bool updateCurve)
{
	mKeys.push_back(Key(time, value));

	if (updateCurve)
	{
		computeControlPoints();
		cacheCurve();
	}
}

int ASplineVec3::insertKey(double time, const vec3& value, bool updateCurve)
{
	if (mKeys.size() == 0)
	{
		appendKey(time, value, updateCurve);
		return 0;
	}

	for (int i = 0; i < mKeys.size(); ++i)
	{
		assert(time != mKeys[i].first);
		if (time < mKeys[i].first)
		{
			mKeys.insert(mKeys.begin() + i, Key(time, value));
			if (updateCurve)
			{
				computeControlPoints();
				cacheCurve();
			}
			return i;
		}
	}

	// Append at the end of the curve
	appendKey(time, value, updateCurve);
	return mKeys.size() - 1;
}

void ASplineVec3::appendKey(const vec3& value, bool updateCurve)
{
	if (mKeys.size() == 0)
	{
		appendKey(0, value, updateCurve);
	}
	else
	{
		double lastT = mKeys[mKeys.size() - 1].first;
		appendKey(lastT + 1, value, updateCurve);
	}
}

void ASplineVec3::deleteKey(int keyID)
{
	assert(keyID >= 0 && keyID < mKeys.size());
	mKeys.erase(mKeys.begin() + keyID);
	computeControlPoints();
	cacheCurve();
}

vec3 ASplineVec3::getKey(int keyID) const
{
	assert(keyID >= 0 && keyID < mKeys.size());
	return mKeys[keyID].second;
}

int ASplineVec3::getNumKeys() const
{
	return mKeys.size();
}

vec3 ASplineVec3::getControlPoint(int ID) const
{
	assert(ID >= 0 && ID < mCtrlPoints.size() + 2);
	if (ID == 0) return mStartPoint;
	else if (ID == mCtrlPoints.size() + 1) return mEndPoint;
	else return mCtrlPoints[ID - 1];
}

int ASplineVec3::getNumControlPoints() const
{
	return mCtrlPoints.size() + 2; // include endpoints
}

void ASplineVec3::clear()
{
	mKeys.clear();
}

double ASplineVec3::getDuration() const
{
	return mKeys.size() == 0 ? 0 : mKeys[mKeys.size() - 1].first;
}

double ASplineVec3::getNormalizedTime(double t) const
{
	return (t / getDuration());
}

double ASplineVec3::getKeyTime(int keyID) const
{
	assert(keyID >= 0 && keyID < mKeys.size());
	return mKeys[keyID].first;
}

vec3 ASplineVec3::getValue(double t) const
{
	if (mCachedCurve.size() == 0 || mKeys.size() == 0) return vec3();
	if (t < mKeys[0].first)
		return mCachedCurve[0];
	else
		t -= mKeys[0].first;

	double dt = mInterpolator->getDeltaTime();
	int rawi = (int)(t / dt); // assumes uniform spacing
	double frac = (t - rawi * dt) / dt;

	int i = mLooping ? rawi % mCachedCurve.size() : std::min<int>(rawi, mCachedCurve.size() - 1);
	int inext = mLooping ? (i + 1) % mCachedCurve.size() : std::min<int>(i + 1, mCachedCurve.size() - 1);

	vec3 v1 = mCachedCurve[i];
	vec3 v2 = mCachedCurve[inext];
	vec3 v = v1 * (1 - frac) + v2 * frac;
	return v;
}

void ASplineVec3::cacheCurve()
{
	mInterpolator->interpolate(mKeys, mCtrlPoints, mCachedCurve);
}

void ASplineVec3::computeControlPoints(bool updateEndPoints)
{
	if (mKeys.size() >= 2 && updateEndPoints)
	{
		int totalPoints = mKeys.size();

		//If there are more than 1 interpolation point, set up the 2 end points to help determine the curve.
		//They lie on the tangent of the first and last interpolation points.
		vec3 tmp = mKeys[0].second - mKeys[1].second;
		double n = tmp.Length();
		mStartPoint = mKeys[0].second + (tmp / n) * n * 0.25; // distance to endpoint is 25% of distance between first 2 points

		tmp = mKeys[totalPoints - 1].second - mKeys[totalPoints - 2].second;
		n = tmp.Length();
		mEndPoint = mKeys[totalPoints - 1].second + (tmp / n) * n * 0.25;
	}
	mInterpolator->computeControlPoints(mKeys, mCtrlPoints, mStartPoint, mEndPoint);
}

vec3* ASplineVec3::getCachedCurveData()
{
	return mCachedCurve.data();
}

vec3* ASplineVec3::getControlPointsData()
{
	return mCtrlPoints.data();
}

int ASplineVec3::getNumCurveSegments() const
{
	return mCachedCurve.size();
}

vec3 ASplineVec3::getCurvePoint(int i) const
{
	return mCachedCurve[i];
}

//---------------------------------------------------------------------
AInterpolatorVec3::AInterpolatorVec3(ASplineVec3::InterpolationType t) : mDt(1.0 / 120.0), mType(t)
{
}

void AInterpolatorVec3::setFramerate(double fps)
{
	mDt = 1.0 / fps;
}

double AInterpolatorVec3::getFramerate() const
{
	return 1.0 / mDt;
}

double AInterpolatorVec3::getDeltaTime() const
{
	return mDt;
}

void AInterpolatorVec3::interpolate(const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints, std::vector<vec3>& curve)
{
	vec3 val = 0.0;
	double u = 0.0;

	curve.clear();

	int numSegments = keys.size() - 1;
	for (int segment = 0; segment < numSegments; segment++)
	{
		for (double t = keys[segment].first; t < keys[segment + 1].first - FLT_EPSILON; t += mDt)
		{
			// TODO: Compute u, fraction of duration between segment and segmentnext, for example,
			// u = 0.0 when t = keys[segment-1].first  
			// u = 1.0 when t = keys[segment].first
			double wholeDuration = keys[segment + 1].first - keys[segment].first;
			double timeDuration = t - keys[segment].first;
			u = timeDuration / wholeDuration;

			val = interpolateSegment(keys, ctrlPoints, segment, u);
			curve.push_back(val);
		}
	}
	// add last point
	if (keys.size() > 1)
	{
		u = 1.0;
		val = interpolateSegment(keys, ctrlPoints, numSegments - 1, u);
		curve.push_back(val);
	}
}


// Interpolate p0 and p1 so that t = 0 returns p0 and t = 1 returns p1
vec3 ALinearInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 curveValue(0, 0, 0);
	vec3 key0 = keys[segment].second;
	vec3 key1 = keys[segment + 1].second;

	// TODO: Linear interpolate between key0 and key1 so that u = 0 returns key0 and u = 1 returns key1
	curveValue = key0 * (1 - u) + key1 * u;
	return curveValue;
}

vec3 ABernsteinInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);
	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  Bernstein polynomials
	vec3 b0_poly;
	vec3 b1_poly;
	vec3 b2_poly;
	vec3 b3_poly;

	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	b0_poly = b0 * (1 - u) * (1 - u) * (1 - u);
	b1_poly = 3 * b1 * u * (1 - u) * (1 - u);
	b2_poly = 3 * b2 * u * u * (1 - u);
	b3_poly = b3 * u * u * u;

	curveValue = b0_poly + b1_poly + b2_poly + b3_poly;
	return curveValue;
}

vec3 ACasteljauInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);

	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  deCsteljau alogithm
	vec3 b01;
	vec3 b12;
	vec3 b23;
	vec3 b0112;
	vec3 b1223;

	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	b01 = b0 * (1 - u) + b1 * u;
	b12 = b1 * (1 - u) + b2 * u;
	b23 = b2 * (1 - u) + b3 * u;
	b0112 = b01 * (1 - u) + b12 * u;
	b1223 = b12 * (1 - u) + b23 * u;
	curveValue = b0112 * (1 - u) + b1223 * u;

	return curveValue;
}

vec3 AMatrixInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);

	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  matrix method f(u) = GMU
	// Hint: Using Eigen::MatrixXd data representations for a matrix operations
	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	Eigen::MatrixXd G(4, 4);
	Eigen::MatrixXd M(4, 4);
	Eigen::Vector4d U;
	Eigen::Vector4d f_u;

	G(0, 0) = b0[0], G(0, 1) = b1[0], G(0, 2) = b2[0], G(0, 3) = b3[0];
	G(1, 0) = b0[1], G(1, 1) = b1[1], G(1, 2) = b2[1], G(1, 3) = b3[1];
	G(2, 0) = b0[2], G(2, 1) = b1[2], G(2, 2) = b2[2], G(2, 3) = b3[2];
	G(3, 0) = 1, G(3, 1) = 1, G(3, 2) = 1, G(3, 3) = 1;

	M(0, 0) = 1; M(0, 1) = -3; M(0, 2) = 3; M(0, 3) = -1;
	M(1, 0) = 0; M(1, 1) = 3; M(1, 2) = -6; M(1, 3) = 3;
	M(2, 0) = 0; M(2, 1) = 0; M(2, 2) = 3; M(2, 3) = -3;
	M(3, 0) = 0; M(3, 1) = 0; M(3, 2) = 0; M(3, 3) = 1;

	U[0] = 1, U[1] = u, U[2] = u * u, U[3] = u * u * u;

	f_u = G * (M * U);
	curveValue = vec3(f_u[0] / f_u[3], f_u[1] / f_u[3], f_u[2] / f_u[3]);

	return curveValue;
}

vec3 AHermiteInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 p0;
	vec3 p1;
	vec3 q0; // slope at p0
	vec3 q1; // slope at p1
	vec3 curveValue(0, 0, 0);

	// TODO: Compute the interpolated value h(u) using a cubic Hermite polynomial  
	p0 = keys[segment].second;
	p1 = keys[segment + 1].second;
	q0 = ctrlPoints[segment];
	q1 = ctrlPoints[segment + 1];

	double p0_cubic = 2 * u * u * u - 3 * u * u + 1;
	double p1_cubic = -2 * u * u * u + 3 * u * u;
	double q0_cubic = u * u * u - 2 * u * u + u;
	double q1_cubic = u * u * u - u * u;

	vec3 h_u = p0 * p0_cubic + p1 * p1_cubic + q0 * q0_cubic + q1 * q1_cubic;
	curveValue = h_u;

	return curveValue;
}

vec3 ABSplineInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 curveValue(0, 0, 0);

	// Hint: Create a recursive helper function N(knots,n,j,t) to calculate BSpline basis function values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n =3 for cubic)
	//     j = curve interval on knot vector in which to interpolate
	//     t = time value	

	// Step 1: determine the index j
	// Step 2: compute the n nonzero Bspline Basis functions N given j
	// Step 3: get the corresponding control points from the ctrlPoints vector
	// Step 4: compute the Bspline curveValue at time t

	return curveValue;
}

void ACubicInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPoint, vec3& endPoint)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;

	for (int i = 1; i < keys.size(); i++)
	{
		vec3 b0, b1, b2, b3;
		// TODO: compute b0, b1, b2, b3
		b0 = keys[i - 1].second;
		b3 = keys[i].second;

		vec3 slope_left, slope_right;
		if (i == 1) {
			slope_left = 3 * (keys[i - 1].second - startPoint);
		}
		else {
			slope_left = (keys[i].second - keys[i - 2].second) / 2;
		}

		if (i == keys.size() - 1) {
			slope_right = 3 * (endPoint - keys[i].second);
		}
		else {
			slope_right = (keys[i + 1].second - keys[i - 1].second) / 2;
		}

		b1 = b0 + slope_left / 3.f;
		b2 = b3 - slope_right / 3.f;

		ctrlPoints.push_back(b0);
		ctrlPoints.push_back(b1);
		ctrlPoints.push_back(b2);
		ctrlPoints.push_back(b3);
	}
}

void AHermiteInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPoint, vec3& endPoint)
{
	ctrlPoints.clear();
	ctrlPoints.resize(keys.size(), vec3(0, 0, 0));
	if (keys.size() <= 1) return;

	// TODO: 
	// For each key point pi, compute the corresonding value of the slope pi_prime.
	// Hints: Using Eigen::MatrixXd for a matrix data structures, 
	// this can be accomplished by solving the system of equations AC=D for C.
	// Don't forget to save the values computed for C in ctrlPoints
	// For clamped endpoint conditions, set 1st derivative at first and last points (p0 and pm) to s0 and s1, respectively
	// For natural endpoints, set 2nd derivative at first and last points (p0 and pm) equal to 0

	// Step 1: Initialize A
	// Step 2: Initialize D
	// Step 3: Solve AC=D for C
	// Step 4: Save control points in ctrlPoints
	int keysSize;
	keysSize = keys.size();

	Eigen::MatrixXd A(keysSize, keysSize);
	Eigen::MatrixXd D(keysSize, 3);
	Eigen::MatrixXd C(keysSize, 3);

	for (int i = 0; i < keysSize; i++) {
		for (int j = 0; j < keysSize; j++) {
			double entry = 0;
			if (i - j == 1 || j - i == 1) {
				entry = 1;
			}
			else if (i == j) {
				if ((i == 0 && j == 0) || (i == keysSize - 1 && j == keysSize - 1)) {
					entry = 2;
				}
				else {
					entry = 4;
				}
			}
			A(i, j) = entry;
		}
	}

	for (int i = 0; i < keysSize; i++) {
		if (i == 0) {
			D(i, 0) = keys[i + 1].second[0] - keys[i].second[0];
			D(i, 1) = keys[i + 1].second[1] - keys[i].second[1];
			D(i, 2) = keys[i + 1].second[2] - keys[i].second[2];
		}
		else if (i == keysSize - 1) {
			D(i, 0) = keys[i].second[0] - keys[i - 1].second[0];
			D(i, 1) = keys[i].second[1] - keys[i - 1].second[1];
			D(i, 2) = keys[i].second[2] - keys[i - 1].second[2];
		}
		else {
			D(i, 0) = keys[i + 1].second[0] - keys[i - 1].second[0];
			D(i, 1) = keys[i + 1].second[1] - keys[i - 1].second[1];
			D(i, 2) = keys[i + 1].second[2] - keys[i - 1].second[2];
		}
	}

	C = A.inverse() * (3 * D);

	for (int i = 0; i < keysSize; i++) {
		ctrlPoints[i] = vec3(C(i, 0), C(i, 1), C(i, 2));
	}

}

void ABSplineInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPt, vec3& endPt)
{
	ctrlPoints.clear();
	ctrlPoints.resize(keys.size() + 2, vec3(0, 0, 0));
	if (keys.size() <= 1) return;

	// TODO:
	// Hints: 
	// 1. use Eigen::MatrixXd to calculate the control points by solving the system of equations AC=D for C

	// 2. Create a recursive helper function dN(knots,n,t,l) to calculate derivative BSpline values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n =3 for cubic)
	//     j = interval on knot vector in which to interpolate
	//     t = time value
	//     l = derivative (l = 1 => 1st derivative)

	// Step 1: Calculate knot vector using a uniform BSpline
	//         (assune knots are evenly spaced 1 apart and the start knot is at time = 0.0)

	// Step 2: Calculate A matrix  for a natural BSpline
	//         (Set 2nd derivative at t0 and tm to zero, where tm is the last point knot; m = #segments)

	// Step 3: Calculate  D matrix composed of our target points to interpolate

	// Step 4: Solve AC=D for C 

	// Step 5: save control points in ctrlPoints
}


vec3 AEulerLinearInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 curveValue(0, 0, 0);
	vec3 key0 = keys[segment].second;
	vec3 key1 = keys[segment + 1].second;
	int quotient;

	// TODO:
	// Linear interpolate between key0 and key1
	// You should convert the angles to find the shortest path for interpolation

	for (int i = 0; i < 3; i++) {
		if (key0[i] < 0) {
			quotient = (- key0[i]) / 360;
			key0[i] = key0[i] + 360.f * quotient;
			}
		if (key0[i] > 0) {
			quotient = (key0[i]) / 360;
			key0[i] = key0[i] - 360.f * quotient;
		}
		if (key1[i] < 0) {
			quotient = (-key1[i]) / 360;
			key1[i] = key1[i] + 360.f * quotient;
		}
		if (key1[i] > 0) {
			quotient = (key1[i]) / 360;
			key1[i] = key1[i] - 360.f * quotient;
		}
    }

	for (int i = 0; i < 3; i++) {
		if (key1[i] - key0[i] > 180.f) {
			key0[i] = key0[i] + 360.f;
		}
		if (key0[i] - key1[i] > 180.f) {
			key1[i] = key1[i] + 360.f;
		}
	}

	curveValue = key0 * (1 - u) + key1 * u;

	return curveValue;
}

vec3 AEulerCubicInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints, int segment, double t)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);
	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  Bernstein polynomials
	// You should convert the angles to find the shortest path for interpolation
	vec3 b0_poly;
	vec3 b1_poly;
	vec3 b2_poly;
	vec3 b3_poly;
	int quotient;
	double u = t;

	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	for (int i = 0; i < 3; i++) {
		if (b0[i] < 0) {
			quotient = (-b0[i]) / 360;
			b0[i] = b0[i] + 360.f * quotient;
		}
		if (b0[i] > 0) {
			quotient = (b0[i]) / 360;
			b0[i] = b0[i] - 360.f * quotient;
		}
		if (b1[i] < 0) {
			quotient = (-b1[i]) / 360;
			b1[i] = b1[i] + 360.f * quotient;
		}
		if (b1[i] > 0) {
			quotient = (b1[i]) / 360;
			b1[i] = b1[i] - 360.f * quotient;
		}
		if (b2[i] < 0) {
			quotient = (-b2[i]) / 360;
			b2[i] = b2[i] + 360.f * quotient;
		}
		if (b2[i] > 0) {
			quotient = (b2[i]) / 360;
			b2[i] = b2[i] - 360.f * quotient;
		}
		if (b3[i] < 0) {
			quotient = (-b3[i]) / 360;
			b3[i] = b3[i] + 360.f * quotient;
		}
		if (b3[i] > 0) {
			quotient = (b3[i]) / 360;
			b3[i] = b3[i] - 360.f * quotient;
		}
	}

	for (int i = 0; i < 3; i++) {
		if (b1[i] - b0[i] > 180.f) {
			b0[i] = b0[i] + 360.f;
		}
		if (b0[i] - b1[i] > 180.f) {
			b1[i] = b1[i] + 360.f;
		}

		if (b2[i] - b1[i] > 180.f) {
			b1[i] = b1[i] + 360.f;
		}
		if (b1[i] - b2[i] > 180.f) {
			b2[i] = b2[i] + 360.f;
		}
		if (b3[i] - b2[i] > 180.f) {
			b2[i] = b2[i] + 360.f;
		}
		if (b2[i] - b3[i] > 180.f) {
			b3[i] = b3[i] + 360.f;
		}
	}

	b0_poly = b0 * (1 - u) * (1 - u) * (1 - u);
	b1_poly = 3 * b1 * u * (1 - u) * (1 - u);
	b2_poly = 3 * b2 * u * u * (1 - u);
	b3_poly = b3 * u * u * u;

	curveValue = b0_poly + b1_poly + b2_poly + b3_poly;
	return curveValue;
}

void AEulerCubicInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints, vec3& startPoint, vec3& endPoint)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;

	// Hint: One naive way is to first convert the keys such that the differences of the x, y, z Euler angles 
	//		 between every two adjacent keys are less than 180 degrees respectively 

	for (int i = 1; i < keys.size(); i++)
	{
		vec3 b0, b1, b2, b3;

		// TODO: compute b0, b1, b2, b3
		b0 = keys[i - 1].second;
		b3 = keys[i].second;
		vec3 slope_left, slope_right;
		int quotient;

		for (int i = 0; i < 3; i++) {
			if (b0[i] < 0) {
				quotient = (-b0[i]) / 360;
				b0[i] = b0[i] + 360.f * quotient;
			}
			if (b0[i] > 0) {
				quotient = (b0[i]) / 360;
				b0[i] = b0[i] - 360.f * quotient;
			}
			if (b3[i] < 0) {
				quotient = (-b3[i]) / 360;
				b3[i] = b3[i] + 360.f * quotient;
			}
			if (b3[i] > 0) {
				quotient = (b3[i]) / 360;
				b3[i] = b3[i] - 360.f * quotient;
			}
		}

		for (int i = 0; i < 3; i++) {
			if (b3[i] - b0[i] > 180.f) {
				b0[i] = b0[i] + 360.f;
			}
			if (b0[i] - b3[i] > 180.f) {
				b3[i] = b3[i] + 360.f;
			}
		}

		if (i == 1) {
			slope_left = 3 * (keys[i - 1].second - startPoint);
		}
		else {
			slope_left = (keys[i].second - keys[i - 2].second) / 2;
		}

		if (i == keys.size() - 1) {
			slope_right = 3 * (endPoint - keys[i].second);
		}
		else {
			slope_right = (keys[i + 1].second - keys[i - 1].second) / 2;
		}

		b1 = b0 + slope_left / 3.f;
		b2 = b3 - slope_right / 3.f;

		ctrlPoints.push_back(b0);
		ctrlPoints.push_back(b1);
		ctrlPoints.push_back(b2);
		ctrlPoints.push_back(b3);
	}
}
