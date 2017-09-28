/*
Copyright (c) 2017, Fabian Prada
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef INTERSECTABLE_OBJECT_INCLUDED
#define INTERSECTABLE_OBJECT_INCLUDED

#include <stdlib.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <Eigen/Dense>
#include <cfloat>

class IntersectableObject{
public:
	IntersectableObject(){
		worldToObjectTranslation = Eigen::Vector3f(0, 0, 0);
		worldToObjectRotation(0, 0) = worldToObjectRotation(1, 1) = worldToObjectRotation(2, 2) = 1.0;
		worldToObjectRotation(0, 1) = worldToObjectRotation(0, 2) = worldToObjectRotation(1, 0) = worldToObjectRotation(1, 2) = worldToObjectRotation(2, 0) = worldToObjectRotation(2,1) = 0.0;

		objectToWorldTranslation = Eigen::Vector3f(0, 0, 0);
		objectToWorldRotation(0, 0) = objectToWorldRotation(1, 1) = objectToWorldRotation(2, 2) = 1.0;
		objectToWorldRotation(0, 1) = objectToWorldRotation(0, 2) = objectToWorldRotation(1, 0) = objectToWorldRotation(1, 2) = objectToWorldRotation(2, 0) = objectToWorldRotation(2, 1) = 0.0;
	}

	std::vector<Eigen::Vector3f> samples;
	virtual void _InitializeSamples(const int numSamples) = 0;
	void InitializeSamples(const int numSamples){
		_InitializeSamples(numSamples);
	}

	void GetSamples(std::vector<Eigen::Vector3f> & samplesPositions){
		samplesPositions.clear(); 
		samplesPositions.resize(samples.size());
		for (int i = 0; i < samples.size(); i++){
			samplesPositions[i] = objectToWorldRotation*samples[i] + objectToWorldTranslation;
		}
	}

	Eigen::Vector3f worldToObjectTranslation;
	Eigen::Matrix3f worldToObjectRotation;
	Eigen::Vector3f objectToWorldTranslation;
	Eigen::Matrix3f objectToWorldRotation;

	virtual Eigen::Vector3f _GetCenterOfMass() const = 0;

	Eigen::Vector3f GetCenterOfMass() const{
		return objectToWorldRotation*_GetCenterOfMass() + objectToWorldTranslation;
	}

	virtual double _GetArea() const = 0;
	double GetArea() const{
		return _GetArea();
	}

	Eigen::Vector3f GetWorldPosition(const Eigen::Vector3f & objectPos) const{
		return objectToWorldRotation*objectPos + objectToWorldTranslation;
	}

	Eigen::Vector3f GetWorldDirection(const Eigen::Vector3f & objectDirection) const{
		return objectToWorldRotation*objectDirection;
	}

	void UpdateTransformations(){
		worldToObjectRotation = objectToWorldRotation.transpose();
		worldToObjectTranslation = -worldToObjectRotation*objectToWorldTranslation;
	}

	virtual bool _isInterior(const Eigen::Vector3f & p) const = 0;
	bool isInterior(const Eigen::Vector3f & p) const{
		Eigen::Vector3f tp = worldToObjectRotation*p + worldToObjectTranslation;
		return _isInterior(tp);
	}

	virtual Eigen::Vector3f _closestSurfacePoint(const Eigen::Vector3f & p) const = 0;
	
	//TODO: Add normal consistency

	Eigen::Vector3f closestSurfacePoint(const Eigen::Vector3f & p) const{
		Eigen::Vector3f tp = worldToObjectRotation*p + worldToObjectTranslation;
		Eigen::Vector3f cp = _closestSurfacePoint(tp);
		Eigen::Vector3f wcp = objectToWorldRotation*cp + objectToWorldTranslation;
		return wcp;
	}

	virtual bool _edgeIntersect(const Eigen::Vector3f & p0, const Eigen::Vector3f & p1, Eigen::Vector3f & isect_p, Eigen::Vector3f & isect_n) const = 0;

	bool edgeIntersect(const Eigen::Vector3f & p0, const Eigen::Vector3f & p1, Eigen::Vector3f & isect_p, Eigen::Vector3f & isect_n){
		Eigen::Vector3f tp0 = worldToObjectRotation*p0 + worldToObjectTranslation;
		Eigen::Vector3f tp1 = worldToObjectRotation*p1 + worldToObjectTranslation;
		bool intersect = _edgeIntersect(tp0, tp1, isect_p,isect_n);
		if (intersect){
			isect_p = objectToWorldRotation*isect_p + objectToWorldTranslation;
			isect_n = objectToWorldRotation*isect_n;
			return true;
		}
		else
			return false;
	}

	virtual bool _rayIntersect(const Eigen::Vector3f & p, const Eigen::Vector3f & d, Eigen::Vector3f & isect_p) const = 0;
	bool rayIntersect(const Eigen::Vector3f & p, const Eigen::Vector3f & d, Eigen::Vector3f & isect_p){
		Eigen::Vector3f tp = worldToObjectRotation*p + worldToObjectTranslation;
		Eigen::Vector3f td = worldToObjectRotation*d;
		bool intersect = _rayIntersect(tp, td, isect_p);
		if (intersect){
			isect_p = objectToWorldRotation*isect_p + objectToWorldTranslation;
			return true;
		}
		else
			return false;
	}
	void InteriorPointQuery(const std::vector<Eigen::Vector3f> & queryPoints, std::vector<bool> & isInteriorPoint){
		isInteriorPoint.resize(queryPoints.size());
		for (int i = 0; i < queryPoints.size(); i++){
			isInteriorPoint[i] = isInterior(queryPoints[i]);
		}
	}
	
	virtual void _Draw() = 0;
	void Draw(){
		glPushMatrix();
		glTranslatef(objectToWorldTranslation[0], objectToWorldTranslation[1], objectToWorldTranslation[2]);
		GLdouble m[4][4];
		for (int i = 0; i < 4; i++)for (int j = 0; j < 4; j++){
			if (i < 3 && j < 3){
				m[i][j] = objectToWorldRotation(j, i);
			}
			else if (i == 3 && j == 3){
				m[i][j] = 1.0;
			}
			else{
				m[i][j] = 0.0;
			}
		}
		glMultMatrixd(&m[0][0]);
		_Draw();
		glPopMatrix();
	}


	virtual Eigen::Vector3f _MinCorner() const = 0;
	Eigen::Vector3f MinCorner() const{
		return objectToWorldRotation*_MinCorner() + objectToWorldTranslation;
	}

	virtual Eigen::Vector3f _MaxCorner() const = 0;
	Eigen::Vector3f MaxCorner() const{
		return objectToWorldRotation*_MaxCorner() + objectToWorldTranslation;
	}

};

class IntersectableSphere : virtual public IntersectableObject{
public:
	IntersectableSphere(double p_radius){
		radius = p_radius;
	}
	double radius;
	bool _isInterior(const Eigen::Vector3f & p) const{
		//printf("%f %f %f \n", p[0], p[1], p[2]);
		return ((p).squaredNorm() < radius*radius);
	}
	Eigen::Vector3f _closestSurfacePoint(const Eigen::Vector3f & p) const{
		Eigen::Vector3f dir = p;
		dir *= (radius/(dir).norm());
		return  dir;
	}
	bool _edgeIntersect(const Eigen::Vector3f & p0, const Eigen::Vector3f & p1, Eigen::Vector3f & isect_p, Eigen::Vector3f & isect_n) const{
		printf("Unimplemented! \n");
		return false;
	}
	Eigen::Vector3f _GetCenterOfMass() const{
		return Eigen::Vector3f(0, 0, 0);
	}

	bool _rayIntersect(const Eigen::Vector3f & p, const Eigen::Vector3f & d, Eigen::Vector3f & isect_p) const{
		double a = d.dot(d);
		double b = 2.0 * d.dot(p);
		double c = p.dot(p) - radius*radius;
		double q = b*b - 4 * a*c;
		if (q > 0){
			double r0 = (-b - sqrt(q)) / (2.0 *a);
			if (r0 > 0){
				isect_p = p + d*r0;
				return true;
			}
			else{
				double r1 = (-b + sqrt(q)) / (2.0 *a);
				if (r1 > 0){
					isect_p = p + d*r1;
					return true;
				}
				else{
					return false;
				}
			}
		}
		else{
			return false;
		}
	}
	void _Draw(){
		glColor3f(0.f, 1.f, 0.f);
		GLUquadric* quad = gluNewQuadric();
		gluSphere(quad,radius, 100, 100);
		gluDeleteQuadric(quad);
	}
	double _GetArea() const{
		return  4 * M_PI*radius*radius;
	}
	Eigen::Vector3f _GetFirstMoment() const{
		return Eigen::Vector3f(0, 0, 0);
	}
	Eigen::Matrix3f  _GetSecondMoment() const{
		Eigen::Matrix3f secondMoment;
		secondMoment(0, 0) = secondMoment(1, 1) = secondMoment(2, 2) = radius*radius / 3.0;
		secondMoment(0, 1) = secondMoment(0, 2) = secondMoment(1, 0) = secondMoment(1, 2) = secondMoment(2, 0) = secondMoment(2, 1) = 0.0;
		return secondMoment;
	}

	Eigen::Vector3f _MinCorner() const{
		return Eigen::Vector3f(-radius, -radius, -radius);
	}
	Eigen::Vector3f _MaxCorner() const{
		return Eigen::Vector3f(radius,radius,radius);
	}

	void _InitializeSamples(const int numSamples){
		samples.resize(numSamples);
		for (int i = 0; i < numSamples; i++){
			samples[i] = Eigen::Vector3f(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
			samples[i] -= Eigen::Vector3f(0.5, 0.5, 0.5);
			samples[i] *= (radius / (samples[i]).norm());
		}
	}
};

//Cylinder aligned with the z-axis. base at z = 0, cap at z = height > 0
class IntersectableCylinder : virtual public IntersectableObject{
public:
	IntersectableCylinder(double p_radius_base, double p_radius_cap, double p_height){
		radius_base = p_radius_base;
		radius_cap = p_radius_cap;
		height = p_height;

		maxRadius2 = std::max<double>(radius_base, radius_cap)*std::max<double>(radius_base, radius_cap);

		
		slopeDirection = Eigen::Vector2f(height, radius_cap - radius_base);
		double slopeDirectionNorm2 = (slopeDirection).norm();
		normalized2SlopeDirection = slopeDirection / slopeDirectionNorm2;
		originValue = normalized2SlopeDirection[1] * radius_base;
		//originValue = Eigen::Vector2f::Dot(normalized2SlopeDirection, Eigen::Vector2f(0, radius_base));

	}
	double radius_base;
	double radius_cap;
	double height;

	double maxRadius2;
	Eigen::Vector2f slopeDirection;
	Eigen::Vector2f normalized2SlopeDirection;
	double originValue;

	Eigen::Vector3f _GetCenterOfMass() const{
		return Eigen::Vector3f(0, 0, height/2.0);
	}

	bool _isInterior(const Eigen::Vector3f & p) const{
		if (p[2] < 0 || p[2] > height){
			return false;
		}
		else{
			double distance_xy_2 = p[0] * p[0] + p[1] * p[1];
			if (distance_xy_2 > maxRadius2){
				return false;
			}
			else{
				double t = (p[2] / height);
				double z_radius = (1 - t)*radius_base + t* radius_cap;
				if (distance_xy_2 > z_radius*z_radius){
					return false;
				}
				else{
					return true;
				}
			}
		}
	}
	
	Eigen::Vector3f _closestSurfacePoint(const Eigen::Vector3f & p) const{
		double distance_xy = sqrt( p[0] * p[0] + p[1] * p[1] );
		Eigen::Vector2f planePoint(p[2], distance_xy);
		double t = planePoint.dot(normalized2SlopeDirection) - originValue;
		t = std::min<double>(std::max<double>(t, 0), 1);
		Eigen::Vector2f closestPlanePoint = slopeDirection *t;
		closestPlanePoint[1] += radius_base;
		return Eigen::Vector3f(p[0] * closestPlanePoint[1] / distance_xy, p[1] * closestPlanePoint[1] / distance_xy, closestPlanePoint[0]);
	}
	bool _edgeIntersect(const Eigen::Vector3f & p0, const Eigen::Vector3f & p1, Eigen::Vector3f & isect_p, Eigen::Vector3f & isect_n) const{
		printf("Unimplemented! \n");
		return false;
	}
	bool _rayIntersect(const Eigen::Vector3f & p, const Eigen::Vector3f & d, Eigen::Vector3f & isect_p) const{
		if (radius_base != radius_cap){
	
			//Solve for x^2 + y^2 = (m + sz)^2 , where  m = radius_base and s = (radius_cap -radius_base) / height, and (x,y,z) = p + td
			double c_lhs = p[0] * p[0] + p[1] * p[1];
			double b_lhs = 2.0*(p[0] * d[0] + p[1] * d[1]);
			double a_lhs = d[0] * d[0] + d[1] * d[1];

			double m = radius_base;
			double s = (radius_cap - radius_base) / height;

			double c_rhs = m*m + 2.0*m*s*p[2] + s * s * p[2] * p[2];
			double b_rhs = 2.0 * m * s * d[2] + s * s * 2.0 *p[2] *d[2];
			double a_rhs = s * s * d[2] * d[2];

			double a = a_lhs - a_rhs;
			double b = b_lhs - b_rhs;
			double c = c_lhs - c_rhs;

			double q = b*b - 4 * a*c;

			if (q > 0){
				double r0 = (-b - sqrt(q)) / (2.0 *a);
				if (r0 > 0){
					isect_p = p + d*r0;
					if (isect_p[2] > 0 && isect_p[2] < height) return true;
					else{
						double r1 = (-b + sqrt(q)) / (2.0 *a);
						if (r1 > 0){
							isect_p = p + d*r1;
							if (isect_p[2] > 0 && isect_p[2] < height) return true;
							else return false;
						}
						else{
							return false;
						}
					}
				}
				else{
					double r1 = (-b + sqrt(q)) / (2.0 *a);
					if (r1 > 0){
						isect_p = p + d*r1;
						if (isect_p[2] > 0 && isect_p[2] < height) return true;
						else return false;
					}
					else{
						return false;
					}
				}
			}
			else{
				return false;
			}


			//printf("Unimplemented! \n");
			//return false;
		}

		Eigen::Vector2f _p(p[0], p[1]);
		Eigen::Vector2f _d(d[0], d[1]);
		double a = _d.dot(_d);
		double b = 2.0 * _d.dot(_p);
		double c = _p.dot(_p) - radius_base*radius_base;
		double q = b*b - 4 * a*c;
		if (q > 0){
			double r0 = (-b - sqrt(q)) / (2.0 *a);
			if (r0 > 0){
				isect_p = p + d*r0;
				if (isect_p[2] > 0 && isect_p[2] < height) return true;
				else{
					double r1 = (-b + sqrt(q)) / (2.0 *a);
					if (r1 > 0){
						isect_p = p + d*r1;
						if (isect_p[2] > 0 && isect_p[2] < height) return true;
						else return false;
					}
					else{
						return false;
					}
				}
			}
			else{
				double r1 = (-b + sqrt(q)) / (2.0 *a);
				if (r1 > 0){
					isect_p = p + d*r1;
					if (isect_p[2] > 0 && isect_p[2] < height) return true;
					else return false;
				}
				else{
					return false;
				}
			}
		}
		else{
			return false;
		}

	}
	void _Draw(){
		glColor3f(0.f, 1.f, 0.f);
		GLUquadric* quad = gluNewQuadric();
		gluCylinder(quad,radius_base,radius_cap,height, 100, 100);
		gluDeleteQuadric(quad);
	}
	double _GetArea() const{
		double radius_mean = (radius_base + radius_cap) / 2.0;
		double radius_diff = radius_base - radius_cap;
		return  2 * M_PI*radius_mean*sqrt(radius_diff*radius_diff + height * height);
	}
	Eigen::Vector3f _GetFirstMoment() const{
		return Eigen::Vector3f(0, 0, height*(radius_base + 2.0 * radius_cap) / (3.0*(radius_base + radius_cap)));
	}
	Eigen::Matrix3f  _GetSecondMoment() const{
		Eigen::Matrix3f secondMoment;
		secondMoment(0, 0) = secondMoment(1, 1) = (radius_base*radius_base*radius_base + radius_base*radius_base*radius_cap + radius_base*radius_cap*radius_cap + radius_cap*radius_cap*radius_cap) / (4.0*(radius_base + radius_cap));
		secondMoment(2, 2) = height*height*(radius_base + 3.0 * radius_cap) / (6.0 *(radius_base + radius_cap));
		secondMoment(0, 1) = secondMoment(0, 2) = secondMoment(1, 0) = secondMoment(1, 2) = secondMoment(2, 0) = secondMoment(2, 1) = 0.0;
		return secondMoment;
	}

	Eigen::Vector3f _MinCorner() const{
		return Eigen::Vector3f(-radius_base, -radius_base, 0);
	}
	Eigen::Vector3f _MaxCorner() const{
		return Eigen::Vector3f(radius_cap, radius_cap, height);
	}

	void _InitializeSamples(const int numSamples){
		samples.resize(numSamples);
		for (int i = 0; i < numSamples; i++){
			double sampleHeight = double(rand())*height / double(RAND_MAX);
			double heighFactor = sampleHeight / height;
			double sampleRadius = (1.0 - heighFactor)*radius_base + heighFactor*radius_cap;
			double sampleAngle = double(rand())*2.0*M_PI/ double(RAND_MAX);
			samples[i] = Eigen::Vector3f(sampleRadius*cos(sampleAngle),sampleRadius*sin(sampleAngle),sampleHeight);
		}
	}
};

class IntersectableMultiObject : virtual public IntersectableObject{
public:
	std::vector<IntersectableObject *> objects;

	Eigen::Vector3f _GetCenterOfMass() const{
		Eigen::Vector3f cumCenter(0, 0, 0);
		double cumArea = 0;
		for (int i = 0; i < objects.size(); i++){
			double area = objects[i]->GetArea();
			cumCenter += (objects[i]->GetCenterOfMass()*area);
			cumArea += area;
		}
		return (cumCenter / cumArea);
	}

	bool _isInterior(const Eigen::Vector3f & p) const{
		for (int i = 0; i < objects.size(); i++){
			if (objects[i]->isInterior(p)){
				return true;
			}
		}
		return false;
	}
	Eigen::Vector3f _closestSurfacePoint(const Eigen::Vector3f & p) const{
		double minDistance2 = DBL_MAX;
		Eigen::Vector3f cumClosest;
		for (int i = 0; i < objects.size(); i++){
			Eigen::Vector3f candidate = objects[i]->closestSurfacePoint(p);
			double distance2 = (p - candidate).squaredNorm();
			if (distance2 < minDistance2){
				minDistance2 = distance2;
				cumClosest = candidate;
			}
		}
		return cumClosest;
	}
	bool _rayIntersect(const Eigen::Vector3f & p, const Eigen::Vector3f & d, Eigen::Vector3f & isect_p) const{
		double minDistance2 = DBL_MAX;
		bool foundIntersection = false;
		for (int i = 0; i < objects.size(); i++){
			Eigen::Vector3f candidate;
			bool isIntersection = objects[i]->rayIntersect(p, d, candidate);
			if (isIntersection){
				foundIntersection = true;
				double distance2 = (p - candidate).squaredNorm();
				if (distance2 < minDistance2){
					minDistance2 = distance2;
					isect_p = candidate;
				}
			}
		}
		return foundIntersection;
	}

	bool _edgeIntersect(const Eigen::Vector3f & p0, const Eigen::Vector3f & p1, Eigen::Vector3f & isect_p, Eigen::Vector3f & isect_n) const{
		printf("Unimplemented! \n");
		return false;
	}
	void _Draw(){
		for (int i = 0; i < objects.size(); i++){
			objects[i]->Draw();
		}
	}
	double _GetArea() const{
		double cumArea = 0.0;
		for (int i = 0; i < objects.size(); i++){
			cumArea += objects[i]->GetArea();
		}
		return cumArea;
	}
	Eigen::Vector3f _MinCorner() const{
		Eigen::Vector3f cumMin = objects[0]->MinCorner();
		for (int i = 0; i < objects.size(); i++){
			Eigen::Vector3f minCorner = objects[i]->MinCorner();
			for (int c = 0; c < 3; c++) cumMin[c] = std::min<double>(cumMin[c], minCorner[c]);
		}
		return cumMin;
	}
	Eigen::Vector3f _MaxCorner() const{
		Eigen::Vector3f cumMax = objects[0]->MaxCorner();
		for (int i = 0; i < objects.size(); i++){
			Eigen::Vector3f maxCorner = objects[i]->MaxCorner();
			for (int c = 0; c < 3; c++) cumMax[c] = std::max<double>(cumMax[c], maxCorner[c]);
		}
		return cumMax;
	}

	void _InitializeSamples(const int numSamples){
		double cumArea = 0.0;
		std::vector<double> normalizedAreas(objects.size());
		for (int i = 0; i < objects.size(); i++){
			double area = objects[i]->GetArea();
			normalizedAreas[i] = area;
			cumArea += area;
		}
		for (int i = 0; i < objects.size(); i++) normalizedAreas[i] /= cumArea;
		for (int i = 0; i < objects.size(); i++){
			//int objectSamples = normalizedAreas[i] * numSamples;
			int objectSamples = (int)(normalizedAreas[i] * double(numSamples));
			double residual = normalizedAreas[i] * double(numSamples) - double(objectSamples);
			objects[i]->_InitializeSamples(objectSamples + 1);
			for (int j = 0; j < objectSamples + 1; j++){
				if (j < objectSamples) samples.push_back(objects[i]->objectToWorldRotation*objects[i]->samples[j] + objects[i]->objectToWorldTranslation);
				else{
					double value = rand() / double(RAND_MAX);
					if (value < residual) samples.push_back(objects[i]->objectToWorldRotation*objects[i]->samples[j] + objects[i]->objectToWorldTranslation);
				}
			}
		}
	}
};

class  IntersectableCylindricalCurve : virtual public IntersectableMultiObject{
public:
	IntersectableCylindricalCurve(const std::vector<Eigen::Vector3f> & nodes, std::vector<double> & radius){
		for (int i = 0; i < nodes.size(); i++){
			IntersectableSphere * sphere = new IntersectableSphere(radius[i]);
			sphere->objectToWorldTranslation = nodes[i];
			sphere->UpdateTransformations();
			objects.push_back(sphere);
		}
		for (int i = 0; i < nodes.size() - 1; i++){
			IntersectableCylinder * cylinder = new IntersectableCylinder(radius[i], radius[i + 1], (nodes[i + 1] - nodes[i]).norm());
			
			Eigen::Vector3f d[3];
			d[0] = Eigen::Vector3f(0,0,1);
			d[1] = d[0].cross(Eigen::Vector3f(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX)));
			d[1] /= (d[1]).norm();
			d[2] = d[0].cross(d[1]);
			d[2] /= (d[2]).norm();

			Eigen::Vector3f _d[3];
			_d[0] = (nodes[i + 1] - nodes[i]) / (nodes[i + 1] - nodes[i]).norm();
			_d[1] = _d[0].cross(Eigen::Vector3f(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX)));
			_d[1] /= (_d[1]).norm();
			_d[2] = _d[0].cross(_d[1]);
			_d[2] /= (_d[2]).norm();

			Eigen::Matrix3f projection;
			for (int ki = 0; ki < 3; ki++)for (int kj = 0; kj < 3; kj++)projection(ki,kj) = d[ki][kj];
			Eigen::Matrix3f rotation;
			for (int ki = 0; ki < 3; ki++) for (int kj = 0; kj < 3; kj++)rotation(kj,ki) = _d[ki][kj];

			cylinder->objectToWorldTranslation = nodes[i];
			cylinder->objectToWorldRotation = rotation*projection;
			cylinder->UpdateTransformations();
			objects.push_back(cylinder);
		}
	}
};
//
//class LinearCilinder{
//public:
//	Eigen::Vector2f normal;
//	double crossing;
//	double GetDistance(const Eigen::Vector2f & p){
//		return 
//	}
//};
//
//class IntersectableRadialLine : virtual public IntersectableObject{
//public:
//	std::vector<Eigen::Vector2f> nodes;
//	std::vector<double> radius;
//
//	//Intersection acceleration
//	double maxRadius;
//	Eigen::Vector2f bboxMin;
//	Eigen::Vector2f bboxMax;
//	SquareMatrix<double, 2> bboxRotation;
//
//	bool _isInterior(const Eigen::Vector3f & p) const{
//		if (abs(p[2]) > maxRadius){
//			return false;
//		}
//		else{
//			Eigen::Vector2f planeProjection(p[0], p[1]);
//			Eigen::Vector2f rotatedProjection = bboxRotation*planeProjection;
//			if (rotatedProjection[0] < bboxMin[0] || rotatedProjection[1] < bboxMin[1] || rotatedProjection[0] > bboxMax[0] || rotatedProjection[1] > bboxMax[1]){
//				return false;
//			}
//			else{
//
//			}
//		}
//
//	}
//
//};


unsigned long SetVoxelKey(int i[3]){
	return	(((static_cast<unsigned long>(i[0]) << 20) & 0x3FF00000) | ((static_cast<unsigned long>(i[1]) << 10) & 0x000FFC00) | (static_cast<unsigned long>(i[2]) & 0x000003FF));
}


void GetVoxelIndices(unsigned long  key, int i[3]){
	i[0] = int((key >> 20) & 0x000003FF);
	i[1] = int((key >> 10) & 0x000003FF);
	i[2] = int((key & 0x000003FF));
}

unsigned long long SetEdgeKey(const unsigned long i0, const unsigned long i1){
	return ((((unsigned long long)i0 << 32) & 0xFFFFFFFF00000000) | ((unsigned long long)i1 & 0x00000000FFFFFFFF));
}

void GetEdgeIndices(unsigned long long key, unsigned long & i0, unsigned long & i1){
	i1 = static_cast<unsigned long>(key & 0x00000000FFFFFFFF);
	i0 = static_cast<unsigned long>((key >> 32) & 0x00000000FFFFFFFF);
}



#endif//INTERSECTABLE_OBJECT_INCLUDED