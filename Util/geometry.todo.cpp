/*
Copyright (c) 2019, Michael Kazhdan
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
#include <cmath>
#include <SVD/SVDFit.h>
#include <SVD/MatrixMNTC.h>
#include <Util/exceptions.h>
#include "geometry.h"

namespace Util
{
	////////////////////////////
	// EulerRotationParameter //
	////////////////////////////
	Matrix3D EulerRotationParameter::operator() ( void ) const
	{
		Point3D angles = this->parameter;

		Matrix3D mx = Matrix3D();
		mx(0, 0) = 1.0;
		mx(1, 1) = cos(angles[0]);
		mx(1, 2) = -sin(angles[0]);
		mx(2, 1) = sin(angles[0]);
		mx(2, 2) = cos(angles[0]);

		Matrix3D my = Matrix3D();
		my(0, 0) = cos(angles[1]);
		my(0, 2) = sin(angles[1]);
		my(1, 1) = 1.0;
		my(2, 0) = -sin(angles[1]);
		my(2, 2) = cos(angles[1]);

		Matrix3D mz = Matrix3D();
		mz(0, 0) = cos(angles[2]);
		mz(0, 1) = -sin(angles[2]);
		mz(1, 0) = sin(angles[2]);
		mz(1, 1) = cos(angles[2]);
		mz(2, 2) = 1.0;

		return (mz * my * mx);
	}

	/////////////////////////////////
	// QuaternionRotationParameter //
	/////////////////////////////////
	Matrix3D QuaternionRotationParameter::operator()( void ) const
	{
		Quaternion q = this->parameter;
		Matrix3D result = Matrix3D();

		result(0, 0) = 1.0 - 2 * q.imag[1] * q.imag[1] - 2 * q.imag[2] * q.imag[2];
		result(0, 1) = 2 * q.imag[0] * q.imag[1] - 2 * q.real * q.imag[2];
		result(0, 2) = 2 * q.imag[0] * q.imag[2] + 2 * q.real * q.imag[1];
		result(1, 0) = 2 * q.imag[0] * q.imag[1] + 2 * q.real * q.imag[2];
		result(1, 1) = 1.0 - 2 * q.imag[0] * q.imag[0] - 2 * q.imag[2] * q.imag[2];
		result(1, 2) = 2 * q.imag[1] * q.imag[2] - 2 * q.real * q.imag[0];
		result(2, 0) = 2 * q.imag[0] * q.imag[2] - 2 * q.real * q.imag[1];
		result(2, 1) = 2 * q.imag[1] * q.imag[2] + 2 * q.real * q.imag[0];
		result(2, 2) = 1.0 - 2 * q.imag[0] * q.imag[0] - 2 * q.imag[1] * q.imag[1];

		return result;
	}
}
