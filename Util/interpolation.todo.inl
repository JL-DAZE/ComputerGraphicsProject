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

#include <math.h>
#include <Util/exceptions.h>

namespace Util
{
	///////////////////
	// Interpolation //
	///////////////////
	template< typename SampleType >
	SampleType Interpolation::Sample( const std::vector< SampleType > &samples , double t , int interpolationType )
	{
		switch (interpolationType)
		{
		case NEAREST:
		{
			t *= samples.size();
			int it1 = (int)floor(t);
			int it2 = (it1 + 1) % samples.size();
			t -= it1;
			if (t < 0.5) return samples[it1];
			else        return samples[it2];
			break;
		}
		case LINEAR:
		{
			t = t * samples.size();
			int it1 = static_cast<int>(floor(t));
			int it2 = (it1 + 1) % samples.size();
			t = t - it1;
			return (1 - t) * samples[it1] + t * samples[it2];
			break;
		}
		case CATMULL_ROM:
		{
			t = t * samples.size();
			int it1 = static_cast<int>(floor(t));
			int it2 = (it1 + 1) % samples.size();
			int it3 = (it1 + 2) % samples.size();
			int it0 = (it1 - 1) % samples.size();
			t = t - it1;

			double s = 0.5;
			double h0 = 2 * pow(t, 3) - 3 * pow(t, 2) + 1;
			double h1 = -2 * pow(t, 3) + 3 * pow(t, 2);
			double h2 = pow(t, 3) - 2 * pow(t, 2) + t;
			double h3 = pow(t, 3) - pow(t, 2);
			SampleType p0 = samples[it1];
			SampleType p1 = samples[it2];
			SampleType t0 = s * (samples[it2] - samples[it0]);
			SampleType t1 = s * (samples[it3] - samples[it1]);

			return (h0 * p0 + h1 * p1 + h2 * t0 + h3 * t1);
			break;
		}
		case UNIFORM_CUBIC_B_SPLINE:
		{
			t = t * samples.size();
			int it1 = static_cast<int>(floor(t));
			int it2 = (it1 + 1) % samples.size();
			int it3 = (it1 + 2) % samples.size();
			int it0 = (it1 - 1) % samples.size();
			t = t - it1;

			double s = 0.5;
			double h0 = 2 * pow(t, 3) - 3 * pow(t, 2) + 1;
			double h1 = -2 * pow(t, 3) + 3 * pow(t, 2);
			double h2 = pow(t, 3) - 2 * pow(t, 2) + t;
			double h3 = pow(t, 3) - pow(t, 2);
			SampleType p0 = (samples[it0] + 4 * samples[it1] + samples[it2]) / 6;
			SampleType p1 = (samples[it1] + 4 * samples[it2] + samples[it3]) / 6;
			SampleType t0 = s * (samples[it2] - samples[it0]);
			SampleType t1 = s * (samples[it3] - samples[it1]);

			return (h0 * p0 + h1 * p1 + h2 * t0 + h3 * t1);
			break;
		}
		}
	}
}