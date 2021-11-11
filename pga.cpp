
#include <vector>
#include <ranges>
#include <iostream>
#include <algorithm>

namespace pga
{
	// +-------------+-----------------------+-----------------+
	// |     Type    |        Values         | Grade/Antigrade |
	// +-------------+-----------------------+-----------------+
	// |    Scalar   |          1            |   0/4    0000   |
	// +-------------+-----------------------+-----------------+
	// |             |          e1           |          1000   |
	// |   Vectors   |          e2           |   1/3    0100   |
	// |             |          e3           |          0010   |
	// |             |          e4           |          0001   |
	// +-------------+-----------------------+-----------------+
	// |             |     e23 = e2 ^ e3     |          0110   |
	// |             |     e31 = e3 ^ e1     |          1010   |
	// |  Bivectors  |     e12 = e1 ^ e2     |   2/2    1100   |
	// |             |     e43 = e4 ^ e3     |          0011   |
	// |             |     e42 = e4 ^ e2     |          0101   |
	// |             |     e41 = e4 ^ e1     |          1001   |
	// +-------------+-----------------------+-----------------+
	// |             |  e321 = e3 ^ e2 ^ e1  |          1110   |
	// | Trivectors/ |  e124 = e1 ^ e2 ^ e4  |   3/1    1101   |
	// | Antivectors |  e314 = e3 ^ e1 ^ e4  |          1011   |
	// |             |  e234 = e2 ^ e3 ^ e4  |          0111   |
	// +-------------+-----------------------+-----------------+
	// | Antiscalar  | 1 = e1 ^ e2 ^ e3 ^ e4 |   4/0    1111   |
	// +-------------+-----------------------+-----------------+


	struct point {
		union { float x, e1; };
		union { float y, e2; };
		union { float z, e3; };
		union { float w, e4; };
	};


	struct line {
		struct {
			union { float x, e41; };
			union { float y, e42; };
			union { float z, e43; };
		} v;

		struct {
			union { float x, e23; };
			union { float y, e31; };
			union { float z, e12; };
		} m;

		point weight() const
		{
			return point{
				.e1 = v.e41,
				.e2 = v.e42,
				.e3 = v.e43,
				.e4 = 0,
			};
		}
	};


	struct plane {
		union { float x, e234; };
		union { float y, e314; };
		union { float z, e124; };
		union { float w, e321; };

		point weight() const
		{
			return point{
				.e1 = e234,
				.e2 = e314,
				.e3 = e124,
				.e4 = 1,
			};
		}
	};


	struct motor {
		motor(const line& l, float phi, float d)
		{
			const auto& v = l.v;
			const auto& m = l.m;

			this->r = {
				.e41 = v.x * sin(phi),
				.e42 = v.y * sin(phi),
				.e43 = v.z * sin(phi),
				.e1234 = cos(phi),
			};

			this->u = {
				.e23 = d * v.x * cos(phi) + m.x * sin(phi),
				.e31 = d * v.y * cos(phi) + m.y * sin(phi),
				.e12 = d * v.z * cos(phi) + m.z * sin(phi),
				.scalar = -d * sin(phi),
			};
		}

		struct {
			union { float x, e41; };
			union { float y, e42; };
			union { float z, e43; };
			union { float w, e1234, antiscalar; };
		} r;

		struct {
			union { float x, e23; };
			union { float y, e31; };
			union { float z, e12; };
			union { float w, scalar; };
		} u;
	};


	// line containing points p and q.
	// zero if p and q are coincident.
	line operator^(const point& p, const point& q)
	{
		return line{
			.v = {
				.e41 = (q.x * p.w - p.x * q.w),
				.e42 = (q.y * p.w - p.y * q.w),
				.e43 = (q.z * p.w - p.z * q.w),
			},
			.m = {
				.e23 = (p.y * q.z - p.z * q.y),
				.e31 = (p.z * q.x - p.x * q.z),
				.e12 = (p.x * q.y - p.y * q.x),
			},
		};
	}


	// plane containing line l and point p.
	// normal is zero if p lies in l
	plane operator^(const line& l, const point& p)
	{
		const auto& v = l.v;
		const auto& m = l.m;

		return {
			.e234 = +(v.y * p.z - v.z * p.y + m.x * p.w),
			.e314 = +(v.z * p.x - v.x * p.z + m.y * p.w),
			.e124 = +(v.x * p.y - v.y * p.x + m.z * p.w),
			.e321 = -(m.x * p.x + m.y * p.y + m.z * p.z),
		};
	}


	// line where planes f and g intersect.
	// direction is zero if f and g are parallel.
	line operator&(const plane& f, const plane& g)
	{
		return line{
			.v = {
				.e41 = (f.z * g.y - f.y * g.z),
				.e42 = (f.x * g.z - f.z * g.x),
				.e43 = (f.y * g.x - f.x * g.y),
			},
			.m = {
				.e23 = (f.x * g.w - g.x * f.w),
				.e31 = (f.y * g.w - g.y * f.w),
				.e12 = (f.z * g.w - g.z * f.w),
			},
		};
	}


	// point where line l intersects plane f.
	// weight is zero if l and f are parallel.
	point operator&(const line& l, const plane& f)
	{
		const auto& v = l.v;
		const auto& m = l.m;

		return point{
			.e1 = +(m.y * f.z - m.z * f.y + v.x * f.w),
			.e2 = +(m.z * f.x - m.x * f.z + v.y * f.w),
			.e3 = +(m.x * f.y - m.y * f.x + v.z * f.w),
			.e4 = -(v.x * f.x + v.y * f.y + v.z * f.z),
		};
	}
}


bool line_from_points(pga::point a, pga::point b, pga::line line)
{
	return [&] {
		auto result = a ^ b;
		return true;
	}();
}


bool plane_from_points()
{
	return [&] {
		pga::point a{ 0, 0, 0, 1 };
		pga::point b{ 0, 1, 0, 1 };
		pga::point c{ 1, 1, 0, 1 };

		auto plane = a ^ b ^ c;
		auto normal = plane.weight();

		return true;
	}();
}


int main(int argc, char** argv)
{
	std::vector<bool> tests{
		line_from_points({ 2, 3, 7, 1 }, { 2, 1, 0, 1 }, { 1, 0, 0, 0, 0, 0 }),
		plane_from_points(),
	};

	std::cout << tests.size() << " tests executed." << std::endl;
	std::cout << std::ranges::count(tests, true) << " tests passed." << std::endl;
	std::cout << std::ranges::count(tests, false) << " tests failed." << std::endl;

	return EXIT_SUCCESS;
}
