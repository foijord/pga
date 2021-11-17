
#include <vector>
#include <ranges>
#include <numbers>
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
		union { float mx, e23; };
		union { float my, e31; };
		union { float mz, e12; };
		union { float vx, e41; };
		union { float vy, e42; };
		union { float vz, e43; };
	};


	struct plane {
		union { float x, e234; };
		union { float y, e314; };
		union { float z, e124; };
		union { float w, e321; };
	};


	template <typename T>
	inline bool operator==(const T& a, const T& b)
	{
		return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
	}


	inline bool operator==(const line& a, const line& b)
	{
		return
			a.mx == b.mx && a.my == b.my && a.mz == b.mz &&
			a.vx == b.vx && a.vy == b.vy && a.vz == b.vz;
	}


	// ===============================================================================
	// Join and Meet
	// ===============================================================================


	// line containing points p and q.
	// zero if p and q are coincident
	inline line operator^(const point& p, const point& q)
	{
		return {
			.e23 = p.y * q.z - p.z * q.y,
			.e31 = p.z * q.x - p.x * q.z,
			.e12 = p.x * q.y - p.y * q.x,
			.e41 = p.w * q.x - p.x * q.w,
			.e42 = p.w * q.y - p.y * q.w,
			.e43 = p.w * q.z - p.z * q.w,
		};
	}


	// plane containing line l and point p.
	// normal is zero if p lies in l.
	inline plane operator^(const line& l, const point& p)
	{
		return {
			.e234 = l.vy * p.z - l.vz * p.y + l.mx * p.w,
			.e314 = l.vz * p.x - l.vx * p.z + l.my * p.w,
			.e124 = l.vx * p.y - l.vy * p.x + l.mz * p.w,
			.e321 = -(l.mx * p.x + l.my * p.y + l.mz * p.z),
		};
	}


	// line where planes f and g intersect.
	// direction is zero if f and g are parallel.
	inline line operator&(const plane& f, const plane& g)
	{
		return {
			.e23 = f.x * g.w - f.w * g.x,
			.e31 = f.y * g.w - f.w * g.y,
			.e12 = f.z * g.w - f.w * g.z,
			.e41 = f.z * g.y - f.y * g.z,
			.e42 = f.x * g.z - f.z * g.x,
			.e43 = f.y * g.x - f.x * g.y,
		};
	}


	// point where line l intersects plane f.
	// weight is zero if l and f are parallel.
	inline point operator&(const line& l, const plane& f)
	{
		return {
			.e1 = l.my * f.z - l.mz * f.y + l.vx * f.w,
			.e2 = l.mz * f.x - l.mx * f.z + l.vy * f.w,
			.e3 = l.mx * f.y - l.my * f.x + l.vz * f.w,
			.e4 = -(l.vx * f.x + l.vy * f.y + l.vz * f.z),
		};
	}


	inline point operator&(const plane& f, const line& l)
	{
		return l & f;
	}

	// ===============================================================================
	// Dualization
	// ===============================================================================

	// plane at infinity
	inline plane operator!(const point& p)
	{
		return {
			.e234 = 0,
			.e314 = 0,
			.e124 = 0,
			.e321 = -p.w,
		};
	}

	// line at infinity perpendicular to line l.
	inline line operator!(const line& l)
	{
		return {
			.e23 = -l.vx,
			.e31 = -l.vy,
			.e12 = -l.vz,
			.e41 = 0,
			.e42 = 0,
			.e43 = 0,
		};
	}

	// point at infinity perpendicular to plane f.
	inline point operator!(const plane& f)
	{
		return {
			.e1 = f.x,
			.e2 = f.y,
			.e3 = f.z,
			.e4 = 0,
		};
	}
}


bool test_line_from_plane_and_point()
{
	return [&] {
		pga::point a{ 1, 0, 0, 1 };
		pga::point b{ 0, 1, 0, 1 };
		pga::point c{ 0, 0, 1, 1 };

		pga::plane f = a ^ b ^ c;
		pga::point p{ 1, 1, 1, 1 };

		// line perpendicular to plane f and passing through point p.
		// direct formula, equivalent to !f ^ p
		pga::line l = {
			.e23 = f.y * p.z - f.z * p.y,
			.e31 = f.z * p.x - f.x * p.z,
			.e12 = f.x * p.y - f.z * p.x,
			.e41 = -f.x * p.w,
			.e42 = -f.y * p.w,
			.e43 = -f.z * p.w,
		};
		// compare result of direct formula with result obtained using dualization
		return l == (!f ^ p);
	}();
}


bool test_plane_from_line_and_point()
{
	return [&] {
		pga::point a{ 1, 0, 0, 1 };
		pga::point b{ 0, 1, 1, 1 };

		pga::line l = a ^ b;
		pga::point p{ 1, 1, 1, 1 };

		// plane perpendicular to line l and containing point p.
		// direct formula, equivalent to !l ^ p
		pga::plane f = {
			.e234 = -l.vx * p.w,
			.e314 = -l.vy * p.w,
			.e124 = -l.vz * p.w,
			.e321 = l.vx * p.x + l.vy * p.y + l.vz * p.z,
		};
		// compare result of direct formula with result obtained using dualization
		return f == (!l ^ p);
	}();
}


bool test_plane_from_plane_and_line()
{
	return [&] {
		// xy plane
		pga::point a{ 1, 0, 0, 1 };
		pga::point b{ 0, 1, 0, 1 };
		pga::point c{ 0, 0, 0, 1 };
		pga::plane f = a ^ b ^ c;

		// line in xz plane
		pga::point p{ 1, 0, 0, 1 };
		pga::point q{ 0, 0, 1, 1 };
		pga::line l = p ^ q;

		// plane perpendicular to plane f and containing line l.
		// direct formula, equivalent to l ^ !f
		pga::plane g = {
			.e234 = l.vy * f.z - l.vz * f.y,
			.e314 = l.vz * f.x - l.vx * f.z,
			.e124 = l.vx * f.y - l.vy * f.x,
			.e321 = -(l.mx * f.x + l.my * f.y + l.mz * f.z),
		};
		// compare result of direct formula with result obtained using dualization
		return g == (l ^ !f);
	}();
}


bool test_project_point_onto_plane()
{
	return [&] {
		// xy plane
		pga::point a{ 1, 0, 0, 1 };
		pga::point b{ 0, 1, 0, 1 };
		pga::point c{ 0, 0, 0, 1 };
		pga::plane f = a ^ b ^ c;

		// point above xy plane
		const pga::point p{ 1, -1, 1, 1 };

		float f2 = f.x * f.x + f.y * f.y + f.z * f.z;
		float fp = f.x * p.x + f.y * p.y + f.z * p.z + f.w * p.w;
		pga::point pf = {
			 .e1 = f2 * p.x - fp * f.x,
			 .e2 = f2 * p.y - fp * f.y,
			 .e3 = f2 * p.z - fp * f.z,
			 .e4 = p.w,
		};
		// compare result of direct formula with result obtained using dualization
		return pf == (((!f) ^ p) & f);
	}();
}


bool test_project_point_onto_line()
{
	return [&] {
		pga::point a{ 1, 0, 0, 1 };
		pga::point b{ 0, 1, 0, 1 };
		pga::line l = a ^ b;

		// point above xy plane
		const pga::point p{ 1, 1, 1, 1 };

		pga::point pl = {
			 .e1 = (l.vy * l.mz - l.vz * l.my) * p.w,
			 .e2 = (l.vz * l.mx - l.vx * l.mz) * p.w,
			 .e3 = (l.vx * l.my - l.vy * l.mx) * p.w,
			 .e4 = l.vx * p.x + l.vy * p.y + l.vz * p.z + (l.vx * l.vx + l.vy * l.vy + l.vz * l.vz) * p.w,
		};
		// compare result of direct formula with result obtained using dualization
		return pl == (((!l) ^ p) & l);
	}();
}


int main(int argc, char** argv)
{
	std::vector<bool> tests{
		test_line_from_plane_and_point(),
		test_plane_from_line_and_point(),
		test_plane_from_plane_and_line(),
		test_project_point_onto_plane(),
		test_project_point_onto_line(),
	};

	std::cout << tests.size() << " tests executed." << std::endl;
	std::cout << std::ranges::count(tests, true) << " tests passed." << std::endl;
	std::cout << std::ranges::count(tests, false) << " tests failed." << std::endl;

	return EXIT_SUCCESS;
}
