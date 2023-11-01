
#ifndef _VEC_H
#define _VEC_H

#include <sycl/sycl.hpp>

struct Vector;
struct Point;

struct Point
{
    SYCL_EXTERNAL Point( ) : x(0), y(0), z(0) {}
    SYCL_EXTERNAL explicit Point( const float _x, const float _y, const float _z ) : x(_x), y(_y), z(_z) {}

    float x, y, z;
};

struct Vector
{
    Vector( ) : x(0), y(0), z(0) {}
    explicit Vector( const float _x, const float _y, const float _z ) : x(_x), y(_y), z(_z) {}
    
    explicit Vector( const Point& a, const Point& b ) : x(b.x - a.x), y(b.y - a.y), z(b.z - a.z) {}

    float x, y, z;
};

inline Vector operator- ( const Point& a, const Point& b )
{
    return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Point operator* ( const float k, const Point& a )
{
    return Point(k * a.x, k * a.y, k * a.z);
}

inline Point operator* ( const Point& a, const float k )
{
    return k * a;
}

inline Point operator/ ( const Point& a, const float k )
{
    float kk= 1.f / k;
    return kk * a;
}

inline Point operator+ ( const Point& a, const Point& b )
{
    return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vector operator- ( const Vector& v )
{
    return Vector(-v.x, -v.y, -v.z);
}

inline Point operator+ ( const Point& a, const Vector& v )
{
    return Point(a.x + v.x, a.y + v.y, a.z + v.z);
}

inline Point operator+ ( const Vector& v, const Point& a )
{
    return a + v;
}

inline Point operator- ( const Vector& v, const Point& a )
{
    return a + (-v);
}

inline Point operator- ( const Point& a, const Vector& v )
{
    return a + (-v);
}

inline Vector operator+ ( const Vector& u, const Vector& v )
{
    return Vector(u.x + v.x, u.y + v.y, u.z + v.z);
}

inline Vector operator- ( const Vector& u, const Vector& v )
{
    return Vector(u.x - v.x, u.y - v.y, u.z - v.z);
}

inline Vector operator* ( const float k, const Vector& v )
{
    return Vector(k * v.x, k * v.y, k * v.z);
}

inline Vector operator* ( const Vector& v, const float k )
{
    return k * v;
}

inline Vector operator* ( const Vector& a, const Vector& b )
{
    return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline Vector operator/ (const Vector& v, const float k)
{
    float kk = 1 / k;
    return kk * v;
}

inline float length2( const Vector& v )
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline float length( const Vector& v )
{
    return sycl::sqrt(length2(v));
}

inline Vector abs(const Vector& v)
{
    return Vector(sycl::abs(v.x), sycl::abs(v.y), sycl::abs(v.z));
}

inline Vector normalize( const Vector& v )
{
    float kk= 1.0f / length(v);
    return kk * v;
}

inline Vector cross( const Vector& u, const Vector& v )
{
    return Vector(
        (u.y * v.z) - (u.z * v.y),
        (u.z * v.x) - (u.x * v.z),
        (u.x * v.y) - (u.y * v.x));
}

inline float dot( const Vector& u, const Vector& v )
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

#endif
