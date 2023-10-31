
#ifndef _VEC_H
#define _VEC_H

#include <sycl/sycl.hpp>

//! \addtogroup math
///@{

//! \file
//! operations sur points et vecteurs

//! declarations anticipees.
struct vec2;
struct vec3;
struct vec4;
struct Vector;
struct Point;

//! representation d'un point 3d.
struct Point
{
    //! constructeur par defaut.
    SYCL_EXTERNAL Point( ) : x(0), y(0), z(0) {}
    SYCL_EXTERNAL explicit Point( const float _x, const float _y, const float _z ) : x(_x), y(_y), z(_z) {}

    //! cree un point a partir des coordonnees du vecteur generique (v.x, v.y, v.z).
    SYCL_EXTERNAL Point( const vec2& v, const float z );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    SYCL_EXTERNAL Point( const vec3& v );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    SYCL_EXTERNAL Point( const vec4& v );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    //! cree un point a partir des coordonnes du vecteur (v.x, v.y, v.z).
    SYCL_EXTERNAL explicit Point( const Vector& v );   // l'implementation se trouve en fin de fichier, la structure vector n'est pas encore connue.
    
    //! renvoie la ieme composante du point.
    SYCL_EXTERNAL float operator() ( const unsigned int i ) const; // l'implementation se trouve en fin de fichier
    SYCL_EXTERNAL float& operator() ( const unsigned int i ); // l'implementation se trouve en fin de fichier
    
    Point& operator=(const Point& other)
    {
        x = other.x;
        y = other.y;
        z = other.z;

        return *this;
    }

    float x, y, z;
};

//! renvoie le point origine (0, 0, 0)
SYCL_EXTERNAL Point Origin( );

//! renvoie la distance etre 2 points.
SYCL_EXTERNAL float distance( const Point& a, const Point& b );
//! renvoie le carre de la distance etre 2 points.
SYCL_EXTERNAL float distance2( const Point& a, const Point& b );

//! renvoie le milieu du segment ab.
SYCL_EXTERNAL Point center( const Point& a, const Point& b );

//! renvoie la plus petite composante de chaque point. x, y, z= min(a.x, b.x), min(a.y, b.y), min(a.z, b.z).
SYCL_EXTERNAL Point min( const Point& a, const Point& b );
//! renvoie la plus grande composante de chaque point. x, y, z= max(a.x, b.x), max(a.y, b.y), max(a.z, b.z).
SYCL_EXTERNAL Point max( const Point& a, const Point& b );


//! representation d'un vecteur 3d.
struct Vector
{
    //! constructeur par defaut.
    Vector( ) : x(0), y(0), z(0) {}
    explicit Vector( const float _x, const float _y, const float _z ) : x(_x), y(_y), z(_z) {}
    
    //! cree le vecteur ab.
    explicit Vector( const Point& a, const Point& b ) : x(b.x - a.x), y(b.y - a.y), z(b.z - a.z) {}

    //! cree un vecteur a partir des coordonnees du vecteur generique (v.x, v.y, v.z).
    Vector( const vec3& v );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    Vector( const vec4& v );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    //! cree un vecteur a partir des coordonnes du vecteur (v.x, v.y, v.z).
    explicit Vector( const Point& a );   // l'implementation se trouve en fin de fichier.
    
    //! renvoie la ieme composante du vecteur.
    float operator() ( const unsigned int i ) const; // l'implementation se trouve en fin de fichier
    float& operator() ( const unsigned int i ); // l'implementation se trouve en fin de fichier
    
    Vector& operator=(const Vector& other)
    {
        x = other.x;
        y = other.y;
        z = other.z;

        return *this;
    }

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

//! vecteur generique, utilitaire.
struct vec2
{
    //! constructeur par defaut.
    vec2( ) : x(0), y(0) {}
    SYCL_EXTERNAL explicit vec2( const float _x, const float _y ) : x(_x), y(_y) {}
    
    //! renvoie la ieme composante du vecteur.
    SYCL_EXTERNAL float operator() ( const unsigned int i ) const { return (&x)[i]; }
    SYCL_EXTERNAL float& operator() ( const unsigned int i ) { return (&x)[i]; }

    float x, y;
};


//! vecteur generique, utilitaire.
struct vec3
{
    //! constructeur par defaut.
    SYCL_EXTERNAL vec3( ) : x(0), y(0), z(0) {}
    SYCL_EXTERNAL explicit vec3( const float _x, const float _y, const float _z ) : x(_x), y(_y), z(_z) {}
    //! constructeur par defaut.
    SYCL_EXTERNAL vec3( const vec2& a, const float _z ) : x(a.x), y(a.y), z(_z) {}

    //! cree un vecteur generique a partir des coordonnees du point a.
    SYCL_EXTERNAL vec3( const Point& a );    // l'implementation se trouve en fin de fichier.
    //! cree un vecteur generique a partir des coordonnees du vecteur v.
    SYCL_EXTERNAL vec3( const Vector& v );    // l'implementation se trouve en fin de fichier.

    //! renvoie la ieme composante du vecteur.
    SYCL_EXTERNAL float operator() ( const unsigned int i ) const { return (&x)[i]; }
    SYCL_EXTERNAL float& operator() ( const unsigned int i ) { return (&x)[i]; }
    
    float x, y, z;
};


//! vecteur generique 4d, ou 3d homogene, utilitaire.
struct vec4
{
    //! constructeur par defaut.
    SYCL_EXTERNAL vec4( ) : x(0), y(0), z(0), w(0) {}
    SYCL_EXTERNAL explicit vec4( const float _x, const float _y, const float _z, const float _w ) : x(_x), y(_y), z(_z), w(_w) {}
    //! constructeur par defaut.
    SYCL_EXTERNAL vec4( const vec2& v, const float _z= 0, const float _w= 0 ) : x(v.x), y(v.y), z(_z), w(_w) {}
    //! constructeur par defaut.
    SYCL_EXTERNAL vec4( const vec3& v, const float _w= 0 ) : x(v.x), y(v.y), z(v.z), w(_w) {}

    //! cree un vecteur generique a partir des coordonnees du point a, (a.x, a.y, a.z, 1).
    SYCL_EXTERNAL vec4( const Point& a );    // l'implementation se trouve en fin de fichier.
    //! cree un vecteur generique a partir des coordonnees du vecteur v, (v.x, v.y, v.z, 0).
    SYCL_EXTERNAL vec4( const Vector& v );    // l'implementation se trouve en fin de fichier.
    
    //! renvoie la ieme composante du vecteur.
    SYCL_EXTERNAL float operator() ( const unsigned int i ) const { return (&x)[i]; }
    SYCL_EXTERNAL float& operator() ( const unsigned int i ) { return (&x)[i]; }

    float x, y, z, w;
};


// implementation des constructeurs explicites.
inline Point::Point( const vec2& v, const float z ) : x(v.x), y(v.y), z(z) {}
inline Point::Point( const vec3& v ) : x(v.x), y(v.y), z(v.z) {}
inline Point::Point( const vec4& v ) : x(v.x), y(v.y), z(v.z) {}
inline Point::Point( const Vector& v ) : x(v.x), y(v.y), z(v.z) {}

inline Vector::Vector( const vec3& v ) : x(v.x), y(v.y), z(v.z) {}
inline Vector::Vector( const vec4& v ) : x(v.x), y(v.y), z(v.z) {}
inline Vector::Vector( const Point& a ) : x(a.x), y(a.y), z(a.z) {}

inline vec3::vec3( const Point& a ) : x(a.x), y(a.y), z(a.z) {}
inline vec3::vec3( const Vector& v ) : x(v.x), y(v.y), z(v.z) {}

inline vec4::vec4( const Point& a ) : x(a.x), y(a.y), z(a.z), w(1.f) {}
inline vec4::vec4( const Vector& v ) : x(v.x), y(v.y), z(v.z), w(0.f) {}

//
inline float Point::operator( ) ( const unsigned int i ) const { return (&x)[i]; }
inline float Vector::operator( ) ( const unsigned int i ) const { return (&x)[i]; }

inline float& Point::operator( ) ( const unsigned int i ) { return (&x)[i]; }
inline float& Vector::operator( ) ( const unsigned int i ) { return (&x)[i]; }

//
#include <iostream>

inline std::ostream& operator<<(std::ostream& o, const Point& p)
{
    o<<"p("<<p.x<<","<<p.y<<","<<p.z<<")";
    return o;
}

inline std::ostream& operator<<(std::ostream& o, const Vector& v)
{
    o<<"v("<<v.x<<","<<v.y<<","<<v.z<<")";
    return o;
}

inline const sycl::stream &operator<<(const sycl::stream& stream, const Point& p)
{
    stream << "(" << p.x << ", " << p.y << ", " << p.z << ")";

    return stream;
}

inline const sycl::stream &operator<<(const sycl::stream& stream, const Vector& v)
{
    stream << "(" << v.x << ", " << v.y << ", " << v.z << ")";

    return stream;
}

///@}
#endif
