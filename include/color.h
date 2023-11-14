#ifndef _COLOR_H
#define _COLOR_H

#include "vec.h"

#include <iostream>

struct Color
{
    Color( ) : r(0.f), g(0.f), b(0.f), a(1.f) {}
    explicit Color(const Vector& vec) : r(vec.x), g(vec.y), b(vec.z), a(1.0f) {}
    explicit Color( const float _r, const float _g, const float _y, const float _x= 1.f ) : r(_r), g(_g), b(_y), a(_x) {}
    explicit Color( const float _value ) : r(_value), g(_value), b(_value), a(1.f) {}
    
    inline Color& operator=(const Vector& vec)
    {
        r = vec.x;
        g = vec.y;
        b = vec.z;

        return *this;
    }

    inline Color& operator+=(const Color &other)
    {
        r += other.r;
        g += other.g;
        b += other.b;

        return *this;
    }

    inline Color& operator*=(const Color &other)
    {
        r *= other.r;
        g *= other.g;
        b *= other.b;

        return *this;
    }

    inline Color& operator*=(float k)
    {
        r *= k;
        g *= k;
        b *= k;

        return *this;
    }

    inline Color& operator/=(const float k)
    {
        r /= k;
        g /= k;
        b /= k;

        return *this;
    }

    float r, g, b, a;
};

inline Color operator+ ( const Color& a, const Color& b )
{
    return Color(a.r + b.r, a.g + b.g, a.b + b.b, a.a + b.a);
}

inline Color operator- ( const Color& c )
{
    return Color(-c.r, -c.g, -c.b, -c.a);
}

inline Color operator- ( const Color& a, const Color& b )
{
    return a + (-b);
}

inline Color operator* ( const Color& a, const Color& b )
{
    return Color(a.r * b.r, a.g * b.g, a.b * b.b, a.a * b.a);
}

inline Color operator* ( const float k, const Color& c )
{
    return Color(c.r * k, c.g * k, c.b * k, c.a * k);
}

inline Color operator* ( const Color& c, const float k )
{
    return k * c;
}

inline Color operator/ ( const Color& a, const Color& b )
{
    return Color(a.r / b.r, a.g / b.g, a.b / b.b, a.a / b.a);
}

inline Color operator/ ( const float k, const Color& c )
{
    return Color(k / c.r, k / c.g, k / c.b, k / c.a);
}

inline Color operator/ ( const Color& c, const float k )
{
    float kk= 1 / k;
    return kk * c;
}

#endif
