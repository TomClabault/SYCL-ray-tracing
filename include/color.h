
#ifndef _COLOR_H
#define _COLOR_H

#include "vec.h"

#include <sycl/sycl.hpp>

#include <iostream>

//! \addtogroup image
///@{

//! \file
//! manipulation de couleurs

//! representation d'une couleur (rgba) transparente ou opaque.
struct Color
{
    //! constructeur par defaut.
    Color( ) : r(0.f), g(0.f), b(0.f), a(1.f) {}
    explicit Color( const std::array<float, 3>& rgb) : r(rgb[0]), g(rgb[1]), b(rgb[2]), a(1.0f) {}
    explicit Color( const float _r, const float _g, const float _y, const float _x= 1.f ) : r(_r), g(_g), b(_y), a(_x) {}
    explicit Color( const float _value ) : r(_value), g(_value), b(_value), a(1.f) {}
    
    //! cree une couleur avec les memes composantes que color, mais remplace sa composante alpha (color.r, color.g, color.b, alpha).
    Color( const Color& color, const float alpha ) : r(color.r), g(color.g), b(color.b), a(alpha) {}  // remplace alpha.

    SYCL_EXTERNAL Color& operator=(const Vector& vec);
    SYCL_EXTERNAL Color& operator+=(const Color& other);
    SYCL_EXTERNAL Color& operator*=(const Color& other);
    SYCL_EXTERNAL Color& operator*=(float k);
    SYCL_EXTERNAL Color& operator/=(const float k);

    friend std::ostream& operator << (std::ostream& os, const Color& color);

    float power( ) const;
    float max( ) const;
    
    float r, g, b, a;
};

//! utilitaire. renvoie une couleur noire.
Color Black( );
//! utilitaire. renvoie une couleur blanche.
Color White( );
//! utilitaire. renvoie une couleur rouge.
Color Red( );
//! utilitaire. renvoie une couleur verte.
Color Green( );
//! utilitaire. renvoie une couleur bleue.
Color Blue( );
//! utilitaire. renvoie une couleur jaune.
Color Yellow( );

SYCL_EXTERNAL Color exp(const Color& col);
SYCL_EXTERNAL Color pow(const Color& col, float k);

SYCL_EXTERNAL Color operator+ ( const Color& a, const Color& b );
SYCL_EXTERNAL Color operator- ( const Color& a, const Color& b );
SYCL_EXTERNAL Color operator- ( const Color& c );
SYCL_EXTERNAL Color operator* ( const Color& a, const Color& b );
SYCL_EXTERNAL Color operator* ( const Color& c, const float k );
SYCL_EXTERNAL Color operator* ( const float k, const Color& c );
SYCL_EXTERNAL Color operator/ ( const Color& a, const Color& b );
SYCL_EXTERNAL Color operator/ ( const float k, const Color& c );
SYCL_EXTERNAL Color operator/ ( const Color& c, const float k );

///@}
#endif
