#ifndef _INCLUDE_LIGHTSOURCE_SPHERICAL_COMMON
#define _INCLUDE_LIGHTSOURCE_SPHERICAL_COMMON
#include "lookup.glsl"

uniform LightParams {
    vec3 position;
    uint mediumIdx;
    
    uvec2 timeProfileTable;

    float contribFwd;
    float contribBwd;
} lightParams;

#endif
