#ifndef _INCLUDE_LIGHTSOURCE_SPHERICAL_BACKWARD
#define _INCLUDE_LIGHTSOURCE_SPHERICAL_BACKWARD

#include "util/jacobian.glsl"
#include "lightsource/spherical/common.glsl"

ForwardRay sampleLight(
    vec3 observer, vec3 normal,
    float wavelength,
    uint mediumIdx,
    uint idx, inout uint dim
) {
    //get direction
    vec3 direction = normalize(observer - lightParams.position);

    //sample startTime
    float time = lookUp(
        Table1D(lightParams.timeProfileTable),
        random(idx, dim)
    );
    //calculate contribution
    float contrib = lightParams.contribBwd * dw_dA(lightParams.position, observer);

    //assemble ray
    return createForwardRay(
        lightParams.position,
        direction,
        wavelength,
        lightParams.mediumIdx, //do not use function parameter here. We might be in a different medium!
        time,
        contrib
    );
}

#endif
