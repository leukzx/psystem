float4 get_acc(__global float4 *pos, size_t id, size_t pnum)
{
    float4 r = pos[id]; // Position of the particle
    float4 a = (float4) 0.0f; // Acceleration

    float rm = 1.0f; // Distance at which the potential reaches its minimum (F=0)
    float rm6 = rm * rm * rm * rm * rm * rm;
    float rm12 = rm6 * rm6;

    float epsilon = 1.0f;

    float4 relr = (float4) 0.0f; // Relative radius-vector of p2 to p1
    float dist2 = 0.0f; // Squared distance from p2 to p1
    float dist7 = 0.0f;
    float dist13 = 0.0f;

    // Accleration of the particle
    for (int j = 0; j < id; ++j) {
        relr.xyz = pos[j].xyz - r.xyz;
        dist2 = relr.x * relr.x + relr.y * relr.y + relr.z * relr.z;
        dist7 = dist2 * dist2 * dist2 * sqrt(dist2);
        dist13 = dist7 * dist7 / sqrt(dist2);
        a.xyz += 12
               * epsilon
               * (rm6/dist7 - rm12/dist13)
               * relr.xyz / sqrt(dist2) / r.w;
    }

    for (int j = id + 1; j < pnum; ++j) {
        relr.xyz = pos[j].xyz - r.xyz;
        dist2 = relr.x * relr.x + relr.y * relr.y + relr.z * relr.z;
        dist7 = dist2 * dist2 * dist2 * sqrt(dist2);
        dist13 = dist7 * dist7 / sqrt(dist2);
        a.xyz += 12
               * epsilon
               * (rm6/dist7 - rm12/dist13)
               * relr.xyz / sqrt(dist2) / r.w;
    }

    a.w = 0;

    return a;
}

__kernel void tsForwardEulerCL(__global float4 *pos,
                         __global float4 *vel,
                         __global float4 *acc,
                         float dt,
                         __global float4 *newPos)
{
    // This function moves particles in new position.
    size_t gid = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r = pos[gid]; // Position
    float4 v = vel[gid]; // Velocity
    //float4 a = (float4) 0.0f ; // Acceleration
    float4 a = acc[gid];

    // Particle's acceleration
    //a = get_acc(pos, gid, pnum);

    // New position and velocity
    r.xyz = r.xyz + v.xyz * dt;
    v.xyz = v.xyz + a.xyz * dt;

    newPos[gid] = r;
    vel[gid] = v;
}

__kernel void calcAccCL(__global float4 *pos, __global float4 *acc)
{
    // This function calculates acceleration of particles in current position.

    size_t gid = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r = pos[gid]; // Position of the particle

    float4 a = (float4) 0.0f; // Acceleration

    a = get_acc(pos, gid, pnum);

    acc[gid] = a;
}

__kernel void estimateDtCL(__global float4 *pos,
                           __global float4 *vel,
                           __global float4 *acc,
                           __global float *est)
{
    // Function estimates collision time of particle gid0 with
    // particle gid1.
    size_t gid0 = get_global_id(0);
    size_t gid1 = get_global_id(1);

    float4 relPos;
    float4 relVel;
    float4 relAcc;

    float d, v, a, t, est0, est1;


    if (gid0 != gid1) {
        est0 = FLT_MAX;
        est1 = FLT_MAX;

        relPos = pos[gid0] - pos[gid1];
        relPos.w = 0;
        d =  sqrt(relPos.x * relPos.x
                + relPos.y * relPos.y
                + relPos.z * relPos.z);

        relVel = vel[gid0] - vel[gid1];
        if (dot(relPos, relVel) < 0) {
            v =  sqrt(relVel.x * relVel.x
                    + relVel.y * relVel.y
                    + relVel.z * relVel.z);
            est0 = d / v;
            if (isnan(est0)) {est0 = FLT_MAX;}
        }

        relAcc = acc[gid0] - acc[gid1];
        if (dot(relPos, relAcc) < 0) {
            a =  sqrt(relAcc.x * relAcc.x
                    + relAcc.y * relAcc.y
                    + relAcc.z * relAcc.z);
            est1 = sqrt(2 * d / a);
            if (isnan(est1)) {est1 = FLT_MAX;}
        }
        //est0 = (v != 0)?(d / v):FLT_MAX;
        //est1 = (a != 0)?(sqrt(2 * d / a)):FLT_MAX;

        t = (est0 < est1)?est0:est1;
    } else {
        t = FLT_MAX;
        }

    size_t width = get_global_size(0);

    est[gid0 * width + gid1] = t;
}

__kernel void checkBoundariesCL(__global float4 *pos,
                                __global float4 *vel,
                                __global float4 *newPos,
                                __global float4 *verts,
                                            int vNum)
{
//
//https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
//
    size_t gid = get_global_id(0);
    size_t vnum = vNum;

    float4 r1 = newPos[gid];
    float4 r0 = pos[gid];

    float4 v = vel[gid];
    // Border vertices
    float4 vertex[3] = {(float4) 0.0f, (float4) 0.0f, (float4) 0.0f};

    // Some variables to get intersection point of displacement and boundary
    float4 tuv = (float4) 0.0f;
    float  det = 0.0f;
    float4 dir, edge1, edge2 , pvec, tvec, qvec;
    dir = edge1 = edge2 = pvec = tvec = qvec = (float4) 0.0f;
    float epsilon = 0.00001f;//10.0f * FLT_MIN; // Some very small number to compare derterminant with
    float4 n0 = (float4) 0.0f; // Normalized plane normal
    float4 xpoint = (float4) 0.0f; // Intersection point
    float4 dr_out = (float4) 0.0f; // Out of space displacement of particle

    // Boundary crossing flag
    // Each time there is a crossing, particle moves to new position
    // Boundaries are checked several times until there are no new crossings
    bool flag = 0;

    do{
    flag = 0;

    // Particle displacement or direction vector, or segment
    dir = r1 - r0;
    dir.w = 0;

    if (length(dir) == 0) {break;}

    for (size_t i = 0; i < vnum; i += 3) {
        // Get another 3 vertices from vertices array
        for (size_t k = 0; k < 3; ++k) {vertex[k] = verts[i+k];}

        // ?Maybe it is better to pre process edges
        edge1 = vertex[1] - vertex[0];
        edge1.w = 0;
        edge2 = vertex[2] - vertex[0];
        edge2.w = 0;

        pvec = cross(dir, edge2);
        det = dot(edge1, pvec);
        // If determinant is 0, displacement lies in plane of triangle
        // If determinat is negative, triangle is crossed from outside.
        // (Vertices of boundaries should be listed in counter clockwise order
        // from inside of area.)
        // In both cases there is no boundary crossing. Goes to next triangle
        if (det < epsilon) {continue;}

        // Calculate vector from triangle's zero vertice to displacement origin
        tvec = r0 - vertex[0];
        tvec.w = 0;

        // Calculate U parameter and test bounds
        tuv.y = dot(tvec, pvec);
        if (tuv.y < 0.0f || tuv.y > det) {continue;}

        // Prepare to test V parameter
        qvec = cross(tvec, edge1);

        // Calculate V parameter and test bounds
        tuv.z = dot(dir, qvec);
        if (tuv.z < 0.0f || (tuv.y + tuv.z) > det) {continue;}

        // Calculate T parameter and test bounds
        tuv.x  = dot (edge2, qvec);
        if (tuv.x < 0.0f || tuv.x > det) {continue;}

        // Segment intersects triangle
        // Scale parameters
        tuv /= det;

        // Set intersection flag
        flag = 1;
        // Calculate intersection point
        xpoint = r0 + dir * tuv.x;
        // Calculate plane normal of triangle
        n0 = normalize(cross(edge1, edge2));

        // New velocity
        v = v - 2 * dot(v, n0) * n0; // Mirrors v against boundary

        // Out of plane displacement
        dr_out = r1 - xpoint;
        dr_out = dr_out - 2 * dot(dr_out, n0) * n0; // Mirrors dr against boundary
        // Save old position. Calculate new position
        r0 = r1;
        r1 = xpoint + dr_out;

        // Save new values to array
        newPos[gid] = r1;
        vel[gid] = v;
        // Stop check other triangles. Do all checks again. (?)
        //break;
    }

    } while (flag); // Do loop until there are no new intersections
}
