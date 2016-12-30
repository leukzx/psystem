__kernel void stepInTime(__global float4 *pos, __global float4 *vel,
                                  float dt)
{
	size_t gid = get_global_id(0);
        float4 r = pos[gid];
        float4 v = vel[gid];

        r.xyz =  r.xyz + v.xyz * dt;
        //r.x = 5.0f;
        //r.y = 2.0f;
        //r.z = 3.0f;

        pos[gid] = r;
        vel[gid] = v;
}
